// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/representations/aerosol_property.hpp>
#include <miam/util/error.hpp>
#include <miam/util/miam_exception.hpp>

#include <micm/system/phase.hpp>

#include <cmath>
#include <map>
#include <numbers>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace miam
{
  /// @brief Two moment log-normal particle size distribution representation
  /// @details Represents a two moment log-normal distribution for aerosol or cloud particle size distributions.
  ///          Characterized by number concentration, geometric mean radius, and geometric standard deviation.
  class TwoMomentMode
  {
   public:
    TwoMomentMode() = delete;

    TwoMomentMode(const std::string& prefix, const std::vector<micm::Phase>& phases)
        : prefix_(prefix),
          phases_(phases),
          default_geometric_standard_deviation_(1.0)
    {
    }

    TwoMomentMode(
        const std::string& prefix,
        const std::vector<micm::Phase>& phases,
        const double geometric_standard_deviation)
        : prefix_(prefix),
          phases_(phases),
          default_geometric_standard_deviation_(geometric_standard_deviation)
    {
    }

    std::tuple<std::size_t, std::size_t> StateSize() const
    {
      std::size_t size = 0;
      for (const auto& phase : phases_)
      {
        size += phase.StateSize();
      }
      size++;              // Number concentration
      return { size, 1 };  // One parameter: geometric standard deviation
    }

    std::set<std::string> StateVariableNames() const
    {
      std::set<std::string> names;
      for (const auto& phase : phases_)
      {
        for (const auto& species : phase.UniqueNames())
        {
          names.insert(prefix_ + "." + species);
        }
      }
      names.insert(prefix_ + ".NUMBER_CONCENTRATION");
      return names;
    }

    std::set<std::string> StateParameterNames() const
    {
      std::set<std::string> names;
      names.insert(prefix_ + ".GEOMETRIC_STANDARD_DEVIATION");
      return names;
    }

    std::string Species(const micm::Phase& phase, const micm::Species& species) const
    {
      return prefix_ + "." + phase.name_ + "." + species.name_;
    }

    std::map<std::string, double> DefaultParameters() const
    {
      return { { prefix_ + ".GEOMETRIC_STANDARD_DEVIATION", default_geometric_standard_deviation_ } };
    }

    std::string NumberConcentration() const
    {
      return prefix_ + ".NUMBER_CONCENTRATION";
    }

    std::string GeometricStandardDeviation() const
    {
      return prefix_ + ".GEOMETRIC_STANDARD_DEVIATION";
    }

    void SetDefaultParameters(auto& state) const
    {
      auto gsd_it = state.custom_rate_parameter_map_.find(GeometricStandardDeviation());
      if (gsd_it == state.custom_rate_parameter_map_.end())
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_STATE_PARAMETER,
            "TwoMomentMode::SetDefaultParameters: Geometric standard deviation parameter not found in state.");
      }
      for (std::size_t i_cell = 0; i_cell < state.variables_.NumRows(); ++i_cell)
      {
        state.custom_rate_parameters_[i_cell][gsd_it->second] = default_geometric_standard_deviation_;
      }
    }

    std::map<std::string, std::size_t> NumPhaseInstances() const
    {
      std::map<std::string, std::size_t> num_instances;
      for (const auto& phase : phases_)
      {
        num_instances[phase.name_] = 1;  // Two moment representation has one instance per phase
      }
      return num_instances;
    }

    /// @brief Returns a map of phase names to sets of state variable prefixes associated with that phase
    ///        The prefix does not include the phase or species names, and each prefix must be unique across all
    ///        representations.
    /// @return Map of phase names to sets of state variable prefixes
    std::map<std::string, std::set<std::string>> PhaseStatePrefixes() const
    {
      std::map<std::string, std::set<std::string>> phase_prefixes;
      for (const auto& phase : phases_)
      {
        phase_prefixes[phase.name_].insert(prefix_);
      }
      return phase_prefixes;
    }

    /// @brief Returns a provider for the requested aerosol property
    template<typename DenseMatrixPolicy>
    AerosolPropertyProvider<DenseMatrixPolicy> GetPropertyProvider(
        AerosolProperty property,
        const auto& state_parameter_indices,
        const auto& state_variable_indices,
        const std::string& target_phase_name = "") const
    {
      AerosolPropertyProvider<DenseMatrixPolicy> provider;

      switch (property)
      {
        case AerosolProperty::EffectiveRadius:
        {
          // r_eff = (3*V_total/(4pi*N))^(1/3) * exp(2.5*ln^2(GSD))
          // Depends on all species variables and NUMBER_CONCENTRATION
          std::size_t gsd_idx = state_parameter_indices.at(GeometricStandardDeviation());
          std::size_t nc_var_idx = state_variable_indices.at(NumberConcentration());
          std::vector<std::size_t> species_indices;
          std::vector<double> molar_volumes;
          for (const auto& phase : phases_)
            for (const auto& ps : phase.phase_species_)
              if (!ps.species_.IsParameterized())
              {
                species_indices.push_back(state_variable_indices.at(prefix_ + "." + phase.name_ + "." + ps.species_.name_));
                molar_volumes.push_back(
                    ps.species_.GetProperty<double>("molecular weight [kg mol-1]") /
                    ps.species_.GetProperty<double>("density [kg m-3]"));
              }
          // dependent_variable_indices: species first, N last
          provider.dependent_variable_indices = species_indices;
          provider.dependent_variable_indices.push_back(nc_var_idx);
          auto storage = std::make_shared<std::pair<
              typename DenseMatrixPolicy::template VectorType<std::size_t>,
              typename DenseMatrixPolicy::template VectorType<double>>>(species_indices, molar_volumes);
          storage->first.CopyToDevice();
          storage->second.CopyToDevice();
          const auto species_view = storage->first.GetView();
          const auto volumes_view = storage->second.GetView();
          DenseMatrixPolicy example_params{ 1, state_parameter_indices.size(), 0.0 };
          DenseMatrixPolicy example_vars{ 1, state_variable_indices.size(), 0.0 };
          DenseMatrixPolicy example_result{ 1, 1, 0.0 };
          DenseMatrixPolicy example_partials{ 1, provider.dependent_variable_indices.size(), 0.0 };
          auto value_function =
            DenseMatrixPolicy::Function(
                MICM_LAMBDA(
                    const typename DenseMatrixPolicy::ConstViewType& params_view,
                    const typename DenseMatrixPolicy::ConstViewType& vars_view,
                    const typename DenseMatrixPolicy::ViewType& result_view) {
                  auto r = result_view.GetColumnView(0);
                  params_view.ForEachRowStrict([](double& v) { v = 0.0; }, r);
                  const std::size_t n = species_view.size();
                  for (std::size_t k = 0; k < n; ++k)
                  {
                    const double mv = volumes_view[k];
                    params_view.ForEachRowStrict(
                        [mv](const double& c, double& V) { V += c * mv; }, vars_view.GetConstColumnView(species_view[k]), r);
                  }
                  params_view.ForEachRowStrict(
                      [](const double& gsd, const double& nc, double& r_eff)
                      {
                        const double ln_gsd = std::log(gsd);
                        const double V_mean = r_eff / nc;  // r_eff held V_total
                        const double r_mean = std::cbrt(3.0 * V_mean / (4.0 * std::numbers::pi));
                        r_eff = r_mean * std::exp(2.5 * ln_gsd * ln_gsd);
                      },
                      params_view.GetConstColumnView(gsd_idx),
                      vars_view.GetConstColumnView(nc_var_idx),
                      r);
                },
                example_params,
                example_vars,
                example_result);
          provider.ComputeValue = [storage, value_function](
                                                    const DenseMatrixPolicy& params,
                                                    const DenseMatrixPolicy& vars,
                                      DenseMatrixPolicy& result) mutable { value_function(params, vars, result); };
          auto value_and_derivatives_function =
            DenseMatrixPolicy::Function(
                MICM_LAMBDA(
                    const typename DenseMatrixPolicy::ConstViewType& params_view,
                    const typename DenseMatrixPolicy::ConstViewType& vars_view,
                    const typename DenseMatrixPolicy::ViewType& result_view,
                    const typename DenseMatrixPolicy::ViewType& partials_view) {
                  // Accumulate V_total
                  auto V_total = result_view.GetRowVariable();
                  params_view.ForEachRowStrict([](double& v) { v = 0.0; }, V_total);
                  const std::size_t n = species_view.size();
                  for (std::size_t k = 0; k < n; ++k)
                  {
                    const double mv = volumes_view[k];
                    params_view.ForEachRowStrict(
                        [mv](const double& c, double& V) { V += c * mv; },
                        vars_view.GetConstColumnView(species_view[k]),
                        V_total);
                  }
                  // Partials w.r.t. species: dr_eff/d[species_k] = r_eff / (3 * V_total) * molar_volume_k
                  for (std::size_t k = 0; k < n; ++k)
                  {
                    const double mv = volumes_view[k];
                    params_view.ForEachRowStrict(
                        [mv](const double& gsd, const double& nc, const double& V_total_row, double& dr)
                        {
                          const double ln_gsd = std::log(gsd);
                          const double V_mean = V_total_row / nc;
                          const double r_mean = std::cbrt(3.0 * V_mean / (4.0 * std::numbers::pi));
                          const double r_eff = r_mean * std::exp(2.5 * ln_gsd * ln_gsd);
                          dr = r_eff * mv / (3.0 * V_total_row);
                        },
                        params_view.GetConstColumnView(gsd_idx),
                        vars_view.GetConstColumnView(nc_var_idx),
                        V_total,
                        partials_view.GetColumnView(k));
                  }
                  // Partial w.r.t. N: dr_eff/dN = -r_eff / (3 * N)
                  params_view.ForEachRowStrict(
                      [](const double& gsd, const double& nc, const double& V_total_row, double& dr_dN)
                      {
                        const double ln_gsd = std::log(gsd);
                        const double V_mean = V_total_row / nc;
                        const double r_mean = std::cbrt(3.0 * V_mean / (4.0 * std::numbers::pi));
                        const double r_eff = r_mean * std::exp(2.5 * ln_gsd * ln_gsd);
                        dr_dN = -r_eff / (3.0 * nc);
                      },
                      params_view.GetConstColumnView(gsd_idx),
                      vars_view.GetConstColumnView(nc_var_idx),
                      V_total,
                      partials_view.GetColumnView(n));
                  params_view.ForEachRowStrict(
                      [](const double& gsd, const double& nc, const double& V_total_row, double& r_out)
                      {
                        const double ln_gsd = std::log(gsd);
                        const double V_mean = V_total_row / nc;
                        const double r_mean = std::cbrt(3.0 * V_mean / (4.0 * std::numbers::pi));
                        r_out = r_mean * std::exp(2.5 * ln_gsd * ln_gsd);
                      },
                      params_view.GetConstColumnView(gsd_idx),
                      vars_view.GetConstColumnView(nc_var_idx),
                      V_total,
                      result_view.GetColumnView(0));
                },
                example_params,
                example_vars,
                example_result,
                example_partials);
          provider.ComputeValueAndDerivatives = [storage, value_and_derivatives_function](
                                                    const DenseMatrixPolicy& params,
                                                    const DenseMatrixPolicy& vars,
                                                    DenseMatrixPolicy& result,
                                                    DenseMatrixPolicy& partials) mutable
          { value_and_derivatives_function(params, vars, result, partials); };
          break;
        }
        case AerosolProperty::NumberConcentration:
        {
          // N = state_variables[NUMBER_CONCENTRATION], d(N)/d(N) = 1
          std::size_t nc_var_idx = state_variable_indices.at(NumberConcentration());
          provider.dependent_variable_indices = { nc_var_idx };
          DenseMatrixPolicy example_params{ 1, state_parameter_indices.size(), 0.0 };
          DenseMatrixPolicy example_vars{ 1, state_variable_indices.size(), 0.0 };
          DenseMatrixPolicy example_result{ 1, 1, 0.0 };
          DenseMatrixPolicy example_partials{ 1, 1, 0.0 };
          auto value_function =
            DenseMatrixPolicy::Function(
                MICM_LAMBDA(
                    const typename DenseMatrixPolicy::ConstViewType& params_view,
                    const typename DenseMatrixPolicy::ConstViewType& vars_view,
                    const typename DenseMatrixPolicy::ViewType& result_view) {
                  params_view.ForEachRowStrict(
                      [](const double& nc, double& N) { N = nc; },
                      vars_view.GetConstColumnView(nc_var_idx),
                      result_view.GetColumnView(0));
                },
                example_params,
                example_vars,
                example_result);
          provider.ComputeValue = [value_function](
                                                    const DenseMatrixPolicy& params,
                                                    const DenseMatrixPolicy& vars,
                                      DenseMatrixPolicy& result) mutable { value_function(params, vars, result); };
          auto value_and_derivatives_function =
            DenseMatrixPolicy::Function(
                MICM_LAMBDA(
                    const typename DenseMatrixPolicy::ConstViewType& params_view,
                    const typename DenseMatrixPolicy::ConstViewType& vars_view,
                    const typename DenseMatrixPolicy::ViewType& result_view,
                    const typename DenseMatrixPolicy::ViewType& partials_view) {
                  params_view.ForEachRowStrict(
                      [](const double& nc, double& N, double& dN_dN)
                      {
                        N = nc;
                        dN_dN = 1.0;
                      },
                      vars_view.GetConstColumnView(nc_var_idx),
                      result_view.GetColumnView(0),
                      partials_view.GetColumnView(0));
                },
                example_params,
                example_vars,
                example_result,
                example_partials);
          provider.ComputeValueAndDerivatives = [value_and_derivatives_function](
                                                    const DenseMatrixPolicy& params,
                                                    const DenseMatrixPolicy& vars,
                                                    DenseMatrixPolicy& result,
                                                    DenseMatrixPolicy& partials) mutable
          { value_and_derivatives_function(params, vars, result, partials); };
          break;
        }
        case AerosolProperty::PhaseVolumeFraction:
        {
          if (phases_.size() == 1)
          {
            provider = MakePhaseVolumeFractionProvider<DenseMatrixPolicy>(
                {}, {}, 0, state_parameter_indices.size(), state_variable_indices.size());
            break;
          }
          if (target_phase_name.empty())
            throw MiamException(
                MIAM_ERROR_CATEGORY_CONFIGURATION,
                MIAM_CONFIGURATION_PHASE_NAME_REQUIRED,
                "TwoMomentMode::GetPropertyProvider: target_phase_name required for PhaseVolumeFraction "
                "with multiple phases");
          std::vector<std::size_t> all_species;
          std::vector<double> all_molar_volumes;
          std::size_t phase_count = 0;
          for (const auto& phase : phases_)
            if (phase.name_ == target_phase_name)
            {
              for (const auto& ps : phase.phase_species_)
                if (!ps.species_.IsParameterized())
                {
                  all_species.push_back(state_variable_indices.at(prefix_ + "." + phase.name_ + "." + ps.species_.name_));
                  all_molar_volumes.push_back(
                      ps.species_.GetProperty<double>("molecular weight [kg mol-1]") /
                      ps.species_.GetProperty<double>("density [kg m-3]"));
                }
              phase_count = all_species.size();
              break;
            }
          for (const auto& phase : phases_)
          {
            if (phase.name_ == target_phase_name)
              continue;
            for (const auto& ps : phase.phase_species_)
              if (!ps.species_.IsParameterized())
              {
                all_species.push_back(state_variable_indices.at(prefix_ + "." + phase.name_ + "." + ps.species_.name_));
                all_molar_volumes.push_back(
                    ps.species_.GetProperty<double>("molecular weight [kg mol-1]") /
                    ps.species_.GetProperty<double>("density [kg m-3]"));
              }
          }
          provider = MakePhaseVolumeFractionProvider<DenseMatrixPolicy>(
              all_species, all_molar_volumes, phase_count, state_parameter_indices.size(), state_variable_indices.size());
          break;
        }
        default:
          throw MiamException(
              MIAM_ERROR_CATEGORY_CONFIGURATION,
              MIAM_CONFIGURATION_UNSUPPORTED_PROPERTY,
              "TwoMomentMode: unsupported AerosolProperty");
      }
      return provider;
    }

   private:
    std::string prefix_;                           // State name prefix to apply to shape properties
    std::vector<micm::Phase> phases_;              // Phases associated with the mode
    double default_geometric_standard_deviation_;  // Default geometric standard deviation
  };
}  // namespace miam