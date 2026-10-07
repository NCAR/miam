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
  /// @brief Sectional particle size distribution representation with uniform sections
  /// @details Represents a sectional distribution with uniform sections for aerosol or cloud particle size distributions.
  ///          Each section is characterized by a fixed size range and variable total volume. Number concentrations are
  ///          derived from the total volume and section size.
  class UniformSection
  {
   public:
    UniformSection() = delete;

    UniformSection(const std::string& prefix, const std::vector<micm::Phase>& phases)
        : prefix_(prefix),
          phases_(phases),
          default_min_radius_(0.0),
          default_max_radius_(0.0)
    {
    }

    UniformSection(
        const std::string& prefix,
        const std::vector<micm::Phase>& phases,
        const double minimum_radius,
        const double maximum_radius)
        : prefix_(prefix),
          phases_(phases),
          default_min_radius_(minimum_radius),
          default_max_radius_(maximum_radius)
    {
    }

    std::tuple<std::size_t, std::size_t> StateSize() const
    {
      std::size_t size = 0;
      for (const auto& phase : phases_)
      {
        size += phase.StateSize();
      }
      return { size, 2 };  // Two parameters: min and max radius
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
      return names;
    }

    std::set<std::string> StateParameterNames() const
    {
      std::set<std::string> names;
      names.insert(prefix_ + ".MIN_RADIUS");
      names.insert(prefix_ + ".MAX_RADIUS");
      return names;
    }

    std::string Species(const micm::Phase& phase, const micm::Species& species) const
    {
      return prefix_ + "." + phase.name_ + "." + species.name_;
    }

    std::map<std::string, double> DefaultParameters() const
    {
      return { { prefix_ + ".MIN_RADIUS", default_min_radius_ }, { prefix_ + ".MAX_RADIUS", default_max_radius_ } };
    }

    std::string MinRadius() const
    {
      return prefix_ + ".MIN_RADIUS";
    }

    std::string MaxRadius() const
    {
      return prefix_ + ".MAX_RADIUS";
    }

    void SetDefaultParameters(auto& state) const
    {
      auto min_radius_it = state.custom_rate_parameter_map_.find(MinRadius());
      if (min_radius_it == state.custom_rate_parameter_map_.end())
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_STATE_PARAMETER,
            "UniformSection::SetDefaultParameters: MIN_RADIUS parameter not found in state for " + prefix_);
      }
      auto max_radius_it = state.custom_rate_parameter_map_.find(MaxRadius());
      if (max_radius_it == state.custom_rate_parameter_map_.end())
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_STATE_PARAMETER,
            "UniformSection::SetDefaultParameters: MAX_RADIUS parameter not found in state for " + prefix_);
      }
      for (std::size_t cell = 0; cell < state.variables_.NumRows(); ++cell)
      {
        state.custom_rate_parameters_[cell][min_radius_it->second] = default_min_radius_;
        state.custom_rate_parameters_[cell][max_radius_it->second] = default_max_radius_;
      }
    }

    std::map<std::string, std::size_t> NumPhaseInstances() const
    {
      std::map<std::string, std::size_t> num_instances;
      for (const auto& phase : phases_)
      {
        num_instances[phase.name_] = 1;  // Uniform section representation has one instance per phase
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
          // r_eff = (r_min + r_max) / 2 - no state variable dependencies
          std::size_t rmin_idx = state_parameter_indices.at(MinRadius());
          std::size_t rmax_idx = state_parameter_indices.at(MaxRadius());
          provider.dependent_variable_indices = {};
          DenseMatrixPolicy example_params{ 1, state_parameter_indices.size(), 0.0 };
          DenseMatrixPolicy example_vars{ 1, state_variable_indices.size(), 0.0 };
          DenseMatrixPolicy example_result{ 1, 1, 0.0 };
          auto value_function =
            DenseMatrixPolicy::Function(
                MICM_LAMBDA(
                    const typename DenseMatrixPolicy::ConstViewType& params_view,
                    const typename DenseMatrixPolicy::ConstViewType& /*vars_view*/,
                    const typename DenseMatrixPolicy::ViewType& result_view) {
                  params_view.ForEachRowStrict(
                      [](const double& r_min, const double& r_max, double& r_eff) { r_eff = 0.5 * (r_min + r_max); },
                      params_view.GetConstColumnView(rmin_idx),
                      params_view.GetConstColumnView(rmax_idx),
                      result_view.GetColumnView(0));
                },
                example_params,
                example_vars,
                example_result);
          provider.ComputeValue = [value_function](
                                      const DenseMatrixPolicy& params,
                                      const DenseMatrixPolicy& vars,
                                      DenseMatrixPolicy& result) mutable { value_function(params, vars, result); };
          provider.ComputeValueAndDerivatives = [compute_value = provider.ComputeValue](
                                                    const DenseMatrixPolicy& params,
                                                    const DenseMatrixPolicy& vars,
                                                    DenseMatrixPolicy& result,
                                                    DenseMatrixPolicy& /*partials*/)
          { compute_value(params, vars, result); };
          break;
        }
        case AerosolProperty::NumberConcentration:
        {
          // N = V_total / V_single, V_single = (4/3)pi*r_eff^3
          std::size_t rmin_idx = state_parameter_indices.at(MinRadius());
          std::size_t rmax_idx = state_parameter_indices.at(MaxRadius());
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
          provider.dependent_variable_indices = species_indices;
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
          DenseMatrixPolicy example_partials{ 1, species_indices.size(), 0.0 };
          auto value_function =
            DenseMatrixPolicy::Function(
                MICM_LAMBDA(
                    const typename DenseMatrixPolicy::ConstViewType& params_view,
                    const typename DenseMatrixPolicy::ConstViewType& vars_view,
                    const typename DenseMatrixPolicy::ViewType& result_view) {
                  auto N = result_view.GetColumnView(0);
                  params_view.ForEachRowStrict([](double& v) { v = 0.0; }, N);
                  const std::size_t n = species_view.size();
                  for (std::size_t k = 0; k < n; ++k)
                  {
                    const double mv = volumes_view[k];
                    params_view.ForEachRowStrict(
                        [mv](const double& c, double& V) { V += c * mv; }, vars_view.GetConstColumnView(species_view[k]), N);
                  }
                  params_view.ForEachRowStrict(
                      [](const double& r_min, const double& r_max, double& N_out)
                      {
                        const double r_eff = 0.5 * (r_min + r_max);
                        const double V_s = (4.0 / 3.0) * std::numbers::pi * r_eff * r_eff * r_eff;
                        N_out /= V_s;
                      },
                      params_view.GetConstColumnView(rmin_idx),
                      params_view.GetConstColumnView(rmax_idx),
                      N);
                },
                example_params,
                example_vars,
                example_result);
          provider.ComputeValue = [storage, value_function](
                  const DenseMatrixPolicy& params,
                  const DenseMatrixPolicy& vars,
                                      DenseMatrixPolicy& result) mutable { value_function(params, vars, result); };
          const std::size_t n = species_indices.size();
          auto partials_function =
            DenseMatrixPolicy::Function(
                MICM_LAMBDA(
                    const typename DenseMatrixPolicy::ConstViewType& params_view,
                    const typename DenseMatrixPolicy::ConstViewType& /*vars_view*/,
                    const typename DenseMatrixPolicy::ViewType& /*result_view*/,
                    const typename DenseMatrixPolicy::ViewType& partials_view) {
                  for (std::size_t k = 0; k < n; ++k)
                  {
                    const double mv = volumes_view[k];
                    params_view.ForEachRowStrict(
                        [mv](const double& r_min, const double& r_max, double& dN)
                        {
                          const double r_eff = 0.5 * (r_min + r_max);
                          const double V_s = (4.0 / 3.0) * std::numbers::pi * r_eff * r_eff * r_eff;
                          dN = mv / V_s;
                        },
                        params_view.GetConstColumnView(rmin_idx),
                        params_view.GetConstColumnView(rmax_idx),
                        partials_view.GetColumnView(k));
                  }
                },
                example_params,
                example_vars,
                example_result,
                example_partials);
          provider.ComputeValueAndDerivatives =
              [storage, compute_value = provider.ComputeValue, partials_function, n](
                  const DenseMatrixPolicy& params,
                  const DenseMatrixPolicy& vars,
                  DenseMatrixPolicy& result,
                  DenseMatrixPolicy& partials) mutable
          {
            compute_value(params, vars, result);
            if (n == 0)
              return;
            partials_function(params, vars, result, partials);
          };
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
                "UniformSection::GetPropertyProvider: target_phase_name required for PhaseVolumeFraction "
                "with multiple phases");
          std::vector<std::size_t> all_species;
          std::vector<double> all_mw_over_rho;
          std::size_t phase_count = 0;
          for (const auto& phase : phases_)
            if (phase.name_ == target_phase_name)
            {
              for (const auto& ps : phase.phase_species_)
                if (!ps.species_.IsParameterized())
                {
                  all_species.push_back(state_variable_indices.at(prefix_ + "." + phase.name_ + "." + ps.species_.name_));
                  all_mw_over_rho.push_back(
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
                all_mw_over_rho.push_back(
                    ps.species_.GetProperty<double>("molecular weight [kg mol-1]") /
                    ps.species_.GetProperty<double>("density [kg m-3]"));
              }
          }
          provider = MakePhaseVolumeFractionProvider<DenseMatrixPolicy>(
              all_species, all_mw_over_rho, phase_count, state_parameter_indices.size(), state_variable_indices.size());
          break;
        }
        default:
          throw MiamException(
              MIAM_ERROR_CATEGORY_CONFIGURATION,
              MIAM_CONFIGURATION_UNSUPPORTED_PROPERTY,
              "UniformSection: unsupported AerosolProperty");
      }
      return provider;
    }

   private:
    std::string prefix_;               // State name prefix to apply to section properties
    std::vector<micm::Phase> phases_;  // Phases associated with the section
    double default_min_radius_;        // Minimum radius of the section
    double default_max_radius_;        // Maximum radius of the section
  };
}  // namespace miam