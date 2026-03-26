// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/aerosol_property.hpp>

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
  namespace representation
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
          throw std::runtime_error(
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
        DenseMatrixPolicy dummy_params{ 1, state_parameter_indices.size(), 0.0 };
        DenseMatrixPolicy dummy_vars{ 1, state_variable_indices.size(), 0.0 };
        DenseMatrixPolicy dummy_result{ 1, 1, 0.0 };

        switch (property)
        {
          case AerosolProperty::EffectiveRadius:
          {
            // r_eff = (3·V_total/(4π·N))^(1/3) · exp(2.5·ln²(GSD))
            // Depends on all species variables and NUMBER_CONCENTRATION
            std::size_t gsd_idx = state_parameter_indices.at(GeometricStandardDeviation());
            std::size_t nc_var_idx = state_variable_indices.at(NumberConcentration());
            std::vector<std::size_t> species_indices;
            std::vector<double> mw_over_rho;
            for (const auto& phase : phases_)
              for (const auto& ps : phase.phase_species_)
                if (!ps.species_.IsParameterized())
                {
                  species_indices.push_back(
                      state_variable_indices.at(prefix_ + "." + phase.name_ + "." + ps.species_.name_));
                  mw_over_rho.push_back(
                      ps.species_.GetProperty<double>("molecular weight [kg mol-1]") /
                      ps.species_.GetProperty<double>("density [kg m-3]"));
                }
            // dependent_variable_indices: species first, N last
            provider.dependent_variable_indices = species_indices;
            provider.dependent_variable_indices.push_back(nc_var_idx);
            DenseMatrixPolicy dummy_partials{ 1, provider.dependent_variable_indices.size(), 0.0 };
            provider.ComputeValue = DenseMatrixPolicy::Function(
                [gsd_idx, nc_var_idx, species_indices, mw_over_rho](auto&& params, auto&& vars, auto&& result)
                {
                  auto r = result.GetColumnView(0);
                  params.ForEachRow([](double& v) { v = 0.0; }, r);
                  for (std::size_t k = 0; k < species_indices.size(); ++k)
                    params.ForEachRow(
                        [mwr = mw_over_rho[k]](const double& c, double& V) { V += c * mwr; },
                        vars.GetConstColumnView(species_indices[k]),
                        r);
                  params.ForEachRow(
                      [](const double& gsd, const double& nc, double& r_eff)
                      {
                        double ln_gsd = std::log(gsd);
                        double V_mean = r_eff / nc;  // r_eff held V_total
                        double r_mean = std::cbrt(3.0 * V_mean / (4.0 * std::numbers::pi));
                        r_eff = r_mean * std::exp(2.5 * ln_gsd * ln_gsd);
                      },
                      params.GetConstColumnView(gsd_idx),
                      vars.GetConstColumnView(nc_var_idx),
                      r);
                },
                dummy_params,
                dummy_vars,
                dummy_result);
            provider.ComputeValueAndDerivatives = DenseMatrixPolicy::Function(
                [gsd_idx, nc_var_idx, species_indices, mw_over_rho](
                    auto&& params, auto&& vars, auto&& result, auto&& partials)
                {
                  auto r = result.GetColumnView(0);
                  // Accumulate V_total into result
                  params.ForEachRow([](double& v) { v = 0.0; }, r);
                  for (std::size_t k = 0; k < species_indices.size(); ++k)
                    params.ForEachRow(
                        [mwr = mw_over_rho[k]](const double& c, double& V) { V += c * mwr; },
                        vars.GetConstColumnView(species_indices[k]),
                        r);
                  // Partials w.r.t. species (while r still holds V_total)
                  // ∂r_eff/∂[species_k] = r_eff / (3·V_total) · (Mw_k/ρ_k)
                  for (std::size_t k = 0; k < species_indices.size(); ++k)
                    params.ForEachRow(
                        [mwr = mw_over_rho[k]](
                            const double& gsd, const double& nc, const double& V_total, double& dr)
                        {
                          double ln_gsd = std::log(gsd);
                          double V_mean = V_total / nc;
                          double r_mean = std::cbrt(3.0 * V_mean / (4.0 * std::numbers::pi));
                          double r_eff = r_mean * std::exp(2.5 * ln_gsd * ln_gsd);
                          dr = r_eff * mwr / (3.0 * V_total);
                        },
                        params.GetConstColumnView(gsd_idx),
                        vars.GetConstColumnView(nc_var_idx),
                        r,
                        partials.GetColumnView(k));
                  // Partial w.r.t. N: ∂r_eff/∂N = -r_eff / (3·N)
                  std::size_t N_col = species_indices.size();
                  params.ForEachRow(
                      [](const double& gsd, const double& nc, const double& V_total, double& dr_dN)
                      {
                        double ln_gsd = std::log(gsd);
                        double V_mean = V_total / nc;
                        double r_mean = std::cbrt(3.0 * V_mean / (4.0 * std::numbers::pi));
                        double r_eff = r_mean * std::exp(2.5 * ln_gsd * ln_gsd);
                        dr_dN = -r_eff / (3.0 * nc);
                      },
                      params.GetConstColumnView(gsd_idx),
                      vars.GetConstColumnView(nc_var_idx),
                      r,
                      partials.GetColumnView(N_col));
                  // Now overwrite result with r_eff
                  params.ForEachRow(
                      [](const double& gsd, const double& nc, double& r_eff)
                      {
                        double ln_gsd = std::log(gsd);
                        double V_mean = r_eff / nc;  // r_eff still holds V_total
                        double r_mean = std::cbrt(3.0 * V_mean / (4.0 * std::numbers::pi));
                        r_eff = r_mean * std::exp(2.5 * ln_gsd * ln_gsd);
                      },
                      params.GetConstColumnView(gsd_idx),
                      vars.GetConstColumnView(nc_var_idx),
                      r);
                },
                dummy_params,
                dummy_vars,
                dummy_result,
                dummy_partials);
            break;
          }
          case AerosolProperty::NumberConcentration:
          {
            // N = state_variables[NUMBER_CONCENTRATION], d(N)/d(N) = 1
            std::size_t nc_var_idx = state_variable_indices.at(NumberConcentration());
            provider.dependent_variable_indices = { nc_var_idx };
            DenseMatrixPolicy dummy_partials{ 1, 1, 0.0 };
            provider.ComputeValue = DenseMatrixPolicy::Function(
                [nc_var_idx](auto&& params, auto&& vars, auto&& result)
                {
                  vars.ForEachRow(
                      [](const double& nc, double& N) { N = nc; },
                      vars.GetConstColumnView(nc_var_idx),
                      result.GetColumnView(0));
                },
                dummy_params,
                dummy_vars,
                dummy_result);
            provider.ComputeValueAndDerivatives = DenseMatrixPolicy::Function(
                [nc_var_idx](auto&& params, auto&& vars, auto&& result, auto&& partials)
                {
                  vars.ForEachRow(
                      [](const double& nc, double& N, double& dN_dN)
                      {
                        N = nc;
                        dN_dN = 1.0;
                      },
                      vars.GetConstColumnView(nc_var_idx),
                      result.GetColumnView(0),
                      partials.GetColumnView(0));
                },
                dummy_params,
                dummy_vars,
                dummy_result,
                dummy_partials);
            break;
          }
          case AerosolProperty::PhaseVolumeFraction:
          {
            if (phases_.size() == 1)
            {
              provider.dependent_variable_indices = {};
              DenseMatrixPolicy dummy_partials{ 1, 0, 0.0 };
              provider.ComputeValue = DenseMatrixPolicy::Function(
                  [](auto&& params, auto&& vars, auto&& result)
                  { params.ForEachRow([](double& phi) { phi = 1.0; }, result.GetColumnView(0)); },
                  dummy_params,
                  dummy_vars,
                  dummy_result);
              provider.ComputeValueAndDerivatives = DenseMatrixPolicy::Function(
                  [](auto&& params, auto&& vars, auto&& result, auto&& partials)
                  { params.ForEachRow([](double& phi) { phi = 1.0; }, result.GetColumnView(0)); },
                  dummy_params,
                  dummy_vars,
                  dummy_result,
                  dummy_partials);
              break;
            }
            if (target_phase_name.empty())
              throw std::runtime_error(
                  "TwoMomentMode::GetPropertyProvider: target_phase_name required for PhaseVolumeFraction "
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
                    all_species.push_back(
                        state_variable_indices.at(prefix_ + "." + phase.name_ + "." + ps.species_.name_));
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
                  all_species.push_back(
                      state_variable_indices.at(prefix_ + "." + phase.name_ + "." + ps.species_.name_));
                  all_mw_over_rho.push_back(
                      ps.species_.GetProperty<double>("molecular weight [kg mol-1]") /
                      ps.species_.GetProperty<double>("density [kg m-3]"));
                }
            }
            provider.dependent_variable_indices = all_species;
            DenseMatrixPolicy dummy_partials{ 1, all_species.size(), 0.0 };
            provider.ComputeValue = DenseMatrixPolicy::Function(
                [all_species, all_mw_over_rho, phase_count](auto&& params, auto&& vars, auto&& result)
                {
                  auto phi = result.GetColumnView(0);
                  auto V_phase = result.GetRowVariable();
                  params.ForEachRow([](double& vt, double& vp) { vt = 0.0; vp = 0.0; }, phi, V_phase);
                  for (std::size_t k = 0; k < all_species.size(); ++k)
                  {
                    if (k < phase_count)
                      params.ForEachRow(
                          [mwr = all_mw_over_rho[k]](const double& c, double& vt, double& vp)
                          {
                            double vol = c * mwr;
                            vt += vol;
                            vp += vol;
                          },
                          vars.GetConstColumnView(all_species[k]),
                          phi,
                          V_phase);
                    else
                      params.ForEachRow(
                          [mwr = all_mw_over_rho[k]](const double& c, double& vt) { vt += c * mwr; },
                          vars.GetConstColumnView(all_species[k]),
                          phi);
                  }
                  params.ForEachRow(
                      [](double& phi, const double& vp) { phi = (phi > 0.0) ? vp / phi : 1.0; }, phi, V_phase);
                },
                dummy_params,
                dummy_vars,
                dummy_result);
            provider.ComputeValueAndDerivatives = DenseMatrixPolicy::Function(
                [all_species, all_mw_over_rho, phase_count](
                    auto&& params, auto&& vars, auto&& result, auto&& partials)
                {
                  auto result_col = result.GetColumnView(0);
                  auto V_phase = result.GetRowVariable();
                  auto V_total = result.GetRowVariable();
                  params.ForEachRow([](double& vt, double& vp) { vt = 0.0; vp = 0.0; }, V_total, V_phase);
                  for (std::size_t k = 0; k < all_species.size(); ++k)
                  {
                    if (k < phase_count)
                      params.ForEachRow(
                          [mwr = all_mw_over_rho[k]](const double& c, double& vt, double& vp)
                          {
                            double vol = c * mwr;
                            vt += vol;
                            vp += vol;
                          },
                          vars.GetConstColumnView(all_species[k]),
                          V_total,
                          V_phase);
                    else
                      params.ForEachRow(
                          [mwr = all_mw_over_rho[k]](const double& c, double& vt) { vt += c * mwr; },
                          vars.GetConstColumnView(all_species[k]),
                          V_total);
                  }
                  params.ForEachRow(
                      [](const double& vp, const double& vt, double& phi)
                      { phi = (vt > 0.0) ? vp / vt : 1.0; },
                      V_phase,
                      V_total,
                      result_col);
                  for (std::size_t k = 0; k < all_species.size(); ++k)
                  {
                    if (k < phase_count)
                      params.ForEachRow(
                          [mwr = all_mw_over_rho[k]](const double& phi, const double& vt, double& dphi)
                          { dphi = (vt > 0.0) ? mwr * (1.0 - phi) / vt : 0.0; },
                          result_col,
                          V_total,
                          partials.GetColumnView(k));
                    else
                      params.ForEachRow(
                          [mwr = all_mw_over_rho[k]](const double& phi, const double& vt, double& dphi)
                          { dphi = (vt > 0.0) ? -mwr * phi / vt : 0.0; },
                          result_col,
                          V_total,
                          partials.GetColumnView(k));
                  }
                },
                dummy_params,
                dummy_vars,
                dummy_result,
                dummy_partials);
            break;
          }
          default: throw std::runtime_error("TwoMomentMode: unsupported AerosolProperty");
        }
        return provider;
      }

     private:
      std::string prefix_;                           // State name prefix to apply to shape properties
      std::vector<micm::Phase> phases_;              // Phases associated with the mode
      double default_geometric_standard_deviation_;  // Default geometric standard deviation
    };
  }  // namespace representation
}  // namespace miam