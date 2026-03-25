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
    /// @brief Single moment log-normal particle size distribution representation
    /// @details Represents a single moment log-normal distribution for aerosol or cloud particle size distributions.
    ///          Characterized by a geometric mean radius and geometric standard deviation.
    class SingleMomentMode
    {
     public:
      SingleMomentMode() = delete;

      SingleMomentMode(const std::string& prefix, const std::vector<micm::Phase>& phases)
          : prefix_(prefix),
            phases_(phases),
            default_geometric_mean_radius_(0.0),
            default_geometric_standard_deviation_(1.0)
      {
      }

      SingleMomentMode(
          const std::string& prefix,
          const std::vector<micm::Phase>& phases,
          const double geometric_mean_radius,
          const double geometric_standard_deviation)
          : prefix_(prefix),
            phases_(phases),
            default_geometric_mean_radius_(geometric_mean_radius),
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
        return { size, 2 };  // Two parameters: geometric mean radius and geometric standard deviation
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
        names.insert(prefix_ + ".GEOMETRIC_MEAN_RADIUS");
        names.insert(prefix_ + ".GEOMETRIC_STANDARD_DEVIATION");
        return names;
      }

      std::string Species(const micm::Phase& phase, const micm::Species& species) const
      {
        return prefix_ + "." + phase.name_ + "." + species.name_;
      }

      std::map<std::string, double> DefaultParameters() const
      {
        return { { prefix_ + ".GEOMETRIC_MEAN_RADIUS", default_geometric_mean_radius_ },
                 { prefix_ + ".GEOMETRIC_STANDARD_DEVIATION", default_geometric_standard_deviation_ } };
      }

      std::string GeometricMeanRadius() const
      {
        return prefix_ + ".GEOMETRIC_MEAN_RADIUS";
      }

      std::string GeometricStandardDeviation() const
      {
        return prefix_ + ".GEOMETRIC_STANDARD_DEVIATION";
      }

      void SetDefaultParameters(auto& state) const
      {
        auto gmd_it = state.custom_rate_parameter_map_.find(GeometricMeanRadius());
        if (gmd_it == state.custom_rate_parameter_map_.end())
        {
          throw std::runtime_error("Geometric mean radius parameter not found in state for " + prefix_);
        }
        auto gsd_it = state.custom_rate_parameter_map_.find(GeometricStandardDeviation());
        if (gsd_it == state.custom_rate_parameter_map_.end())
        {
          throw std::runtime_error("Geometric standard deviation parameter not found in state for " + prefix_);
        }
        for (std::size_t i_cell = 0; i_cell < state.variables_.NumRows(); ++i_cell)
        {
          state.custom_rate_parameters_[i_cell][gmd_it->second] = default_geometric_mean_radius_;
          state.custom_rate_parameters_[i_cell][gsd_it->second] = default_geometric_standard_deviation_;
        }
      }

      std::map<std::string, std::size_t> NumPhaseInstances() const
      {
        std::map<std::string, std::size_t> num_instances;
        for (const auto& phase : phases_)
        {
          num_instances[phase.name_] = 1;  // Single moment representation has one instance per phase
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
      /// @tparam DenseMatrixPolicy The dense matrix type used for state data
      /// @param property The aerosol property to provide
      /// @param state_parameter_indices Map of parameter names to column indices
      /// @param state_variable_indices Map of variable names to column indices
      /// @param target_phase_name Phase name for PhaseVolumeFraction (required if multi-phase)
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
            std::size_t gmd_idx = state_parameter_indices.at(GeometricMeanRadius());
            std::size_t gsd_idx = state_parameter_indices.at(GeometricStandardDeviation());
            provider.dependent_variable_indices = {};
            DenseMatrixPolicy dummy_partials{ 1, 0, 0.0 };
            provider.ComputeValue = DenseMatrixPolicy::Function(
                [gmd_idx, gsd_idx](auto&& params, auto&& vars, auto&& result)
                {
                  params.ForEachRow(
                      [](const double& gmd, const double& gsd, double& r_eff)
                      {
                        double ln_gsd = std::log(gsd);
                        r_eff = gmd * std::exp(2.5 * ln_gsd * ln_gsd);
                      },
                      params.GetConstColumnView(gmd_idx),
                      params.GetConstColumnView(gsd_idx),
                      result.GetColumnView(0));
                },
                dummy_params,
                dummy_vars,
                dummy_result);
            provider.ComputeValueAndDerivatives = DenseMatrixPolicy::Function(
                [gmd_idx, gsd_idx](auto&& params, auto&& vars, auto&& result, auto&& partials)
                {
                  params.ForEachRow(
                      [](const double& gmd, const double& gsd, double& r_eff)
                      {
                        double ln_gsd = std::log(gsd);
                        r_eff = gmd * std::exp(2.5 * ln_gsd * ln_gsd);
                      },
                      params.GetConstColumnView(gmd_idx),
                      params.GetConstColumnView(gsd_idx),
                      result.GetColumnView(0));
                },
                dummy_params,
                dummy_vars,
                dummy_result,
                dummy_partials);
            break;
          }
          case AerosolProperty::NumberConcentration:
          {
            std::size_t gmd_idx = state_parameter_indices.at(GeometricMeanRadius());
            std::size_t gsd_idx = state_parameter_indices.at(GeometricStandardDeviation());
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
            provider.dependent_variable_indices = species_indices;
            DenseMatrixPolicy dummy_partials{ 1, species_indices.size(), 0.0 };
            provider.ComputeValue = DenseMatrixPolicy::Function(
                [gmd_idx, gsd_idx, species_indices, mw_over_rho](auto&& params, auto&& vars, auto&& result)
                {
                  auto N = result.GetColumnView(0);
                  params.ForEachRow([](double& v) { v = 0.0; }, N);
                  for (std::size_t k = 0; k < species_indices.size(); ++k)
                    params.ForEachRow(
                        [mwr = mw_over_rho[k]](const double& c, double& V) { V += c * mwr; },
                        vars.GetConstColumnView(species_indices[k]),
                        N);
                  params.ForEachRow(
                      [](const double& gmd, const double& gsd, double& N)
                      {
                        double ln_gsd = std::log(gsd);
                        double V_s = (4.0 / 3.0) * std::numbers::pi * gmd * gmd * gmd *
                                     std::exp(4.5 * ln_gsd * ln_gsd);
                        N /= V_s;
                      },
                      params.GetConstColumnView(gmd_idx),
                      params.GetConstColumnView(gsd_idx),
                      N);
                },
                dummy_params,
                dummy_vars,
                dummy_result);
            provider.ComputeValueAndDerivatives = DenseMatrixPolicy::Function(
                [gmd_idx, gsd_idx, species_indices, mw_over_rho](
                    auto&& params, auto&& vars, auto&& result, auto&& partials)
                {
                  auto N = result.GetColumnView(0);
                  params.ForEachRow([](double& v) { v = 0.0; }, N);
                  for (std::size_t k = 0; k < species_indices.size(); ++k)
                    params.ForEachRow(
                        [mwr = mw_over_rho[k]](const double& c, double& V) { V += c * mwr; },
                        vars.GetConstColumnView(species_indices[k]),
                        N);
                  params.ForEachRow(
                      [](const double& gmd, const double& gsd, double& N)
                      {
                        double ln_gsd = std::log(gsd);
                        double V_s = (4.0 / 3.0) * std::numbers::pi * gmd * gmd * gmd *
                                     std::exp(4.5 * ln_gsd * ln_gsd);
                        N /= V_s;
                      },
                      params.GetConstColumnView(gmd_idx),
                      params.GetConstColumnView(gsd_idx),
                      N);
                  for (std::size_t k = 0; k < species_indices.size(); ++k)
                    params.ForEachRow(
                        [mwr = mw_over_rho[k]](const double& gmd, const double& gsd, double& dN)
                        {
                          double ln_gsd = std::log(gsd);
                          double V_s = (4.0 / 3.0) * std::numbers::pi * gmd * gmd * gmd *
                                       std::exp(4.5 * ln_gsd * ln_gsd);
                          dN = mwr / V_s;
                        },
                        params.GetConstColumnView(gmd_idx),
                        params.GetConstColumnView(gsd_idx),
                        partials.GetColumnView(k));
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
                  "SingleMomentMode::GetPropertyProvider: target_phase_name required for PhaseVolumeFraction "
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
          default: throw std::runtime_error("SingleMomentMode: unsupported AerosolProperty");
        }
        return provider;
      }

     private:
      std::string prefix_;                           // State name prefix to apply to mode properties
      std::vector<micm::Phase> phases_;              // Phases associated with the mode
      double default_geometric_mean_radius_;         // Default geometric mean radius
      double default_geometric_standard_deviation_;  // Default geometric standard deviation
    };
  }  // namespace representation
}  // namespace miam