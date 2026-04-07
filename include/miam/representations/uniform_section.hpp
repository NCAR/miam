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
          throw std::runtime_error("Custom parameter map missing MIN_RADIUS for " + prefix_);
        }
        auto max_radius_it = state.custom_rate_parameter_map_.find(MaxRadius());
        if (max_radius_it == state.custom_rate_parameter_map_.end())
        {
          throw std::runtime_error("Custom parameter map missing MAX_RADIUS for " + prefix_);
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
        DenseMatrixPolicy dummy_params{ 1, state_parameter_indices.size(), 0.0 };
        DenseMatrixPolicy dummy_vars{ 1, state_variable_indices.size(), 0.0 };
        DenseMatrixPolicy dummy_result{ 1, 1, 0.0 };

        switch (property)
        {
          case AerosolProperty::EffectiveRadius:
          {
            // r_eff = (r_min + r_max) / 2 — no state variable dependencies
            std::size_t rmin_idx = state_parameter_indices.at(MinRadius());
            std::size_t rmax_idx = state_parameter_indices.at(MaxRadius());
            provider.dependent_variable_indices = {};
            DenseMatrixPolicy dummy_partials{ 1, 0, 0.0 };
            provider.ComputeValue = DenseMatrixPolicy::Function(
                [rmin_idx, rmax_idx](auto&& params, auto&& vars, auto&& result)
                {
                  params.ForEachRow(
                      [](const double& r_min, const double& r_max, double& r_eff) { r_eff = 0.5 * (r_min + r_max); },
                      params.GetConstColumnView(rmin_idx),
                      params.GetConstColumnView(rmax_idx),
                      result.GetColumnView(0));
                },
                dummy_params,
                dummy_vars,
                dummy_result);
            provider.ComputeValueAndDerivatives = DenseMatrixPolicy::Function(
                [rmin_idx, rmax_idx](auto&& params, auto&& vars, auto&& result, auto&& partials)
                {
                  params.ForEachRow(
                      [](const double& r_min, const double& r_max, double& r_eff) { r_eff = 0.5 * (r_min + r_max); },
                      params.GetConstColumnView(rmin_idx),
                      params.GetConstColumnView(rmax_idx),
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
            // N = V_total / V_single, V_single = (4/3)π·r_eff³
            std::size_t rmin_idx = state_parameter_indices.at(MinRadius());
            std::size_t rmax_idx = state_parameter_indices.at(MaxRadius());
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
                [rmin_idx, rmax_idx, species_indices, mw_over_rho](auto&& params, auto&& vars, auto&& result)
                {
                  auto N = result.GetColumnView(0);
                  params.ForEachRow([](double& v) { v = 0.0; }, N);
                  for (std::size_t k = 0; k < species_indices.size(); ++k)
                    params.ForEachRow(
                        [mwr = mw_over_rho[k]](const double& c, double& V) { V += c * mwr; },
                        vars.GetConstColumnView(species_indices[k]),
                        N);
                  params.ForEachRow(
                      [](const double& r_min, const double& r_max, double& N)
                      {
                        double r_eff = 0.5 * (r_min + r_max);
                        double V_s = (4.0 / 3.0) * std::numbers::pi * r_eff * r_eff * r_eff;
                        N /= V_s;
                      },
                      params.GetConstColumnView(rmin_idx),
                      params.GetConstColumnView(rmax_idx),
                      N);
                },
                dummy_params,
                dummy_vars,
                dummy_result);
            provider.ComputeValueAndDerivatives = DenseMatrixPolicy::Function(
                [rmin_idx, rmax_idx, species_indices, mw_over_rho](
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
                      [](const double& r_min, const double& r_max, double& N)
                      {
                        double r_eff = 0.5 * (r_min + r_max);
                        double V_s = (4.0 / 3.0) * std::numbers::pi * r_eff * r_eff * r_eff;
                        N /= V_s;
                      },
                      params.GetConstColumnView(rmin_idx),
                      params.GetConstColumnView(rmax_idx),
                      N);
                  for (std::size_t k = 0; k < species_indices.size(); ++k)
                    params.ForEachRow(
                        [mwr = mw_over_rho[k]](const double& r_min, const double& r_max, double& dN)
                        {
                          double r_eff = 0.5 * (r_min + r_max);
                          double V_s = (4.0 / 3.0) * std::numbers::pi * r_eff * r_eff * r_eff;
                          dN = mwr / V_s;
                        },
                        params.GetConstColumnView(rmin_idx),
                        params.GetConstColumnView(rmax_idx),
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
            provider.dependent_variable_indices = all_species;
            DenseMatrixPolicy dummy_partials{ 1, all_species.size(), 0.0 };
            provider.ComputeValue = DenseMatrixPolicy::Function(
                [all_species, all_mw_over_rho, phase_count](auto&& params, auto&& vars, auto&& result)
                {
                  auto phi = result.GetColumnView(0);
                  auto V_phase = result.GetRowVariable();
                  params.ForEachRow(
                      [](double& vt, double& vp)
                      {
                        vt = 0.0;
                        vp = 0.0;
                      },
                      phi,
                      V_phase);
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
                  params.ForEachRow([](double& phi, const double& vp) { phi = (phi > 0.0) ? vp / phi : 1.0; }, phi, V_phase);
                },
                dummy_params,
                dummy_vars,
                dummy_result);
            provider.ComputeValueAndDerivatives = DenseMatrixPolicy::Function(
                [all_species, all_mw_over_rho, phase_count](auto&& params, auto&& vars, auto&& result, auto&& partials)
                {
                  auto result_col = result.GetColumnView(0);
                  auto V_phase = result.GetRowVariable();
                  auto V_total = result.GetRowVariable();
                  params.ForEachRow(
                      [](double& vt, double& vp)
                      {
                        vt = 0.0;
                        vp = 0.0;
                      },
                      V_total,
                      V_phase);
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
                      [](const double& vp, const double& vt, double& phi) { phi = (vt > 0.0) ? vp / vt : 1.0; },
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
          default: throw std::runtime_error("UniformSection: unsupported AerosolProperty");
        }
        return provider;
      }

     private:
      std::string prefix_;               // State name prefix to apply to section properties
      std::vector<micm::Phase> phases_;  // Phases associated with the section
      double default_min_radius_;        // Minimum radius of the section
      double default_max_radius_;        // Maximum radius of the section
    };
  }  // namespace representation
}  // namespace miam