// Copyright (C) 2026 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/model/aerosol_scheme.hpp>
#include <miam/util/utils.hpp>

#include <micm/system/phase.hpp>

#include <cmath>
#include <string>
#include <vector>

namespace miam
{

  class Mode : public AerosolScheme
  {
   public:
    /// @brief Phases associated with this mode (e.g., aqueous, organic)
    std::vector<micm::Phase> phases_;
    /// @brief Type of distribution (single moment or two moment)
    DistributionType distribution_;
    /// @brief Geometric mean diameter [m] - the center of the log-normal size distribution
    double geometric_mean_diameter_;
    /// @brief Geometric standard deviation (unitless) - width of the log-normal size distribution
    double geometric_standard_deviation_;

    /// @brief Construct a Mode with specified physical properties
    Mode(
        std::string name,
        std::vector<micm::Phase> phases,
        DistributionType distribution,
        double geometric_mean_diameter,
        double geometric_standard_deviation)
        : AerosolScheme(std::move(name)),
          phases_(std::move(phases)),
          distribution_(distribution),
          geometric_mean_diameter_(geometric_mean_diameter),
          geometric_standard_deviation_(geometric_standard_deviation)
    {
    }

    /// @brief Get effective radius for the single-moment mode (uses fixed geometric mean diameter)
    /// @return Effective radius [m]
    inline double GetEffectiveRadius() const
    {
      const double r_g = 0.5 * geometric_mean_diameter_;
      const double ln_sig = std::log(geometric_standard_deviation_);
      return r_g * std::exp(2.5 * ln_sig * ln_sig);
    }

    // TODO - GetEffectiveRadius() assumes SetStateIndices() has already been called
    /// @brief Get effective radius for the two-moment mode (calculates from mass and number concentration)
    /// @tparam StateType Type of the state object (e.g., micm::State)
    /// @param state The state object containing variable data
    /// @param cell The index of the grid cell
    /// @return Effective radius [m]
    template<typename StateType>
    double GetEffectiveRadius(const StateType& state, std::size_t cell)
    {
      double total_mass = 0.0;
      for (const auto& [species_key, id] : state_idx_.state_id_map)
      {
        total_mass += state.variables_[cell][id];
      }

      double number_concentration = state.variables_[cell][state_idx_.number_id];
      double density = state.variables_[cell][state_idx_.density_id];

      return CalculateEffectiveRadius(total_mass, number_concentration, density, geometric_standard_deviation_);
    }

    /// @brief Set the state indices for accessing mode variables in the state vector
    /// @tparam StateType Type of the state object (e.g., micm::State)
    /// @param state The state object containing variable map
    template<typename StateType>
    void SetStateIndices(const StateType& state)
    {
      // Find mass indices for all species in all phases
      for (const auto& phase : phases_)
      {
        for (const auto& phase_species : phase.phase_species_)
        {
          // NAME: MODE.PHASE.SPECIES
          std::string key = JoinStrings({ name_, phase.name_, phase_species.species_.name_ });
          state_idx_.state_id_map[key] = state.variable_map_.at(key);
        }
      }

      // Find number concentration index
      state_idx_.number_id = state.variable_map_.at(JoinStrings({ name_, AerosolScheme::AEROSOL_MOMENTS_[0] }));

      // Find density index
      state_idx_.density_id = state.variable_map_.at(JoinStrings({ name_, AerosolScheme::AEROSOL_MOMENTS_[1] }));

      has_initialized_state_idx_ = true;
    }

   private:
    /// @brief Calculate the effective radius for a log-normal aerosol mode
    /// @param mass Total mass concentration of particles in the mode [kg/m³]
    /// @param N Total number concentration of particles in the mode [#/m³]
    /// @param density Particle density [kg/m³]
    /// @param sig_g Geometric standard deviation (unitless)
    /// @return Effective radius [m]
    double CalculateEffectiveRadius(double mass, double N, double density, double sig_g);
  };

}  // namespace miam
