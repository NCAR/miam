// Copyright (C) 2026 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/model/aerosol_scheme.hpp>
#include <miam/util/utils.hpp>

#include <micm/system/phase.hpp>

#include <cmath>
#include <stdexcept>
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
    double GetEffectiveRadius() const;

    /// @brief Get effective radius for the two-moment mode (calculates from mass and number concentration)
    /// @tparam StateType Type of the state object (e.g., micm::State)
    /// @param state The state object containing variable data
    /// @param cell The index of the grid cell (default 0)
    /// @return Effective radius [m]
    /// @throws std::runtime_error if SetStateIndices() has not been called
    template<typename StateType>
    double GetEffectiveRadius(const StateType& state, std::size_t cell = 0)
    {
      if (!has_initialized_state_idx_)
      {
        throw std::runtime_error("State indices not initialized. Call SetStateIndices() before GetEffectiveRadius().");
      }

      double total_mass = 0.0;
      for (const auto& [species_key, id] : state_idx_.state_id_map)
      {
        total_mass += state.variables_[cell][id];
      }
      double number_concentration = state.variables_[cell][state_idx_.number_id];
      double density = state.variables_[cell][state_idx_.density_id];

      return CalculateEffectiveRadius(total_mass, number_concentration, density, geometric_standard_deviation_);
    }

    /// @brief Set the state indices for accessing section variables in the state vector
    /// @tparam StateType Type of the state object (e.g., micm::State)
    /// @param state The state object containing variable map
    /// @throws std::runtime_error If keys are not found in state
    template<typename StateType>
    void SetStateIndices(const StateType& state)
    {
      // Find mass indices for all species in all phases
      for (const auto& phase : phases_)
      {
        for (const auto& phase_species : phase.phase_species_)
        {
          // NAME: SECTION.PHASE.SPECIES
          std::string species_key = JoinStrings({ name_, phase.name_, phase_species.species_.name_ });
          auto species_it = state.variable_map_.find(species_key);
          if (species_it == state.variable_map_.end())
          {
            throw std::runtime_error(std::format(("Species '{}' not found in state for '{}'", species_key, name_)));
          }
          state_idx_.state_id_map[species_key] = species_it->second;
        }
      }

      // Find number concentration index
      std::string number_key = JoinStrings({ name_, AerosolScheme::AEROSOL_MOMENTS_[0] });
      auto number_it = state.variable_map_.find(number_key);
      if (number_it == state.variable_map_.end())
      {
        throw std::runtime_error(std::format(("Variable '{}' not found in state for '{}'", number_key, name_)));
      }
      state_idx_.number_id = number_it->second;

      // Find density index
      std::string density_key = JoinStrings({ name_, AerosolScheme::AEROSOL_MOMENTS_[1] });
      auto density_it = state.variable_map_.find(density_key);
      if (density_it == state.variable_map_.end())
      {
        throw std::runtime_error(std::format(("Variable '{}' not found in state for '{}'", density_key, name_)));
      }
      state_idx_.density_id = density_it->second;

      has_initialized_state_idx_ = true;
    }

   private:
    /// @brief Calculate the effective radius for a log-normal aerosol mode
    /// @param mass Total mass concentration of particles in the mode [kg m-3]
    /// @param N Total number concentration of particles in the mode [# m-3]
    /// @param density Particle density [kg m-3]
    /// @param sig_g Geometric standard deviation (unitless)
    /// @return Effective radius [m]
    double CalculateEffectiveRadius(double mass, double N, double density, double sig_g);
  };

}  // namespace miam
