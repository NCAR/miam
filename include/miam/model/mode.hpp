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
#include <format>

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

    std::unordered_map<std::string, std::size_t> map_state_id_;  // species - index pairs in state
    std::size_t number_id_;                                      // Index for number concentration
    std::size_t density_id_;                                     // Index for density
    std::size_t radius_id_;                                      // Index for radius
  
   private:
    double fixed_radius_;

   public:
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
      fixed_radius_ = GetEffectiveRadius();
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
    double GetEffectiveRadius(StateType& state, std::size_t cell = 0)
    {
      double total_mass = 0.0;
      for (const auto& [species_key, id] : map_state_id_)
      {
        total_mass += state.variables_[cell][id];
      }
      double number_concentration = state.variables_[cell][number_id_];
      double density = state.variables_[cell][density_id_];

      return CalculateEffectiveRadius(total_mass, number_concentration, density, geometric_standard_deviation_);
    }

    /// @brief Set the state indices for accessing mode variables in the state vector
    /// @tparam StateType Type of the state object (e.g., micm::State)
    /// @param state The state object containing variable map
    /// @throws std::runtime_error If keys are not found in state
    // template<typename StateType>
    // void SetStateIndices(const StateType& state)
    // {
    //   std::cout <<  "void SetStateIndices(const StateType& state) MODE" << std::endl;
    //   // Find mass indices for all species in all phases
    //   for (const auto& phase : phases_)
    //   {
    //     for (const auto& phase_species : phase.phase_species_)
    //     {
    //       // NAME: SECTION.PHASE.SPECIES
    //       std::string species_key = JoinStrings({ name_, phase.name_, phase_species.species_.name_ });
    //       auto species_it = state.variable_map_.find(species_key);
    //       if (species_it == state.variable_map_.end())
    //       {
    //         throw std::runtime_error(std::format("Species '{}' not found in state for '{}'", species_key, name_));
    //       }
    //       map_state_id_[species_key] = species_it->second;
    //     }
    //   }

    //   // Find number concentration index
    //   std::string number_key = JoinStrings({ name_, AerosolModel::AEROSOL_MOMENTS_[0] });
    //   auto number_it = state.variable_map_.find(number_key);
    //   if (number_it == state.variable_map_.end())
      // {
      //   throw std::runtime_error(std::format("Variable '{}' not found in state for '{}'", number_key, name_));
      // }
      // number_id = number_it->second;

      // // Find density index
      // std::string density_key = JoinStrings({ name_, AerosolModel::AEROSOL_MOMENTS_[1] });
      // auto density_it = state.variable_map_.find(density_key);
      // if (density_it == state.variable_map_.end())
      // {
      //   throw std::runtime_error(std::format("Variable '{}' not found in state for '{}'", density_key, name_));
      // }
      // density_id = density_it->second;

      // // Find radius index
      // std::string radius_key = JoinStrings({ name_, AerosolModel::AEROSOL_MOMENTS_[2] });
      // auto radius_it = state.variable_map_.find(radius_key);
      // if (radius_it == state.variable_map_.end())
      // {
    //     throw std::runtime_error(std::format("Variable '{}' not found in state for '{}'", radius_key, name_));
    //   }
    //   radius_id = radius_it->second;

    //       for (auto& [key, value] : map_state_id_)
    // {
    //   std::cout << "key: " << key << " value: " << value << std::endl;
    // }
    // }

    /// @brief Set the effective radius in the state for this mode
    /// @tparam StateType Type of the state object (e.g., micm::State)
    /// @param state The state object containing variable map
    /// @param cell The index of the grid cell (default 0)
    /// @throws std::runtime_error If radius key is not found in state
    template<typename StateType>
    void SetRadius(StateType& state, std::size_t cell = 0)
    {
      std::cout << "Set Radius: " << name_ << std::endl;
      if (distribution_ == DistributionType::SingleMoment)
      {
        std::cout << "fixed_radius_;" << fixed_radius_ << std::endl;
        state.variables_[cell][radius_id_] = fixed_radius_;
      }
        else
      {
        std::cout << "Not fixed" << "GetEffectiveRadius(state, cell);"<< std::endl;
        state.variables_[cell][radius_id_] = GetEffectiveRadius(state, cell);
      }
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
