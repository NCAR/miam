// Copyright (C) 2026 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/util/utils.hpp>

#include <micm/system/phase.hpp>

#include <iostream>
#include <string>
#include <unordered_map>
#include <format>

namespace miam
{

  enum class DistributionType
  {
    SingleMoment,  // tracks mass, fixed radius, number calculated
    TwoMoment      // tracks mass and number concentration, radius calculated
  };

  /// @brief Base AerosolScheme class (common to both modal and sectional)
  class AerosolScheme
  {
   public:
    // Indices from the state
    struct StateIndices
    {
      std::unordered_map<std::string, std::size_t> state_id_map;  // species - index pairs in state
      std::size_t number_id;                                      // Index for number concentration
      std::size_t density_id;                                     // Index for density
      std::size_t radius_id;                                      // Index for radius
    };

    inline static const std::vector<std::string> AEROSOL_MOMENTS_ = { "NUMBER_CONCENTRATION", "DENSITY", "RADIUS" };

    /// @brief Name of the mode or bin (e.g., "aitken", "accumulation", "large_drop")
    std::string name_;

    AerosolScheme(std::string name)
        : name_(std::move(name))
    {
    }

    virtual ~AerosolScheme() = default;

    inline std::string GetScope() const
    {
      return name_;
    }

    /// @brief Set concentration of a species
    /// @tparam StateType Type of the state object (e.g., micm::State)
    /// @param state The state object to modify
    /// @param phase The phase where the species belongs
    /// @param species The species to set concentration for
    /// @param concentration The concentration value to set [mol m-3]
    /// @param cell The grid cell index (default 0)
    /// @throws std::runtime_error If the species is not found in state
    template<typename StateType>
    void SetConcentration(
        StateType& state,
        const micm::Phase& phase,
        const micm::Species& species,
        double concentration,
        std::size_t cell = 0)
    {
      const std::string species_key = JoinStrings({ name_, phase.name_, species.name_ });
      std::size_t index;

      if (has_initialized_state_idx_)
      {
        auto it = state_idx_.state_id_map.find(species_key);
        if (it == state_idx_.state_id_map.end())
        {
          throw std::runtime_error(std::format("Species '{}' not found in state_id_map for '{}'", species_key, name_));
        }
        index = it->second;
      }
      else
      {
        auto it = state.variable_map_.find(species_key);
        if (it == state.variable_map_.end())
        {
          throw std::runtime_error(std::format("Species '{}' not found in state for '{}'", species_key, name_));
        }
        index = it->second;
      }

      state.variables_[cell][index] = concentration;
    }

    /// @brief Set number concentration
    /// @param state StateType Type of the state object
    /// @param concentration The number concentration value to set [# m-3]
    /// @param cell The grid cell index (default 0)
    /// @throws std::runtime_error If the number concentration variable is not found in state
    template<typename StateType>
    void SetNumberConcentration(StateType& state, double concentration, std::size_t cell = 0)
    {
      std::size_t index;

      if (has_initialized_state_idx_)
        index = state_idx_.number_id;
      else
      {
        std::string number_key = JoinStrings({ name_, AerosolScheme::AEROSOL_MOMENTS_[0] });
        auto it = state.variable_map_.find(number_key);
        if (it == state.variable_map_.end())
        {
          throw std::runtime_error(std::format("Variable '{}' not found in state for '{};", number_key, name_));
        }
        index = it->second;
      }

      state.variables_[cell][index] = concentration;
    }

   protected:
    /// @brief Indices for accessing data in the state
    StateIndices state_idx_;

    bool has_initialized_state_idx_ = false;
    bool is_radius_fixed_ = false;

    double fixed_radius_ = 0.0;
  };

}  // namespace miam