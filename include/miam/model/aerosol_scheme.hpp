// Copyright (C) 2026 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/util/utils.hpp>

#include <micm/system/phase.hpp>

#include <string>
#include <unordered_map>

#include <iostream>

namespace miam
{

enum class DistributionType
{
  SingleMoment, // tracks mass, fixed radius, number calculated
  TwoMoment     // tracks mass and number concentration, radius calculated
};

/// @brief Base AerosolScheme class (common to both modal and sectional)
class AerosolScheme
{
public:

  // Indices from the state
  struct StateIndices
  {
    std::unordered_map<std::string, std::size_t> state_id_map; // species - index pairs in state
    std::size_t number_id;                                     // Index for number concentration
    std::size_t density_id;                                    // Index for density
  };

  inline static const std::vector<std::string> AEROSOL_MOMENTS_ = {"NUMBER_CONCENTRATION", "DENSITY"};

  /// @brief Name of the mode or bin (e.g., "aitken", "accumulation", "large_drop")
  std::string name_;

  AerosolScheme(std::string name)
    : name_(std::move(name)) {}

  virtual ~AerosolScheme() = default;

  inline std::string GetScope() const {return name_; }

  /// @brief Set concentration of a species in this mode
  /// @tparam StateType Type of the state object (e.g., micm::State)
  /// @param state The state object to modify
  /// @param species The species to set concentration for
  /// @param concentration The concentration value to set [mol m-3]
  /// @param cell The grid cell index (default 0)
  /// @throws std::runtime_error If the species is not found in any phase of this mode
  template<typename StateType>
  void SetConcentration(StateType& state, 
                        const micm::Phase& phase, 
                        const micm::Species& species, 
                        double concentration, 
                        std::size_t cell = 0)
  {
    const std::string species_key = JoinStrings({ name_, phase.name_, species.name_ });
    std::size_t index;

    if (has_initialized_state_idx_)
    {
      try {
        index = state_idx_.state_id_map.at(species_key);
      } catch (const std::out_of_range& e) {
        throw std::runtime_error("Species " + species_key + " not found in state_id_map for '" + name_ + "'");
      }
    }
    else
    {
      auto it = state.variable_map_.find(species_key);
      if (it == state.variable_map_.end()) 
      {
        throw std::runtime_error("Species " + species_key + " not found in state");
      }
      index = it->second;
    }

    state.variables_[cell][index] = concentration;
  }

  /// @brief Set number concentration
  /// @param state StateType type of the state object
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
      std::string number_key = JoinStrings({ name_, AerosolScheme::AEROSOL_MOMENTS_[0]});
      auto it = state.variable_map_.find(number_key);
      if (it == state.variable_map_.end()) 
      {
        throw std::runtime_error("Variable " + number_key + " not found in state");
      }
      index = it->second;
    }

    state.variables_[cell][index] = concentration;
  }

protected:

  /// @brief Indices for accessing mode data in the state vector
  StateIndices state_idx_;
  
  bool has_initialized_state_idx_ = false;
};

} // namespace miam