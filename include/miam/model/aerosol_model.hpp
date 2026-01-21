// Copyright (C) 2026 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/model/mode.hpp>
#include <miam/model/section.hpp>

#include <string>
#include <vector>

#include <iostream>

namespace miam
{

  /// @brief Represents an aerosol or cloud containing modes and/or sections
  class AerosolModel
  {
   public:
    std::string name_;
    std::vector<Mode> modes_;
    std::vector<Section> sections_;

    inline static const std::vector<std::string> AEROSOL_MOMENTS_ = { "NUMBER_CONCENTRATION", "DENSITY", "RADIUS" };

   protected:
    std::unordered_map<std::string, std::size_t> map_state_id_;  // species - index pairs in state

   public:
    AerosolModel(std::string name, std::vector<Mode> modes = {}, std::vector<Section> sections = {})
        : name_(std::move(name)),
          modes_(std::move(modes)),
          sections_(std::move(sections))
    {
    }


    template<typename SchemeType>
    std::string GetScope(SchemeType& scheme)
    {
      return name_ + "." + scheme.name_;
    }

    /// @brief Wrapper for initializing state indices for modes and sections
    /// @tparam StateType Type of the state object (e.g., micm::State)
    /// @param state The state object containing variable map
    template<typename StateType>
    void InitializeStateIndices(StateType& state)
    {
      for (auto& mode : modes_)
        SetStateIndices(state, mode);
      for (auto& section : sections_)
        SetStateIndices(state, section);
    }

    /// @brief Set the state indices for accessing mode variables in the state vector
    /// @tparam StateType Type of the state object (e.g., micm::State)
    /// @param state The state object containing variable map
    /// @throws std::runtime_error If keys are not found in state
    template<typename StateType, typename SchemeType>
    void SetStateIndices(StateType& state, SchemeType& scheme)
    {
      // Find mass indices for all species in all phases
      for (const auto& phase : scheme.phases_)
      {
        for (const auto& phase_species : phase.phase_species_)
        {
          // NAME: SECTION.PHASE.SPECIES
          const std::string species_key = JoinStrings({ name_, scheme.name_, phase.name_, phase_species.species_.name_ });
          auto it = state.variable_map_.find(species_key);
          if (it == state.variable_map_.end())
          {
            throw std::runtime_error(std::format("Species '{}' not found in state for '{}'.", species_key, scheme.name_));
          }
          map_state_id_[species_key] = it->second;
          // Initialize each schemes state map
          scheme.map_state_id_[species_key] = it->second;
        }
      }

      // Find number concentration index
      const std::string number_key = JoinStrings({ name_, scheme.name_, AerosolModel::AEROSOL_MOMENTS_[0] });
      auto number_it = state.variable_map_.find(number_key);
      if (number_it == state.variable_map_.end())
      {
        throw std::runtime_error(std::format("Variable '{}' not found in state for '{}'", number_key, scheme.name_));
      }
      map_state_id_[number_key] = number_it->second;
      scheme.number_id_ = number_it->second;

      // Find density index
      const std::string density_key = JoinStrings({ name_, scheme.name_, AerosolModel::AEROSOL_MOMENTS_[1] });
      auto density_it = state.variable_map_.find(density_key);
      if (density_it == state.variable_map_.end())
      {
        throw std::runtime_error(std::format("Variable '{}' not found in state for '{}'", density_key, scheme.name_));
      }
      map_state_id_[density_key] = density_it->second;
      scheme.density_id_ = density_it->second;

      // Find radius index
      const std::string radius_key = JoinStrings({ name_, scheme.name_, AerosolModel::AEROSOL_MOMENTS_[2] });
      auto radius_it = state.variable_map_.find(radius_key);
      if (radius_it == state.variable_map_.end())
      {
        throw std::runtime_error(std::format("Variable '{}' not found in state for '{}'", radius_key, scheme.name_));
      }
      map_state_id_[radius_key] = radius_it->second;
      scheme.radius_id_ = radius_it->second;

      // debug
      for (auto& [key, value] : map_state_id_)
      {
        std::cout << "key: " << key << " value: " << value << std::endl;
      }
    }

    /// @brief Set concentration for a species in a specific mode/section and phase
    /// @tparam StateType Type of the state object (e.g., micm::State)
    /// @tparam SchemeType Type of the aerosol scheme (Mode or Section)
    /// @param state The state object containing variable data
    /// @param scheme The mode or section to set concentration for
    /// @param phase The phase containing the species
    /// @param species The species to set concentration for
    /// @param concentration The concentration value to set [mol m-3]
    /// @param cell The index of the grid cell (default 0)
    template<typename StateType, typename SchemeType>
    void SetConcentration(
        StateType& state,
        const SchemeType& scheme,
        const micm::Phase& phase,
        const micm::Species& species,
        double concentration,
        std::size_t cell = 0)
  {
    if (map_state_id_.empty())
      throw std::runtime_error(std::format("State indices for '{}' not initialized.", name_));

    const std::string species_key = JoinStrings({ name_, scheme.name_, phase.name_, species.name_ });

    auto it = map_state_id_.find(species_key);
    if (it == map_state_id_.end())
    {
      throw std::runtime_error(std::format("Species '{}' not found in state index map for '{}'.", species_key, scheme.name_));
    }

    state.variables_[cell][it->second] = concentration;
  }
    
    /// @brief Set number concentration
    /// @param state StateType Type of the state object
    /// @param concentration The number concentration value to set [# m-3]
    /// @param cell The grid cell index (default 0)
    /// @throws std::runtime_error If the number concentration variable is not found in state
    template<typename StateType, typename SchemeType>
    void SetNumberConcentration(StateType& state, SchemeType& scheme, double concentration, std::size_t cell = 0)
    {
      if (map_state_id_.empty())
        throw std::runtime_error(std::format("State indices for '{}' not initialized.", name_));

      const std::string number_key = JoinStrings({ name_, scheme.name_, AEROSOL_MOMENTS_[0] });

      auto it = map_state_id_.find(number_key);
      if (it == map_state_id_.end())
      {
        throw std::runtime_error(std::format("Variable '{}' not found in state for '{}'", number_key, scheme.name_));
      }

      state.variables_[cell][it->second] = concentration;
    }

    /// @brief Set the effective radius in the state for this mode
    /// @tparam StateType Type of the state object (e.g., micm::State)
    /// @param state The state object containing variable map
    /// @param cell The index of the grid cell (default 0)
    /// @throws std::runtime_error If radius key is not found in state
    template<typename StateType, typename SchemeType>
    void SetRadius(StateType& state, SchemeType& scheme, std::size_t cell = 0)
    {
      if (map_state_id_.empty())
        throw std::runtime_error(std::format("State indices for '{}' not initialized.", name_));
      
      const std::string radius_key = JoinStrings({ name_, scheme.name_, AEROSOL_MOMENTS_[2] });

      auto it = map_state_id_.find(radius_key);
      if (it == map_state_id_.end())
      {
        throw std::runtime_error(std::format("Variable '{}' not found in state for '{}'", radius_key, scheme.name_));
      }

      scheme.SetRadius(state, cell);
    } 

  };
}  // namespace miam