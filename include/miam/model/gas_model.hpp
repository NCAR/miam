// Copyright (C) 2026 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/util/utils.hpp>

#include <micm/system/phase.hpp>

#include <stdexcept>
#include <string>
#include <unordered_map>
#include <format>

namespace miam
{

  /// @brief Represents a gas phase model
  class GasModel
  {
   public:
    /// @brief The gas phase containing species
    micm::Phase phase_;

    GasModel(micm::Phase phase)
        : phase_(std::move(phase))
    {
    }

    /// @brief Set concentration of a gas species
    /// @tparam StateType Type of the state object (e.g., micm::State)
    /// @param state The state object to modify
    /// @param species The species to set concentration for
    /// @param concentration The concentration value to set [mol m-3]
    /// @param cell The grid cell index (default 0)
    /// @throws std::runtime_error If the species is not found in state
    template<typename StateType>
    void SetConcentration(
        StateType& state,
        const micm::Species& species,
        double concentration,
        std::size_t cell = 0)
    {
      if (state_idx_.empty())
        throw std::runtime_error(std::format("State indices for '{}' not initialized. Call SetStateIndices().", phase_.name_));

      auto it = state_idx_.find(species.name_);
      if (it == state_idx_.end())
      {
        throw std::runtime_error(std::format(
          "Species '{}' not found in state index map for gas phase '{}'.", species.name_, phase_.name_));
      }

      state.variables_[cell][it->second] = concentration;
    }

    /// @brief Get concentration of a gas species
    /// @tparam StateType Type of the state object (e.g., micm::State)
    /// @param state The state object to read from
    /// @param species The species to get concentration for
    /// @param cell The grid cell index (default 0)
    /// @return Concentration value [mol m-3]
    /// @throws std::runtime_error If the species is not found in state
    template<typename StateType>
    double GetConcentration(
        const StateType& state,
        const micm::Species& species,
        std::size_t cell = 0) const
    {
      if (state_idx_.empty())
        throw std::runtime_error(std::format("State indices for '{}' not initialized. Call SetStateIndices().", phase_.name_));

      auto it = state_idx_.find(species.name_);
      if (it == state_idx_.end())
      {
        throw std::runtime_error(std::format(
          "Species '{}' not found in state index map for gas phase '{}'.", species.name_, phase_.name_));
      }

      return state.variables_[cell][it->second];
    }

    /// @brief Set the state indices for accessing gas species in the state vector
    /// @tparam StateType Type of the state object (e.g., micm::State)
    /// @param state The state object containing variable map
    template<typename StateType>
    void InitializeStateIndices(const StateType& state)
    {
      // Find indices for all gas species
      for (const auto& phase_species : phase_.phase_species_)
      {
        const std::string& species_key = phase_species.species_.name_;
        auto it = state.variable_map_.find(species_key);
        if (it == state.variable_map_.end())
        {
          throw std::runtime_error(std::format(
            "Species '{}' not found in state for gas phase '{}'.", species_key, phase_.name_));
        }
        state_idx_[species_key] = it->second;
      }

    }

   private:
    /// @brief Map of species keys to their indices in the state
    std::unordered_map<std::string, std::size_t> state_idx_;
  };

}  // namespace miam
