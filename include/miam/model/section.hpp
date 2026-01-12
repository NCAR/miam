// Copyright (C) 2026 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/model/aerosol_scheme.hpp>
#include <miam/util/utils.hpp>

#include <vector>
#include <string>

namespace miam
{

class Section : public AerosolScheme
{
public:

  /// @brief Phases associated with this section (e.g., aqueous, organic)
  std::vector<micm::Phase> phases_;
  /// @brief Type of distribution (single moment or two moment)
  DistributionType distribution_;
  /// @brief Minimum diameter [m]
  double min_diameter_;
  /// @brief Maximum diameter [m]
  double max_diameter_;

  /// @brief Construct a Section with specified physical properties
  Section(std::string name,
          std::vector<micm::Phase> phases,
          DistributionType distribution,
          double min_diameter,
          double max_diameter)
  : AerosolScheme(std::move(name)),
    phases_(std::move(phases)),
    distribution_(distribution),
    min_diameter_(min_diameter),
    max_diameter_(max_diameter)
  {}

  /// @brief Set the state indices for accessing section variables in the state vector
  /// @tparam StateType Type of the state object (e.g., micm::State)
  /// @param state The state object containing variable map
  template<typename StateType>
  void SetStateIndices(const StateType& state) 
  { 
    // Find mass indices for all species in all phases
    for (const auto& phase : phases_) {
      for (const auto& phase_species : phase.phase_species_) {
        // NAME: SECTION.PHASE.SPECIES
        std::string key = JoinStrings({ name_, phase.name_, phase_species.species_.name_ });
        state_idx_.state_id_map[key] = state.variable_map_.at(key);
      }
    }

    // Find number concentration index
    std::string number_key = JoinStrings({ name_, AerosolScheme::AEROSOL_MOMENTS_[0]});
    auto number_it = state.variable_map_.find(number_key);
    if (number_it == state.variable_map_.end()) {
      throw std::runtime_error("Variable " + number_key + " not found in state");
    }
    state_idx_.number_id = number_it->second;

    // Find density index
    std::string density_key = JoinStrings({ name_, AerosolScheme::AEROSOL_MOMENTS_[1]});
    auto density_it = state.variable_map_.find(density_key);
    if (density_it == state.variable_map_.end()) {
      throw std::runtime_error("Variable " + density_key + " not found in state");
    }
    state_idx_.density_id = density_it->second;

    has_initialized_state_idx_ = true; 
  }

};

} // namespace miam