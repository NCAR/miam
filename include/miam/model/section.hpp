// Copyright (C) 2026 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/model/aerosol_scheme.hpp>
#include <miam/util/utils.hpp>

#include <string>
#include <vector>
#include <format>

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


    std::unordered_map<std::string, std::size_t> map_state_id_;  // species - index pairs in state
    std::size_t number_id_;                                      // Index for number concentration
    std::size_t density_id_;                                     // Index for density
    std::size_t radius_id_;                                      // Index for radius

    /// @brief Construct a Section with specified physical properties
    Section(
        std::string name,
        std::vector<micm::Phase> phases,
        DistributionType distribution,
        double min_diameter,
        double max_diameter)
        : AerosolScheme(std::move(name)),
          phases_(std::move(phases)),
          distribution_(distribution),
          min_diameter_(min_diameter),
          max_diameter_(max_diameter)
    {
      // if (distribution == DistributionType::SingleMoment)
      //   is_radius_fixed_ = true;
    }

    /// @brief Set the state indices for accessing section variables in the state vector
    /// @tparam StateType Type of the state object (e.g., micm::State)
    /// @param state The state object containing variable map
    /// @throws std::runtime_error If keys are not found in state
    template<typename StateType>
    void SetStateIndices(StateType& state)
    {

            std::cout <<  "void SetStateIndices(const StateType& state) SECTION" << std::endl;


      // // Find mass indices for all species in all phases
      // for (const auto& phase : phases_)
      // {
      //   for (const auto& phase_species : phase.phase_species_)
      //   {
      //     // NAME: SECTION.PHASE.SPECIES
      //     std::string species_key = JoinStrings({ name_, phase.name_, phase_species.species_.name_ });
      //     auto species_it = state.variable_map_.find(species_key);
      //     if (species_it == state.variable_map_.end())
      //     {
      //       throw std::runtime_error(std::format("Species '{}' not found in state for '{}'", species_key, name_));
      //     }
      //     state_idx_.state_id_map[species_key] = species_it->second;
      //   }
      // }

      // // Find number concentration index
      // std::string number_key = JoinStrings({ name_, AerosolModel::AEROSOL_MOMENTS_[0] });
      // auto number_it = state.variable_map_.find(number_key);
      // if (number_it == state.variable_map_.end())
      // {
      //   throw std::runtime_error(std::format("Variable '{}' not found in state for '{}'", number_key, name_));
      // }
      // state_idx_.number_id = number_it->second;

      // // Find density index
      // std::string density_key = JoinStrings({ name_, AerosolModel::AEROSOL_MOMENTS_[1] });
      // auto density_it = state.variable_map_.find(density_key);
      // if (density_it == state.variable_map_.end())
  //     {
  //       throw std::runtime_error(std::format("Variable '{}' not found in state for '{}'", density_key, name_));
  //     }
  //     state_idx_.density_id = density_it->second;
  //   // }


  //   ///
  //   for (auto& [key, value] : state_idx_.state_id_map)
  //   {
  //     std::cout << "key: " << key << " value: " << value << std::endl;
  //   }
  }
  };

}  // namespace miam