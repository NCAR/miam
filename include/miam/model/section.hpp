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

   private:
    std::unordered_map<std::string, std::size_t> map_state_id_;  // species - index pairs in state
    std::size_t number_id_;                                      // Index for number concentration
    std::size_t density_id_;                                     // Index for density
    std::size_t radius_id_;                                      // Index for radius

   public:
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
    }

  };

}  // namespace miam