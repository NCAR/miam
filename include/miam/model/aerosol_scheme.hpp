// Copyright (C) 2026 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/util/utils.hpp>

#include <micm/system/phase.hpp>

#include <format>
#include <iostream>
#include <string>
#include <unordered_map>

namespace miam
{

  namespace AerosolMoment
  {
    inline const std::string NUMBER_CONCENTRATION = "NUMBER_CONCENTRATION";
    inline const std::string RADIUS = "RADIUS";
    inline const std::string DENSITY = "DENSITY";
  }  // namespace AerosolMoment

  enum class DistributionType
  {
    SingleMoment,  // tracks mass, fixed radius, number calculated
    TwoMoment      // tracks mass and number concentration, radius calculated
  };

  /// @brief Base AerosolScheme class (common to both modal and sectional)
  class AerosolScheme
  {
   public:
    /// @brief Name of the mode or bin (e.g., "aitken", "accumulation", "large_drop")
    std::string name_;

    AerosolScheme(std::string name)
        : name_(std::move(name))
    {
    }

    virtual ~AerosolScheme() = default;

    std::string GetScope() const;
    std::string Species(const micm::Phase& phase, const micm::Species& species) const;
    std::string NumberConcentration() const;
    std::string Radius() const;
    std::string Density() const;
  };

}  // namespace miam