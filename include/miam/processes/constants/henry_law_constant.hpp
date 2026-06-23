// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/processes/constants/vant_hoff.hpp>

#include <micm/system/conditions.hpp>

#include <cmath>

namespace miam
{
  /// @brief Parameters for a Henry's Law constant
  struct HenryLawConstantParameters
  {
    /// @brief Henry's Law constant at reference temperature [mol m⁻³ Pa⁻¹]
    double HLC_ref_{ 1.0 };
    /// @brief Temperature dependence parameter [K]
    double C_{ 0.0 };
    /// @brief Reference temperature [K]
    double T0_{ 298.15 };
  };

  /// @brief A Henry's Law constant dependent on temperature
  /// @details Calculates Henry's Law constant as:
  ///          HLC(T) = HLC_ref * exp( C * ( 1 / T - 1 / T0 ) )
  class HenryLawConstant
  {
   public:
    const HenryLawConstantParameters parameters_;

    /// @brief Default constructor
    HenryLawConstant()
        : parameters_()
    {
    }

    /// @brief Constructor with parameters
    /// @param parameters A set of Henry's Law constant parameters
    HenryLawConstant(const HenryLawConstantParameters& parameters)
        : parameters_(parameters)
    {
    }

    /// @brief Calculate the Henry's Law constant
    /// @param conditions The current environmental conditions of the chemical system
    /// @return A Henry's Law constant based off of the conditions in the system [mol m⁻³ Pa⁻¹]
    double Calculate(const micm::Conditions& conditions) const
    {
      return Calculate(conditions.temperature_);
    }

    /// @brief Calculate the Henry's Law constant
    /// @param temperature The temperature [K]
    /// @return A Henry's Law constant [mol m⁻³ Pa⁻¹]
    double Calculate(const double& temperature) const
    {
      // HLC = HLC_ref * exp( C * (1/T - 1/T0) ) = HLC_ref * exp( -C * (1/T0 - 1/T) ),
      // i.e. the same van 't Hoff kernel with C negated (solubility rises as T falls).
      return CalculateVantHoff({ parameters_.HLC_ref_, -parameters_.C_, parameters_.T0_ }, temperature);
    }
  };
}  // namespace miam