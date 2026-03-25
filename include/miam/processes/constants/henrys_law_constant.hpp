// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/system/conditions.hpp>

#include <cmath>

namespace miam
{
  namespace process
  {
    namespace constant
    {
      /// @brief Parameters for a Henry's Law constant
      struct HenrysLawConstantParameters
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
      class HenrysLawConstant
      {
       public:
        const HenrysLawConstantParameters parameters_;

        /// @brief Default constructor
        HenrysLawConstant()
            : parameters_()
        {
        }

        /// @brief Constructor with parameters
        /// @param parameters A set of Henry's Law constant parameters
        HenrysLawConstant(const HenrysLawConstantParameters& parameters)
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
          return parameters_.HLC_ref_ * std::exp(parameters_.C_ * (1.0 / temperature - 1.0 / parameters_.T0_));
        }
      };
    }  // namespace constant
  }  // namespace process
}  // namespace miam
