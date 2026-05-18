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
      /// @brief Parameters for an equilibrium constant
      /// @details The pre-exponential factor A_ is the equilibrium constant at the
      ///          reference temperature T0_ in MIAM's solvent-normalized units.
      ///          To convert from a literature value K_lit in molar units:
      ///          A_ = K_lit / c_H2O^(n_p - n_r), where c_H2O = 55.556 mol/L.
      ///          After this conversion, A_ is always dimensionless.
      struct EquilibriumConstantParameters
      {
        /// @brief Pre-exponential factor [dimensionless]
        /// @details K_eq at T0_ in MIAM's solvent-normalized units.
        ///          Conversion from literature: A_ = K_lit / c_H2O^(n_p - n_r)
        double A_{ 1.0 };
        /// @brief Temperature dependence parameter [K]
        double C_{ 0.0 };
        /// @brief Reference temperature [K]
        double T0_{ 298.15 };
      };

      /// @brief An equilibrium constant dependent on temperature
      /// @details Calculates equilibrium constant as:
      ///          K_eq = A * exp( C * ( 1 / T0 - 1 / T ) )
      class EquilibriumConstant
      {
       public:
        const EquilibriumConstantParameters parameters_;

        /// @brief Default constructor
        EquilibriumConstant()
            : parameters_()
        {
        }

        /// @brief Constructor with parameters
        /// @param parameters A set of equilibrium constant parameters
        EquilibriumConstant(const EquilibriumConstantParameters& parameters)
            : parameters_(parameters)
        {
        }

        /// @brief Calculate the equilibrium constant
        /// @param conditions The current environmental conditions of the chemical system
        /// @return An equilibrium constant based off of the conditions in the system
        double Calculate(const micm::Conditions& conditions) const
        {
          return Calculate(conditions.temperature_);
        }

        /// @brief Calculate the equilibrium constant
        /// @param temperature The temperature [K]
        /// @return An equilibrium constant
        double Calculate(const double& temperature) const
        {
          return parameters_.A_ * std::exp(parameters_.C_ * (1.0 / parameters_.T0_ - 1.0 / temperature));
        }
      };
    }  // namespace constant
  }  // namespace process
}  // namespace miam
