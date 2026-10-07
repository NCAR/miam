// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/system/conditions.hpp>
#include <micm/util/types.hpp>

#include <cmath>

namespace miam
{
  /// @brief Parameters for the van 't Hoff temperature dependence
  ///        \f$ f(T) = \mathrm{pre\_factor} \cdot \exp\!\left( C \left( \frac{1}{T_0} - \frac{1}{T} \right) \right) \f$
  ///        Arrhenius reaction rate constants and equilibrium constants use this form directly
  ///        Henry's law uses the same form with the opposite temperature trend, expressed by negating \f$ C \f$.
  struct VantHoffParameters
  {
    double A_{ 1.0 };      ///< Value at the reference temperature T0_ [units vary by use]
    double C_{ 0.0 };      ///< Temperature-dependence parameter [K] (e.g. Ea/R)
    double T0_{ 298.15 };  ///< Reference temperature [K]
  };

  /// @brief The van 't Hoff form: A * exp( C * (1/T0 - 1/T) )
  /// @param p Parameters
  /// @param temperature Temperature [K]
  MICM_INLINE_DEVICE_FUNCTION double CalculateVantHoff(const VantHoffParameters& p, double temperature)
  {
    return p.A_ * std::exp(p.C_ * (1.0 / p.T0_ - 1.0 / temperature));
  }

  MICM_INLINE_DEVICE_FUNCTION double Calculate(const VantHoffParameters& p, const micm::Conditions& conditions)
  {
    return CalculateVantHoff(p, conditions.temperature_);
  }
}  // namespace miam
