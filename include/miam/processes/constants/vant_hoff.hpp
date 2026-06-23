// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cmath>

namespace miam
{
  /// @brief Parameters for the van 't Hoff temperature dependence
  /// @details Shared functional form behind MIAM's temperature-dependent constants:
  ///
  ///          \f$ f(T) = \mathrm{pre\_factor} \cdot \exp\!\left( C \left( \frac{1}{T_0} - \frac{1}{T} \right) \right) \f$
  ///
  ///          Reaction rate constants and equilibrium constants use this form directly
  ///          (positive \f$ C \f$ ⇒ increases with temperature). Henry's law uses the same
  ///          kernel with the opposite temperature trend, expressed by negating \f$ C \f$
  ///          (see HenryLawConstant).
  struct VantHoffParameters
  {
    double A_{ 1.0 };      ///< Value at the reference temperature T0_ [units vary by use]
    double C_{ 0.0 };      ///< Temperature-dependence parameter [K] (e.g. Ea/R or ΔH/R)
    double T0_{ 298.15 };  ///< Reference temperature [K]
  };

  /// @brief Evaluate the van 't Hoff form: A * exp( C * (1/T0 - 1/T) )
  /// @param p   Parameters
  /// @param temperature Temperature [K]
  inline double CalculateVantHoff(const VantHoffParameters& p, double temperature)
  {
    return p.A_ * std::exp(p.C_ * (1.0 / p.T0_ - 1.0 / temperature));
  }
}  // namespace miam
