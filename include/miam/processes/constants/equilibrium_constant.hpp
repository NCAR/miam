// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/processes/constants/vant_hoff.hpp>

#include <micm/system/conditions.hpp>

#include <cmath>

namespace miam
{
  /// @brief An equilibrium constant: K_eq = A * exp( C * ( 1 / T0 - 1 / T ) )
  /// @details A_ is K_eq at T0_ in MIAM's solvent-normalized units.
  ///          To convert from a literature value K_lit in molar units:
  ///          A_ = K_lit / c_H2O^(n_p - n_r), where c_H2O = 55.51 mol/L
  ///          (1000 g/L ÷ 18.015 g/mol).
  ///          After this conversion, A_ is always dimensionless.
  struct EquilibriumConstant
  {
    double A_{ 1.0 };      ///< K_eq at T0_ [dimensionless]
    double C_{ 0.0 };      ///< Temperature dependence parameter [K]
    double T0_{ 298.15 };  ///< Reference temperature [K]
  };

  MICM_INLINE_DEVICE_FUNCTION double Calculate(const EquilibriumConstant& k_eq, const micm::Conditions& conditions)
  {
    return CalculateVantHoff({ k_eq.A_, k_eq.C_, k_eq.T0_ }, conditions.temperature_);
  }
}  // namespace miam
