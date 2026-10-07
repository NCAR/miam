// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/processes/constants/vant_hoff.hpp>

#include <micm/system/conditions.hpp>

#include <cmath>

namespace miam
{
  /// @brief A Henry's Law constant: HLC(T) = HLC_ref * exp( C * ( 1 / T - 1 / T0 ) )
  struct HenrysLawConstant
  {
    double HLC_ref_{ 1.0 };  ///< Henry's Law constant at T0_ [mol m⁻³ Pa⁻¹]
    double C_{ 0.0 };        ///< Temperature dependence parameter [K]
    double T0_{ 298.15 };    ///< Reference temperature [K]
  };

  MICM_INLINE_DEVICE_FUNCTION double Calculate(const HenrysLawConstant& hlc, const micm::Conditions& conditions)
  {
    // The van 't Hoff form with C negated
    return CalculateVantHoff({ hlc.HLC_ref_, -hlc.C_, hlc.T0_ }, conditions.temperature_);
  }
}  // namespace miam
