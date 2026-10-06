// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

namespace miam
{
  /// @brief Aerosol/cloud particle properties that representations can provide
  enum class AerosolProperty
  {
    EffectiveRadius,      // [m]
    NumberConcentration,  // [# m^-3]
    PhaseVolumeFraction   // [dimensionless, 0-1]
  };
}  // namespace miam
