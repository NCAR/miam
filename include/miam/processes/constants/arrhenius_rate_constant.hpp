// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>

namespace miam
{
  namespace process
  {
    namespace constant
    {
      /// @brief Parameters for an Arrhenius rate constant
      /// @details Alias for micm::ArrheniusRateConstantParameters.
      ///          k = A * exp(C / T) * (T / D)^B * (1 + E * P)
      using ArrheniusRateConstantParameters = micm::ArrheniusRateConstantParameters;
    }  // namespace constant
  }  // namespace process
}  // namespace miam
