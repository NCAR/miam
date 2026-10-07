// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/processes/constants/vant_hoff.hpp>

#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>
#include <micm/process/rate_constant/rate_constant_functions.hpp>
#include <micm/system/conditions.hpp>
#include <micm/util/types.hpp>

#include <variant>

namespace miam
{
  /// @brief A rate constant: a fixed value, the van 't Hoff form, or the MICM Arrhenius form
  using RateConstant = std::variant<double, VantHoffParameters, micm::ArrheniusRateConstantParameters>;

  MICM_INLINE_DEVICE_FUNCTION double Calculate(double rate_constant, const micm::Conditions& /*conditions*/)
  {
    return rate_constant;
  }

  MICM_INLINE_DEVICE_FUNCTION double Calculate(
      const micm::ArrheniusRateConstantParameters& parameters,
      const micm::Conditions& conditions)
  {
    return micm::CalculateArrhenius(parameters, conditions.temperature_, conditions.pressure_);
  }

  inline double Calculate(const RateConstant& rate_constant, const micm::Conditions& conditions)
  {
    return std::visit([&](const auto& form) { return Calculate(form, conditions); }, rate_constant);
  }
}  // namespace miam
