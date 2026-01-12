// Copyright (C) 2026 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
#include <vector>

namespace miam
{

  /// @brief Join multiple strings with a delimiter (.')
  /// @param names Vector of strings to join
  /// @return Joined string (e.g., "MODE.PHASE.SPECIES")
  std::string JoinStrings(const std::vector<std::string>& names);

}  // namespace miam
