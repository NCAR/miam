// Copyright (C) 2026 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <miam/util/solver_utils.hpp>

#include <string>
#include <vector>

namespace miam
{

std::string JoinStrings(const std::vector<std::string>& names)
{
  std::string result;
  for (size_t i = 0; i < names.size(); ++i)
  {
    if (!names[i].empty())
    {
      if (!result.empty())
        result += ".";
      result += names[i];
    }
  }
  return result;
}

} // namespace miam