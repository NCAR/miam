// Copyright (C) 2026 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <miam/util/solver_utils.hpp>

#include <string>

namespace miam
{
  micm::System ConfigureSystem(
    const micm::Phase& gas, 
    const std::vector<Mode>& modes,
    const std::vector<Section>& sections)
  {
    micm::SystemParameters params{ .gas_phase_ = gas};

    for (const auto& mode : modes)
    {
      for (const auto& phase : mode.phases_)
      {
        std::string key = Join({ mode.name_, phase.name_ });
        params.phases_[key] = phase;
      }
      params.others_.push_back(Join({ mode.name_, AerosolMoment::NUMBER_CONCENTRATION }));
      params.others_.push_back(Join({ mode.name_, AerosolMoment::DENSITY }));
      params.others_.push_back(Join({ mode.name_, AerosolMoment::RADIUS }));
    }

    // Process sections
    for (const auto& section : sections)
    {
      for (const auto& phase : section.phases_)
      {
        std::string key = Join({ section.name_, phase.name_ });
        params.phases_[key] = phase;
      }
      params.others_.push_back(Join({ section.name_, AerosolMoment::NUMBER_CONCENTRATION }));
      params.others_.push_back(Join({ section.name_, AerosolMoment::DENSITY }));
      params.others_.push_back(Join({ section.name_, AerosolMoment::RADIUS }));
    }

    return micm::System(std::move(params));
  }

}  // namespace miam
