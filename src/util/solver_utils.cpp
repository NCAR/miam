// Copyright (C) 2026 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <miam/model/aerosol_scheme.hpp>
#include <miam/util/solver_utils.hpp>

#include <string>

namespace miam
{

  micm::System ConfigureSystem(const GasModel& gas_model, const std::vector<AerosolModel>& aerosol_models)
  {
    micm::SystemParameters params{ .gas_phase_ = gas_model.phase_ };

    for (const auto& model : aerosol_models)
    {
      for (const auto& mode : model.modes_)
      {
        for (const auto& phase : mode.phases_)
        {
          std::string key = JoinStrings({ mode.name_, phase.name_ });
          params.phases_[key] = phase;
        }

        for (const auto& moment : AerosolScheme::AEROSOL_MOMENTS_)
        {
          params.others_.push_back(JoinStrings({ mode.name_, moment }));
        }
      }
      for (const auto& section : model.sections_)
      {
        for (const auto& phase : section.phases_)
        {
          std::string key = JoinStrings({ section.name_, phase.name_ });
          params.phases_[key] = phase;
        }

        for (const auto& moment : AerosolScheme::AEROSOL_MOMENTS_)
        {
          params.others_.push_back(JoinStrings({ section.name_, moment }));
        }
      }
    }

    return micm::System(std::move(params));
  }

}  // namespace miam
