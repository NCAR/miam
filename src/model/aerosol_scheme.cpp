// Copyright (C) 2026 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <miam/model/aerosol_scheme.hpp>
#include <miam/util/utils.hpp>

#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>

#include <string>

namespace miam
{
  std::string AerosolScheme::GetScope() const
  {
    return name_;
  }

  std::string AerosolScheme::Species(const micm::Phase& phase, const micm::Species& species) const
  {
    return Join({ name_, phase.name_, species.name_ });
  }

  std::string AerosolScheme::NumberConcentration() const
  {
    return Join({ name_, AerosolMoment::NUMBER_CONCENTRATION });
  }

  std::string AerosolScheme::Radius() const
  {
    return Join({ name_, AerosolMoment::RADIUS });
  }

  std::string AerosolScheme::Density() const
  {
    return Join({ name_, AerosolMoment::DENSITY });
  }
}  // namespace miam