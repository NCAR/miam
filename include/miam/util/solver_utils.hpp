#pragma once

#include "aerosol_model.hpp"

#include <micm/system/phase.hpp>
#include <micm/system/system.hpp>

#include <vector>

namespace miam
{

micm::System ConfigureSystem(
  const micm::Phase& gas_phase, 
  const std::vector<AerosolModel>& aerosol_models);

} // namespace miam
