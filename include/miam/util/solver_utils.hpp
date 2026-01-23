// Copyright (C) 2026 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/model/mode.hpp>
#include <miam/model/section.hpp>

#include <micm/system/phase.hpp>
#include <micm/system/system.hpp>

#include <vector>

namespace miam
{
  micm::System
  ConfigureSystem(const micm::Phase& gas, const std::vector<Mode>& modes = {}, const std::vector<Section>& sections = {});

}  // namespace miam