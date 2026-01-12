// Copyright (C) 2026 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/model/mode.hpp>
#include <miam/model/section.hpp>

#include <string>
#include <vector>

namespace miam
{

  /// @brief Represents an aerosol or cloud containing modes and/or sections
  class AerosolModel
  {
   public:
    std::string name_;
    std::vector<Mode> modes_;
    std::vector<Section> sections_;

    AerosolModel(std::string name, std::vector<Mode> modes = {}, std::vector<Section> sections = {})
        : name_(std::move(name)),
          modes_(std::move(modes)),
          sections_(std::move(sections))
    {
    }
  };

}  // namespace miam