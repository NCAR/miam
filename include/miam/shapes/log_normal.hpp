// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

namespace miam
{
  namespace shape
  {

    /// @brief Log-normal size distribution shape
    /// @details Represents a log-normal distribution shape for aerosol or cloud particle size distributions.
    ///          Characterized by a geometric mean diameter and geometric standard deviation.
    class LogNormal
    {
     public:
      LogNormal() = default;
      ~LogNormal() = default;
    };
  }  // namespace shape
}  // namespace miam