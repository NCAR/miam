// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

namespace miam
{
  namespace shape
  {

    /// @brief Delta-function size distribution shape
    /// @details Represents a delta-function distribution shape for aerosol or cloud particle size distributions.
    ///          All particles are assumed to have the same size.
    class DeltaFunction
    {
     public:
      DeltaFunction() = default;
      ~DeltaFunction() = default;
    };
  }  // namespace shape
}  // namespace miam