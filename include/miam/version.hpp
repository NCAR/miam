// Copyright (C) 2025 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#pragma once

#ifdef __cplusplus
namespace miam {
extern "C" {
#endif

  inline const char* GetMiamVersion()
  {
    return "0.0.0";
  }
  inline unsigned GetMiamVersionMajor()
  {
    return 0;
  }
  inline unsigned GetMiamVersionMinor()
  {
    return 0+0;
  }
  inline unsigned GetMiamVersionPatch()
  {
    return 0+0;
  }
  inline unsigned GetMiamVersionTweak()
  {
    return +0;
  }

#ifdef __cplusplus
}  // extern "C"
}  // namespace miam
#endif
