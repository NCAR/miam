// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <miam/util/error.hpp>

#include <stdexcept>
#include <string>

namespace miam
{
  struct MiamException : public std::runtime_error
  {
    const char* category_;  // Must point to string literal
    int code_;

    MiamException(const char* category, int code, const std::string& message)
        : std::runtime_error(message),
          category_(category),
          code_(code)
    {
    }

    const char* Category() const noexcept
    {
      return category_;
    }

    int Code() const noexcept
    {
      return code_;
    }
  };

}  // namespace miam
