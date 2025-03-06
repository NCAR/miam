// Copyright (C) 2025 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <ostream>
#include <string>

class Message {
public:
  Message(const std::string& message)
  : message_(message) 
  {}

  friend std::ostream &operator<<(std::ostream &ostream, Message& message)
  {
    return message.print(ostream);
  }

private:
  std::string message_;
  std::ostream& print(std::ostream &ostream);
};
