// Copyright (C) 2025 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <miam/util/print.hpp>

std::ostream& Message::print(std::ostream &ostream)
std::ostream &Message::print(std::ostream &os) 
{
  os << "This is MIAM repository." << std::endl;
  os << message_;

  return os;
}
