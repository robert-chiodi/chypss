// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef CHYPS_STRING_MANIPULATION_HPP_
#define CHYPS_STRING_MANIPULATION_HPP_

#include <string>

#include "chyps/debug_assert.hpp"

namespace chyps {

inline std::string ZeroFill(const uint32_t a_number,
                            const uint32_t a_total_length) {
  std::string number = std::to_string(a_number);
  DEBUG_ASSERT(number.size() <= a_total_length, global_assert{},
               DebugLevel::CHEAP{},
               "Total length not large enough to store number: " + number);
  return std::string(a_total_length - number.size(), '0') + number;
}

}  // namespace chyps

#endif  // CHYPS_STRING_MANIPULATION_HPP_
