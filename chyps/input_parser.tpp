// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef CHYPS_INPUT_PARSER_TPP_
#define CHYPS_INPUT_PARSER_TPP_

#include "chyps/debug_assert.hpp"

namespace chyps {

template <class Type>
void InputParser::DirectSet(const std::string& a_name, const Type& a_value) {
  parsed_input_m[a_name] = a_value;
}

template <class ValueType>
void InputParser::AddOption(const std::string& a_name,
                            const std::string& a_description,
                            const ValueType& a_default_value) {
  DEBUG_ASSERT(
      !option_description_m.Contains(a_name), global_assert{},
      DebugLevel::CHEAP{},
      "Option \"" + a_name + "\" already exists in InputParser object");
  option_description_m[a_name] = a_description;
  default_values_m[a_name] = a_default_value;
}

}  // namespace chyps

#endif  // CHYPS_INPUT_PARSER_TPP_
