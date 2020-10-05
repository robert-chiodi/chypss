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

#include <cassert>

#include <iostream>

namespace chyps {

template <class ValueType>
void InputParser::AddOption(const std::string& a_name,
                            const std::string& a_short_flag,
                            const std::string& a_long_flag,
                            const std::string& a_description,
                            const ValueType& a_default_value,
                            const OptionType a_option_type) {
  option_description_m.push_back(std::array<std::string, 4>{
      {a_name, a_short_flag, a_long_flag, a_description}});
  default_value_m.push_back(CommonType(a_default_value));
  input_storage_m.push_back(CommonType());
  parsed_input_m[a_name] = input_storage_m.size() - 1;
  // FIXME: Make an exception
  assert(static_cast<int>(a_option_type) < 10);  // Not Required
  option_type_m.push_back(a_option_type);
}

template <class ValueType>
void InputParser::AddOption(const std::string& a_name,
                            const std::string& a_short_flag,
                            const std::string& a_long_flag,
                            const std::string& a_description,
                            const OptionType a_option_type) {
  option_description_m.push_back(std::array<std::string, 4>{
      {a_name, a_short_flag, a_long_flag, a_description}});
  default_value_m.push_back(CommonType());
  input_storage_m.push_back(CommonType());
  parsed_input_m[a_name] = input_storage_m.size() - 1;
  // FIXME: Make an exception
  assert(static_cast<int>(a_option_type) >= 10);  // Required
  option_type_m.push_back(a_option_type);
}

}  // namespace chyps

#endif  // CHYPS_INPUT_PARSER_TPP_
