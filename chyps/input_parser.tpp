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

// template <class ValueType>
// InputType TypeToInputType(void) {
//   static_assert(
//       false,
//       "Invalid type supplied to TypeToInputType. Verify type has existing "
//       "implementation in input_parser.tpp\n");
//   return InputType::INVALID;
// }

// template <class ValueType>
// ValueType CommonType::GetPointer(void) {
//   static_assert(false,
//                 "Invalid type supplied to CommonType::GetPointer. Verify type
//                 " "has existing " "implementation in input_parser.tpp\n");
//   return nullptr;
// }

template <class ValueType>
void InputParser::AddOption(const std::string& a_name,
                            const std::string& a_short_flag,
                            const std::string& a_long_flag,
                            const std::string& a_description,
                            const ValueType& a_default_value,
                            const OptionType a_option_type) {
  option_description_m.push_back(std::array<std::string, 4>{
      {a_name, a_short_flag, a_long_flag, a_description}});
  input_storage_m.push_back(CommonType(a_default_value));
  type_m.push_back(TypeToInputType<ValueType>());
  // FIXME: Make an exception
  assert(parsed_input_m.find(a_name) == parsed_input_m.end());
  parsed_input_m[a_name] = input_storage_m.size() - 1;
  option_type_m.push_back(a_option_type);
  negative_bool_statement_m.push_back(std::array<std::string, 2>());
}

template <class ValueType>
void InputParser::AddOption(const std::string& a_name,
                            const std::string& a_short_flag,
                            const std::string& a_long_flag,
                            const std::string& a_description,
                            const OptionType a_option_type) {
  option_description_m.push_back(std::array<std::string, 4>{
      {a_name, a_short_flag, a_long_flag, a_description}});
  input_storage_m.push_back(CommonType());
  type_m.push_back(TypeToInputType<ValueType>());
  // FIXME: Make an exception
  assert(parsed_input_m.find(a_name) == parsed_input_m.end());
  parsed_input_m[a_name] = input_storage_m.size() - 1;
  option_type_m.push_back(MakeOptionRequired(a_option_type));
  negative_bool_statement_m.push_back(std::array<std::string, 2>());
}

}  // namespace chyps

#endif  // CHYPS_INPUT_PARSER_TPP_
