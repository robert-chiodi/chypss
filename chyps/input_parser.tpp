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

template <class Type>
void InputParser::DirectSet(const std::string& a_name, const Type& a_value) {
  parsed_input_m[a_name] = a_value;
}

template <class ValueType>
void InputParser::AddOption(const std::string& a_name,
                            const std::string& a_short_flag,
                            const std::string& a_long_flag,
                            const std::string& a_description,
                            const ValueType& a_default_value,
                            const OptionType a_option_type) {
  option_description_m.push_back(std::array<std::string, 4>{
      {a_name, a_short_flag, a_long_flag, a_description}});
  type_m.push_back(TypeToInputType<ValueType>());
  // FIXME: Make an exception
  assert(!parsed_input_m.contains(a_name));
  parsed_input_m[a_name] = a_default_value;
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
  type_m.push_back(TypeToInputType<ValueType>());
  // FIXME: Make an exception
  assert(!parsed_input_m.contains(a_name));
  parsed_input_m[a_name] = nlohmann::json::object();
  option_type_m.push_back(MakeOptionRequired(a_option_type));
}

template <class ValueType>
void InputParser::AddOption(const std::string& a_name,
                            const std::string& a_description,
                            const ValueType& a_default_value) {
  option_description_m.push_back(
      std::array<std::string, 4>{{a_name, "", "", a_description}});
  // FIXME: Make an exception
  assert(!parsed_input_m.contains(a_name));
  parsed_input_m[a_name] = a_default_value;
  option_type_m.push_back(OptionType::INPUT_FILE);
}

template <class ValueType>
void InputParser::AddOption(const std::string& a_name,
                            const std::string& a_description) {
  option_description_m.push_back(
      std::array<std::string, 4>{{a_name, "", "", a_description}});
  // FIXME: Make an exception
  assert(!parsed_input_m.contains(a_name));
  parsed_input_m[a_name] = nlohmann::json::object();
  option_type_m.push_back(MakeOptionRequired(OptionType::INPUT_FILE));
}

}  // namespace chyps

#endif  // CHYPS_INPUT_PARSER_TPP_
