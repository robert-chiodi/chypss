// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chyps/input_parser.hpp"

#include <mfem/mfem.hpp>

namespace chyps {

CommonType::CommonType(void) : type_id_m{InputType::INVALID}, value_m(false) {}
CommonType::CommonType(const bool a_bool)
    : type_id_m{InputType::BOOL}, value_m(a_bool) {}
CommonType::CommonType(const int a_int)
    : type_id_m{InputType::INT}, value_m(a_int) {}
CommonType::CommonType(const double a_double)
    : type_id_m{InputType::DOUBLE}, value_m(a_double) {}
CommonType::CommonType(const char* a_string)
    : type_id_m{InputType::STRING}, value_m(std::string(a_string)) {}
CommonType::CommonType(const std::string& a_string)
    : type_id_m{InputType::STRING}, value_m(a_string) {}

CommonType::CommonType(const CommonType& a_other) {
  type_id_m = a_other.type_id_m;
  switch (type_id_m) {
    case InputType::BOOL: {
      value_m.the_bool = a_other.value_m.the_bool;
      break;
    }
    case InputType::INT: {
      value_m.the_int = a_other.value_m.the_int;
      break;
    }
    case InputType::DOUBLE: {
      value_m.the_double = a_other.value_m.the_double;
      break;
    }
    case InputType::STRING: {
      ::new (&value_m.the_string) decltype(a_other.value_m.the_string)(
          a_other.value_m.the_string);
      break;
    }
  }
}
CommonType::CommonType(CommonType&& a_other) {
  type_id_m = a_other.type_id_m;
  switch (type_id_m) {
    case InputType::BOOL: {
      value_m.the_bool = a_other.value_m.the_bool;
      break;
    }
    case InputType::INT: {
      value_m.the_int = a_other.value_m.the_int;
      break;
    }
    case InputType::DOUBLE: {
      value_m.the_double = a_other.value_m.the_double;
      break;
    }
    case InputType::STRING: {
      ::new (&value_m.the_string) decltype(a_other.value_m.the_string)(
          std::move(a_other.value_m.the_string));
      break;
    }
  }
  a_other.type_id_m = InputType::INVALID;
}

CommonType& CommonType::operator=(const CommonType& a_other) {
  if (&a_other != this) {
    type_id_m = a_other.type_id_m;
    switch (type_id_m) {
      case InputType::BOOL: {
        value_m.the_bool = a_other.value_m.the_bool;
        break;
      }
      case InputType::INT: {
        value_m.the_int = a_other.value_m.the_int;
        break;
      }
      case InputType::DOUBLE: {
        value_m.the_double = a_other.value_m.the_double;
        break;
      }
      case InputType::STRING: {
        ::new (&value_m.the_string) decltype(a_other.value_m.the_string)(
            a_other.value_m.the_string);
        break;
      }
    }
  }
  return *this;
}
CommonType& CommonType::operator=(CommonType&& a_other) {
  if (&a_other != this) {
    type_id_m = a_other.type_id_m;
    switch (type_id_m) {
      case InputType::BOOL: {
        value_m.the_bool = a_other.value_m.the_bool;
        break;
      }
      case InputType::INT: {
        value_m.the_int = a_other.value_m.the_int;
        break;
      }
      case InputType::DOUBLE: {
        value_m.the_double = a_other.value_m.the_double;
        break;
      }
      case InputType::STRING: {
        ::new (&value_m.the_string) decltype(a_other.value_m.the_string)(
            std::move(a_other.value_m.the_string));
        break;
      }
    }
    a_other.type_id_m = InputType::INVALID;
  }
  return *this;
}

CommonType& CommonType::operator=(const bool a_value) {
  type_id_m = InputType::BOOL;
  value_m.the_bool = a_value;
  return *this;
}
CommonType& CommonType::operator=(const int a_value) {
  type_id_m = InputType::INT;
  value_m.the_int = a_value;
  return *this;
}
CommonType& CommonType::operator=(const double a_value) {
  type_id_m = InputType::DOUBLE;
  value_m.the_double = a_value;
  return *this;
}
CommonType& CommonType::operator=(const char* a_value) {
  type_id_m = InputType::STRING;
  ::new (&value_m.the_string) std::string(a_value);
  return *this;
}
CommonType& CommonType::operator=(const std::string& a_value) {
  type_id_m = InputType::STRING;
  ::new (&value_m.the_string) std::string(a_value);
  return *this;
}

CommonType::operator bool(void) const {
  // FIXME : Make an exception
  assert(this->GetType() == InputType::BOOL);
  return value_m.the_bool;
}
CommonType::operator int(void) const {
  // FIXME : Make an exception
  assert(this->GetType() == InputType::INT);
  return value_m.the_int;
}
CommonType::operator double(void) const {
  // FIXME : Make an exception
  assert(this->GetType() == InputType::DOUBLE);
  return value_m.the_double;
}
CommonType::operator std::string(void) const {
  // FIXME : Make an exception
  assert(this->GetType() == InputType::STRING);
  return value_m.the_string;
}

InputType CommonType::GetType(void) const { return type_id_m; }

template <>
bool* CommonType::GetPointer<bool*>(void) {
  assert(this->GetType() == InputType::BOOL);
  return &value_m.the_bool;
}

template <>
int* CommonType::GetPointer<int*>(void) {
  assert(this->GetType() == InputType::INT);
  return &value_m.the_int;
}

template <>
double* CommonType::GetPointer<double*>(void) {
  assert(this->GetType() == InputType::DOUBLE);
  return &value_m.the_double;
}

template <>
const char* CommonType::GetPointer<const char*>(void) {
  assert(this->GetType() == InputType::STRING);
  return value_m.the_string.data();
}

template <>
std::string* CommonType::GetPointer<std::string*>(void) {
  assert(this->GetType() == InputType::STRING);
  return &value_m.the_string;
}

CommonType::~CommonType(void) {
  if (this->GetType() == InputType::STRING) {
    value_m.the_string.std::string::~string();
  }
}

bool InputParser::CommandLineOption(const OptionType a_type) const {
  return a_type == OptionType::ANY || a_type == OptionType::COMMAND_LINE ||
         a_type == OptionType::ANY_REQUIRED ||
         a_type == OptionType::COMMAND_LINE_REQUIRED;
}

bool InputParser::InputFileOption(const OptionType a_type) const {
  return a_type == OptionType::ANY || a_type == OptionType::INPUT_FILE ||
         a_type == OptionType::ANY_REQUIRED ||
         a_type == OptionType::INPUT_FILE_REQUIRED;
}

void InputParser::ParseCL(int argc, char** argv) {
  mfem::OptionsParser mfem_parser(argc, argv);
  const std::size_t number_of_options = option_description_m.size();
  for (std::size_t n = 0; n < number_of_options; ++n) {
    if (!this->CommandLineOption(option_type_m[n])) {
      continue;
    }
    const bool required = static_cast<int>(option_type_m[n]) >= 10;
    // FIXME : Throw exception here.
    const auto& string_set = option_description_m[n];
    switch (input_storage_m[n].GetType()) {
      case InputType::BOOL: {
        std::string bool_no_short = string_set[1];
        assert(bool_no_short[0] == '-');
        bool_no_short.insert(1, "no-");
        std::string bool_no_long = string_set[2];
        assert(bool_no_long[0] == '-');
        assert(bool_no_long[1] == '-');
        bool_no_long.insert(2, "no-");
        auto value = input_storage_m[n].GetPointer<bool*>();
        mfem_parser.AddOption(value, string_set[1].c_str(),
                              string_set[2].c_str(), bool_no_short.c_str(),
                              bool_no_long.c_str(), string_set[3].c_str(),
                              required);
        break;
      }

      case InputType::INT: {
        auto value = input_storage_m[n].GetPointer<int*>();
        mfem_parser.AddOption(value, string_set[1].c_str(),
                              string_set[2].c_str(), string_set[3].c_str(),
                              required);
        break;
      }

      case InputType::DOUBLE: {
        auto value = input_storage_m[n].GetPointer<double*>();
        mfem_parser.AddOption(value, string_set[1].c_str(),
                              string_set[2].c_str(), string_set[3].c_str(),
                              required);
        break;
      }

      case InputType::STRING: {
        auto value = input_storage_m[n].GetPointer<std::string*>();
        value->reserve(CommonType::MAX_STRING_LENGTH);
        const char* ptr = value->data();
        mfem_parser.AddOption(&ptr, string_set[1].c_str(),
                              string_set[2].c_str(), string_set[3].c_str(),
                              required);
        break;
      }

      default:
        std::cout << "Cannot receive value for type "
                  << static_cast<int>(input_storage_m[n].GetType())
                  << std::endl;
        break;
    }
  }
  mfem_parser.Parse();  // This will also populate the stored parse data for us.

  // Need to handle default values too and fill those into the parse
}

const CommonType& InputParser::operator[](const std::string& a_name) const {
  // FIXME : Make an exception.
  auto location = parsed_input_m.find(a_name);
  assert(location != parsed_input_m.end());
  assert(input_storage_m[location->second].GetType() != InputType::INVALID);
  return input_storage_m[location->second];
}

void InputParser::ClearOptions(void) {
  option_description_m.clear();
  option_type_m.clear();
}

}  // namespace chyps
