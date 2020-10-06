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

template <>
InputType TypeToInputType<bool>(void) {
  return InputType::BOOL;
}

template <>
InputType TypeToInputType<int>(void) {
  return InputType::INT;
}

template <>
InputType TypeToInputType<double>(void) {
  return InputType::DOUBLE;
}

template <>
InputType TypeToInputType<const char*>(void) {
  return InputType::STRING;
}

template <>
InputType TypeToInputType<std::string>(void) {
  return InputType::STRING;
}

OptionType MakeOptionRequired(const OptionType a_option) {
  switch (a_option) {
    case OptionType::ANY: {
      return OptionType::ANY_REQUIRED;
    }

    case OptionType::COMMAND_LINE: {
      return OptionType::COMMAND_LINE_REQUIRED;
    }

    case OptionType::INPUT_FILE: {
      return OptionType::INPUT_FILE_REQUIRED;
    }

    default:
      // FIXME : Make actual exit function and use.
      std::cout << "Unknown option type: " << static_cast<int>(a_option);
      std::exit(-1);
  }
}

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
    case InputType::INVALID: {
      break;
    }
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
    case InputType::INVALID: {
      break;
    }
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
      case InputType::INVALID: {
        break;
      }
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
      case InputType::INVALID: {
        break;
      }
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
  type_id_m = InputType::BOOL;
  return &value_m.the_bool;
}

template <>
int* CommonType::GetPointer<int*>(void) {
  type_id_m = InputType::INT;
  return &value_m.the_int;
}

template <>
double* CommonType::GetPointer<double*>(void) {
  type_id_m = InputType::DOUBLE;
  return &value_m.the_double;
}

template <>
const char* CommonType::GetPointer<const char*>(void) {
  if (this->GetType() != InputType::STRING) {
    ::new (&value_m.the_string) std::string();
  }
  type_id_m = InputType::STRING;
  return value_m.the_string.c_str();
}

template <>
std::string* CommonType::GetPointer<std::string*>(void) {
  if (this->GetType() != InputType::STRING) {
    ::new (&value_m.the_string) std::string();
  }
  type_id_m = InputType::STRING;
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

// TODO Command lines specified as ANY but required should not fail
// if not found here but potentially found later in an input file (or separate
// manner of supply)
//
// MFEM simply moves the const char* you provide to the const char* entry in
// argv. This will not effect the std::string we are actually storing, so we
// need to do some special things here to handle string input. Namely, we create
// the vector of string pointers and point them to our initial strings in
// input_storage_m. After parsing the command line with mfem, we then compare
// and deep copy any who have new pointer locations.
void InputParser::ParseCL(int argc, char** argv) {
  mfem::OptionsParser mfem_parser(argc, argv);
  const std::size_t number_of_options = option_description_m.size();
  std::vector<const char*> string_pointers(number_of_options, nullptr);
  for (std::size_t n = 0; n < number_of_options; ++n) {
    if (!this->CommandLineOption(option_type_m[n])) {
      continue;
    }
    const bool required = static_cast<int>(option_type_m[n]) >= 10;
    const auto& string_set = option_description_m[n];
    switch (type_m[n]) {
      case InputType::INVALID: {
        std::cout << "Option added to parser with InputType::INVALID type. All "
                     "options must use a valid type (see InputType enum for "
                     "other options."
                  << std::endl;
        // FIXME : Make an exception.
        std::exit(-1);
      }
      case InputType::BOOL: {
        negative_bool_statement_m[n][0] = string_set[1];
        assert(negative_bool_statement_m[n][0][0] == '-');
        negative_bool_statement_m[n][0].insert(1, "no-");
        negative_bool_statement_m[n][1] = string_set[2];
        assert(negative_bool_statement_m[n][1][0] == '-');
        assert(negative_bool_statement_m[n][1][1] == '-');
        negative_bool_statement_m[n][1].insert(2, "no-");
        auto value = input_storage_m[n].GetPointer<bool*>();
        mfem_parser.AddOption(value, string_set[1].c_str(),
                              string_set[2].c_str(),
                              negative_bool_statement_m[n][0].c_str(),
                              negative_bool_statement_m[n][1].c_str(),
                              string_set[3].c_str(), required);
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
        string_pointers[n] = value->data();
        mfem_parser.AddOption(&(string_pointers[n]), string_set[1].c_str(),
                              string_set[2].c_str(), string_set[3].c_str(),
                              required);
        break;
      }

      default:
        std::cout
            << "Cannot transmit InputType to MFEM::OptionsParser. Type is:  "
            << static_cast<int>(input_storage_m[n].GetType()) << std::endl;
        break;
    }
  }
  mfem_parser.Parse();  // This will also populate the stored parse data for us.

  // Perform a deep copy of strings if they changed from reading in the command
  // line.
  for (std::size_t n = 0; n < number_of_options; ++n) {
    if (!this->CommandLineOption(option_type_m[n]) ||
        type_m[n] != InputType::STRING) {
      continue;
    }
    auto value = input_storage_m[n].GetPointer<std::string*>();
    if (string_pointers[n] != nullptr && string_pointers[n] != value->data()) {
      *value = std::string(string_pointers[n]);
    }
  }
}

const CommonType& InputParser::operator[](const std::string& a_name) const {
  auto location = parsed_input_m.find(a_name);
  // FIXME : Make an exception.
  assert(location != parsed_input_m.end());
  assert(input_storage_m[location->second].GetType() != InputType::INVALID);
  return input_storage_m[location->second];
}

bool InputParser::AllOptionsSet(void) const {
  // TODO Add test for this.
  for (const auto& elem : input_storage_m) {
    if (elem.GetType() == InputType::INVALID) {
      return false;
    }
  }
  return true;
}

void InputParser::ClearOptions(void) {
  option_description_m.clear();
  option_type_m.clear();
  type_m.clear();
  negative_bool_statement_m.clear();
}

}  // namespace chyps
