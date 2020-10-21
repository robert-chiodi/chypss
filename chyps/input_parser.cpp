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

#include <fstream>

#include <mfem/mfem.hpp>

#include "chyps/logger.hpp"

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

void CommonType::SetType(const InputType a_type) { type_id_m = a_type; }

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

bool InputParser::RequiredOption(const OptionType a_type) const {
  return a_type == OptionType::ANY_REQUIRED ||
         a_type == OptionType::COMMAND_LINE_REQUIRED ||
         a_type == OptionType::INPUT_FILE_REQUIRED;
}

bool InputParser::OptionalOption(const OptionType a_type) const {
  return a_type == OptionType::ANY || a_type == OptionType::COMMAND_LINE ||
         a_type == OptionType::INPUT_FILE;
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
  std::vector<std::array<std::string, 2>> negative_bool_statement(
      number_of_options);
  std::vector<CommonType> common_type(number_of_options);

  for (std::size_t n = 0; n < number_of_options; ++n) {
    const auto& string_set = option_description_m[n];
    const bool required = static_cast<int>(option_type_m[string_set[0]]) >= 10;
    if (!this->CommandLineOption(option_type_m[string_set[0]])) {
      continue;
    }
    switch (type_m[string_set[0]]) {
      case InputType::INVALID: {
        std::cout << "Option added to parser with InputType::INVALID type. All "
                     "options must use a valid type (see InputType enum for "
                     "other options."
                  << std::endl;
        // FIXME : Make an exception.
        std::exit(-1);
      }
      case InputType::BOOL: {
        negative_bool_statement[n][0] = string_set[1];
        assert(negative_bool_statement[n][0][0] == '-');
        negative_bool_statement[n][0].insert(1, "no-");
        negative_bool_statement[n][1] = string_set[2];
        assert(negative_bool_statement[n][1][0] == '-');
        assert(negative_bool_statement[n][1][1] == '-');
        negative_bool_statement[n][1].insert(2, "no-");
        auto value = common_type[n].GetPointer<bool*>();
        if (this->OptionalOption(option_type_m[string_set[0]])) {
          *value = static_cast<bool>(parsed_input_m[string_set[0]]);
        }
        mfem_parser.AddOption(value, string_set[1].c_str(),
                              string_set[2].c_str(),
                              negative_bool_statement[n][0].c_str(),
                              negative_bool_statement[n][1].c_str(),
                              string_set[3].c_str(), required);
        break;
      }

      case InputType::INT: {
        auto value = common_type[n].GetPointer<int*>();
        if (this->OptionalOption(option_type_m[string_set[0]])) {
          *value = static_cast<int>(parsed_input_m[string_set[0]]);
        }
        mfem_parser.AddOption(value, string_set[1].c_str(),
                              string_set[2].c_str(), string_set[3].c_str(),
                              required);
        break;
      }

      case InputType::DOUBLE: {
        auto value = common_type[n].GetPointer<double*>();
        if (this->OptionalOption(option_type_m[string_set[0]])) {
          *value = static_cast<double>(parsed_input_m[string_set[0]]);
        }
        mfem_parser.AddOption(value, string_set[1].c_str(),
                              string_set[2].c_str(), string_set[3].c_str(),
                              required);
        break;
      }

      case InputType::STRING: {
        auto value = common_type[n].GetPointer<std::string*>();
        if (this->OptionalOption(option_type_m[string_set[0]])) {
          *value = static_cast<std::string>(parsed_input_m[string_set[0]]);
        }
        string_pointers[n] = value->data();
        mfem_parser.AddOption(&(string_pointers[n]), string_set[1].c_str(),
                              string_set[2].c_str(), string_set[3].c_str(),
                              required);
        break;
      }

      default:
        std::cout << "Cannot transmit InputType to MFEM::OptionsParser. "
                     "Type is:  "
                  << static_cast<int>(type_m[string_set[0]]) << std::endl;
        break;
    }
  }
  mfem_parser.Parse();  // This will also populate the stored parse data for us.
  if (!mfem_parser.Good()) {
    std::cout << "problem with input arguments" << std::endl;
  }

  // Perform a deep copy of strings if they changed from reading in the
  // command line.
  for (std::size_t n = 0; n < number_of_options; ++n) {
    const auto& string_set = option_description_m[n];
    if (!this->CommandLineOption(option_type_m[string_set[0]])) {
      continue;
    }
    switch (type_m[string_set[0]]) {
      case InputType::INVALID: {
        std::cout << "Option added to parser with InputType::INVALID type. All "
                     "options must use a valid type (see InputType enum for "
                     "other options."
                  << std::endl;
        // FIXME : Make an exception.
        std::exit(-1);
      }
      case InputType::BOOL: {
        parsed_input_m[option_description_m[n][0]] =
            static_cast<bool>(common_type[n]);
        break;
      }
      case InputType::INT: {
        parsed_input_m[option_description_m[n][0]] =
            static_cast<int>(common_type[n]);
        break;
      }
      case InputType::DOUBLE: {
        parsed_input_m[option_description_m[n][0]] =
            static_cast<double>(common_type[n]);
        break;
      }
      case InputType::STRING: {
        if ((string_pointers[n] != nullptr &&
             string_pointers[n] != common_type[n].GetPointer<const char*>()) ||
            this->RequiredOption(option_type_m[string_set[0]])) {
          parsed_input_m[option_description_m[n][0]] =
              std::string(string_pointers[n]);
        }
        break;
      }
    }
  }
}

void InputParser::ParseFromFile(const std::string& a_file_name) {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Parsing inputs from file {}", a_file_name);
  std::ifstream myfile(a_file_name.c_str());
  assert(myfile.good());
  // Parse file while stripping comments given by // or /* ...  */
  // Allows comments to be used and parse to still be of proper
  // (standard conforming) JSON
  nlohmann::json read_input =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();
  for (const auto& element : read_input.items()) {
    // FIXME: Make an exception.
    assert(parsed_input_m.contains(element.key()));
    parsed_input_m[element.key()] = element.value();
  }
  SPDLOG_LOGGER_INFO(
      MAIN_LOG, "Finished parsing inputs from {}. Total of {} inputs found.",
      a_file_name, parsed_input_m.size());
}

void InputParser::WriteToFile(const std::string& a_file_name) const {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Writing parser to JSON with filename {}",
                     a_file_name);
  std::ofstream myfile(a_file_name);
  myfile << std::setw(4) << parsed_input_m << std::endl;
  myfile.close();
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Finished writing JSON to file {}", a_file_name);
}

std::vector<std::uint8_t> InputParser::ToBSON(void) const {
  return nlohmann::json::to_bson(parsed_input_m);
}

void InputParser::SetFromBSON(const std::vector<std::uint8_t>& a_bson) {
  parsed_input_m = nlohmann::json::from_bson(a_bson);
}

const nlohmann::json& InputParser::operator[](const std::string& a_name) const {
  // FIXME : Make an exception.
  assert(parsed_input_m.contains(a_name));
  return parsed_input_m[a_name];
}

void InputParser::AddOption(const std::string& a_name,
                            const std::string& a_description) {
  option_description_m.push_back(
      std::array<std::string, 4>{{a_name, "", "", a_description}});
  // FIXME: Make an exception
  assert(!parsed_input_m.contains(a_name));
  parsed_input_m[a_name] = nlohmann::json::object();
  option_type_m[a_name] = MakeOptionRequired(OptionType::INPUT_FILE);
}

bool InputParser::AllOptionsSet(void) const {
  // FIXME: Add test for this.
  // FIXME: Update for using json
  for (const auto& elem : parsed_input_m) {
    if (elem.empty()) {
      return false;
    }
  }
  return true;
}

void InputParser::ClearOptions(void) {
  option_description_m.clear();
  option_type_m.clear();
  type_m.clear();
}

void InputParser::PrintOptions(void) const {
  for (const auto& elem : option_description_m) {
    std::cout << "Option name: " << elem[0] << '\n';
    std::cout << "\tDescription: " << elem[3] << '\n';
    if (this->RequiredOption(option_type_m.at(elem[0]))) {
      std::cout << "\tREQUIRED \n";
    } else {
      std::cout << "\tDefault value: " << parsed_input_m[elem[0]] << '\n';
    }

    if (this->CommandLineOption(option_type_m.at(elem[0]))) {
      std::cout << "\tCommand line flags: " << elem[1] << " " << elem[2];
      if (type_m.at(elem[0]) == InputType::BOOL) {
        std::array<std::string, 2> negative_statements;
        negative_statements[0] = elem[1];
        assert(negative_statements[0][0] == '-');
        negative_statements[0].insert(1, "no-");
        negative_statements[1] = elem[2];
        assert(negative_statements[1][0] == '-');
        assert(negative_statements[1][1] == '-');
        negative_statements[1].insert(2, "no-");

        std::cout << " " << negative_statements[0] << " "
                  << negative_statements[1];
      }
      std::cout << '\n';

    } else {
      std::cout << "\tInput file only \n";
    }
    std::cout << '\n';
  }
  std::cout << std::endl;
}

}  // namespace chyps
