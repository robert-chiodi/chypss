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

#include <algorithm>
#include <fstream>

#include <mfem/mfem.hpp>

#include "chyps/logger.hpp"

namespace chyps {

nlohmann::json& DirectoryJSON::operator[](const std::string& a_name) {
  nlohmann::json* lowest_dir;
  std::string lowest_name;
  std::tie(lowest_dir, lowest_name) = this->LowestObject(a_name);
  return (*lowest_dir)[lowest_name];
}

const nlohmann::json& DirectoryJSON::operator[](
    const std::string& a_name) const {
  const nlohmann::json* lowest_dir;
  std::string lowest_name;
  std::tie(lowest_dir, lowest_name) = this->LowestObject(a_name);
  assert(lowest_dir->contains(lowest_name));
  return (*lowest_dir)[lowest_name];
}

bool DirectoryJSON::Contains(const std::string& a_name) const {
  const nlohmann::json* lowest_dir;
  std::string lowest_name;
  std::tie(lowest_dir, lowest_name) = this->LowestObject(a_name);
  return lowest_dir->contains(lowest_name);
}

std::pair<nlohmann::json*, std::string> DirectoryJSON::LowestObject(
    const std::string a_name) {
  nlohmann::json* lowest_dir = &json_m;
  std::size_t name_start_location = 0;
  std::size_t slash_location = a_name.find('/');
  while (slash_location != a_name.npos) {
    lowest_dir = &((*lowest_dir)[a_name.substr(
        name_start_location, slash_location - name_start_location)]);
    name_start_location = slash_location + 1;
    slash_location = a_name.find('/', name_start_location);
  }
  return std::make_pair(lowest_dir, a_name.substr(name_start_location));
}

std::pair<const nlohmann::json*, std::string> DirectoryJSON::LowestObject(
    const std::string a_name) const {
  const nlohmann::json* lowest_dir = &json_m;
  std::size_t name_start_location = 0;
  std::size_t slash_location = a_name.find('/');
  while (slash_location != a_name.npos) {
    lowest_dir = &((*lowest_dir)[a_name.substr(
        name_start_location, slash_location - name_start_location)]);
    name_start_location = slash_location + 1;
    slash_location = a_name.find('/', name_start_location);
  }
  return std::make_pair(lowest_dir, a_name.substr(name_start_location));
}

void InputParser::ParseCL(int argc, char** argv) {
  assert(argc % 2 ==
         0);  // Always expect flag, value in ordering, therefore must be even

  for (int n = 0; n < argc; n += 2) {
    std::string argument(argv[n]);
    std::string value(argv[n + 1]);

    assert(argument[0] == '-');
    const std::string name = argument.substr(1);

    assert(option_description_m.Contains(name));
    if (isdigit(value[0]) || value[0] == '-' || value[0] == '+') {
      std::size_t decimal_location = 0;
      decimal_location = value.find('.');
      if (decimal_location == value.npos) {
        // Was an int
        parsed_input_m[name] = std::stoi(value);
      } else {
        parsed_input_m[name] = std::stod(value);
      }
    } else {
      parsed_input_m[name] = value;
    }
  }
}

bool InputParser::ParseFromFile(const std::string& a_file_name) {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Parsing inputs from file {}", a_file_name);
  std::ifstream myfile(a_file_name.c_str());
  if (!myfile.good()) {
    return false;
  }
  // Parse file while stripping comments given by // or /* ...  */
  // Allows comments to be used and parse to still be of proper
  // (standard conforming) JSON
  const auto read_input = nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  this->RecursiveInsert(read_input, parsed_input_m.json_m);

  SPDLOG_LOGGER_INFO(
      MAIN_LOG, "Finished parsing inputs from {}. Total of {} inputs found.",
      a_file_name, read_input.size());
  return true;
}

void InputParser::AddOptionNoDefault(const std::string& a_name,
                                     const std::string& a_description,
                                     const bool a_required,
                                     const std::string& a_dependency) {
  option_description_m[a_name] = a_description;
  parsed_input_m[a_name] = nlohmann::json::object();
  assert(option_required_status_m.find(a_name) ==
         option_required_status_m.end());
  if (a_dependency != "") {
    option_required_status_m[a_name] = static_cast<int>(dependencies_m.size());
    dependencies_m.push_back(a_dependency);
  } else {
    option_required_status_m[a_name] = a_required ? -1 : -2;
  }
}

void InputParser::RecursiveInsert(const nlohmann::json& a_input,
                                  nlohmann::json& a_parsed_nest) {
  // FIXME : Need way to make sure there is a corresponding option for this
  // input.
  for (nlohmann::json::const_iterator it = a_input.begin(); it != a_input.end();
       ++it) {
    if (it->is_object()) {
      // Is a nested object object
      this->RecursiveInsert(*it, a_parsed_nest[it.key()]);
    } else {
      a_parsed_nest[it.key()] = it.value();
    }
  }
}

void InputParser::WriteToFile(const std::string& a_file_name) const {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Writing parser to JSON with filename {}",
                     a_file_name);
  std::ofstream myfile(a_file_name);
  myfile << std::setw(4) << parsed_input_m.json_m << std::endl;
  myfile.close();
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Finished writing JSON to file {}", a_file_name);
}

std::string InputParser::WriteToString(void) const {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Writing JSON to string");
  return parsed_input_m.json_m.dump(4);
}

std::vector<std::uint8_t> InputParser::ToBSON(void) const {
  return nlohmann::json::to_bson(parsed_input_m.json_m);
}

void InputParser::SetFromBSON(const std::vector<std::uint8_t>& a_bson) {
  parsed_input_m.json_m = nlohmann::json::from_bson(a_bson);
}

const nlohmann::json& InputParser::operator[](const std::string& a_name) const {
  // FIXME : Make an exception.
  return parsed_input_m[a_name];
}

bool InputParser::AllOptionsSet(void) const {
  // FIXME: Add test for this.
  for (const auto& elem : parsed_input_m.json_m) {
    if (elem.empty()) {
      return false;
    }
  }
  return true;
}

bool InputParser::OptionSet(const std::string& a_name) const {
  return parsed_input_m.Contains(a_name);
}

void InputParser::ClearOptions(void) {
  option_description_m.json_m.clear();
  option_required_status_m.clear();
  dependencies_m.clear();
}

void InputParser::PrintOptions(void) const {
  std::cout << option_description_m.json_m.dump(4);
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Finished writing help to std::cout");
}

}  // namespace chyps
