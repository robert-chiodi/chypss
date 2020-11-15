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

#include <debug_assert/debug_assert.hpp>
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
  DEBUG_ASSERT(lowest_dir->contains(lowest_name), global_assert{},
               DebugLevel::CHEAP{},
               "DirectoryJSON object does not contain \"" + a_name + "\"");
  return (*lowest_dir)[lowest_name];
}

const nlohmann::json* DirectoryJSON::find(const std::string& a_name) const {
  const nlohmann::json* lowest_dir;
  std::string lowest_name;
  std::tie(lowest_dir, lowest_name) = this->LowestObject(a_name);
  return lowest_dir != nullptr && lowest_dir->contains(lowest_name)
             ? &((*lowest_dir)[lowest_name])
             : nullptr;
}

bool DirectoryJSON::Contains(const std::string& a_name) const {
  const nlohmann::json* lowest_dir;
  std::string lowest_name;
  std::tie(lowest_dir, lowest_name) = this->LowestObject(a_name);
  return lowest_dir != nullptr ? lowest_dir->contains(lowest_name) : false;
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
    const std::string& next_object_name = a_name.substr(
        name_start_location, slash_location - name_start_location);
    if (!(lowest_dir->contains(next_object_name))) {
      return std::make_pair(nullptr, "");
    }
    lowest_dir = &((*lowest_dir)[next_object_name]);
    name_start_location = slash_location + 1;
    slash_location = a_name.find('/', name_start_location);
  }
  return std::make_pair(lowest_dir, a_name.substr(name_start_location));
}

void InputParser::ParseCL(int argc, char** argv) {
  DEBUG_ASSERT(argc % 2 == 0, global_assert{}, DebugLevel::ALWAYS{},
               "Even number of command line arguments required due to use of "
               "flag, value paring. Number of arguments supplied:" +
                   std::to_string(argc));

  for (int n = 0; n < argc; n += 2) {
    std::string argument(argv[n]);
    std::string value(argv[n + 1]);

    DEBUG_ASSERT(
        argument[0] == '-', global_assert{}, DebugLevel::ALWAYS{},
        "Flag requires first character of -. Flag passed: " + argument);
    const std::string name = argument.substr(1);

    DEBUG_ASSERT(option_description_m.Contains(name), global_assert{},
                 DebugLevel::ALWAYS{},
                 "Option with key " + name + " is unknown.");

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
      if (value == "true") {
        parsed_input_m[name] = true;
      } else if (value == "false") {
        parsed_input_m[name] = false;
      } else {
        parsed_input_m[name] = value;
      }
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

bool InputParser::ParseFromFile(const std::string& a_file_name,
                                const MPIParallel& a_mpi_session) {
  std::vector<std::uint8_t> v_bson;
  std::size_t size = 0;
  if (a_mpi_session.IAmRoot()) {
    const bool good_read = this->ParseFromFile(a_file_name);
    if (!good_read) {
      DEBUG_ASSERT(
          good_read, global_assert{}, DebugLevel::ALWAYS{},
          "Trouble reading input file: " + a_file_name + '\n' +
              "First argument should be name of input file.\n" +
              "Use the input file name \"ignore\" to use no input file.\n" +
              "Use the input file name \"help\" to export available "
              "options.");
    }
    v_bson = this->ToBSON();
    size = v_bson.size();
  }
  MPI_Bcast(&size, 1, my_MPI_SIZE_T, 0, a_mpi_session.GetComm());
  v_bson.resize(size);  // Will do nothing for rank 0
  MPI_Bcast(v_bson.data(), size, MPI_BYTE, 0, a_mpi_session.GetComm());
  if (a_mpi_session.IAmNotRoot()) {
    this->SetFromBSON(v_bson);
  }
  return true;
}

void InputParser::AddOption(const std::string& a_name,
                            const std::string& a_description) {
  DEBUG_ASSERT(
      !option_description_m.Contains(a_name), global_assert{},
      DebugLevel::CHEAP{},
      "Option \"" + a_name + "\" already exists in InputParser object");
  option_description_m[a_name] = a_description;
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
  const auto given_object = parsed_input_m.find(a_name);
  if (given_object == nullptr) {
    // Check default options
    const auto default_object = default_values_m.find(a_name);
    DEBUG_ASSERT(
        default_object != nullptr, global_assert{}, DebugLevel::CHEAP{},
        "Cannot find \"" + a_name + "\" as a given or default option.");
    return *default_object;
  }
  return *given_object;
}

bool InputParser::OptionSet(const std::string& a_name) const {
  return parsed_input_m.Contains(a_name);
}

void InputParser::ClearOptions(void) { option_description_m.json_m.clear(); }

void InputParser::PrintOptions(void) const {
  std::cout << "\n\n// NOTE: This is not a valid JSON File.\n\n";
  this->RecursiveOptionPrint(option_description_m.json_m, "", 0);
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Finished writing help to std::cout");
}

void InputParser::RecursiveOptionPrint(const nlohmann::json& a_input,
                                       const std::string& a_path_name,
                                       const int a_nest_level) const {
  for (nlohmann::json::const_iterator it = a_input.begin(); it != a_input.end();
       ++it) {
    if (it->is_object()) {
      // Is a nested object object
      std::cout << InputParser::AddTabs(a_nest_level, 4) << '"' << it.key()
                << "\":"
                << "{\n\n";
      this->RecursiveOptionPrint(*it, a_path_name + it.key() + '/',
                                 a_nest_level + 1);
      std::cout << InputParser::AddTabs(a_nest_level, 4) << "} // End " << '"'
                << it.key() << "\"\n\n ";
    } else {
      const std::string full_name = a_path_name + it.key();
      std::cout << InputParser::AddTabs(a_nest_level, 4) << '"' << it.key()
                << "\":"
                << "\n";
      std::cout << InputParser::AddTabs(a_nest_level + 2, 4) << "DESCRIPTION: "
                << InputParser::ReTab(InputParser::BreakupString(
                                          a_input[it.key()].get<std::string>()),
                                      a_nest_level + 2, 4)
                << '\n';

      std::cout << InputParser::AddTabs(a_nest_level + 2, 4)
                << "DEFAULT VALUE: ";
      if (!default_values_m.Contains(full_name)) {
        std::cout << "NONE";
      } else {
        std::cout << default_values_m[full_name].dump();
      }
      std::cout << "\n\n";
    }
  }
}

std::string InputParser::AddTabs(const int a_number_of_tabs,
                                 const int a_tab_size) {
  return std::string(a_tab_size * a_number_of_tabs, ' ');
}

std::string InputParser::ReTab(const std::string& a_string,
                               const int a_number_of_tabs,
                               const int a_tab_size) {
  std::string return_string = a_string;
  std::size_t newline_location = return_string.find('\n');
  while (newline_location != return_string.npos) {
    return_string.replace(
        newline_location, 1,
        '\n' + InputParser::AddTabs(a_number_of_tabs, a_tab_size));
    newline_location = return_string.find(
        '\n', newline_location + a_number_of_tabs * a_tab_size + 1);
  }
  return return_string;
}

std::string InputParser::BreakupString(const std::string& a_string) {
  static constexpr std::size_t line_length = 80;
  if (a_string.size() < line_length) {
    return a_string;
  } else {
    std::string return_string = a_string;
    std::size_t string_end = line_length;
    std::size_t space_location = 0;
    while (return_string.size() > string_end) {
      space_location = return_string.find(" ", space_location + line_length);
      if (space_location == return_string.npos) {
        return return_string;
      }
      return_string[space_location] = '\n';
      string_end += line_length;
    }
    return return_string;
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

}  // namespace chyps
