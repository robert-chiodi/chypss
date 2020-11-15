// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// \file input_parser.hpp
/// \brief File containing the InputParser class and helper classes.

#ifndef CHYPS_INPUT_PARSER_HPP_
#define CHYPS_INPUT_PARSER_HPP_

#include <array>
#include <string>
#include <unordered_map>
#include <vector>

#include <nlohmann_json/json.hpp>

#include "chyps/mpi_parallel.hpp"

namespace chyps {

/// \class DirectoryJSON input_parser.hpp chyps/input_parser.hpp
/// \brief Light wrapper for nlohmann::json that treats
/// the character '/' in keys as indicating nesting
/// in the JSON file.
class DirectoryJSON {
 public:
  DirectoryJSON(void) = default;

  nlohmann::json& operator[](const std::string& a_name);
  const nlohmann::json& operator[](const std::string& a_name) const;
  bool Contains(const std::string& a_name) const;

  /// \brief Finds the directory-formatted name `a_name`. If it does not
  /// exist, returns nullptr, otherwise returns pointer to the object.
  const nlohmann::json* find(const std::string& a_name) const;

  ~DirectoryJSON(void) = default;

  nlohmann::json json_m;

 private:
  std::pair<nlohmann::json*, std::string> LowestObject(
      const std::string a_name);
  std::pair<const nlohmann::json*, std::string> LowestObject(
      const std::string a_name) const;
};

/// \class InputParser input_parser.hpp chyps/input_parser.hpp
/// \brief Handles options available and setting of them through
/// command line, an input file, or direct setting.
class InputParser {
 public:
  /// \brief Initializes empty InputParser object. Options should be added
  /// to this object before calling a (or several) parsing methods.
  InputParser(void) = default;

  /// \brief Look for command-line arguments (supplied via argc and argv)
  /// that match added options. Allows options to after be checked with their
  /// given name.
  void ParseCL(int argc, char** argv);

  /// \brief Parse input from a supplied input file in JSON format.
  bool ParseFromFile(const std::string& a_file_name);

  /// \brief Parse and communicate file. Reads file in on root processor only
  /// and broadcasts to all others.
  bool ParseFromFile(const std::string& a_file_name,
                     const MPIParallel& a_mpi_session);

  /// \brief Write parsed input to file at a_file_name.
  void WriteToFile(const std::string& a_file_name) const;

  /// \brief Write parsed input to string in JSON format.
  std::string WriteToString(void) const;

  /// \brief Convert parsed input to BSON format stored in the returned vector.
  std::vector<std::uint8_t> ToBSON(void) const;

  /// \brief Convert the provided vector of bytes in BSON format to the parsed
  /// JSON input.
  void SetFromBSON(const std::vector<std::uint8_t>& a_bson);

  /// \brief Lookup of option value through given name.
  ///
  /// As suggested by nlohmann::json, should explicitly
  /// case value to appropriate type via .get<Type>()
  /// called on the returned reference.
  const nlohmann::json& operator[](const std::string& a_name) const;

  /// \brief Directly set values for parsed_input_m that can then be
  /// used through operator[].
  template <class Type>
  void DirectSet(const std::string& a_name, const Type& a_value);

  /// \brief Add option with a_name and the given description. The default value
  /// will be used if the option is not provided in the input file or on the
  /// command line.
  ///
  /// In order to indicate nesting inside the input file (JSON), use '/'
  /// in a_name.
  ///
  /// NOTE: Use of '/' is not valid in a_name and will be assumed to represent
  /// nesting.
  template <class ValueType>
  void AddOption(const std::string& a_name, const std::string& a_description,
                 const ValueType& a_default_value);

  /// \brief Add option with a_name and the given description. A value must be
  /// supplied via the input file or command line since no default exists.
  ///
  /// In order to indicate nesting inside the input file (JSON), use '/'
  /// in a_name.
  ///
  /// NOTE: Use of '/' is not valid in a_name and will be assumed to represent
  /// nesting.
  void AddOption(const std::string& a_name, const std::string& a_description);

  /// \brief Returns whether the option is set during parsing, including
  /// whether the option was added with a default value.
  bool OptionSet(const std::string& a_name) const;

  /// \brief Remove all option metadata used when initially setting up available
  /// options. Any parsed options will still be available through the
  /// operator[].
  void ClearOptions(void);

  /// \brief Print out the options and their default value (if provided).
  void PrintOptions(void) const;

 private:
  void RecursiveInsert(const nlohmann::json& a_input,
                       nlohmann::json& a_parsed_nest);
  bool RecursiveOptionSetCheck(const nlohmann::json& a_input,
                               const std::string& a_path_name) const;
  void RecursiveOptionPrint(const nlohmann::json& a_input,
                            const std::string& a_path_name,
                            const int a_nest_level) const;
  static std::string AddTabs(const int a_number_of_tabs, const int a_tab_size);
  static std::string ReTab(const std::string& a_string,
                           const int a_number_of_tabs, const int a_tab_size);
  static std::string BreakupString(const std::string& a_string);

  DirectoryJSON parsed_input_m;
  DirectoryJSON default_values_m;
  DirectoryJSON option_description_m;
};

}  // namespace chyps

#include "chyps/input_parser.tpp"

#endif  // CHYPS_INPUT_PARSER_HPP_
