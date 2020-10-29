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

namespace chyps {

/// \enum InputType input_parser.hpp chyps/input_parser.hpp
/// \brief Enum storing the type of input to set flag for CommonType union.
enum class InputType { INVALID = -1, BOOL = 0, INT, DOUBLE, STRING };

/// \function TypeToInputType input_parser.hpp chyps/input_parser.hpp
/// \brief Given a value type, responds the correct InputType enum value.
///
/// Only ValueTypes valid are bool, int, double, and std::string, which
/// correspond to the types in InputType.
template <class ValueType>
InputType TypeToInputType(void);

/// \enum OptionType input_parser.hpp chyps/input_parser.hpp
/// \brief Enum storing the type of option added to the parser.
///
/// Users should only use the options ANY, COMMAND_LINE, and INPUT_FILE.
/// Whether or not they are required are determined via supply of a
/// default value in InputParser::AddOption. Including a default value
/// makes the Option NOT required. Lack of default value makes the
/// option required, and the OptionType will be updated to reflect that
/// in the InputParser class object.
///
///  Required options are to have values >= 10, and should be aligned
/// with their non-required equivalent so that ``required_type =
/// non_required_type+10``.
enum class OptionType {
  ANY = 0,
  COMMAND_LINE,
  INPUT_FILE,
  ANY_REQUIRED = 10,
  COMMAND_LINE_REQUIRED,
  INPUT_FILE_REQUIRED
};

/// \function MakeOptionRequired
/// \brief Takes the OptionType stored in a_option and returns the _REQUIRED
/// equivalent.
OptionType MakeOptionRequired(const OptionType a_option);

/// \class CommonType input_parser.hpp chyps/input_parser.hpp
/// \brief A base type capable of storing a boolean, integer, double, or
/// std::string. Intended for use in InputParser only.
class CommonType {
  union CT {
    bool the_bool;
    int the_int;
    double the_double;
    std::string the_string;

    /// \brief Default constructor. Marked as InputType::INVALID type when used
    /// in CommonType.
    CT(void) : the_bool(false) {}

    /// \brief Construct to bool value. CommonType marks as InputType::BOOL.
    CT(const bool a_value) : the_bool(a_value) {}

    /// \brief Construct to int value. CommonType marks as InputType::INT.
    CT(const int a_value) : the_int(a_value) {}

    /// \brief Construct to double value. CommonType marks as InputType::DOUBLE.
    CT(const double a_value) : the_double(a_value) {}

    /// \brief Construct to std::string value. CommonType marks as
    /// InputType::STRING.
    CT(const std::string a_value) : the_string(a_value) {}

    ~CT(void) {}
  };

 public:
  /// \brief Default constructor. Sets type to InputType::INVALID.
  CommonType(void);

  /// \brief Constructs with value a_bool, sets type to InputType::BOOL.
  explicit CommonType(const bool a_bool);

  /// \brief Constructs with value a_int, sets type to InputType::INT.
  explicit CommonType(const int a_int);

  /// \brief Constructs with value a_double, sets type to InputType::DOUBLE.
  explicit CommonType(const double a_double);

  /// \brief Constructs a std::string from  value a_string, sets type to
  /// InputType::STRING.
  explicit CommonType(const char* a_string);

  /// \brief Constructs with value a_string, sets type to InputType::STRING.
  explicit CommonType(const std::string& a_string);

  /// \brief Copy constructor for CommonType.
  CommonType(const CommonType& a_other);

  /// \brief Move constructor for CommonType.
  CommonType(CommonType&& a_other);

  /// \brief Copy assignment for CommonType.
  CommonType& operator=(const CommonType& a_other);

  /// \brief Move assignment for CommonType.
  CommonType& operator=(CommonType&& a_other);

  /// \brief Copy assignment. Sets to InputType::BOOL with value of a_value.
  CommonType& operator=(const bool a_value);

  /// \brief Copy assignment. Sets to InputType::INT with value of a_value.
  CommonType& operator=(const int a_value);

  /// \brief Copy assignment. Sets to InputType::DOUBLE with value of a_value.
  CommonType& operator=(const double a_value);

  /// \brief Copy assignment. Sets to InputType::STRING with value of a_value.
  CommonType& operator=(const char* a_value);

  /// \brief Copy assignment. Sets to InputType::STRING with value of a_value.
  CommonType& operator=(const std::string& a_value);

  /// \brief Conversion operator allowing conversion of CommonType to bool.
  operator bool(void) const;

  /// \brief Conversion operator allowing conversion of CommonType to int.
  operator int(void) const;

  /// \brief Conversion operator allowing conversion of CommonType to double.
  operator double(void) const;

  /// \brief Conversion operator allowing conversion of CommonType to
  /// std::string.
  operator std::string(void) const;

  /// \brief Set type to treat this as. See InputType enum for options.
  void SetType(const InputType a_type);

  /// \brief Returns the type CommonType is currently being used as.
  InputType GetType(void) const;

  /// \brief Returns a pointer to CommonType interpreted as ``Type``.
  ///
  /// This method has template specializations for each member in
  /// the InputType enum, with the exception of InputType::INVALID.
  /// Calling this method will change CommonType to be Type. This means
  /// if GetPointer is called when type_id_m is not for Type, the current
  /// object will be invalid according to the returned pointer.
  template <class Type>
  Type GetPointer(void);

  /// \brief Custom destructor to properly handle and delete std::string in CT.
  ~CommonType(void);

 private:
  InputType type_id_m;
  CT value_m;
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

  /// \brief Write parsed input to file at a_file_name.
  void WriteToFile(const std::string& a_file_name) const;

  /// \brief Write parsed into to string in JSON format.
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

  /// \brief Add an option to be looked for when parsing. If not found,
  /// a_default value will be used. OptionType must NOT be a *_REQUIRED kind.
  /// Type must be a type from InputType enum.
  template <class ValueType>
  void AddOption(const std::string& a_name, const std::string& a_short_flag,
                 const std::string& a_long_flag,
                 const std::string& a_description,
                 const ValueType& a_default_value,
                 const OptionType a_input_type = OptionType::ANY);

  /// \brief Add an option to be looked for when parsing. OptionType must be
  /// a *_REQUIRED kind. Type must be a type from InputType enum.
  template <class ValueType>
  void AddOption(const std::string& a_name, const std::string& a_short_flag,
                 const std::string& a_long_flag,
                 const std::string& a_description,
                 const OptionType a_input_type = OptionType::ANY);

  /// \brief Add an option to be looked for when parsing. If not found,
  /// a_default value will be used. OptionType must NOT be a *_REQUIRED kind.
  /// Since no flag is specified, this must be an InputType::INPUT_FILE option.
  /// These types can be any of those valid in nlohmann::json.
  template <class ValueType>
  void AddOption(const std::string& a_name, const std::string& a_description,
                 const ValueType& a_default_value);

  /// \brief Add an option to be looked for when parsing. OptionType must be
  /// a *_REQUIRED kind. Since no flag is specified, this must be an
  /// OptionType::INPUT_FILE option. These types can be any of those valid
  /// in nlohmann::json. If marked as optional, it does not need to be supplied.
  // FIXME : Add template to check that requested type is going to be valid.
  void AddNoDefaultOption(const std::string& a_name,
                          const std::string& a_description,
                          const bool is_required = true);

  /// \brief Checks that all options added have been specified or have an
  /// available default value.
  bool AllOptionsSet(void) const;

  bool OptionSet(const std::string& a_name) const;

  /// \brief Remove all option storage to save space. Parsed values will still
  /// be available through operator[].
  void ClearOptions(void);

  /// \brief Print out the options and their default value (if provided).
  void PrintOptions(void) const;

 private:
  bool CommandLineOption(const OptionType a_type) const;
  bool InputFileOption(const OptionType a_type) const;
  bool RequiredOption(const OptionType a_type) const;
  bool OptionalOption(const OptionType a_type) const;

  nlohmann::json parsed_input_m;
  std::vector<std::array<std::string, 4>> option_description_m;
  std::unordered_map<std::string, OptionType> option_type_m;
  std::unordered_map<std::string, InputType> type_m;
};

}  // namespace chyps

#include "chyps/input_parser.tpp"

#endif  // CHYPS_INPUT_PARSER_HPP_
