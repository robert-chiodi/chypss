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

#include <array>
#include <fstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "tests/helper/command_line_input.hpp"

namespace {

using namespace chyps;

TEST(CommonType, DefaultInvalid) {
  CommonType ctype;
  ASSERT_TRUE(ctype.GetType() == InputType::INVALID);
}
TEST(CommonType, BoolConstruct) {
  CommonType ctype(false);
  EXPECT_FALSE(ctype);
  ctype = true;
  EXPECT_TRUE(ctype);
}
TEST(CommonType, IntConstruct) {
  CommonType ctype(43);
  EXPECT_EQ(static_cast<int>(ctype), 43);
  ctype = -10;
  EXPECT_EQ(static_cast<int>(ctype), -10);
}
TEST(CommonType, DoubleConstruct) {
  CommonType ctype(-38.4962);
  EXPECT_DOUBLE_EQ(ctype, -38.4962);
  ctype = 43.9682;
  EXPECT_DOUBLE_EQ(ctype, 43.9682);
}
TEST(CommonType, ConstCharConstruct) {
  CommonType ctype("BestOfTimes");
  EXPECT_TRUE(StringEqual(ctype, "BestOfTimes"));
  ctype = "WorstOfTimes";
  EXPECT_TRUE(StringEqual(ctype, "WorstOfTimes"));
}
TEST(CommonType, StringConstruct) {
  CommonType ctype(std::string("BestOfTimes"));
  EXPECT_TRUE(StringEqual(ctype, "BestOfTimes"));
  ctype = std::string("WorstOfTimes");
  EXPECT_TRUE(StringEqual(ctype, "WorstOfTimes"));
}
TEST(CommonType, ChangeAfterConstruct) {
  CommonType ctype;
  EXPECT_TRUE(ctype.GetType() == InputType::INVALID);
  ctype = 12;
  EXPECT_EQ(static_cast<int>(ctype), 12);
  ctype = 108.124694;
  EXPECT_DOUBLE_EQ(ctype, 108.124694);
  ctype = std::string("StillKickin");
  EXPECT_TRUE(StringEqual(ctype, "StillKickin"));
}

TEST(InputParserCL, BoolDefault) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data();
  int argc = input_vec.size();
  InputParser parser;
  parser.AddOption("test", "-t", "--test", "A test value for testing", false,
                   OptionType::COMMAND_LINE);
  parser.ParseCL(argc, argv);
  ASSERT_FALSE(parser["test"].get<bool>());
  DeleteCommandLineInput(input_vec);
}

TEST(InputParserCL, IntDefault) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data();
  int argc = input_vec.size();
  InputParser parser;
  static constexpr int correct = {43};
  parser.AddOption("test", "-t", "--test", "A test value for testing", correct,
                   OptionType::COMMAND_LINE);
  parser.ParseCL(argc, argv);
  int result = parser["test"].get<int>();
  ASSERT_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParserCL, DoubleDefault) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data();
  int argc = input_vec.size();
  InputParser parser;
  static constexpr double correct = {13.7892};
  parser.AddOption("test", "-t", "--test", "A test value for testing", correct,
                   OptionType::COMMAND_LINE);
  parser.ParseCL(argc, argv);
  double result = parser["test"].get<double>();
  ASSERT_DOUBLE_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParserCL, StringDefault) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data();
  int argc = input_vec.size();
  InputParser parser;
  std::string correct = "A_Real_Test";
  parser.AddOption("test", "-t", "--test", "A test value for testing", correct,
                   OptionType::COMMAND_LINE);
  parser.ParseCL(argc, argv);
  std::string result = parser["test"].get<std::string>();
  ASSERT_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParserCL, BoolDefaultOverwrite) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  inputs.emplace_back("-t");
  inputs.emplace_back("--no-second-test");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data();
  int argc = input_vec.size();
  InputParser parser;
  parser.AddOption("test", "-t", "--test", "A test value for testing", false,
                   OptionType::COMMAND_LINE);
  parser.AddOption("test2", "-st", "--second-test", "A test value for testing",
                   true, OptionType::COMMAND_LINE);
  parser.ParseCL(argc, argv);
  ASSERT_TRUE(parser["test"].get<bool>());
  ASSERT_FALSE(parser["test2"].get<bool>());
  DeleteCommandLineInput(input_vec);
}

TEST(InputParserCL, IntDefaultOverwrite) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  inputs.emplace_back("-t");
  inputs.emplace_back("43");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data();
  int argc = input_vec.size();
  InputParser parser;
  static constexpr int incorrect = {22};
  static constexpr int correct = {43};
  parser.AddOption("test", "-t", "--test", "A test value for testing",
                   incorrect, OptionType::COMMAND_LINE);
  parser.ParseCL(argc, argv);
  int result = parser["test"].get<int>();
  ASSERT_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParserCL, DoubleDefaultOverwrite) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  inputs.emplace_back("--test");
  inputs.emplace_back("13.7892");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data();
  int argc = input_vec.size();
  InputParser parser;
  static constexpr double incorrect = {962.00};
  static constexpr double correct = {13.7892};
  parser.AddOption("test", "-t", "--test", "A test value for testing",
                   incorrect, OptionType::COMMAND_LINE);
  parser.ParseCL(argc, argv);
  double result = parser["test"].get<double>();
  ASSERT_DOUBLE_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParserCL, StringDefaultOverwrite) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  inputs.emplace_back("--test");
  inputs.emplace_back("A_Real_Test");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data();
  int argc = input_vec.size();
  InputParser parser;
  std::string incorrect = "Faker";
  std::string correct = "A_Real_Test";
  parser.AddOption("test", "-t", "--test", "A test value for testing",
                   incorrect, OptionType::COMMAND_LINE);
  parser.ParseCL(argc, argv);
  std::string result = parser["test"].get<std::string>();
  ASSERT_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParserCL, BoolRequired) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  inputs.emplace_back("-t");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data();
  int argc = input_vec.size();
  InputParser parser;
  parser.AddOption<bool>("test", "-t", "--test", "A test value for testing",
                         OptionType::COMMAND_LINE);
  parser.ParseCL(argc, argv);
  ASSERT_TRUE(parser["test"].get<bool>());
  DeleteCommandLineInput(input_vec);
}

TEST(InputParserCL, IntRequired) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  inputs.emplace_back("-t");
  inputs.emplace_back("43");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data();
  int argc = input_vec.size();
  InputParser parser;
  static constexpr int correct = {43};
  parser.AddOption<int>("test", "-t", "--test", "A test value for testing",
                        OptionType::COMMAND_LINE);
  parser.ParseCL(argc, argv);
  int result = parser["test"].get<int>();
  ASSERT_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParserCL, DoubleRequired) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  inputs.emplace_back("--test");
  inputs.emplace_back("13.7892");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data();
  int argc = input_vec.size();
  InputParser parser;
  static constexpr double correct = {13.7892};
  parser.AddOption<double>("test", "-t", "--test", "A test value for testing",
                           OptionType::COMMAND_LINE);
  parser.ParseCL(argc, argv);
  double result = parser["test"].get<double>();
  ASSERT_DOUBLE_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParserCL, StringRequired) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  inputs.emplace_back("--test");
  inputs.emplace_back("A_Real_Test");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data();
  int argc = input_vec.size();
  InputParser parser;
  std::string correct = "A_Real_Test";
  parser.AddOption<std::string>("test", "-t", "--test",
                                "A test value for testing",
                                OptionType::COMMAND_LINE);
  parser.ParseCL(argc, argv);
  std::string result = parser["test"].get<std::string>();
  ASSERT_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParserCL, MultiOptions) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  inputs.emplace_back("--string_test");
  inputs.emplace_back("A_Real_Test");
  inputs.emplace_back("--no-on-bool");
  inputs.emplace_back("-d");
  inputs.emplace_back("-482.968");
  inputs.emplace_back("--int-test");
  inputs.emplace_back("1042");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data();
  int argc = input_vec.size();

  InputParser parser;
  // Bool, Default, overwritten
  parser.AddOption("false_bool", "-b", "--on-bool", "A test value for testing",
                   true, OptionType::COMMAND_LINE);
  // Int, Required
  parser.AddOption<int>("int_test", "-i", "--int-test", "Test for int",
                        OptionType::COMMAND_LINE);
  // Double, Required
  parser.AddOption<double>("double_test", "-d", "--double-test",
                           "Test for double", OptionType::COMMAND_LINE);
  // String, Default, overwritten
  parser.AddOption("string_test", "-s", "--string_test", "Test for string",
                   std::string("Faker"), OptionType::COMMAND_LINE);

  // Default bool, not overwritten
  parser.AddOption("on_bool_def", "-v", "--on-bool-def",
                   "A test value for testing (take default)", false,
                   OptionType::COMMAND_LINE);

  // Default string, not overwritten
  parser.AddOption("string test default", "-sd", "--string_test-def",
                   "Test for string (take default)", std::string("Faker"),
                   OptionType::COMMAND_LINE);

  parser.ParseCL(argc, argv);

  EXPECT_FALSE(parser["false_bool"].get<bool>());
  EXPECT_EQ(parser["int_test"].get<int>(), 1042);
  EXPECT_DOUBLE_EQ(parser["double_test"].get<double>(), -482.968);
  EXPECT_TRUE(
      StringEqual(parser["string_test"].get<std::string>(), "A_Real_Test"));
  EXPECT_FALSE(parser["on_bool_def"].get<bool>());
  EXPECT_TRUE(
      StringEqual(parser["string test default"].get<std::string>(), "Faker"));
  DeleteCommandLineInput(input_vec);
}

TEST(InputParserFILE, ParseFile) {
  InputParser parser;
  parser.AddOption("test_bool", "Testing reading of bool.");
  parser.AddOption("test_int", "Testing reading of int.");
  parser.AddOption("test_double", "Testing reading of double.");
  parser.AddOption("test_string1", "Test string without spaces");
  parser.AddOption("test_string2", "Test string with spaces");
  parser.AddOption("test_double_vec", "Vector of doubles.");
  parser.AddOption("test_int_vec", "Vector of ints.");

  parser.ParseFromFile("tests/unit/serial/data/test_input.json");
  EXPECT_FALSE(parser["test_bool"].get<bool>());
  EXPECT_EQ(parser["test_int"].get<int>(), 42);
  EXPECT_DOUBLE_EQ(parser["test_double"].get<double>(), 68.492);
  EXPECT_TRUE(StringEqual(parser["test_string1"], "TheBestString"));
  EXPECT_TRUE(StringEqual(parser["test_string2"], "TheBestString With Spaces"));
  const auto double_vector =
      parser["test_double_vec"].get<std::vector<double>>();
  const std::vector<double> double_correct{{-82.3, 0.0, 15.67, 22.8}};
  for (std::size_t n = 0; n < double_correct.size(); ++n) {
    EXPECT_DOUBLE_EQ(double_vector[n], double_correct[n]);
  }

  const auto int_vector = parser["test_int_vec"].get<std::vector<int>>();
  const std::vector<int> int_correct{{1, 2, 3, -6, 4}};
  for (std::size_t n = 0; n < int_correct.size(); ++n) {
    EXPECT_EQ(int_vector[n], int_correct[n]);
  }
}

TEST(InputParser, CLOverwrite) {
  InputParser parser;
  parser.AddOption("test_bool", "Testing reading of bool.");
  parser.AddOption("test_int", "Testing reading of int.");
  parser.AddOption("test_string1", "Test string without spaces");
  parser.AddOption("test_string2", "Test string with spaces");
  parser.AddOption("test_double_vec", "Vector of doubles.");
  parser.AddOption("test_int_vec", "Vector of ints.");

  parser.AddOption<double>("test_double", "-t", "--test",
                           "A test value for testing",
                           OptionType::COMMAND_LINE);

  parser.ParseFromFile("tests/unit/serial/data/test_input.json");

  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  inputs.emplace_back("--test");
  inputs.emplace_back("13.7892");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data();
  int argc = input_vec.size();
  parser.ParseCL(argc, argv);
  EXPECT_DOUBLE_EQ(parser["test_double"].get<double>(), 13.7892);
}

TEST(InputParser, OptionsToFile) {
  InputParser parser;
  parser.AddOption("test_bool", "Testing reading of bool.");
  parser.AddOption("test_int", "Testing reading of int.");
  parser.AddOption("test_string1", "Test string without spaces");
  parser.AddOption("test_string2", "Test string with spaces");
  parser.AddOption("test_double_vec", "Vector of doubles.");
  parser.AddOption("test_int_vec", "Vector of ints.");

  parser.AddOption<double>("test_double", "-t", "--test",
                           "A test value for testing",
                           OptionType::COMMAND_LINE);

  parser.ParseFromFile("tests/unit/serial/data/test_input.json");

  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  inputs.emplace_back("--test");
  inputs.emplace_back("13.7892");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data();
  int argc = input_vec.size();
  parser.ParseCL(argc, argv);

  std::string write_file_name = "tests/unit/serial/parser_writing_test.json";
  parser.WriteToFile(write_file_name);

  InputParser read_parser;
  read_parser.AddOption("test_bool", "Testing reading of bool.");
  read_parser.AddOption("test_int", "Testing reading of int.");
  read_parser.AddOption("test_double", "Testing reading of double.");
  read_parser.AddOption("test_string1", "Test string without spaces");
  read_parser.AddOption("test_string2", "Test string with spaces");
  read_parser.AddOption("test_double_vec", "Vector of doubles.");
  read_parser.AddOption("test_int_vec", "Vector of ints.");
  read_parser.ParseFromFile(write_file_name);

  EXPECT_FALSE(parser["test_bool"].get<bool>());
  EXPECT_EQ(parser["test_int"].get<int>(), 42);
  EXPECT_DOUBLE_EQ(parser["test_double"].get<double>(), 13.7892);
  EXPECT_TRUE(StringEqual(parser["test_string1"], "TheBestString"));
  EXPECT_TRUE(StringEqual(parser["test_string2"], "TheBestString With Spaces"));
  const auto double_vector =
      parser["test_double_vec"].get<std::vector<double>>();
  const std::vector<double> double_correct{{-82.3, 0.0, 15.67, 22.8}};
  for (std::size_t n = 0; n < double_correct.size(); ++n) {
    EXPECT_DOUBLE_EQ(double_vector[n], double_correct[n]);
  }

  const auto int_vector = parser["test_int_vec"].get<std::vector<int>>();
  const std::vector<int> int_correct{{1, 2, 3, -6, 4}};
  for (std::size_t n = 0; n < int_correct.size(); ++n) {
    EXPECT_EQ(int_vector[n], int_correct[n]);
  }
  remove(write_file_name.c_str());
}
TEST(InputParser, RoundTripBSON) {
  InputParser parser;
  parser.AddOption("test_bool", "Testing reading of bool.");
  parser.AddOption("test_int", "Testing reading of int.");
  parser.AddOption("test_double", "Testing reading of double.");
  parser.AddOption("test_string1", "Test string without spaces");
  parser.AddOption("test_string2", "Test string with spaces");
  parser.AddOption("test_double_vec", "Vector of doubles.");
  parser.AddOption("test_int_vec", "Vector of ints.");

  parser.ParseFromFile("tests/unit/serial/data/test_input.json");

  auto v_bson = parser.ToBSON();

  InputParser recreate_from_bson;
  recreate_from_bson.SetFromBSON(v_bson);

  EXPECT_EQ(parser["test_bool"].get<bool>(),
            recreate_from_bson["test_bool"].get<bool>());
  EXPECT_EQ(parser["test_int"].get<int>(),
            recreate_from_bson["test_int"].get<int>());
  EXPECT_DOUBLE_EQ(parser["test_double"].get<double>(),
                   recreate_from_bson["test_double"].get<double>());
  EXPECT_TRUE(
      StringEqual(parser["test_string1"],
                  recreate_from_bson["test_string1"].get<std::string>()));
  EXPECT_TRUE(
      StringEqual(parser["test_string2"],
                  recreate_from_bson["test_string2"].get<std::string>()));
  const auto double_vector =
      parser["test_double_vec"].get<std::vector<double>>();
  const auto double_correct =
      recreate_from_bson["test_double_vec"].get<std::vector<double>>();
  for (std::size_t n = 0; n < double_correct.size(); ++n) {
    EXPECT_DOUBLE_EQ(double_vector[n], double_correct[n]);
  }

  const auto int_vector = parser["test_int_vec"].get<std::vector<int>>();
  const auto int_correct =
      recreate_from_bson["test_int_vec"].get<std::vector<int>>();

  for (std::size_t n = 0; n < int_correct.size(); ++n) {
    EXPECT_EQ(int_vector[n], int_correct[n]);
  }
}

}  // namespace
