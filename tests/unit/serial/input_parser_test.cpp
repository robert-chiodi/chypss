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

TEST(InputParserCL, BoolDefault) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  inputs.emplace_back("input_name");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data() + 2;
  int argc = input_vec.size() - 2;
  InputParser parser;
  parser.AddOption("test", "A test value for testing", false);
  parser.ParseCL(argc, argv);
  ASSERT_FALSE(parser["test"].get<bool>());
  DeleteCommandLineInput(input_vec);
}

TEST(InputParserCL, IntDefault) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  inputs.emplace_back("input_name");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data() + 2;
  int argc = input_vec.size() - 2;
  InputParser parser;
  static constexpr int correct = {43};
  parser.AddOption("test", "A test value for testing", correct);
  parser.ParseCL(argc, argv);
  int result = parser["test"].get<int>();
  ASSERT_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParserCL, DoubleDefault) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  inputs.emplace_back("input_name");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data() + 2;
  int argc = input_vec.size() - 2;
  InputParser parser;
  static constexpr double correct = {13.7892};
  parser.AddOption("test", "A test value for testing", correct);
  parser.ParseCL(argc, argv);
  double result = parser["test"].get<double>();
  ASSERT_DOUBLE_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParserCL, StringDefault) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  inputs.emplace_back("input_name");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data() + 2;
  int argc = input_vec.size() - 2;
  InputParser parser;
  std::string correct = "A_Real_Test";
  parser.AddOption("test", "A test value for testing", correct);
  parser.ParseCL(argc, argv);
  std::string result = parser["test"].get<std::string>();
  ASSERT_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParserCL, BoolDefaultOverwrite) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  inputs.emplace_back("input_name");
  inputs.emplace_back("-test");
  inputs.emplace_back("true");
  inputs.emplace_back("-test2");
  inputs.emplace_back("false");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data() + 2;
  int argc = input_vec.size() - 2;
  InputParser parser;
  parser.AddOption("test", "A test value for testing", false);
  parser.AddOption("test2", "A test value for testing");
  parser.ParseCL(argc, argv);
  ASSERT_TRUE(parser["test"].get<bool>());
  ASSERT_FALSE(parser["test2"].get<bool>());
  DeleteCommandLineInput(input_vec);
}

TEST(InputParserCL, IntDefaultOverwrite) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  inputs.emplace_back("input_name");
  inputs.emplace_back("-test");
  inputs.emplace_back("43");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data() + 2;
  int argc = input_vec.size() - 2;
  InputParser parser;
  static constexpr int incorrect = {22};
  static constexpr int correct = {43};
  parser.AddOption("test", "A test value for testing", incorrect);
  parser.ParseCL(argc, argv);
  int result = parser["test"].get<int>();
  ASSERT_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParserCL, DoubleDefaultOverwrite) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  inputs.emplace_back("input_name");
  inputs.emplace_back("-test");
  inputs.emplace_back("13.7892");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data() + 2;
  int argc = input_vec.size() - 2;
  InputParser parser;
  static constexpr double incorrect = {962.00};
  static constexpr double correct = {13.7892};
  parser.AddOption("test", "A test value for testing", incorrect);
  parser.ParseCL(argc, argv);
  double result = parser["test"].get<double>();
  ASSERT_DOUBLE_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParserCL, StringDefaultOverwrite) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  inputs.emplace_back("input_name");
  inputs.emplace_back("-test");
  inputs.emplace_back("A_Real_Test");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data() + 2;
  int argc = input_vec.size() - 2;
  InputParser parser;
  std::string incorrect = "Faker";
  std::string correct = "A_Real_Test";
  parser.AddOption("test", "A test value for testing", incorrect);
  parser.ParseCL(argc, argv);
  std::string result = parser["test"].get<std::string>();
  ASSERT_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParserCL, BoolRequired) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  inputs.emplace_back("input_name");
  inputs.emplace_back("-test");
  inputs.emplace_back("true");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data() + 2;
  int argc = input_vec.size() - 2;
  InputParser parser;
  parser.AddOption("test", "A test value for testing");
  parser.ParseCL(argc, argv);
  ASSERT_TRUE(parser["test"].get<bool>());
  DeleteCommandLineInput(input_vec);
}

TEST(InputParserCL, IntRequired) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  inputs.emplace_back("input_name");
  inputs.emplace_back("-test");
  inputs.emplace_back("43");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data() + 2;
  int argc = input_vec.size() - 2;
  InputParser parser;
  static constexpr int correct = {43};
  parser.AddOption("test", "A test value for testing");
  parser.ParseCL(argc, argv);
  int result = parser["test"].get<int>();
  ASSERT_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParserCL, DoubleRequired) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  inputs.emplace_back("input_name");
  inputs.emplace_back("-test");
  inputs.emplace_back("13.7892");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data() + 2;
  int argc = input_vec.size() - 2;
  InputParser parser;
  static constexpr double correct = {13.7892};
  parser.AddOption("test", "A test value for testing");
  parser.ParseCL(argc, argv);
  double result = parser["test"].get<double>();
  ASSERT_DOUBLE_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParserCL, StringRequired) {
  std::vector<std::string> inputs;
  inputs.emplace_back("input_name");
  inputs.emplace_back("exe_name");
  inputs.emplace_back("-test");
  inputs.emplace_back("A_Real_Test");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data() + 2;
  int argc = input_vec.size() - 2;
  InputParser parser;
  std::string correct = "A_Real_Test";
  parser.AddOption("test", "A test value for testing");
  parser.ParseCL(argc, argv);
  std::string result = parser["test"].get<std::string>();
  ASSERT_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParserCL, MultiOptions) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  inputs.emplace_back("input_name");
  inputs.emplace_back("-string_test");
  inputs.emplace_back("A_Real_Test");
  inputs.emplace_back("-false_bool");
  inputs.emplace_back("false");
  inputs.emplace_back("-double_test");
  inputs.emplace_back("-482.968");
  inputs.emplace_back("-int_test");
  inputs.emplace_back("1042");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data() + 2;
  int argc = input_vec.size() - 2;

  InputParser parser;
  // Bool, Default, overwritten
  parser.AddOption("false_bool", "A test value for testing", true);
  // Int, Required
  parser.AddOption("int_test", "--int-test", "Test for int");
  // Double, Required
  parser.AddOption("double_test", "Test for double");
  // String, Default, overwritten
  parser.AddOption("string_test", "Test for string", std::string("Faker"));

  // Default bool, not overwritten
  parser.AddOption("on_bool_def", "A test value for testing (take default)",
                   false);

  // Default string, not overwritten
  parser.AddOption("string test default", "Test for string (take default)",
                   std::string("Faker"));

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

  parser.AddOption("test_double", "A test value for testing");

  parser.ParseFromFile("tests/unit/serial/data/test_input.json");

  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  inputs.emplace_back("input_name");
  inputs.emplace_back("-test_double");
  inputs.emplace_back("13.7892");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data() + 2;
  int argc = input_vec.size() - 2;
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

  parser.AddOption("test_double", "A test value for testing");

  parser.ParseFromFile("tests/unit/serial/data/test_input.json");

  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  inputs.emplace_back("input_name");
  inputs.emplace_back("-test_double");
  inputs.emplace_back("13.7892");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data() + 2;
  int argc = input_vec.size() - 2;
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
