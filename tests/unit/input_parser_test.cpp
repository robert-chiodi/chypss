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
#include <string>

#include <gtest/gtest.h>

namespace {

using namespace chyps;

static std::vector<char*> FakeCommandLineInput(
    const std::vector<std::string>& a_input);

static void DeleteCommandLineInput(std::vector<char*>& a_input);

static std::vector<char*> FakeCommandLineInput(
    const std::vector<std::string>& a_input) {
  std::vector<char*> vec;
  for (const auto& s : a_input) {
    char* c_str = new char[s.size() + 1];
    std::copy(s.begin(), s.end(), c_str);
    c_str[s.size()] = '\0';
    vec.push_back(c_str);
  }
  return vec;
}

static void DeleteCommandLineInput(std::vector<char*>& a_input) {
  for (auto& elem : a_input) {
    delete[] elem;
    elem = nullptr;
  }
}

static bool StringEqual(const std::string& a_s1, const std::string& a_s2);
static bool StringEqual(const std::string& a_s1, const std::string& a_s2) {
  if (a_s1.size() != a_s2.size()) {
    return false;
  } else {
    for (std::size_t n = 0; n < a_s1.size(); ++n) {
      if (a_s1[n] != a_s2[n]) {
        return false;
      }
    }
  }
  return true;
}

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

TEST(InputParser_CL, BoolDefault) {
  std::vector<std::string> inputs;
  inputs.emplace_back("exe_name");
  auto input_vec = FakeCommandLineInput(inputs);
  char** argv = input_vec.data();
  int argc = input_vec.size();
  InputParser parser;
  parser.AddOption("test", "-t", "--test", "A test value for testing", false,
                   OptionType::COMMAND_LINE);
  parser.ParseCL(argc, argv);
  ASSERT_FALSE(parser["test"]);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParser_CL, IntDefault) {
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
  int result = parser["test"];
  ASSERT_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParser_CL, DoubleDefault) {
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
  double result = parser["test"];
  ASSERT_DOUBLE_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParser_CL, StringDefault) {
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
  std::string result = parser["test"];
  ASSERT_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParser_CL, BoolDefaultOverwrite) {
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
  ASSERT_TRUE(parser["test"]);
  ASSERT_FALSE(parser["test2"]);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParser_CL, IntDefaultOverwrite) {
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
  int result = parser["test"];
  ASSERT_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParser_CL, DoubleDefaultOverwrite) {
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
  double result = parser["test"];
  ASSERT_DOUBLE_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParser_CL, StringDefaultOverwrite) {
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
  std::string result = parser["test"];
  ASSERT_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParser_CL, BoolRequired) {
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
  ASSERT_TRUE(parser["test"]);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParser_CL, IntRequired) {
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
  int result = parser["test"];
  ASSERT_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParser_CL, DoubleRequired) {
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
  double result = parser["test"];
  ASSERT_DOUBLE_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParser_CL, StringRequired) {
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
  std::string result = parser["test"];
  ASSERT_EQ(result, correct);
  DeleteCommandLineInput(input_vec);
}

TEST(InputParser_CL, MultiOptions) {
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

  EXPECT_FALSE(parser["false_bool"]);
  EXPECT_EQ(static_cast<int>(parser["int_test"]), 1042);
  EXPECT_DOUBLE_EQ(parser["double_test"], -482.968);
  EXPECT_TRUE(StringEqual(parser["string_test"], "A_Real_Test"));
  EXPECT_FALSE(parser["on_bool_def"]);
  EXPECT_TRUE(StringEqual(parser["string test default"], "Faker"));
  DeleteCommandLineInput(input_vec);
}

}  // namespace
