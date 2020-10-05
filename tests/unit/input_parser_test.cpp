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

// Need to implement below
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

}  // namespace
