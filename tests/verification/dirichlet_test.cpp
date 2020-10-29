// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>
#include "tests/unit/parallel/mpi_session.hpp"

#include "tests/helper/command_line_input.hpp"
#include "tests/helper/compute_error.hpp"

#include "chyps/io.hpp"
#include "chyps/logger.hpp"
#include "chyps/simulation.hpp"

namespace {

using namespace chyps;

TEST(Dirichlet, Left) {
  std::vector<std::string> input_string;
  input_string.push_back("Executable_name");
  input_string.push_back(
      "tests/verification/data/dirichlet_left_test_input.json");
  input_string.push_back("-tf");
  input_string.push_back("0.1");
  input_string.push_back("-od");
  input_string.push_back("tests/verification/data/Dirichlet");
  auto input_char = FakeCommandLineInput(input_string);
  int argc = static_cast<int>(input_char.size());
  char** argv = input_char.data();
  chyps::main(argc, argv, *mpi_session, SpdlogLevel::OFF);

  IO run_file(*mpi_session, "OneRun");
  run_file.SetRead("tests/verification/data/Dirichlet");
  std::vector<double> temperature_field;
  run_file.GetImmediateBlock("HeatSolver/Temperature", temperature_field);
  const double global_max = GlobalMaxValue(temperature_field, *mpi_session);
  EXPECT_DOUBLE_EQ(global_max, 10.0);  // 10.0 is Dirichlet value from input
}

TEST(Dirichlet, Right) {
  std::vector<std::string> input_string;
  input_string.push_back("Executable_name");
  input_string.push_back(
      "tests/verification/data/dirichlet_right_test_input.json");
  input_string.push_back("-tf");
  input_string.push_back("0.1");
  input_string.push_back("-od");
  input_string.push_back("tests/verification/data/Dirichlet");
  auto input_char = FakeCommandLineInput(input_string);
  int argc = static_cast<int>(input_char.size());
  char** argv = input_char.data();
  chyps::main(argc, argv, *mpi_session, SpdlogLevel::OFF);

  IO run_file(*mpi_session, "OneRun");
  run_file.SetRead("tests/verification/data/Dirichlet");
  std::vector<double> temperature_field;
  run_file.GetImmediateBlock("HeatSolver/Temperature", temperature_field);
  const double global_max = GlobalMaxValue(temperature_field, *mpi_session);
  EXPECT_DOUBLE_EQ(global_max, 5.0);  // 5.0 is Dirichlet value from input
}

TEST(Dirichlet, Bottom) {
  std::vector<std::string> input_string;
  input_string.push_back("Executable_name");
  input_string.push_back(
      "tests/verification/data/dirichlet_bottom_test_input.json");
  input_string.push_back("-tf");
  input_string.push_back("0.1");
  input_string.push_back("-od");
  input_string.push_back("tests/verification/data/Dirichlet");
  auto input_char = FakeCommandLineInput(input_string);
  int argc = static_cast<int>(input_char.size());
  char** argv = input_char.data();
  chyps::main(argc, argv, *mpi_session, SpdlogLevel::OFF);

  IO run_file(*mpi_session, "OneRun");
  run_file.SetRead("tests/verification/data/Dirichlet");
  std::vector<double> temperature_field;
  run_file.GetImmediateBlock("HeatSolver/Temperature", temperature_field);
  const double global_max = GlobalMaxValue(temperature_field, *mpi_session);
  EXPECT_DOUBLE_EQ(global_max, 3.0);  // 3.0 is Dirichlet value from input
}

TEST(Dirichlet, Top) {
  std::vector<std::string> input_string;
  input_string.push_back("Executable_name");
  input_string.push_back(
      "tests/verification/data/dirichlet_top_test_input.json");
  input_string.push_back("-tf");
  input_string.push_back("0.1");
  input_string.push_back("-od");
  input_string.push_back("tests/verification/data/Dirichlet");
  auto input_char = FakeCommandLineInput(input_string);
  int argc = static_cast<int>(input_char.size());
  char** argv = input_char.data();
  chyps::main(argc, argv, *mpi_session, SpdlogLevel::OFF);

  IO run_file(*mpi_session, "OneRun");
  run_file.SetRead("tests/verification/data/Dirichlet");
  std::vector<double> temperature_field;
  run_file.GetImmediateBlock("HeatSolver/Temperature", temperature_field);
  const double global_max = GlobalMaxValue(temperature_field, *mpi_session);
  EXPECT_DOUBLE_EQ(global_max, 8.0);  // 8.0 is Dirichlet value from input
}
}  // namespace
