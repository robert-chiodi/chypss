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

TEST(Checkpoint, ReproduceLongRunLowOrder) {
  // Full length run
  std::vector<std::string> input_string;
  input_string.push_back("Executable_name");
  input_string.push_back(
      "tests/verification/data/checkpoint_loworder_test_input.json");
  input_string.push_back("-tf");
  input_string.push_back("0.5");
  input_string.push_back("-od");
  input_string.push_back("tests/verification/data/checkpointLO_OutputFull");
  auto input_char = FakeCommandLineInput(input_string);
  int argc = static_cast<int>(input_char.size());
  char** argv = input_char.data();
  chyps::main(argc, argv, *mpi_session, SpdlogLevel::OFF);

  MPI_Barrier(mpi_session->GetComm());
  DeleteCommandLineInput(input_char);
  input_string[3] = "0.25";
  input_string[5] = "tests/verification/data/checkpointLO_OutputHalf";
  input_char = FakeCommandLineInput(input_string);
  argc = static_cast<int>(input_char.size());
  argv = input_char.data();
  chyps::main(argc, argv, *mpi_session, SpdlogLevel::OFF);

  MPI_Barrier(mpi_session->GetComm());
  DeleteCommandLineInput(input_char);
  input_string[3] = "0.5";
  input_string[5] = "tests/verification/data/checkpointLO_OutputRestartFull";
  input_string.push_back("-id");
  input_string.push_back("tests/verification/data/checkpointLO_OutputHalf");
  input_char = FakeCommandLineInput(input_string);
  argc = static_cast<int>(input_char.size());
  argv = input_char.data();
  chyps::main(argc, argv, *mpi_session, SpdlogLevel::OFF);

  // Now check that the last steps in two files below match
  // "tests/verification/data/checkpointLO_OutputHalf"
  // "tests/verification/data/checkpointLO_OutputFull"
  IO one_run_file(*mpi_session, "OneRun");
  one_run_file.SetRead("tests/verification/data/checkpointLO_OutputFull");
  IO multi_run_file(*mpi_session, "MultiRun");
  multi_run_file.SetRead(
      "tests/verification/data/checkpointLO_OutputRestartFull");

  std::vector<double> full_temperature_data;
  one_run_file.GetImmediateBlock("HeatSolver/Temperature",
                                 full_temperature_data);
  std::vector<double> restart_temperature_data;
  multi_run_file.GetImmediateBlock("HeatSolver/Temperature",
                                   restart_temperature_data);

  const double error = GlobalL2Diff_Normalized(
      full_temperature_data, restart_temperature_data, *mpi_session);

  static constexpr double tolerance = 1.0e-14;
  EXPECT_TRUE(error < tolerance);
}

TEST(Checkpoint, ReproduceLongRunHigherOrder) {
  // Full length run
  std::vector<std::string> input_string;
  input_string.push_back("Executable_name");
  input_string.push_back(
      "tests/verification/data/checkpoint_highorder_test_input.json");
  input_string.push_back("-tf");
  input_string.push_back("0.5");
  input_string.push_back("-od");
  input_string.push_back("tests/verification/data/checkpointLO_OutputFull");
  auto input_char = FakeCommandLineInput(input_string);
  int argc = static_cast<int>(input_char.size());
  char** argv = input_char.data();
  chyps::main(argc, argv, *mpi_session, SpdlogLevel::OFF);

  MPI_Barrier(mpi_session->GetComm());
  DeleteCommandLineInput(input_char);
  input_string[3] = "0.25";
  input_string[5] = "tests/verification/data/checkpointLO_OutputHalf";
  input_char = FakeCommandLineInput(input_string);
  argc = static_cast<int>(input_char.size());
  argv = input_char.data();
  chyps::main(argc, argv, *mpi_session, SpdlogLevel::OFF);

  MPI_Barrier(mpi_session->GetComm());
  DeleteCommandLineInput(input_char);
  input_string[3] = "0.5";
  input_string[5] = "tests/verification/data/checkpointLO_OutputRestartFull";
  input_string.push_back("-id");
  input_string.push_back("tests/verification/data/checkpointLO_OutputHalf");
  input_char = FakeCommandLineInput(input_string);
  argc = static_cast<int>(input_char.size());
  argv = input_char.data();
  chyps::main(argc, argv, *mpi_session, SpdlogLevel::OFF);

  // Now check that the last steps in two files below match
  // "tests/verification/data/checkpointLO_OutputHalf"
  // "tests/verification/data/checkpointLO_OutputFull"
  IO one_run_file(*mpi_session, "OneRun");
  one_run_file.SetRead("tests/verification/data/checkpointLO_OutputFull");
  IO multi_run_file(*mpi_session, "MultiRun");
  multi_run_file.SetRead(
      "tests/verification/data/checkpointLO_OutputRestartFull");

  std::vector<double> full_temperature_data;
  one_run_file.GetImmediateBlock("HeatSolver/Temperature",
                                 full_temperature_data);
  std::vector<double> restart_temperature_data;
  multi_run_file.GetImmediateBlock("HeatSolver/Temperature",
                                   restart_temperature_data);

  const double error = GlobalL2Diff_Normalized(
      full_temperature_data, restart_temperature_data, *mpi_session);

  static constexpr double tolerance = 1.0e-14;
  EXPECT_TRUE(error < tolerance);
}

}  // namespace
