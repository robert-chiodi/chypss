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
#include <gtest-mpi-listener/gtest-mpi-listener.hpp>

#define CHYPS_LOGGER_OFF
#include "chyps/logger.hpp"

using namespace chyps;

// Custom main function needed to turn off Spdlog.
// int main(int argc, char** argv) {
//   ::testing::InitGoogleTest(&argc, argv);
//   MPI_Init(&argc, &argv);
//   StartLogger(0, 1, SpdlogLevel::OFF);
//   auto result = RUN_ALL_TESTS();
//   MPI_Finalize();
//   return result;
// }

int main(int argc, char** argv) {
  // Filter out Google Test arguments
  ::testing::InitGoogleTest(&argc, argv);

  // Initialize MPI
  MPI_Init(&argc, &argv);

  // Turn off logger
  StartLogger(0, 1, SpdlogLevel::OFF);

  // Add object that will finalize MPI on exit; Google Test owns this pointer
  ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);

  // Get the event listener list.
  ::testing::TestEventListeners& listeners =
      ::testing::UnitTest::GetInstance()->listeners();

  // Remove default listener: the default printer and the default XML printer
  ::testing::TestEventListener* l =
      listeners.Release(listeners.default_result_printer());

  // Adds MPI listener; Google Test owns this pointer
  listeners.Append(new GTestMPIListener::MPIWrapperPrinter(l, MPI_COMM_WORLD));
  // Run tests, then clean up and exit. RUN_ALL_TESTS() returns 0 if all tests
  // pass and 1 if some test fails.
  int result = RUN_ALL_TESTS();

  return result;  // Run tests, then clean up and exit
}
