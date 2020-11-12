// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chyps/mpi_parallel.hpp"
#include "chyps/simulation_initializer.hpp"

int main(int argc, char** argv) {
  chyps::MPIParallel mpi_session(&argc, &argv);
  chyps::SimulationInitializer(argc, argv, mpi_session);
  return 0;
}
