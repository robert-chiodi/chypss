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

namespace chyps {

MPIParallel::MPIParallel(void) : mpi_comm_m(MPI_COMM_WORLD) {
  int zero = 0;
  char** no_arg = nullptr;
  MPI_Init(&zero, &no_arg);
  MPI_Comm_rank(mpi_comm_m, &my_rank_m);
  MPI_Comm_size(mpi_comm_m, &number_of_ranks_m);
}

MPIParallel::MPIParallel(int* a_argc, char*** a_argv)
    : mpi_comm_m(MPI_COMM_WORLD) {
  MPI_Init(a_argc, a_argv);
  MPI_Comm_rank(mpi_comm_m, &my_rank_m);
  MPI_Comm_size(mpi_comm_m, &number_of_ranks_m);
}

int MPIParallel::MyRank(void) const { return my_rank_m; }

int MPIParallel::NumberOfRanks(void) const { return number_of_ranks_m; }

bool MPIParallel::IAmRoot(void) const { return my_rank_m == 0; }

bool MPIParallel::IAmNotRoot(void) const { return !this->IAmRoot(); }

const MPI_Comm& MPIParallel::GetComm(void) const { return mpi_comm_m; }

MPIParallel::~MPIParallel(void) { MPI_Finalize(); }

}  // namespace chyps
