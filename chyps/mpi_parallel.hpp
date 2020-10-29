// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// \file mpi_parallel.hpp
/// \brief RAII wrapper for MPI session.

#ifndef CHYPS_MPI_PARALLEL_HPP_
#define CHYPS_MPI_PARALLEL_HPP_

#include <limits.h>
#include <stdint.h>

#include <mpi.h>

// Taken from StackOverflow
// https://stackoverflow.com/questions/40807833/sending-size-t-type-data-with-mpi
#if SIZE_MAX == UCHAR_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
#error "Unknown match in MPI for size_t"
#endif

namespace chyps {

class MPIParallel {
 public:
  /// \brief Initialize MPI_COMM_WORLD without any command line arguments.
  MPIParallel(void);

  /// \brief Initialize MPI_COMM_WORLD with command line arguments.
  MPIParallel(int* a_argc, char*** a_argv);

  int MyRank(void) const;
  int NumberOfRanks(void) const;
  bool IAmRoot(void) const;
  bool IAmNotRoot(void) const;
  const MPI_Comm& GetComm(void) const;

  ~MPIParallel(void);

 private:
  MPI_Comm mpi_comm_m;
  int my_rank_m;
  int number_of_ranks_m;
};

}  // namespace chyps

#endif  // CHYPS_MPI_PARALLEL_HPP_
