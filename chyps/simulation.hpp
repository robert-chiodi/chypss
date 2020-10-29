// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef CHYPS_SIMULATION_HPP_
#define CHYPS_SIMULATION_HPP_

#include "chyps/logger.hpp"
#include "chyps/mpi_parallel.hpp"

namespace chyps {
int main(int argc, char** argv, MPIParallel& mpi_session,
         SpdlogLevel a_log_level = SpdlogLevel::INFO);
}

#endif  // CHYPS_SIMULATION_HPP_
