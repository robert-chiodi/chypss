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

#include <stdio.h>
#include <fstream>
#include <iostream>

#include "chyps/heat_solver.hpp"
#include "chyps/input_parser.hpp"
#include "chyps/io.hpp"
#include "chyps/logger.hpp"
#include "chyps/mesh.hpp"
#include "chyps/mpi_parallel.hpp"
#include "chyps/precice_adapter.hpp"

namespace chyps {
int main(int argc, char** argv, MPIParallel& mpi_session,
         const SpdlogLevel a_log_level = SpdlogLevel::INFO);

struct IterationInfo {
  IterationInfo(const std::size_t a_iteration, const double a_time,
                const double a_dt)
      : iteration(a_iteration), time(a_time), dt(a_dt) {}

  uint64_t iteration;
  double time;
  double dt;
};

struct SimulationRestrictions {
  double max_dt;
};

struct SimulationOutput {
  SimulationOutput(void) : output_screen(nullptr), visualization_steps(-1) {}

  FILE* output_screen;
  int visualization_steps;
};

struct SimulationGoals {
  std::size_t max_iteration;
  double end_time;
};

class Simulation {
 public:
  Simulation(MPIParallel& a_mpi_session,
             const SpdlogLevel a_log_level = SpdlogLevel::INFO);

  void Initialize(int argc, char** argv);

  void RunToEnd(void);

  bool PreciceActive(void) const;

  ~Simulation(void);

 private:
  void GatherOptions(void);
  void ParseOptions(int argc, char** argv);
  void SetupFileIO(void);
  void ActivatePrecice(void);
  void InitializeStepInfo(void);
  void InitializeRestrictions(void);
  void InitializeGoals(void);
  void InitializeOutputs(void);
  void WriteInitialConditions(void);
  void WriteIterationConditions(const IterationInfo& a_iter);
  void WriteProgressHeaderToScreen(void);
  void WriteProgressToScreen(void);

  MPIParallel& mpi_session_m;
  InputParser parser_m;
  IO file_io_m;
  Mesh mesh_m;
  HeatSolver heat_solver_m;
  PreciceAdapter* precice_m;
  IterationInfo step_info_m;
  SimulationRestrictions restrictions_m;
  SimulationGoals goals_m;
  SimulationOutput output_m;
};
}  // namespace chyps

#endif  // CHYPS_SIMULATION_HPP_
