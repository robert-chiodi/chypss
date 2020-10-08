// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <mpi.h>
#include <mfem/mfem.hpp>

#include "chyps/boundary_condition.hpp"
#include "chyps/heat_solver.hpp"
#include "chyps/input_parser.hpp"
#include "chyps/precice_adapter.hpp"

namespace chyps {

int main(int argc, char** argv) {
  int num_procs, myid;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  // PreciceAdapter("Heat_Solver", "TestXML", mpi.WorldRank(), mpi.WorldSize());

  InputParser input_parser;
  input_parser.AddOption("end_time", "-tf", "--t-final",
                         "Final time; start time is 0.", 0.5);
  input_parser.AddOption("time_step", "-dt", "--time-step", "Time step.",
                         1.0e-2);
  input_parser.AddOption(
      "use_visit", "-visit", "--visit-datafiles",
      "Save data files for VisIt (visit.l1lnl.gov) visualization.", false);
  input_parser.AddOption("viz_steps", "-vs", "--visualization-steps",
                         "Visualize every n-th timestep.", 5);

  HeatSolver solver(MPI_COMM_WORLD, input_parser);
  input_parser.ParseCL(argc, argv);

  double time = 0.0;
  bool last_step = false;
  double dt = input_parser["time_step"];
  const double final_time = input_parser["end_time"];
  const int viz_steps = input_parser["viz_steps"];
  solver.Initialize();

  auto condition = BoundaryCondition(BoundaryConditionType::DIRICHLET);
  condition.SetValues(10.0);
  solver.SetBoundaryCondition(1, condition);
  condition.SetValues(10.0);
  solver.SetBoundaryCondition(2, condition);
  condition = BoundaryCondition(BoundaryConditionType::HOMOGENEOUS_NEUMANN);
  solver.SetBoundaryCondition(3, condition);
  solver.SetBoundaryCondition(4, condition);

  solver.CommitBoundaryConditions();
  solver.ExportVisIt(0, 0.0);
  for (int ti = 1; !last_step; ++ti) {
    if (time + dt > final_time) {
      dt = final_time - time;
    }
    if (std::fabs(time + dt - final_time) < 1.0e-14) {
      last_step = true;
    }
    dt = solver.AdjustTimeStep(dt);
    dt = solver.Advance(time, dt);
    int number_of_boundaries_updated = solver.UpdateBoundaryConditions();
    ++number_of_boundaries_updated;  // Just to prevent warning while developing
    time += dt;

    if (ti % viz_steps == 0 || last_step) {
      std::cout << "step " << ti << ", t = " << time << std::endl;
      solver.ExportVisIt(ti, time);
    }
  }

  MPI_Finalize();
  return 0;
}

}  // namespace chyps

int main(int argc, char** argv) { return chyps::main(argc, argv); }
