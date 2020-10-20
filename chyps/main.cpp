// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cmath>
#include <utility>

#include <mpi.h>

#include "chyps/boundary_condition.hpp"
#include "chyps/heat_solver.hpp"
#include "chyps/input_parser.hpp"
#include "chyps/logger.hpp"
#include "chyps/mesh.hpp"
#include "chyps/precice_adapter.hpp"

namespace chyps {

int main(int argc, char** argv) {
  // FIXME : Bring mesh building to outside of the solver class. Solver should
  // just be told the mesh to solve the heat equation on
  int num_procs, myid;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  StartLogger(myid, num_procs, SpdlogLevel::INFO);
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Beginning simulation.");

  InputParser input_parser;
  input_parser.AddOption("use_precice", "-p", "--precice",
                         "If preCICE will be used to couple to another solver.",
                         false);
  input_parser.AddOption("precice_config", "-pc", "--preCICE-config",
                         "XML File holding the configuration for preCICE",
                         std::string("../data/precice-config.xml"));
  input_parser.AddOption("end_time", "-tf", "--t-final",
                         "Final time; start time is 0.", 0.5);
  input_parser.AddOption("time_step", "-dt", "--time-step", "Time step.",
                         1.0e-2);
  input_parser.AddOption(
      "use_visit", "-visit", "--visit-datafiles",
      "Save data files for VisIt (visit.l1lnl.gov) visualization.", false);
  input_parser.AddOption("viz_steps", "-vs", "--visualization-steps",
                         "Visualize every n-th timestep.", 5);
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Added main parser options.");

  Mesh mesh(MPI_COMM_WORLD, input_parser);
  HeatSolver solver(input_parser);
  input_parser.ParseCL(argc, argv);
  assert(input_parser.AllOptionsSet());

  mesh.Initialize();

  PreciceAdapter* precice = nullptr;
  const bool use_precice = input_parser["use_precice"];

  double* temperature_bc = nullptr;
  std::vector<double> vertex_positions;
  std::vector<int> vertex_indices;
  if (use_precice) {
    precice =
        new PreciceAdapter("HeatSolver", "HeatSolverMesh",
                           input_parser["precice_config"], myid, num_procs);
    std::tie(vertex_positions, vertex_indices) = mesh.GetBoundaryVertices(4);
    precice->SetVertexPositions(vertex_positions);
    precice->AddData("Temperature", DataOperation::READ);
    temperature_bc = new double[precice->ScalarDataSize()];
    // Note, BC will not be used prior to being set in
    // solver.UpdateBoundaryConditions.
    std::fill(temperature_bc, temperature_bc + precice->ScalarDataSize(),
              -10.0);
  }

  BoundaryCondition(BoundaryConditionType::HOMOGENEOUS_NEUMANN);
  auto condition =
      BoundaryCondition(BoundaryConditionType::DIRICHLET, true, false);
  // TESTING START
  std::tie(vertex_positions, vertex_indices) = mesh.GetBoundaryVertices(4);
  temperature_bc = new double[vertex_indices.size()];
  for (std::size_t n = 0; n < vertex_indices.size(); ++n) {
    temperature_bc[n] = vertex_positions[2 * n];
  }
  condition.SetValues(vertex_indices.size(), temperature_bc,
                      vertex_indices.data(), false);
  mesh.SetBoundaryCondition(4, condition);
  // TESTING END
  condition = BoundaryCondition(BoundaryConditionType::NEUMANN, false, true);
  condition.SetValues(1.0);
  mesh.SetBoundaryCondition(1, condition);
  condition.SetValues(5.0);
  mesh.SetBoundaryCondition(2, condition);
  condition.SetValues(10.0);
  mesh.SetBoundaryCondition(3, condition);
  // mesh.SetBoundaryCondition(4, condition);

  mesh.CommitBoundaryConditions();
  solver.Initialize(mesh);
  solver.ExportVisIt(0, 0.0);

  double time = 0.0;
  bool last_step = false;
  double dt = input_parser["time_step"];
  double max_dt = dt;
  if (use_precice) {
    max_dt = precice->Initialize();
  }
  const double final_time = input_parser["end_time"];
  const int viz_steps = input_parser["viz_steps"];

  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Performing simulation to time {} with time steps of {}",
                     final_time, dt);

  for (int ti = 1; !last_step; ++ti) {
    SPDLOG_LOGGER_INFO(MAIN_LOG, "Starting time-iteration {} for time {:8.6E}",
                       ti, time);
    if (use_precice) {
      precice->ReadBlockScalarData("Temperature", temperature_bc);
    }

    double dt_adj = solver.AdjustTimeStep(dt);
    dt_adj = std::min(max_dt, dt_adj);

    if (time + dt_adj > final_time) {
      dt_adj = final_time - time;
    }
    if (std::fabs(time + dt_adj - final_time) < 1.0e-14) {
      last_step = true;
      SPDLOG_LOGGER_INFO(MAIN_LOG, "Begun last step.");
    }

    dt = solver.Advance(time, dt_adj);
    time += dt;

    max_dt = use_precice ? precice->Advance(dt) : max_dt;

    if (ti % viz_steps == 0 || last_step) {
      if (myid == 0) {
        std::cout << "step " << ti << ", t = " << time << std::endl;
      }
      solver.ExportVisIt(ti, time);
    }

    if (use_precice && !precice->IsCouplingOngoing()) {
      SPDLOG_LOGGER_INFO(MAIN_LOG,
                         "Precice coupling done at iteration {} and time {}",
                         ti, time);
      break;
    }
  }

  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Final simulation time of {:8.6E} reached. Finalizing "
                     "simulation and cleaning up.",
                     time);

  if (use_precice) {
    precice->Finalize();
  }
  MPI_Finalize();
  return 0;
}

}  // namespace chyps

int main(int argc, char** argv) { return chyps::main(argc, argv); }
