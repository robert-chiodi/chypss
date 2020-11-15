// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chyps/simulation.hpp"

#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <utility>

#include <mpi.h>

#include "chyps/boundary_condition.hpp"
#include "chyps/boundary_condition_manager.hpp"
#include "chyps/debug_assert.hpp"
#include "chyps/git.hpp"
#include "chyps/heat_solver.hpp"
#include "chyps/input_parser.hpp"
#include "chyps/io.hpp"
#include "chyps/logger.hpp"
#include "chyps/mesh.hpp"
#include "chyps/mpi_parallel.hpp"
#include "chyps/precice_adapter.hpp"

namespace chyps {

int main(int argc, char** argv, MPIParallel& mpi_session,
         SpdlogLevel a_log_level) {
  StartLogger(mpi_session, a_log_level);
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Beginning simulation.");

  InputParser input_parser;
  input_parser.AddOption("Simulation/use_precice",
                         "If preCICE will be used to couple to another solver.",
                         false);
  input_parser.AddOption("Simulation/precice_config",
                         "XML File holding the configuration for preCICE",
                         std::string("../data/precice-config.xml"));
  input_parser.AddOption("Simulation/end_time", "Final time; start time is 0.",
                         0.5);
  input_parser.AddOption("Simulation/time_step", "Time step.", 1.0e-2);
  input_parser.AddOption(
      "Simulation/use_visit",
      "Save data files for VisIt (visit.llnl.gov) visualization.", false);

  input_parser.AddOption("Simulation/viz_steps",
                         "Visualize every n-th timestep.", 5);
  input_parser.AddOption(
      "Simulation/in_data",
      "Name of file (or BP4 directory) holding "
      "initial data. Do not include extension. Initial data can be generated "
      "from a previous simulation check point or using the tool "
      "chyps_initializer.");
  input_parser.AddOption(
      "Simulation/out_data",
      "Name of file (or BP4 directory) to write that holds "
      "data to restart from. Do not include extension. Pass ignore if you "
      "do not wish to write files.");
  input_parser.AddOption(
      "Simulation/option_output_name",
      "Name of file to write the options used for the simulation. Provides way "
      "to rerun the same simulation. Should end in the extension .json",
      "simulation_configuration.json");

  SPDLOG_LOGGER_INFO(MAIN_LOG, "Added main parser options.");

  IO file_io(mpi_session, "FILEIO");
  Mesh mesh(mpi_session, input_parser, &file_io);
  HeatSolver solver(input_parser, &file_io);

  const std::string input_file_name = argv[1];
  if (input_file_name == "help") {
    if (mpi_session.IAmRoot()) {
      input_parser.PrintOptions();
    }
    return 0;
  }
  if (input_file_name != "ignore") {
    input_parser.ParseFromFile(input_file_name, mpi_session);
  }
  argc -= 2;  // Skip executable and input file name
  input_parser.ParseCL(argc, argv + 2);
  if (mpi_session.IAmRoot()) {
    std::ofstream options_used(
        input_parser["Simulation/option_output_name"].get<std::string>());
    options_used << "/*\n" << GetGitDescriptionString() << "\n*/\n";
    const auto now = std::chrono::system_clock::now();
    const auto time = std::chrono::system_clock::to_time_t(now);
    options_used << "/*\n"
                 << "Simulation start time:\n"
                 << std::ctime(&time) << "*/\n";
    options_used << input_parser.WriteToString();
    options_used.close();
  }
  const auto in_data_name =
      input_parser["Simulation/in_data"].get<std::string>();
  const auto out_data_name =
      input_parser["Simulation/out_data"].get<std::string>();
  DEBUG_ASSERT(in_data_name != out_data_name, global_assert{},
               DebugLevel::ALWAYS{}, "Cannot read and write from same file.");
  file_io.SetRead(in_data_name);
  if (out_data_name != "ignore") {
    DEBUG_ASSERT(in_data_name != out_data_name, global_assert{},
                 DebugLevel::ALWAYS{}, "Cannot read and write from same file.");
    file_io.SetWrite(out_data_name);
    file_io.RootWriteAttribute("InputFile", input_parser.WriteToString());
  }
  file_io.RootWriteAttribute("GitInfo", GetGitDescriptionString());

  mesh.Initialize();

  PreciceAdapter* precice = nullptr;
  const bool use_precice = input_parser["Simulation/use_precice"].get<bool>();

  double* temperature_bc = nullptr;
  std::vector<double> vertex_positions;
  std::vector<int> vertex_indices;
  if (use_precice) {
    std::cout << "Need to figure new way to handle preCICE coupling within new "
                 "boundary condition manager "
              << std::endl;
    return -1;
    // precice =
    //     new PreciceAdapter("HeatSolver", "HeatSolverMesh",
    //                        input_parser["Simulation/precice_config"].get<std::string>(),
    //                        mpi_session.MyRank(),
    //                        mpi_session.NumberOfRanks());
    // std::tie(vertex_positions, vertex_indices) = mesh.GetBoundaryVertices(4);
    // precice->SetVertexPositions(vertex_positions);
    // precice->AddData("Temperature", DataOperation::READ);
    // temperature_bc = new double[precice->ScalarDataSize()];
    // // Note, BC will not be used prior to being set in
    // // ConductionOperator during solver.Advance().
    // std::fill(temperature_bc, temperature_bc + precice->ScalarDataSize(),
    //           -10.0);
  }

  // BoundaryCondition(BoundaryConditionType::HOMOGENEOUS_NEUMANN);
  // auto condition =
  //     BoundaryCondition(BoundaryConditionType::DIRICHLET, true, true);
  // // TESTING START
  // std::tie(vertex_positions, vertex_indices) = mesh.GetBoundaryVertices(4);
  // temperature_bc = new double[vertex_indices.size()];
  // for (std::size_t n = 0; n < vertex_indices.size(); ++n) {
  //   temperature_bc[n] = vertex_positions[2 * n];
  // }
  // condition.SetValues(vertex_indices.size(), temperature_bc,
  //                     vertex_indices.data(), false);
  // mesh.SetBoundaryCondition(4, condition);
  // // TESTING END
  // condition = BoundaryCondition(BoundaryConditionType::HOMOGENEOUS_NEUMANN,
  //                               false, false);
  // //  condition.SetValues(1.0);
  // mesh.SetBoundaryCondition(1, condition);
  // //  condition.SetValues(1.0);
  // mesh.SetBoundaryCondition(2, condition);
  // //  condition.SetValues(2.5);
  // mesh.SetBoundaryCondition(3, condition);
  // mesh.CommitBoundaryConditions();
  solver.Initialize(mesh);

  // Clear options from parser. Can still parse but
  // can not check all options supplied.
  input_parser.ClearOptions();

  int ti = 0;
  double time = 0.0;
  auto dt = input_parser["Simulation/time_step"].get<double>();
  if (file_io.IsReadModeActive()) {
    file_io.GetImmediateValue("CYCLE", &ti);
    file_io.GetImmediateValue("TIME", &time);
    double read_in_dt;
    file_io.GetImmediateValue("DT", &read_in_dt);
    if (read_in_dt > 0.0) {  // Initializer tool writes out a negative dt
      dt = std::min(dt, read_in_dt);
    }
  }
  double max_dt = dt;
  if (use_precice) {
    max_dt = precice->Initialize();
  }
  const auto final_time = input_parser["Simulation/end_time"].get<double>();
  const auto viz_steps = input_parser["Simulation/viz_steps"].get<int>();

  SPDLOG_LOGGER_INFO(MAIN_LOG, "Starting simulation at time {}", time);
  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Performing simulation to time {} with time steps of {}",
                     final_time, dt);

  // Initial writing of data
  file_io.BeginWriteStep(ti, time, dt);
  mesh.WriteMesh();
  solver.WriteFields(ti, time);
  file_io.EndWriteStep();
  file_io.WriteXMLSchema();

  bool last_step = false;
  while (!last_step) {
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
    ++ti;

    max_dt = use_precice ? precice->Advance(dt) : max_dt;

    if (ti % viz_steps == 0 || last_step) {
      if (mpi_session.IAmRoot()) {
        std::cout << "step " << ti << ", t = " << time << std::endl;
      }

      file_io.BeginWriteStep(ti, time, dt);
      solver.WriteFields(ti, time);
      file_io.EndWriteStep();
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
  return 0;
}

}  // namespace chyps
