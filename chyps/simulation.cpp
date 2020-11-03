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

#include <cmath>
#include <iostream>
#include <utility>

#include <mpi.h>

#include "chyps/boundary_condition.hpp"
#include "chyps/boundary_condition_manager.hpp"
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
  input_parser.AddOptionDefault(
      "Simulation/use_precice",
      "If preCICE will be used to couple to another solver.", false);
  input_parser.AddOptionDefault(
      "Simulation/precice_config",
      "XML File holding the configuration for preCICE",
      std::string("../data/precice-config.xml"));
  input_parser.AddOptionDefault("Simulation/end_time",
                                "Final time; start time is 0.", 0.5);
  input_parser.AddOptionDefault("Simulation/time_step", "Time step.", 1.0e-2);
  input_parser.AddOptionDefault(
      "Simulation/use_visit",
      "Save data files for VisIt (visit.llnl.gov) visualization.", false);

  input_parser.AddOptionDefault("Simulation/viz_steps",
                                "Visualize every n-th timestep.", 5);
  input_parser.AddOptionDefault(
      "Simulation/in_data",
      "Name of file (or BP4 directory) holding "
      "data for restart. Do not include extension. Pass ignore if you "
      "do not wish to restart from a file.",
      std::string("ignore"));
  input_parser.AddOptionDefault(
      "Simulation/out_data",
      "Name of file (or BP4 directory) to write that holds "
      "data to restart from. Do not include extension. Pass ignore if you "
      "do not wish to write files.",
      std::string("CHyPSDataOut"));
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
    // Just have root processor parse file, then communicate
    std::vector<std::uint8_t> v_bson;
    int size = 0;
    if (mpi_session.IAmRoot()) {
      const bool good_read = input_parser.ParseFromFile(input_file_name);
      if (!good_read) {
        std::cout << "Trouble reading input file " << input_file_name
                  << std::endl;
        std::cout << "First argument should be name of input file. \n";
        std::cout
            << "Use the input file name \"ignore\" to use no input file. \n";
        std::cout
            << "Use the input file name \"help\" to export available options."
            << std::endl;
        std::exit(-1);
      }
      v_bson = input_parser.ToBSON();
      size = static_cast<int>(v_bson.size());
    }
    MPI_Bcast(&size, 1, MPI_INT, 0, mpi_session.GetComm());
    v_bson.resize(size);  // Will do nothing for rank 0
    MPI_Bcast(v_bson.data(), size, MPI_BYTE, 0, mpi_session.GetComm());
    if (mpi_session.IAmNotRoot()) {
      input_parser.SetFromBSON(v_bson);
    }
  }
  argc += -2;  // Skip executable and input file name
  input_parser.ParseCL(argc, argv + 2);
  input_parser.WriteToFile("simulation_configuration.json");
  assert(input_parser.AllOptionsSet("Simulation"));

  const auto in_data_name =
      input_parser["Simulation/in_data"].get<std::string>();
  const auto out_data_name =
      input_parser["Simulation/out_data"].get<std::string>();
  if (in_data_name != "ignore") {
    assert(in_data_name != out_data_name);
    file_io.SetRead(in_data_name);
  }
  if (out_data_name != "ignore") {
    assert(in_data_name != out_data_name);
    file_io.SetWrite(out_data_name);
    file_io.RootWriteAttribute("InputFile", input_parser.WriteToString());
  }

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
    dt = std::min(dt, read_in_dt);
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
