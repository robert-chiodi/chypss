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
#include <utility>

#include "chyps/debug_assert.hpp"
#include "chyps/git.hpp"

namespace chyps {

int main(int argc, char** argv, MPIParallel& mpi_session,
         const SpdlogLevel a_log_level) {
  Simulation simulation(mpi_session, a_log_level);
  simulation.Initialize(argc, argv);
  simulation.RunToEnd();
  return 0;
}

Simulation::Simulation(MPIParallel& a_mpi_session,
                       const SpdlogLevel a_log_level)
    : mpi_session_m(a_mpi_session),
      parser_m(),
      file_io_m(mpi_session_m, "FILEIO"),
      mesh_m(parser_m, *this),
      heat_solver_m(parser_m, *this),
      precice_m(nullptr),
      step_info_m(0.0, 0.0, 0.0),
      restrictions_m(),
      goals_m(),
      output_m() {
  StartLogger(mpi_session_m, a_log_level);
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Beginning simulation.");

  this->GatherOptions();

  SPDLOG_LOGGER_INFO(MAIN_LOG, "Added main parser options.");
}

void Simulation::Initialize(int argc, char** argv) {
  this->ParseOptions(argc, argv);
  this->SetupFileIO();
  // Note, order is important here. Mesh must be initialized first.
  mesh_m.Initialize();
  heat_solver_m.Initialize();
  this->ActivatePrecice();
  this->InitializeStepInfo();
  this->InitializeRestrictions();
  this->InitializeGoals();
  this->InitializeOutputs();

  SPDLOG_LOGGER_INFO(MAIN_LOG, "Starting simulation at time {}", time);
  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Performing simulation to time {} with time steps of {}",
                     goals_m.end_time, step_info_m.dt);
  this->WriteInitialConditions();
  this->WriteProgressHeaderToScreen();
}

void Simulation::RunToEnd(void) {
  bool last_step = false;
  while (!last_step) {
    SPDLOG_LOGGER_INFO(MAIN_LOG, "Starting time-iteration {} for time {:8.6E}",
                       step_info_m.iteration, step_info_m.time);
    // if (this->PreciceActive()) {
    //   precice_m->ReadBlockScalarData("Temperature", temperature_bc);
    // }

    double dt_adj = heat_solver_m.AdjustTimeStep(step_info_m.dt);
    dt_adj = std::min(restrictions_m.max_dt, dt_adj);

    if (step_info_m.time + dt_adj > goals_m.end_time) {
      dt_adj = goals_m.end_time - step_info_m.time;
    }
    if (std::fabs(step_info_m.time + dt_adj - goals_m.end_time) < 1.0e-14) {
      last_step = true;
      SPDLOG_LOGGER_INFO(MAIN_LOG, "Begun last step.");
    }

    step_info_m.dt = heat_solver_m.Advance(step_info_m.time, dt_adj);
    step_info_m.time += step_info_m.dt;
    ++step_info_m.iteration;

    restrictions_m.max_dt = this->PreciceActive()
                                ? precice_m->Advance(step_info_m.dt)
                                : restrictions_m.max_dt;

    if (mpi_session_m.IAmRoot()) {
      this->WriteProgressToScreen();
    }

    if (step_info_m.iteration % output_m.visualization_steps == 0 ||
        last_step) {
      this->WriteIterationConditions(step_info_m);
    }

    if (this->PreciceActive() && !precice_m->IsCouplingOngoing()) {
      SPDLOG_LOGGER_INFO(MAIN_LOG,
                         "Precice coupling done at iteration {} and time {}",
                         step_info_m.iteration, step_info_m.time);
      break;
    }
  }

  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Final simulation time of {:8.6E} reached. Finalizing "
                     "simulation and cleaning up.",
                     time);

  if (this->PreciceActive()) {
    precice_m->Finalize();
  }
}

bool Simulation::PreciceActive(void) const { return precice_m != nullptr; }

Mesh& Simulation::GetMesh(void) { return mesh_m; }

const Mesh& Simulation::GetMesh(void) const { return mesh_m; }

IO& Simulation::GetIO(void) { return file_io_m; }
const IO& Simulation::GetIO(void) const { return file_io_m; }

const MPIParallel& Simulation::GetMPI(void) const { return mpi_session_m; }

Simulation::~Simulation(void) {
  delete precice_m;
  precice_m = nullptr;
  if (parser_m["Simulation/output_screen"] != "cout") {
    fclose(output_m.output_screen);
    delete output_m.output_screen;
  }
  output_m.output_screen = nullptr;
}

void Simulation::GatherOptions(void) {
  parser_m.AddOption("Simulation/use_precice",
                     "If preCICE will be used to couple to another solver.",
                     false);
  parser_m.AddOption("Simulation/precice_config",
                     "XML File holding the configuration for preCICE",
                     std::string("../data/precice-config.xml"));
  parser_m.AddOption("Simulation/end_time", "Final time; start time is 0.",
                     0.5);
  parser_m.AddOption("Simulation/time_step", "Time step.", 1.0e-2);
  parser_m.AddOption(
      "Simulation/use_visit",
      "Save data files for VisIt (visit.llnl.gov) visualization.", false);

  parser_m.AddOption("Simulation/viz_steps", "Visualize every n-th timestep.",
                     5);
  parser_m.AddOption("Simulation/output_screen",
                     "Location to write output that would normally go to the "
                     "screen. Supply \"cout\" to write to terminal.",
                     "cout");
  parser_m.AddOption(
      "Simulation/in_data",
      "Name of file (or BP4 directory) holding "
      "initial data. Do not include extension. Initial data can be "
      "generated "
      "from a previous simulation check point or using the tool "
      "chyps_initializer.");
  parser_m.AddOption(
      "Simulation/out_data",
      "Name of file (or BP4 directory) to write that holds "
      "data to restart from. Do not include extension. Pass ignore if you "
      "do not wish to write files.");
  parser_m.AddOption(
      "Simulation/option_output_name",
      "Name of file to write the options used for the simulation. Provides way "
      "to rerun the same simulation. Should end in the extension .json",
      "simulation_configuration.json");
}

void Simulation::ParseOptions(int argc, char** argv) {
  const std::string input_file_name = argv[1];
  if (input_file_name == "help") {
    if (mpi_session_m.IAmRoot()) {
      parser_m.PrintOptions();
    }
    DEBUG_ASSERT(false, global_assert{}, DebugLevel::ALWAYS{});
  }
  if (input_file_name != "ignore") {
    parser_m.ParseFromFile(input_file_name, mpi_session_m);
  }
  argc -= 2;  // Skip executable and input file name
  parser_m.ParseCL(argc, argv + 2);
  if (mpi_session_m.IAmRoot()) {
    std::ofstream options_used(
        parser_m["Simulation/option_output_name"].get<std::string>());
    options_used << "/*\n" << GetGitDescriptionString() << "\n*/\n";
    const auto now = std::chrono::system_clock::now();
    const auto time = std::chrono::system_clock::to_time_t(now);
    options_used << "/*\n"
                 << "Simulation start time:\n"
                 << std::ctime(&time) << "*/\n";
    options_used << parser_m.WriteToString();
    options_used.close();
  }
  parser_m.ClearOptions();
}

void Simulation::SetupFileIO(void) {
  const auto in_data_name = parser_m["Simulation/in_data"].get<std::string>();
  const auto out_data_name = parser_m["Simulation/out_data"].get<std::string>();
  DEBUG_ASSERT(in_data_name != out_data_name, global_assert{},
               DebugLevel::ALWAYS{}, "Cannot read and write from same file.");

  file_io_m.SetRead(in_data_name);
  if (out_data_name != "ignore") {
    DEBUG_ASSERT(in_data_name != out_data_name, global_assert{},
                 DebugLevel::ALWAYS{}, "Cannot read and write from same file.");
    file_io_m.SetWrite(out_data_name);
    file_io_m.RootWriteAttribute("InputFile", parser_m.WriteToString());
  }
  file_io_m.RootWriteAttribute("GitInfo", GetGitDescriptionString());
}

void Simulation::ActivatePrecice(void) {
  precice_m =
      parser_m["Simulation/use_precice"].get<bool>()
          ? new PreciceAdapter(
                "HeatSolver", "HeatSolverMesh",
                parser_m["Simulation/precice_config"].get<std::string>(),
                mpi_session_m.MyRank(), mpi_session_m.NumberOfRanks())
          : nullptr;

  // double* temperature_bc = nullptr;
  std::vector<double> vertex_positions;
  std::vector<int> vertex_indices;
  if (this->PreciceActive()) {
    std::cout << "Need to figure new way to handle preCICE coupling within new "
                 "boundary condition manager "
              << std::endl;
    return;
    // std::tie(vertex_positions, vertex_indices) =
    // mesh_m.GetBoundaryVertices(4);
    // precice_m->SetVertexPositions(vertex_positions);
    // precice_m->AddData("Temperature", DataOperation::READ);
    // temperature_bc = new double[precice_m->ScalarDataSize()];
    // // Note, BC will not be used prior to being set in
    // // ConductionOperator during solver.Advance().
    // std::fill(temperature_bc, temperature_bc + precice_m->ScalarDataSize(),
    //           -10.0);
  }

  // BoundaryCondition(BoundaryConditionType::HOMOGENEOUS_NEUMANN);
  // auto condition =
  //     BoundaryCondition(BoundaryConditionType::DIRICHLET, true, true);
  // // TESTING START
  // std::tie(vertex_positions, vertex_indices) = mesh_m.GetBoundaryVertices(4);
  // temperature_bc = new double[vertex_indices.size()];
  // for (std::size_t n = 0; n < vertex_indices.size(); ++n) {
  //   temperature_bc[n] = vertex_positions[2 * n];
  // }
  // condition.SetValues(vertex_indices.size(), temperature_bc,
  //                     vertex_indices.data(), false);
  // mesh_m.SetBoundaryCondition(4, condition);
  // // TESTING END
  // condition = BoundaryCondition(BoundaryConditionType::HOMOGENEOUS_NEUMANN,
  //                               false, false);
  // //  condition.SetValues(1.0);
  // mesh_m.SetBoundaryCondition(1, condition);
  // //  condition.SetValues(1.0);
  // mesh_m.SetBoundaryCondition(2, condition);
  // //  condition.SetValues(2.5);
  // mesh_m.SetBoundaryCondition(3, condition);
  // mesh_m.CommitBoundaryConditions();
}

void Simulation::InitializeStepInfo(void) {
  step_info_m.dt = parser_m["Simulation/time_step"].get<double>();
  if (file_io_m.IsReadModeActive()) {
    file_io_m.GetImmediateValue("CYCLE", &step_info_m.iteration);
    file_io_m.GetImmediateValue("TIME", &step_info_m.time);
    double read_in_dt;
    file_io_m.GetImmediateValue("DT", &read_in_dt);
    if (read_in_dt > 0.0) {  // Initializer tool writes out a negative dt
      step_info_m.dt = std::min(step_info_m.dt, read_in_dt);
    }
  }
}

void Simulation::InitializeRestrictions(void) {
  restrictions_m.max_dt = step_info_m.dt;
  if (this->PreciceActive()) {
    restrictions_m.max_dt = precice_m->Initialize();
  }
}

void Simulation::InitializeGoals(void) {
  goals_m.end_time = parser_m["Simulation/end_time"].get<double>();
}

void Simulation::InitializeOutputs(void) {
  const std::string& screen_name =
      parser_m["Simulation/output_screen"].get<std::string>();
  if (screen_name == "cout") {
    output_m.output_screen = stdout;
  } else if (screen_name == "off") {
    output_m.output_screen = nullptr;
  } else {
    output_m.output_screen = fopen(screen_name.c_str(), "w");
  }

  output_m.visualization_steps = parser_m["Simulation/viz_steps"].get<int>();
}

void Simulation::WriteInitialConditions(void) {
  file_io_m.BeginWriteStep(step_info_m.iteration, step_info_m.time,
                           step_info_m.dt);
  mesh_m.WriteMesh();
  heat_solver_m.WriteFields(step_info_m.iteration, step_info_m.time);
  file_io_m.EndWriteStep();
  file_io_m.WriteXMLSchema();
}

void Simulation::WriteIterationConditions(const IterationInfo& a_iter) {
  file_io_m.BeginWriteStep(a_iter.iteration, a_iter.time, a_iter.dt);
  heat_solver_m.WriteFields(step_info_m.iteration, step_info_m.time);
  file_io_m.EndWriteStep();
}

void Simulation::WriteProgressHeaderToScreen(void) {
  output_m.printf("%-15s %-17s %-17s\n", "Iteration", "Time [s]", "dt [s]");
}

void Simulation::WriteProgressToScreen(void) {
  output_m.printf("%-15lld %-17.6E %-17.6E\n", step_info_m.iteration,
                  step_info_m.time, step_info_m.dt);
}

}  // namespace chyps
