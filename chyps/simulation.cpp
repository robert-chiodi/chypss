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
  StartLogger(mpi_session, a_log_level);
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
      monitor_manager_m(),
      timer_manager_m(&mpi_session_m),
      mesh_m(parser_m, *this),
      heat_solver_m(parser_m, *this),
      precice_adapter_m(nullptr),
      step_info_m(0.0, 0.0, 0.0),
      restrictions_m(),
      goals_m(),
      output_m() {
  // Log needs to have been started before this so that logging in
  // constructors above do not segfault. This will reset
  // log level to a_log_level though.
  StartLogger(mpi_session_m, a_log_level);
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Beginning simulation.");

  this->GatherOptions();
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Added main parser options.");

  this->GetTimerManager().AddTimer("Total");
}

void Simulation::Initialize(int argc, char** argv) {
  this->GetTimerManager().StartTimer("Total");
  this->ParseOptions(argc, argv);
  if (parser_m["Simulation/write_to_monitor"].get<bool>() &&
      mpi_session_m.IAmRoot()) {
    this->GetMonitorManager().Initialize("monitor");
  }
  this->SetupFileIO();
  // Note, order is important here. Mesh must be initialized first.
  this->ActivatePrecice();
  mesh_m.Initialize();
  heat_solver_m.Initialize();
  this->InitializeStepInfo();
  this->InitializeRestrictions();
  this->InitializeGoals();
  this->InitializeOutputs();

  SPDLOG_LOGGER_INFO(MAIN_LOG, "Starting simulation at time {}",
                     step_info_m.time);
  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Performing simulation to time {} with time steps of {}",
                     goals_m.end_time, step_info_m.dt);
  this->WriteInitialConditions();
  this->GetTimerManager().StopTimer("Total");
  this->GetTimerManager().CreateMonitorFile("timing",
                                            this->GetMonitorManager());
  this->GetTimerManager().PushTimesToMonitor();
  if (mpi_session_m.IAmRoot()) {
    this->WriteProgressHeaderToScreen();
    this->WriteProgressToScreen();
    this->GetMonitorManager().WriteStepToFiles(
        step_info_m.iteration, step_info_m.time, step_info_m.dt);
  }
}

void Simulation::RunToEnd(void) {
  bool last_step = false;
  while (!last_step) {
    this->GetTimerManager().ResetAllTimers();  // Reset all timers to 0
    this->GetTimerManager().StartTimer("Total");

    SPDLOG_LOGGER_INFO(MAIN_LOG, "Starting time-iteration {} for time {:8.6E}",
                       step_info_m.iteration, step_info_m.time);
    // if (this->PreciceActive()) {
    //   precice_adapter_m->ReadBlockScalarData("Temperature", temperature_bc);
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

    if (this->PreciceActive()) {
      heat_solver_m.WriteDataToPrecice();
    }
    step_info_m.dt = heat_solver_m.Advance(step_info_m.time, dt_adj);
    step_info_m.time += step_info_m.dt;
    ++step_info_m.iteration;

    restrictions_m.max_dt = this->PreciceActive()
                                ? precice_adapter_m->Advance(step_info_m.dt)
                                : restrictions_m.max_dt;

    // Write out visualization files and checkpoints
    if (step_info_m.iteration % output_m.visualization_steps == 0 ||
        last_step) {
      this->WriteIterationConditions(step_info_m);
    }

    this->GetTimerManager().StopTimer("Total");
    this->GetTimerManager().PushTimesToMonitor();

    // Write out progress to screen and monitor files.
    if (mpi_session_m.IAmRoot()) {
      this->WriteProgressToScreen();
      this->GetMonitorManager().WriteStepToFiles(
          step_info_m.iteration, step_info_m.time, step_info_m.dt);
    }

    if (this->PreciceActive() && !precice_adapter_m->IsCouplingOngoing()) {
      SPDLOG_LOGGER_INFO(MAIN_LOG,
                         "Precice coupling done at iteration {} and time {}",
                         step_info_m.iteration, step_info_m.time);
      break;
    }
  }

  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Final simulation time of {:8.6E} reached. Finalizing "
                     "simulation and cleaning up.",
                     step_info_m.time);

  if (this->PreciceActive()) {
    precice_adapter_m->Finalize();
  }
}

InputParser& Simulation::GetParser(void) { return parser_m; }

const InputParser& Simulation::GetParser(void) const { return parser_m; }

Mesh& Simulation::GetMesh(void) { return mesh_m; }

const Mesh& Simulation::GetMesh(void) const { return mesh_m; }

IO& Simulation::GetIO(void) { return file_io_m; }
const IO& Simulation::GetIO(void) const { return file_io_m; }

const MPIParallel& Simulation::GetMPI(void) const { return mpi_session_m; }

MonitorManager& Simulation::GetMonitorManager(void) {
  return monitor_manager_m;
}

const MonitorManager& Simulation::GetMonitorManager(void) const {
  return monitor_manager_m;
}

TimerManager& Simulation::GetTimerManager(void) { return timer_manager_m; }
const TimerManager& Simulation::GetTimerManager(void) const {
  return timer_manager_m;
}

bool Simulation::PreciceActive(void) const {
  return precice_adapter_m != nullptr;
}

PreciceAdapter& Simulation::GetPreciceAdapter(void) {
  DEBUG_ASSERT(this->PreciceActive(), global_assert{}, DebugLevel::CHEAP{},
               "Precice must be active to request the adapter.");
  return *precice_adapter_m;
}

const PreciceAdapter& Simulation::GetPreciceAdapter(void) const {
  DEBUG_ASSERT(this->PreciceActive(), global_assert{}, DebugLevel::CHEAP{},
               "Precice must be active to request the adapter.");
  return *precice_adapter_m;
}

Simulation::~Simulation(void) {
  delete precice_adapter_m;
  precice_adapter_m = nullptr;
  if (parser_m["Simulation/output_screen"].get<std::string>() != "cout" &&
      parser_m["Simulation/output_screen"].get<std::string>() != "off") {
    fclose(output_m.output_screen);
  }
  output_m.output_screen = nullptr;
}

void Simulation::GatherOptions(void) {
  parser_m.AddOption("Simulation/use_precice",
                     "If preCICE will be used to couple to another solver.",
                     false);
  parser_m.AddOption("Simulation/precice_solver_name",
                     "Name to be used for solver in preCICE.", "CHyPS");
  parser_m.AddOption("Simulation/precice_config",
                     "XML File holding the configuration for preCICE",
                     std::string("../data/precice-config.xml"));
  parser_m.AddOption("Simulation/end_time", "Final time; start time is 0.",
                     0.5);
  parser_m.AddOption("Simulation/time_step", "Time step.", 1.0e-2);
  parser_m.AddOption(
      "Simulation/use_visit",
      "Save data files for VisIt (visit.llnl.gov) visualization.", false);
  parser_m.AddOption("Simulation/visit_directory",
                     "Directory to export VisIt files to.", "Visit/");

  parser_m.AddOption("Simulation/viz_steps", "Visualize every n-th timestep.",
                     5);
  parser_m.AddOption("Simulation/output_screen",
                     "Location to write output that would normally go to the "
                     "screen. Supply \"cout\" to write to terminal.",
                     "cout");
  parser_m.AddOption(
      "Simulation/write_to_monitor",
      "Determines whether MonitorManager object will be initialized "
      "and monitor files written to.",
      true);
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
  precice_adapter_m =
      parser_m["Simulation/use_precice"].get<bool>()
          ? new PreciceAdapter(
                parser_m["Simulation/precice_solver_name"].get<std::string>(),
                parser_m["Simulation/precice_config"].get<std::string>(),
                mpi_session_m.MyRank(), mpi_session_m.NumberOfRanks())
          : nullptr;
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
    restrictions_m.max_dt = precice_adapter_m->Initialize();
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
