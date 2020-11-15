// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chyps/heat_solver.hpp"

#include <array>

#include "chyps/debug_assert.hpp"
#include "chyps/logger.hpp"

namespace chyps {

HeatSolver::HeatSolver(InputParser& a_parser, IO* a_file_io)
    : parser_m(a_parser),
      file_io_m(a_file_io),
      mesh_m(nullptr),
      boundary_conditions_m(),
      ode_solver_m(nullptr),
      element_collection_m(nullptr),
      element_space_m(nullptr),
      coarse_element_collection_m(nullptr),
      coarse_element_space_m(nullptr),
      operator_m(nullptr),
      temperature_m(),
      visit_collection_m(nullptr) {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Constructing HeatSolver object.");
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Gathering options of HeatSolver class.");
  this->GatherOptions();
  this->CreateBoundaryConditionManagers();
}

void HeatSolver::CreateBoundaryConditionManagers(void) {
  boundary_conditions_m.emplace(std::make_pair(
      "temperature",
      BoundaryConditionManager(parser_m,
                               "HeatSolver/BoundaryConditions/temperature")));
}

void HeatSolver::Initialize(Mesh& a_mesh) {
  // Construct mesh, allocate and construct operators, and perform all setup for
  // time advancement
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Initializing HeatSolver class.");
  mesh_m = &a_mesh;
  this->InitializeBoundaryConditions();
  this->AllocateVariablesAndOperators();
  this->SetODESolver();
  this->RegisterFieldsForIO();
}

double HeatSolver::AdjustTimeStep(const double a_proposed_dt) const {
  const double adjusted_dt = a_proposed_dt;
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Proposed dt of {:8.6E} adjusted to {:8.6E}",
                     a_proposed_dt, adjusted_dt);
  return adjusted_dt;
}

double HeatSolver::Advance(const double a_time, const double a_time_step) {
  double taken_time_step = a_time_step;
  double time = a_time;
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Setting operator parameters at time {:8.6E}",
                     a_time);
  operator_m->UpdateBoundaryConditions(temperature_m);
  operator_m->SetParameters(temperature_m);
  ode_solver_m->Step(temperature_m, time, taken_time_step);
  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Advanced solution from time {:8.6E} to time {:8.6E}",
                     a_time, time);
  return taken_time_step;
}

void HeatSolver::WriteFields(const int a_cycle, const double a_time) {
  if (parser_m["Simulation/use_visit"].get<bool>()) {
    visit_collection_m->UpdateField("temperature", temperature_m);
    visit_collection_m->WriteOutFields(a_cycle, a_time);
  }

  if (this->FileWritingEnabled()) {
    DEBUG_ASSERT(file_io_m->IsWriteModeActive(), global_assert{},
                 DebugLevel::CHEAP{},
                 "IO must be in write mode for writing to fields.");
    DEBUG_ASSERT(file_io_m->OngoingWriteStep(), global_assert{},
                 DebugLevel::CHEAP{},
                 "An ongoing IO step is required for writing.");
    mfem::ParGridFunction temperature_gf(element_space_m);
    temperature_gf.SetFromTrueDofs(temperature_m);
    file_io_m->PutDeferred("HeatSolver/temperature", temperature_gf);
    file_io_m->PerformPuts();
  }
}

void HeatSolver::GatherOptions(void) {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Gathering options to look for in parser.");

  parser_m.AddOption("HeatSolver/order",
                     "Order (degree) of the finite elements.", 2);
  parser_m.AddOption(
      "HeatSolver/time_advancement",
      "ODE solver: 1 - Backward Euler, 2 - SDIRK2, 3 - SDIRK3,\n\t"
      "\t   11 - Forward Euler, 12 - RK2, 13 - RK3 SSP, 14 - RK4.",
      3);
  parser_m.AddOption(
      "HeatSolver/kappa",
      "Array of thermal coefficients representing tensor in column major "
      "ordering. Should be MESH_DIM*MESH_DIM in size.");
  // TODO Add a flag and way to specify which variables we wish to export to
  // VisIt.

  SPDLOG_LOGGER_INFO(MAIN_LOG, "Add all options to parser.");
}

void HeatSolver::SetODESolver(void) {
  const auto ode_solver_type =
      parser_m["HeatSolver/time_advancement"].get<int>();
  switch (ode_solver_type) {
    // Implicit L-stable methods
    case 1:
      ode_solver_m = new mfem::BackwardEulerSolver;
      SPDLOG_LOGGER_INFO(MAIN_LOG,
                         "Implicit time advancer chosen: BackwardEulerSolver");
      break;
    case 2:
      ode_solver_m = new mfem::SDIRK23Solver(2);
      SPDLOG_LOGGER_INFO(MAIN_LOG,
                         "Implicit time advancer chosen: SDIRK23Solver");
      break;
    case 3:
      ode_solver_m = new mfem::SDIRK33Solver;
      SPDLOG_LOGGER_INFO(MAIN_LOG,
                         "Implicit time advancer chosen: SDIRK33Solver");
      break;
    // Explicit methods
    case 11:
      ode_solver_m = new mfem::ForwardEulerSolver;
      SPDLOG_LOGGER_INFO(MAIN_LOG,
                         "Explicit time advancer chosen: ForwardEulerSolver");
      break;
    case 12:
      ode_solver_m = new mfem::RK2Solver(0.5);
      SPDLOG_LOGGER_INFO(MAIN_LOG, "Explicit time advancer chosen: RK2Solver");
      break;  // midpoint method
    case 13:
      ode_solver_m = new mfem::RK3SSPSolver;
      SPDLOG_LOGGER_INFO(MAIN_LOG,
                         "Explicit time advancer chosen: RK3SSPSolver");
      break;
    case 14:
      ode_solver_m = new mfem::RK4Solver;
      SPDLOG_LOGGER_INFO(MAIN_LOG, "Explicit time advancer chosen: RK4Solver");
      break;
    case 15:
      ode_solver_m = new mfem::GeneralizedAlphaSolver(0.5);
      SPDLOG_LOGGER_INFO(
          MAIN_LOG, "Explicit time advancer chosen: GeneralizedAlphaSolver");
      break;
    // Implicit A-stable methods (not L-stable)
    case 22:
      ode_solver_m = new mfem::ImplicitMidpointSolver;
      SPDLOG_LOGGER_INFO(
          MAIN_LOG, "Implicit time advancer chosen: ImplicitMidpointSolver");
      break;
    case 23:
      ode_solver_m = new mfem::SDIRK23Solver;
      SPDLOG_LOGGER_INFO(MAIN_LOG,
                         "Implicit time advancer chosen: SDIRK23Solver");
      break;
    case 24:
      ode_solver_m = new mfem::SDIRK34Solver;
      SPDLOG_LOGGER_INFO(MAIN_LOG,
                         "Implicit time advancer chosen: SDIRK34Solver");
      break;
    default:
      SPDLOG_LOGGER_CRITICAL(MAIN_LOG, "Unknown ODE solver type: {}",
                             ode_solver_type);
      DEBUG_ASSERT(
          false, global_assert{}, DebugLevel::ALWAYS{},
          "Unknown ODE solver type: " + std::to_string(ode_solver_type));
  }

  SPDLOG_LOGGER_INFO(MAIN_LOG, "Initializing ode_solver");
  ode_solver_m->Init(*operator_m);
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Initialized ode_solver");
}

void HeatSolver::AllocateVariablesAndOperators(void) {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Creating H1 collection with {} order elements",
                     parser_m["HeatSolver/order"].get<int>());
  element_collection_m = new mfem::H1_FECollection(
      parser_m["HeatSolver/order"].get<int>(), mesh_m->GetDimension());

  SPDLOG_LOGGER_INFO(MAIN_LOG, "Generating element space from H1 collection");
  element_space_m = new mfem::ParFiniteElementSpace(&(mesh_m->GetMfemMesh()),
                                                    element_collection_m);

  coarse_element_collection_m =
      new mfem::H1_FECollection(1, mesh_m->GetDimension());

  SPDLOG_LOGGER_INFO(MAIN_LOG, "Generating element space from H1 collection");
  coarse_element_space_m = new mfem::ParFiniteElementSpace(
      &(mesh_m->GetMfemMesh()), coarse_element_collection_m);

  // Allocates temperature and sets initial values.
  this->SetInitialConditions();

  const std::vector<double> tensor_kappa =
      parser_m["HeatSolver/kappa"].get<std::vector<double>>();
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Creating new conduction operator with kappa.");
  operator_m = new ConductionOperator(*mesh_m, boundary_conditions_m,
                                      *coarse_element_space_m, *element_space_m,
                                      temperature_m, tensor_kappa);
}

void HeatSolver::RegisterFieldsForIO(void) {
  if (parser_m["Simulation/use_visit"].get<bool>()) {
    SPDLOG_LOGGER_INFO(MAIN_LOG, "Registering fields for export via VisIt");
    visit_collection_m = new MfemVisItCollection(
        mesh_m->GetMPIComm(), "HeatSolver", mesh_m->GetMfemMesh());
    visit_collection_m->RegisterField("temperature", element_space_m);
    SPDLOG_LOGGER_INFO(MAIN_LOG, "All fields registered for export");
  }

  // Below is ADIOS2 io calls. Will replace VisIt once working.
  if (this->FileWritingEnabled()) {
    file_io_m->AddVariableForGridFunction("HeatSolver/temperature",
                                          *element_space_m, true);
  }
}

void HeatSolver::SetInitialConditions(void) {
  mfem::ParGridFunction temperature_grid_function(element_space_m);
  DEBUG_ASSERT(this->RestartFileActive(), global_assert{}, DebugLevel::CHEAP{},
               "Cannot read from restart file.");
  temperature_grid_function = 0.0;
  file_io_m->GetImmediateBlock(
      "HeatSolver/temperature", {0},
      {static_cast<std::size_t>(temperature_grid_function.Size())},
      temperature_grid_function.GetData());
  temperature_grid_function.GetTrueDofs(temperature_m);
  SPDLOG_LOGGER_INFO(MAIN_LOG, "temperature field allocated and set");
}

HeatSolver::~HeatSolver(void) {
  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Destructing HeatSolver and freeing associated memory");
  delete visit_collection_m;
  visit_collection_m = nullptr;
  delete ode_solver_m;
  ode_solver_m = nullptr;
  delete element_collection_m;
  element_collection_m = nullptr;
  delete element_space_m;
  element_space_m = nullptr;
  delete coarse_element_collection_m;
  coarse_element_collection_m = nullptr;
  delete coarse_element_space_m;
  coarse_element_space_m = nullptr;
  delete operator_m;
  operator_m = nullptr;
  SPDLOG_LOGGER_INFO(MAIN_LOG, "HeatSolver successfully destructed");
}

void HeatSolver::InitializeBoundaryConditions(void) {
  for (auto& manager : boundary_conditions_m) {
    manager.second.Initialize(*mesh_m);
  }
}

bool HeatSolver::FileWritingEnabled(void) const {
  return file_io_m != nullptr && file_io_m->IsWriteModeActive();
}

bool HeatSolver::RestartFileActive(void) const {
  return file_io_m != nullptr && file_io_m->IsReadModeActive();
}

}  // namespace chyps
