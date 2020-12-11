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
#include <tuple>

#include "chyps/debug_assert.hpp"
#include "chyps/logger.hpp"
#include "chyps/mesh.hpp"
#include "chyps/simulation.hpp"
#include "chyps/string_manipulation.hpp"

namespace chyps {

HeatSolver::HeatSolver(InputParser& a_parser, Simulation& a_simulation)
    : parser_m(a_parser),
      sim_m(a_simulation),
      boundary_conditions_m(),
      ode_solver_m(nullptr),
      element_collection_m(nullptr),
      element_space_m(nullptr),
      coarse_element_collection_m(nullptr),
      coarse_element_space_m(nullptr),
      operator_m(nullptr),
      temperature_m(),
      rho_m(),
      cp_m(),
      kappa_m(nullptr),
      solution_monitor_m(nullptr),
      monitor_file_active_m(false),
      visit_collection_m(nullptr),
      precice_write_names_m() {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Constructing HeatSolver object.");
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Gathering options of HeatSolver class.");
  this->GatherOptions();
  this->CreateBoundaryConditionManagers();
}

void HeatSolver::Initialize(void) {
  // Construct mesh, allocate and construct operators, and perform all setup for
  // time advancement
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Initializing HeatSolver class.");
  this->CreateMonitorFile();
  this->InitializeBoundaryConditions();
  this->AllocateVariablesAndOperators();
  this->SetODESolver();
  this->RegisterFieldsForIO();
  this->PushDataToMonitorFile();
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
  this->UpdateBoundaryConditions(temperature_m);
  operator_m->SetParameters(temperature_m);
  ode_solver_m->Step(temperature_m, time, taken_time_step);
  this->PushDataToMonitorFile();
  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Advanced solution from time {:8.6E} to time {:8.6E}",
                     a_time, time);
  return taken_time_step;
}

void HeatSolver::WriteDataToPrecice(void) {
  DEBUG_ASSERT(sim_m.PreciceActive(), global_assert{}, DebugLevel::CHEAP{},
               "preCICE is not active in simulation.");

  mfem::ParGridFunction tmp_fine_gf(element_space_m);
  tmp_fine_gf.SetFromTrueDofs(temperature_m);
  mfem::GridFunctionCoefficient map_coeff(&tmp_fine_gf);
  mfem::ParGridFunction tmp_coarse_gf(coarse_element_space_m);

  mfem::Array<int> border(sim_m.GetMesh().GetNumberOfBoundaryTagValues());
  mfem::Vector values;

  std::pair<const std::vector<double>*, const std::vector<int>*> bdr_info_pair;

  const BoundaryConditionManager& temperature_bc =
      boundary_conditions_m.at("temperature");
  const int size = static_cast<int>(precice_write_names_m.size());
  for (int n = 0; n < size; ++n) {
    if (precice_write_names_m[n] == "") {
      // Not a precice condition
      continue;
    }
    // Get mesh information
    bdr_info_pair = temperature_bc.GetBoundaryVertices(n + 1);
    auto vertex_indices = bdr_info_pair.second;

    // Wrap data in mfem objects.
    // CASTING AWAY CONST TO WRAP IN MFEM VALUES. WILL NOT MODIFY
    mfem::Array<int> indices(const_cast<int*>(vertex_indices->data()),
                             static_cast<int>(vertex_indices->size()));
    border = 0;
    border[n] = 1;
    tmp_coarse_gf.ProjectBdrCoefficient(map_coeff, border);
    tmp_coarse_gf.GetSubVector(indices, values);
    sim_m.GetPreciceAdapter().WriteBlockScalarData(precice_write_names_m[n],
                                                   values.GetData());
  }
}

void HeatSolver::WriteFields(const int a_cycle, const double a_time) {
  if (parser_m["Simulation/use_visit"].get<bool>()) {
    visit_collection_m->UpdateField("temperature", temperature_m);
    visit_collection_m->UpdateField("rho", rho_m);
    visit_collection_m->UpdateField("cp", cp_m);
    visit_collection_m->WriteOutFields(a_cycle, a_time);
  }

  if (this->FileWritingEnabled()) {
    DEBUG_ASSERT(sim_m.GetIO().IsWriteModeActive(), global_assert{},
                 DebugLevel::CHEAP{},
                 "IO must be in write mode for writing to fields.");
    DEBUG_ASSERT(sim_m.GetIO().OngoingWriteStep(), global_assert{},
                 DebugLevel::CHEAP{},
                 "An ongoing IO step is required for writing.");
    sim_m.GetIO().PutDeferred("HeatSolver/temperature", temperature_m);
    sim_m.GetIO().PutDeferred("HeatSolver/rho", rho_m);
    sim_m.GetIO().PutDeferred("HeatSolver/cp", cp_m);
    if (kappa_m->CanTimeVary()) {
      if (kappa_m->IsScalarCoefficient()) {
        sim_m.GetIO().PutDeferred(
            "HeatSolver/Conductivity/kappa",
            kappa_m->GetElementVaryingScalarCoefficient().GetData().data());
      } else if (kappa_m->IsMatrixCoefficient()) {
        sim_m.GetIO().PutDeferred(
            "HeatSolver/Conductivity/kappa",
            kappa_m->GetElementVaryingMatrixCoefficient().GetData());
      }
    }

    sim_m.GetIO().PerformPuts();
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
      "HeatSolver/PreciceWrite/temperature/",
      "Handles ability to write temperature on boundaries to preCICE for "
      "coupling. To add a boundary, include boundary as a key value with the "
      "precice name as the value. As an example:\n \"1\": "
      "\"temperature_name_in_precice\"");
  parser_m.AddOption(
      "HeatSolver/Conductivity/type",
      "Conductivity model to use. Types are: \"constant_scalar\", "
      "\"constant_matrix\", \"material_varying_scalar\", "
      "\"material_varying_matrix\", \"element_varying_scalar\", and "
      "\"element_varying_matrix\". See Documentation for or comments in "
      "conduction_operator.cpp for definition of each type. Types may require "
      "additional arguments.");

  // TODO Add a flag and way to specify which variables we wish to export to
  // VisIt.

  // Add options for conduction operator.
  ConductionOperator::GatherOptions(parser_m);

  SPDLOG_LOGGER_INFO(MAIN_LOG, "Add all options to parser.");
}

void HeatSolver::CreateBoundaryConditionManagers(void) {
  boundary_conditions_m.emplace(std::make_pair(
      "temperature",
      BoundaryConditionManager(parser_m,
                               "HeatSolver/BoundaryConditions/temperature")));
}

void HeatSolver::CreateMonitorFile(void) {
  if (sim_m.GetMPI().IAmRoot()) {
    solution_monitor_m = sim_m.GetMonitorManager().CreateMonitorFile(
        "heat_solver",
        {"MIN_T", "MAX_T", "MIN_RHO", "MAX_RHO", "MIN_CP", "MAX_CP"},
        {FieldType::DOUBLE, FieldType::DOUBLE, FieldType::DOUBLE,
         FieldType::DOUBLE, FieldType::DOUBLE, FieldType::DOUBLE});
  }
  // Inform all processors if monitor file is active.
  monitor_file_active_m = solution_monitor_m != nullptr;
  int bool_as_int = monitor_file_active_m ? 1 : 0;
  MPI_Bcast(&bool_as_int, 1, MPI_INT, sim_m.GetMPI().GetRootRank(),
            sim_m.GetMPI().GetComm());
  monitor_file_active_m = bool_as_int == 1;
}

void HeatSolver::PushDataToMonitorFile(void) {
  if (!this->MonitorFileActive()) {
    return;
  }

  // 0 - MIN_TEMPERATURE
  // 1 - MAX_TEMPERATURE
  // 2 - MIN_RHO
  // 3 - MAX_RHO
  // 4 - MIN_CP
  // 5 - MAX_CP
  std::array<double, 6> minmax_values;
  minmax_values[0] = temperature_m.Min();
  minmax_values[1] = rho_m.Min();
  minmax_values[2] = cp_m.Min();
  minmax_values[3] = temperature_m.Max();
  minmax_values[4] = rho_m.Max();
  minmax_values[5] = cp_m.Max();

  std::array<double, 6> global_minmax_values;
  MPI_Reduce(minmax_values.data(), global_minmax_values.data(), 3, MPI_DOUBLE,
             MPI_MIN, sim_m.GetMPI().GetRootRank(), sim_m.GetMPI().GetComm());

  MPI_Reduce(minmax_values.data() + 3, global_minmax_values.data() + 3, 3,
             MPI_DOUBLE, MPI_MAX, sim_m.GetMPI().GetRootRank(),
             sim_m.GetMPI().GetComm());

  if (sim_m.GetMPI().IAmRoot()) {
    solution_monitor_m->SetEntries(
        {global_minmax_values[0], global_minmax_values[3],
         global_minmax_values[1], global_minmax_values[4],
         global_minmax_values[2], global_minmax_values[5]});
  }
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
      parser_m["HeatSolver/order"].get<int>(), sim_m.GetMesh().GetDimension());

  SPDLOG_LOGGER_INFO(MAIN_LOG, "Generating element space from H1 collection");
  element_space_m = new mfem::ParFiniteElementSpace(
      &sim_m.GetMesh().GetMfemMesh(), element_collection_m);

  coarse_element_collection_m =
      new mfem::H1_FECollection(1, sim_m.GetMesh().GetDimension());

  SPDLOG_LOGGER_INFO(MAIN_LOG, "Generating element space from H1 collection");
  coarse_element_space_m = new mfem::ParFiniteElementSpace(
      &sim_m.GetMesh().GetMfemMesh(), coarse_element_collection_m);

  kappa_m = new Conductivity(sim_m, *element_space_m);

  if (parser_m.OptionSet("HeatSolver/PreciceWrite/temperature")) {
    DEBUG_ASSERT(sim_m.PreciceActive(), global_assert{}, DebugLevel::CHEAP{},
                 "Precice coupling must be active to set write-fields.");
    precice_write_names_m.resize(sim_m.GetMesh().GetNumberOfBoundaryTagValues(),
                                 "");
    // Fill with precice write names.
    const auto& write_group = parser_m["HeatSolver/PreciceWrite/temperature"];
    for (const auto& elem : write_group.items()) {
      const int index = std::stoi(elem.key());
      DEBUG_ASSERT(index > 0, global_assert{}, DebugLevel::CHEAP{},
                   "Tag value must be >0");
      DEBUG_ASSERT(
          static_cast<std::size_t>(index - 1) < precice_write_names_m.size(),
          global_assert{}, DebugLevel::CHEAP{},
          "Tag value must exist on mesh. Mesh expects all attribute values are "
          "contiguous and in range [1,number_of_attribute_tags]");
      precice_write_names_m[index - 1] = elem.value().get<std::string>();
      sim_m.GetPreciceAdapter().AddData(
          precice_write_names_m[index - 1],
          sim_m.GetMesh().GetPreciceMeshName(index), DataOperation::WRITE);
    }
  }

  // Allocates temperature and sets initial values.
  this->SetInitialConditions();

  SPDLOG_LOGGER_INFO(MAIN_LOG, "Creating new conduction operator.");
  operator_m = new ConductionOperator(parser_m, sim_m, boundary_conditions_m,
                                      *coarse_element_space_m, *element_space_m,
                                      temperature_m, rho_m, cp_m, *kappa_m);
}

void HeatSolver::RegisterFieldsForIO(void) {
  if (parser_m["Simulation/use_visit"].get<bool>()) {
    SPDLOG_LOGGER_INFO(MAIN_LOG, "Registering fields for export via VisIt");
    visit_collection_m = new MfemVisItCollection(
        sim_m.GetMPI().GetComm(), "HeatSolver", sim_m.GetMesh().GetMfemMesh());
    visit_collection_m->RegisterField("temperature", element_space_m);
    visit_collection_m->RegisterField("rho", element_space_m);
    visit_collection_m->RegisterField("cp", element_space_m);
    SPDLOG_LOGGER_INFO(MAIN_LOG, "All fields registered for export");
  }

  // Below is ADIOS2 io calls. Will replace VisIt once working.
  if (this->FileWritingEnabled()) {
    sim_m.GetIO().AddVariableForTrueDofs("HeatSolver/temperature",
                                         *element_space_m, true);
    sim_m.GetIO().AddVariableForTrueDofs("HeatSolver/rho", *element_space_m,
                                         true);
    sim_m.GetIO().AddVariableForTrueDofs("HeatSolver/cp", *element_space_m,
                                         true);
    if (kappa_m->CanTimeVary()) {
      if (kappa_m->IsScalarCoefficient()) {
        sim_m.GetIO().AddVariableForMesh<double>(
            "HeatSolver/Conductivity/kappa", sim_m.GetMesh(),
            MeshElement::ELEMENT);
      } else if (kappa_m->IsMatrixCoefficient()) {
        sim_m.GetIO().AddMatrixForMesh("HeatSolver/Conductivity/kappa",
                                       sim_m.GetMesh(), MeshElement::ELEMENT,
                                       sim_m.GetMesh().GetDimension(),
                                       sim_m.GetMesh().GetDimension());
      }
    }
  }
}

void HeatSolver::SetInitialConditions(void) {
  DEBUG_ASSERT(this->RestartFileActive(), global_assert{}, DebugLevel::CHEAP{},
               "Cannot read from restart file.");
  DEBUG_ASSERT(element_space_m->GetVDim() == 1, global_assert{},
               DebugLevel::CHEAP{},
               "Finite element space should be for a scalar component.");
  const int truedof_size = element_space_m->GetTrueVSize();

  temperature_m.SetSize(truedof_size);
  sim_m.GetIO().GetImmediateBlock("HeatSolver/temperature",
                                  {static_cast<std::size_t>(truedof_size)},
                                  temperature_m.GetData());

  rho_m.SetSize(truedof_size);
  sim_m.GetIO().GetImmediateBlock("HeatSolver/rho",
                                  {static_cast<std::size_t>(truedof_size)},
                                  rho_m.GetData());

  cp_m.SetSize(truedof_size);
  sim_m.GetIO().GetImmediateBlock("HeatSolver/cp",
                                  {static_cast<std::size_t>(truedof_size)},
                                  cp_m.GetData());

  if (kappa_m->CanTimeVary()) {
    if (kappa_m->IsScalarCoefficient()) {
      sim_m.GetIO().GetImmediateBlock(
          "HeatSolver/Conductivity/kappa",
          kappa_m->GetElementVaryingScalarCoefficient().GetData());
    } else if (kappa_m->IsMatrixCoefficient()) {
      sim_m.GetIO().GetImmediateBlock(
          "HeatSolver/Conductivity/kappa",
          kappa_m->GetElementVaryingMatrixCoefficient().GetData());
    }
  }
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Heat solver fields allocated and set");
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
  delete kappa_m;
  kappa_m = nullptr;
  delete operator_m;
  operator_m = nullptr;
  SPDLOG_LOGGER_INFO(MAIN_LOG, "HeatSolver successfully destructed");
}

void HeatSolver::InitializeBoundaryConditions(void) {
  for (auto& manager : boundary_conditions_m) {
    manager.second.Initialize(sim_m);
  }
  if (sim_m.PreciceActive()) {
    for (auto& manager : boundary_conditions_m) {
      for (int n = 0; n < manager.second.GetNumberOfBoundaryConditions(); ++n) {
        if (manager.second.IsPreciceBoundaryCondition(n + 1)) {
          sim_m.GetPreciceAdapter().AddData(
              manager.second.PreciceName(n + 1),
              sim_m.GetMesh().GetPreciceMeshName(n + 1), DataOperation::READ);
        }
      }
    }
  }
}

void HeatSolver::UpdateBoundaryConditions(mfem::Vector& a_temperature) {
  if (sim_m.PreciceActive()) {
    for (auto& manager : boundary_conditions_m) {
      for (int n = 0; n < manager.second.GetNumberOfBoundaryConditions(); ++n) {
        if (manager.second.IsPreciceBoundaryCondition(n + 1)) {
          sim_m.GetPreciceAdapter().ReadBlockScalarData(
              manager.second.PreciceName(n + 1),
              manager.second.GetDataBuffer(n + 1));
        }
      }
    }
  }

  operator_m->UpdateBoundaryConditions(temperature_m);
}

bool HeatSolver::FileWritingEnabled(void) const {
  return sim_m.GetIO().IsWriteModeActive();
}

bool HeatSolver::RestartFileActive(void) const {
  return sim_m.GetIO().IsReadModeActive();
}

bool HeatSolver::MonitorFileActive(void) const { return monitor_file_active_m; }

}  // namespace chyps
