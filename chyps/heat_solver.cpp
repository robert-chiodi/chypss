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

namespace chyps {

HeatSolver::HeatSolver(const MPI_Comm& a_mpi_comm, InputParser& a_parser)
    : parser_m(a_parser),
      mpi_comm_m(a_mpi_comm),
      visit_collection_m(nullptr),
      dimension_m(-1),
      parallel_mesh_m(nullptr),
      ode_solver_m(nullptr),
      element_collection_m(nullptr),
      element_space_m(nullptr),
      operator_m(nullptr),
      temperature_m() {
  this->GatherOptions();
}

void HeatSolver::Initialize(void) {
  // FIXME : make this an exception.
  if (!this->AllOptionsSupplied()) {
    std::cout << "Not all options needed for HeatSolver supplied" << std::endl;
    std::cout << "Make sure that the InputParser has been parsed before "
                 "calling Initialize and that all required options are "
                 "specified or have a valid default value."
              << std::endl;
    std::exit(-1);
  }

  // Construct mesh, allocate and construct operators, and perform all setup for
  // time advancement
  this->ReadAndRefineMesh();
  this->AllocateVariablesAndOperators();
  this->SetInitialConditions();
  this->SetODESolver();
  this->RegisterVisItFields();
}

double HeatSolver::AdjustTimeStep(const double a_proposed_dt) const {
  return a_proposed_dt;
}

double HeatSolver::Advance(const double a_time, const double a_time_step) {
  double taken_time_step = a_time_step;
  double time = a_time;
  operator_m->SetParameters(temperature_m);
  ode_solver_m->Step(temperature_m, time, taken_time_step);
  return taken_time_step;
}

void HeatSolver::GatherOptions(void) {
  parser_m.AddOption("mesh_file", "-m", "--mesh", "Mesh file to use.",
                     std::string("../data/star.mesh"));
  parser_m.AddOption("serial_refine", "-rs", "--refine-serial",
                     "Number of times to refine the mesh uniformly in serial.",
                     2);
  parser_m.AddOption(
      "parallel_refine", "-rp", "--refine-parallel",
      "Number of times to refine the mesh uniformly in parallel.", 1);
  parser_m.AddOption("order", "-o", "--order",
                     "Order (degree) of the finite elements.", 2);
  parser_m.AddOption(
      "time_advancement", "-s", "--ode-solver",
      "ODE solver: 1 - Backward Euler, 2 - SDIRK2, 3 - SDIRK3,\n\t"
      "\t   11 - Forward Euler, 12 - RK2, 13 - RK3 SSP, 14 - RK4.",
      3);
  parser_m.AddOption("alpha", "-a", "--alpha", "Alpha coefficient.", 1.0e-2);
  parser_m.AddOption("kappa", "-k", "--kappa", "Kappa coefficient offset.",
                     0.5);
  // TODO Add a flag and way to specify which variables we wish to export to
  // VisIt.
}

bool HeatSolver::AllOptionsSupplied(void) const {
  return parser_m.AllOptionsSet();
}

void HeatSolver::ReadAndRefineMesh(void) {
  const std::string file_name = parser_m["mesh_file"];
  mfem::Mesh mesh(file_name.c_str(), 1, 1);
  dimension_m = mesh.Dimension();
  const auto serial_refinement = static_cast<int>(parser_m["serial_refine"]);
  const auto parallel_refinement =
      static_cast<int>(parser_m["parallel_refine"]);
  for (int lev = 0; lev < serial_refinement; ++lev) {
    mesh.UniformRefinement();
  }
  parallel_mesh_m = new mfem::ParMesh(mpi_comm_m, mesh);
  for (int lev = 0; lev < parallel_refinement; ++lev) {
    parallel_mesh_m->UniformRefinement();
  }
}

void HeatSolver::SetODESolver(void) {
  auto ode_solver_type = static_cast<int>(parser_m["time_advancement"]);
  switch (ode_solver_type) {
    // Implicit L-stable methods
    case 1:
      ode_solver_m = new mfem::BackwardEulerSolver;
      break;
    case 2:
      ode_solver_m = new mfem::SDIRK23Solver(2);
      break;
    case 3:
      ode_solver_m = new mfem::SDIRK33Solver;
      break;
    // Explicit methods
    case 11:
      ode_solver_m = new mfem::ForwardEulerSolver;
      break;
    case 12:
      ode_solver_m = new mfem::RK2Solver(0.5);
      break;  // midpoint method
    case 13:
      ode_solver_m = new mfem::RK3SSPSolver;
      break;
    case 14:
      ode_solver_m = new mfem::RK4Solver;
      break;
    case 15:
      ode_solver_m = new mfem::GeneralizedAlphaSolver(0.5);
      break;
    // Implicit A-stable methods (not L-stable)
    case 22:
      ode_solver_m = new mfem::ImplicitMidpointSolver;
      break;
    case 23:
      ode_solver_m = new mfem::SDIRK23Solver;
      break;
    case 24:
      ode_solver_m = new mfem::SDIRK34Solver;
      break;
    default:
      // FIXME : Make actual error handler
      std::cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
      std::exit(-1);
  }
  ode_solver_m->Init(*operator_m);
}

void HeatSolver::AllocateVariablesAndOperators(void) {
  element_collection_m =
      new mfem::H1_FECollection(parser_m["order"], dimension_m);
  element_space_m =
      new mfem::ParFiniteElementSpace(parallel_mesh_m, element_collection_m);

  // FIXME : Check this is the right size
  temperature_m.SetSize(element_space_m->GetTrueVSize());

  operator_m = new ConductionOperator(*element_space_m, parser_m["alpha"],
                                      parser_m["kappa"]);
}

void HeatSolver::RegisterVisItFields(void) {
  if (!parser_m["use_visit"]) {
    return;
  }
  visit_collection_m =
      new MfemVisItCollection(mpi_comm_m, "HeatSolver", *parallel_mesh_m);
  visit_collection_m->RegisterField("temperature", element_space_m);
  visit_collection_m->RegisterField("ThermalCoefficient", element_space_m);
}

static double SetInitialTemperature(const mfem::Vector& x) {
  if (x.Norml2() < 0.5) {
    return 2.0;
  } else {
    return 1.0;
  }
}

void HeatSolver::SetInitialConditions(void) {
  mfem::FunctionCoefficient initial_temperature(SetInitialTemperature);
  mfem::ParGridFunction temperature_grid_function(element_space_m);
  temperature_grid_function.ProjectCoefficient(initial_temperature);
  temperature_grid_function.GetTrueDofs(temperature_m);
}

void HeatSolver::ExportVisIt(const int a_cycle, const double a_time) {
  if (!parser_m["use_visit"]) {
    return;
  }
  visit_collection_m->UpdateField("temperature", temperature_m);
  visit_collection_m->UpdateField("ThermalCoefficient",
                                  operator_m->GetThermalCoefficient());
  visit_collection_m->WriteOutFields(a_cycle, a_time);
}

HeatSolver::~HeatSolver(void) {
  delete visit_collection_m;
  visit_collection_m = nullptr;
  delete parallel_mesh_m;
  parallel_mesh_m = nullptr;
  delete ode_solver_m;
  ode_solver_m = nullptr;
  delete element_collection_m;
  element_collection_m = nullptr;
  delete element_space_m;
  element_space_m = nullptr;
  delete operator_m;
  operator_m = nullptr;
}

}  // namespace chyps
