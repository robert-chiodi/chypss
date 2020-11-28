// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// \file heat_solver.hpp
/// \brief Overarching solver that dictates behavior and controls execution
/// of the solution of a heat equation.

#ifndef CHYPS_HEAT_SOLVER_HPP_
#define CHYPS_HEAT_SOLVER_HPP_

#include <mfem/mfem.hpp>

#include "chyps/boundary_condition.hpp"
#include "chyps/boundary_condition_manager.hpp"
#include "chyps/conduction_operator.hpp"
#include "chyps/conductivity.hpp"
#include "chyps/input_parser.hpp"
#include "chyps/io.hpp"
#include "chyps/mesh.hpp"
#include "chyps/mfem_visit_collection.hpp"
#include "chyps/solver_interface.hpp"

namespace chyps {

// Forward declaration
class Simulation;

/// \class HeatSolver heat_solver.hpp chyps/heat_solver.hpp
/// \brief Executes and controls solution of the heat equation.
///
/// This class uses MFEM to formulate and solve for an FEM solution of the
/// heat equation.
class HeatSolver : public SolverInterface {
 public:
  HeatSolver(void) = delete;

  /// \brief Initialize HeatSolver and collect all options passed to
  /// parser.
  HeatSolver(InputParser& a_parser, Simulation& a_simulation);

  /// \brief Initialize solver, including the construction of operators and MFEM
  /// objects.
  virtual void Initialize(void) override final;

  /// \brief Compute time step the solver would like to take. Will be either
  /// a_proposed_dt or smaller.
  virtual double AdjustTimeStep(
      const double a_proposed_dt) const override final;

  /// \brief Advance the solver by a dt (or less) amount of time.
  virtual double Advance(const double a_time, const double dt) override final;

  /// \brief Write data (if needed to preCICE).
  void WriteDataToPrecice(void);

  /// \brief Write out data to disk for visualization in VisIt.
  virtual void WriteFields(const int a_cycle,
                           const double a_time) override final;

  /// \brief Destructor that free all heap allocated variables.
  ~HeatSolver(void);

 private:
  void GatherOptions(void);

  void CreateBoundaryConditionManagers(void);

  void InitializeBoundaryConditions(void);

  void AllocateVariablesAndOperators(void);

  void SetODESolver(void);

  void SetInitialConditions(void);

  void RegisterFieldsForIO(void);

  void UpdateBoundaryConditions(mfem::Vector& a_temperature);

  bool FileWritingEnabled(void) const;
  bool RestartFileActive(void) const;

  InputParser& parser_m;
  Simulation& sim_m;
  std::unordered_map<std::string, BoundaryConditionManager>
      boundary_conditions_m;
  mfem::ODESolver* ode_solver_m;
  mfem::FiniteElementCollection* element_collection_m;
  mfem::ParFiniteElementSpace* element_space_m;
  mfem::FiniteElementCollection* coarse_element_collection_m;
  mfem::ParFiniteElementSpace* coarse_element_space_m;
  ConductionOperator* operator_m;
  mfem::Vector temperature_m;
  mfem::Vector rho_m;
  mfem::Vector cp_m;
  Conductivity* kappa_m;
  MfemVisItCollection* visit_collection_m;
  std::vector<std::string> precice_write_names_m;
};

}  // namespace chyps

#endif  // CHYPS_HEAT_SOLVER_HPP_
