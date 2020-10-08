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

#include <mpi.h>
#include <mfem/mfem.hpp>

#include "chyps/boundary_condition.hpp"
#include "chyps/conduction_operator.hpp"
#include "chyps/input_parser.hpp"
#include "chyps/mfem_visit_collection.hpp"
#include "chyps/solver_interface.hpp"

namespace chyps {

/// \class HeatSolver heat_solver.hpp chyps/heat_solver.hpp
/// \brief Executes and controls solution of the heat equation.
///
/// This class uses MFEM to formulate and solve for an FEM solution of the
/// heat equation.
class HeatSolver : public SolverInterface {
 public:
  HeatSolver(void) = delete;

  /// \brief Default initialize HeatSolver and collect all options passed to
  /// parser.
  HeatSolver(const MPI_Comm& a_mpi_comm, InputParser& a_parser);

  /// \brief Initialize solver, including the construction of operators and MFEM
  /// objects.
  virtual void Initialize(void) override final;

  /// \brief Compute time step the solver would like to take. Will be either
  /// a_proposed_dt or smaller.
  virtual double AdjustTimeStep(
      const double a_proposed_dt) const override final;

  /// \brief Advance the solver by a dt (or less) amount of time.
  virtual double Advance(const double a_time, const double dt) override final;

  /// \brief Sets boundary condition for boundary elements tagged with a_tag.
  ///
  /// These boundary conditions are stored for calling again later.
  void SetBoundaryCondition(const int a_tag,
                            const BoundaryCondition& a_condition);

  /// \brief Applies boundary conditions set by calls to SetBoundaryCondition.
  ///
  /// Note: Any unset boundary condition will be assumed to be a homogeneous
  /// Dirichlet condition. This method MUST BE CALLED before calls to
  /// AdjustTimeStep or Advance. This is true even if the default boundary
  /// conditions are desired.
  void CommitBoundaryConditions(void);

  /// \brief Reapplies the boundary conditions marked as time varying.
  /// Returns the number of boundary conditions updated.
  int UpdateBoundaryConditions(void);

  /// \brief Write out data to disk for visualization in VisIt.
  virtual void ExportVisIt(const int a_cycle,
                           const double a_time) override final;

  /// \brief Destructor that free all heap allocated variables.
  ~HeatSolver(void);

 private:
  // Confirms all runtime options have been specified or have default value;
  bool AllOptionsSupplied(void) const;

  void GatherOptions(void);

  void ReadAndRefineMesh(void);

  void SetODESolver(void);

  void AllocateVariablesAndOperators(void);

  void SetInitialConditions(void);

  void RegisterVisItFields(void);

  InputParser& parser_m;
  const MPI_Comm& mpi_comm_m;
  MfemVisItCollection* visit_collection_m;
  int dimension_m;
  mfem::ParMesh* parallel_mesh_m;
  mfem::ODESolver* ode_solver_m;
  mfem::FiniteElementCollection* element_collection_m;
  mfem::ParFiniteElementSpace* element_space_m;
  ConductionOperatorBase* operator_m;
  mfem::Vector temperature_m;
  std::vector<BoundaryCondition> boundary_conditions_m;
};

}  // namespace chyps

#endif  // CHYPS_HEAT_SOLVER_HPP_
