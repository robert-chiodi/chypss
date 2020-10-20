// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef CHYPS_CONDUCTION_OPERATOR_HPP_
#define CHYPS_CONDUCTION_OPERATOR_HPP_

#include <vector>

#include <mfem/mfem.hpp>

#include "chyps/conduction_operator_base.hpp"
#include "chyps/mesh.hpp"

namespace chyps {
class ConductionOperator : public ConductionOperatorBase {
 public:
  ConductionOperator(Mesh& a_mesh, mfem::ParFiniteElementSpace& f_linear,
                     mfem::ParFiniteElementSpace& f, mfem::Vector& u,
                     const double a_alpha, const double a_kappa);

  virtual void Mult(const mfem::Vector& u,
                    mfem::Vector& du_dt) const override final;
  /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
      This is the only requirement for high-order SDIRK implicit integration.*/
  virtual void ImplicitSolve(const double dt, const mfem::Vector& u,
                             mfem::Vector& k) override final;

  /// \brief Update the diffusion BilinearForm K using the given true-dof vector
  /// `u`.
  void SetParameters(const mfem::Vector& u) override final;

  /// \brief Update boundary conditions if they are time-varying. Might require
  /// changing boundaries in u.
  virtual void UpdateBoundaryConditions(mfem::Vector& u) override final;

  /// \brief Return a vector of the Thermal Coefficient for each element.
  virtual const mfem::ParGridFunction& GetThermalCoefficient(
      void) const override final;

  virtual ~ConductionOperator(void) override final;

 protected:
  mfem::ParFiniteElementSpace& fespace_linear_m;

  Mesh& mesh_m;
  mfem::ParBilinearForm* M;
  mfem::ParBilinearForm* K;
  mfem::ParBilinearForm* T;  // T = M + dt K
  mfem::ParLinearForm* neumann_m;

  double current_dt;

  mutable mfem::CGSolver
      M_solver;  // Krylov solver for inverting the mass matrix M
  mfem::OperatorJacobiSmoother* M_prec;  // Preconditioner for the mass matrix M

  mfem::CGSolver T_solver;  // Implicit solver for T = M + dt K
  mfem::OperatorJacobiSmoother*
      T_prec;  // Preconditioner for the implicit solver

  mfem::OperatorPtr T_op, M_op, K_op;

  double alpha, kappa;

  mutable mfem::Vector z;  // auxiliary vector
  mutable mfem::OperatorHandle A;
  mutable mfem::Vector B;

  std::vector<mfem::Coefficient*> neumann_coefficient_m;
  std::vector<mfem::Array<int>> boundary_marker_m;

 private:
  void SetTrueDofsFromVertexData(const std::size_t a_size,
                                 const int* a_vertex_list,
                                 const double* a_vertex_data,
                                 mfem::Array<int>& a_boundary,
                                 mfem::ParGridFunction& a_temperature_gf);

  void ResetNeumannCondition(void);

  void AddNeumannCondition(const int a_tag,
                           const mfem::ParGridFunction& a_grid_function);
  void AddNeumannCondition(const int a_tag, const double a_value);
  void FinalizeNeumannCondition(void);

  mfem::ParGridFunction thermal_coefficient_m;
  bool neumann_active_m;
};
}  // namespace chyps

#endif  // CHYPS_CONDUCTION_OPERATOR_HPP_
