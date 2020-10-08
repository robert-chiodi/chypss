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

#include <mfem/mfem.hpp>

#include "chyps/conduction_operator_base.hpp"

namespace chyps {
class ConductionOperator : public ConductionOperatorBase {
 public:
  ConductionOperator(mfem::ParFiniteElementSpace& f, const double a_alpha,
                     const double a_kappa);

  /// \brief Build operators that will not change over course of the simulation.
  virtual void BuildStaticOperators(void) override final;

  virtual void Mult(const mfem::Vector& u,
                    mfem::Vector& du_dt) const override final;
  /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
      This is the only requirement for high-order SDIRK implicit integration.*/
  virtual void ImplicitSolve(const double dt, const mfem::Vector& u,
                             mfem::Vector& k) override final;

  /// Update the diffusion BilinearForm K using the given true-dof vector `u`.
  void SetParameters(const mfem::Vector& u) override final;

  /// \brief Return a vector of the Thermal Coefficient for each element.
  virtual const mfem::ParGridFunction& GetThermalCoefficient(
      void) const override final;

  virtual ~ConductionOperator(void) override final;

 protected:
  mfem::ParBilinearForm* M;
  mfem::ParBilinearForm* K;

  mfem::ParBilinearForm* T;  // T = M + dt K
  double current_dt;

  mfem::CGSolver M_solver;  // Krylov solver for inverting the mass matrix M
  mfem::OperatorJacobiSmoother* M_prec;  // Preconditioner for the mass matrix M

  mfem::CGSolver T_solver;  // Implicit solver for T = M + dt K
  mfem::OperatorJacobiSmoother*
      T_prec;  // Preconditioner for the implicit solver

  mfem::OperatorPtr T_op, M_op, K_op;

  double alpha, kappa;

  mutable mfem::Vector z;  // auxiliary vector
  mutable mfem::OperatorHandle A;
  mutable mfem::Vector X, B;
  // mutable mfem::Vector tmp_u;

 private:
  mfem::ParGridFunction thermal_coefficient_m;
};
}  // namespace chyps

#endif  // CHYPS_CONDUCTION_OPERATOR_HPP_
