// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef CHYPSS_CONDUCTION_OPERATOR_H_
#define CHYPSS_CONDUCTION_OPERATOR_H_

#include <mfem.hpp>

namespace chyps {
class ConductionOperator : public TimeDependentOperator {
 protected:
  mfem::ParFiniteElementSpace& fespace;
  mfem::Array<int>
      ess_tdof_list;  // this list remains empty for pure Neumann b.c.

  mfem::ParBilinearForm* M;
  mfem::ParBilinearForm* K;

  mfem::HypreParMatrix Mmat;
  mfem::HypreParMatrix Kmat;
  mfem::HypreParMatrix* T;  // T = M + dt K
  double current_dt;

  mfem::CGSolver M_solver;     // Krylov solver for inverting the mass matrix M
  mfem::HypreSmoother M_prec;  // Preconditioner for the mass matrix M

  mfem::CGSolver T_solver;     // Implicit solver for T = M + dt K
  mfem::HypreSmoother T_prec;  // Preconditioner for the implicit solver

  double alpha, kappa;

  mutable mfem::Vector z;  // auxiliary vector

 public:
  ConductionOperator(mfem::ParFiniteElementSpace& f, double alpha, double kappa,
                     const mfem::Vector& u);

  virtual void Mult(const mfem::Vector& u, mfem::Vector& du_dt) const;
  /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
      This is the only requirement for high-order SDIRK implicit integration.*/
  virtual void ImplicitSolve(const double dt, const mfem::Vector& u,
                             mfem::Vector& k);

  /// Update the diffusion BilinearForm K using the given true-dof vector `u`.
  void SetParameters(const mfem::Vector& u);

  virtual ~ConductionOperator();
};
}  // namespace chyps

#endif  // CHYPSS_CONDUCTION_OPERATOR_H_
