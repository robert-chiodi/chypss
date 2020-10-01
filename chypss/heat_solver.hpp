// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef CHYPSS_HEAT_SOLVER_H_
#define CHYPSS_HEAT_SOLVER_H_

namespace chyps {

class HeatSolver {
 public:
  HeatSolver(void) = default;

 private:
};

class ConductionOperator : public TimeDependentOperator {
 protected:
  ParFiniteElementSpace& fespace;
  Array<int> ess_tdof_list;  // this list remains empty for pure Neumann b.c.

  ParBilinearForm* M;
  ParBilinearForm* K;

  HypreParMatrix Mmat;
  HypreParMatrix Kmat;
  HypreParMatrix* T;  // T = M + dt K
  double current_dt;

  CGSolver M_solver;     // Krylov solver for inverting the mass matrix M
  HypreSmoother M_prec;  // Preconditioner for the mass matrix M

  CGSolver T_solver;     // Implicit solver for T = M + dt K
  HypreSmoother T_prec;  // Preconditioner for the implicit solver

  double alpha, kappa;

  mutable Vector z;  // auxiliary vector

 public:
  ConductionOperator(ParFiniteElementSpace& f, double alpha, double kappa,
                     const Vector& u);

  virtual void Mult(const Vector& u, Vector& du_dt) const;
  /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
      This is the only requirement for high-order SDIRK implicit integration.*/
  virtual void ImplicitSolve(const double dt, const Vector& u, Vector& k);

  /// Update the diffusion BilinearForm K using the given true-dof vector `u`.
  void SetParameters(const Vector& u);

  virtual ~ConductionOperator();
};

}  // namespace chyps

#endif  // CHYPSS_HEAT_SOLVER_H_
