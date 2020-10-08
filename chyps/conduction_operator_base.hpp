// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef CHYPS_CONDUCTION_OPERATOR_BASE_HPP_
#define CHYPS_CONDUCTION_OPERATOR_BASE_HPP_

#include <mfem/mfem.hpp>

namespace chyps {
class ConductionOperatorBase : public mfem::TimeDependentOperator {
 protected:
  mfem::ParFiniteElementSpace& fespace;
  mfem::Array<int>
      ess_tdof_list;  // this list remains empty for pure Neumann b.c.

 public:
  ConductionOperatorBase(mfem::ParFiniteElementSpace& f);

  virtual void BuildStaticOperators(void) = 0;

  virtual void Mult(const mfem::Vector& u, mfem::Vector& du_dt) const = 0;
  /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
      This is the only requirement for high-order SDIRK implicit integration.*/
  virtual void ImplicitSolve(const double dt, const mfem::Vector& u,
                             mfem::Vector& k) = 0;

  /// Update the diffusion BilinearForm K using the given true-dof vector
  /// `u`.
  virtual void SetParameters(const mfem::Vector& u) = 0;

  virtual const mfem::Vector& GetThermalCoefficient(void) const = 0;

  void AddToEssentialDOF(const mfem::Array<int>& a_list);

  void ClearEssentialDOF(void);

  virtual ~ConductionOperatorBase(void) = default;
};
}  // namespace chyps

#endif  // CHYPS_CONDUCTION_OPERATOR_BASE_H_
