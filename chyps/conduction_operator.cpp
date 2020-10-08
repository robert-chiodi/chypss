// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chyps/conduction_operator.hpp"

namespace chyps {

ConductionOperator::ConductionOperator(mfem::ParFiniteElementSpace& f,
                                       const double a_alpha,
                                       const double a_kappa)
    : ConductionOperatorBase(f),
      M(NULL),
      K(NULL),
      T(NULL),
      current_dt(0.0),
      M_solver(f.GetComm()),
      T_solver(f.GetComm()),
      alpha(a_alpha),
      kappa(a_kappa),
      z(height),
      thermal_coefficient_m(&f) {
  const double rel_tol = 1e-8;
  M_solver.iterative_mode = false;
  M_solver.SetRelTol(rel_tol);
  M_solver.SetAbsTol(0.0);
  M_solver.SetMaxIter(100);
  M_solver.SetPrintLevel(0);
  M_prec.SetType(mfem::HypreSmoother::Jacobi);
  M_solver.SetPreconditioner(M_prec);

  T_solver.iterative_mode = false;
  T_solver.SetRelTol(rel_tol);
  T_solver.SetAbsTol(0.0);
  T_solver.SetMaxIter(100);
  T_solver.SetPrintLevel(0);
  T_solver.SetPreconditioner(T_prec);
}

void ConductionOperator::BuildStaticOperators(void) {
  M = new mfem::ParBilinearForm(&fespace);
  M->AddDomainIntegrator(new mfem::MassIntegrator());
  M->Assemble(0);  // keep sparsity pattern of M and K the same
  M->FormSystemMatrix(ess_tdof_list, Mmat);

  M_solver.SetOperator(Mmat);
}

void ConductionOperator::Mult(const mfem::Vector& u,
                              mfem::Vector& du_dt) const {
  // Compute:
  //    du_dt = M^{-1}*-K(u)
  // for du_dt
  Kmat.Mult(u, z);
  z.Neg();  // z = -z
  M_solver.Mult(z, du_dt);
}

void ConductionOperator::ImplicitSolve(const double dt, const mfem::Vector& u,
                                       mfem::Vector& du_dt) {
  // Solve the equation:
  //    du_dt = M^{-1}*[-K(u + dt*du_dt)]
  // for du_dt
  if (!T) {
    T = Add(1.0, Mmat, dt, Kmat);
    current_dt = dt;
    T_solver.SetOperator(*T);
  }
  MFEM_VERIFY(dt == current_dt, "");  // SDIRK methods use the same dt
  Kmat.Mult(u, z);
  z.Neg();
  T_solver.Mult(z, du_dt);
}

void ConductionOperator::SetParameters(const mfem::Vector& u) {
  thermal_coefficient_m.SetFromTrueDofs(u);
  for (int i = 0; i < thermal_coefficient_m.Size(); ++i) {
    thermal_coefficient_m(i) = kappa + alpha * thermal_coefficient_m(i);
  }

  delete K;
  K = new mfem::ParBilinearForm(&fespace);

  mfem::GridFunctionCoefficient u_coeff(&thermal_coefficient_m);

  K->AddDomainIntegrator(new mfem::DiffusionIntegrator(u_coeff));
  K->Assemble(0);  // keep sparsity pattern of M and K the same
  K->FormSystemMatrix(ess_tdof_list, Kmat);
  delete T;
  T = NULL;  // re-compute T on the next ImplicitSolve
}

const mfem::Vector& ConductionOperator::GetThermalCoefficient(void) const {
  return thermal_coefficient_m;
}

ConductionOperator::~ConductionOperator() {
  delete T;
  delete M;
  delete K;
}

}  // namespace chyps
