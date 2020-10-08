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

#include <iostream>

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
      M_prec(nullptr),
      T_solver(f.GetComm()),
      T_prec(nullptr),
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

  T_solver.iterative_mode = false;
  T_solver.SetRelTol(rel_tol);
  T_solver.SetAbsTol(0.0);
  T_solver.SetMaxIter(100);
  T_solver.SetPrintLevel(0);
}

void ConductionOperator::BuildStaticOperators(void) {
  M = new mfem::ParBilinearForm(&fespace);
  M->SetAssemblyLevel(mfem::AssemblyLevel::PARTIAL);
  M->AddDomainIntegrator(new mfem::MassIntegrator());
  M->Assemble();  // keep sparsity pattern of M and K the same
  M->FormSystemMatrix(ess_tdof_list, M_op);

  M_prec = new mfem::OperatorJacobiSmoother(*M, ess_tdof_list);
  M_solver.SetPreconditioner(*M_prec);
  M_solver.SetOperator(*M_op);
}

void ConductionOperator::Mult(const mfem::Vector& u,
                              mfem::Vector& du_dt) const {
  // Compute:
  //    du_dt = M^{-1}*-K(u)
  // for du_dt
  K_op->Mult(u, z);
  z.Neg();  // z = -z

  mfem::ParGridFunction tmp_u(&fespace);
  tmp_u.SetFromTrueDofs(u);

  mfem::ParGridFunction tmp_z(&fespace);
  tmp_z.SetFromTrueDofs(z);
  mfem::GridFunctionCoefficient z_coeff(&tmp_z);

  K->FormLinearSystem(ess_tdof_list, tmp_u, tmp_z, A, X, B);
  M_solver.Mult(B, du_dt);
  // Force no change in temperature along Dirichlet boundary
  du_dt.SetSubVector(ess_tdof_list, 0.0);
}

void ConductionOperator::ImplicitSolve(const double dt, const mfem::Vector& u,
                                       mfem::Vector& du_dt) {
  // Solve the equation:
  //    du_dt = M^{-1}*[-K(u + dt*du_dt)]
  // for du_dt
  if (T == nullptr) {
    T = new mfem::ParBilinearForm(&fespace);
    T->SetAssemblyLevel(mfem::AssemblyLevel::PARTIAL);
    T->AddDomainIntegrator(new mfem::MassIntegrator());

    mfem::ConstantCoefficient dt_coeff(dt);
    mfem::GridFunctionCoefficient k_coeff(&thermal_coefficient_m);
    mfem::ProductCoefficient prod(dt_coeff, k_coeff);
    T->AddDomainIntegrator(new mfem::DiffusionIntegrator(prod));
    T->Assemble();
    T->FormSystemMatrix(ess_tdof_list, T_op);
    if (UsesTensorBasis(fespace)) {
      delete T_prec;
      T_prec = new mfem::OperatorJacobiSmoother(*T, ess_tdof_list);
      T_solver.SetPreconditioner(*T_prec);
    }
    T_solver.SetOperator(*T_op);
    current_dt = dt;
  }
  MFEM_VERIFY(dt == current_dt, "");  // SDIRK methods use the same dt

  K_op->Mult(u, z);
  z.Neg();

  mfem::ParGridFunction tmp_u(&fespace);
  tmp_u.SetFromTrueDofs(u);

  mfem::ParGridFunction tmp_z(&fespace);
  tmp_z.Distribute(z);
  mfem::GridFunctionCoefficient z_coeff(&tmp_z);

  K->FormLinearSystem(ess_tdof_list, tmp_u, tmp_z, A, X, B);
  T_solver.Mult(B, du_dt);
  du_dt.SetSubVector(ess_tdof_list, 0.0);
}

void ConductionOperator::SetParameters(const mfem::Vector& u) {
  thermal_coefficient_m.SetFromTrueDofs(u);
  for (int i = 0; i < thermal_coefficient_m.Size(); ++i) {
    thermal_coefficient_m(i) = kappa + alpha * thermal_coefficient_m(i);
  }

  delete K;
  K = new mfem::ParBilinearForm(&fespace);
  K->SetAssemblyLevel(mfem::AssemblyLevel::PARTIAL);
  mfem::GridFunctionCoefficient grid_coeff(&thermal_coefficient_m);
  K->AddDomainIntegrator(new mfem::DiffusionIntegrator(grid_coeff));
  K->Assemble();  // keep sparsity pattern of M and K the same
  K->FormSystemMatrix(ess_tdof_list, K_op);

  delete T;
  T = nullptr;  // Delete here to be reset on first entrance to SolveImplicit
}

const mfem::ParGridFunction& ConductionOperator::GetThermalCoefficient(
    void) const {
  return thermal_coefficient_m;
}

ConductionOperator::~ConductionOperator(void) {
  delete M;
  delete K;
  delete T;
  delete T_prec;
  delete M_prec;
}

}  // namespace chyps
