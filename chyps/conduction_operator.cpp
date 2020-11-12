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

#include "chyps/debug_assert.hpp"
#include "chyps/logger.hpp"

namespace chyps {

ConductionOperator::ConductionOperator(
    Mesh& a_mesh,
    const std::unordered_map<std::string, BoundaryConditionManager>&
        a_boundary_conditions,
    mfem::ParFiniteElementSpace& f_linear, mfem::ParFiniteElementSpace& f,
    mfem::Vector& u, const std::vector<double>& a_tensor_kappa)
    : ConductionOperatorBase(f),
      fespace_linear_m(f_linear),
      mesh_m(a_mesh),
      boundary_conditions_m(a_boundary_conditions),
      M(NULL),
      K(NULL),
      T(NULL),
      neumann_m(nullptr),
      current_dt(0.0),
      M_solver(f.GetComm()),
      M_prec(nullptr),
      T_solver(f.GetComm()),
      T_prec(nullptr),
      tensor_kappa_m(a_tensor_kappa),
      z(height),  // Note, height inherited from mfem::TimeDependentOperator
      neumann_coefficient_m(a_mesh.GetMfemMesh().bdr_attributes.Max(), nullptr),
      boundary_marker_m(a_mesh.GetMfemMesh().bdr_attributes.Max()),
      tensor_thermal_coeff_m(nullptr),
      dt_tensor_thermal_coeff_m(nullptr),
      inhomogeneous_neumann_active_m(false) {
  DEBUG_ASSERT(
      tensor_kappa_m.size() == static_cast<std::size_t>(mesh_m.GetDimension() *
                                                        mesh_m.GetDimension()),
      global_assert{}, DebugLevel::CHEAP{},
      "Thermal coefficient tensor (kappa) of incorrect size. Provide as a "
      "column-major array of size MESH_DIM*MESH_DIM.");

  const double rel_tol = 1e-12;
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

  // Preset boundary arrays
  for (int n = 0; n < static_cast<int>(boundary_marker_m.size()); ++n) {
    boundary_marker_m[n].SetSize(boundary_marker_m.size());
    boundary_marker_m[n] = 0;
    boundary_marker_m[n][n] = 1;
  }

  // Set initial boundary conditions.
  mfem::Array<int> tdof_list;

  DEBUG_ASSERT(
      boundary_conditions_m.find("temperature") != boundary_conditions_m.end(),
      global_assert{}, DebugLevel::CHEAP{},
      "\"temperature\" not found in supplied boundary conditions.");
  const auto& bc_manager = boundary_conditions_m.at("temperature");
  if (bc_manager.GetNumberOfNeumannConditions() > 0) {
    SPDLOG_LOGGER_INFO(
        MAIN_LOG,
        "Neumann boundary conditions exist. Turning on Neumann ParLinearForm");
    inhomogeneous_neumann_active_m = true;
    this->ResetNeumannCondition();
  }

  const int number_of_boundary_conditions =
      bc_manager.GetNumberOfBoundaryConditions();
  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Committing {} boundary conditions in ConductionOperator",
                     number_of_boundary_conditions);
  mfem::Vector boundary_temp;
  mfem::ParGridFunction temperature_gf(&fespace);
  for (int n = 0; n < number_of_boundary_conditions; ++n) {
    const auto& condition = bc_manager.GetBoundaryCondition(n + 1);

    mfem::Array<int>& boundary = boundary_marker_m[n];
    fespace.GetEssentialTrueDofs(boundary, tdof_list);

    SPDLOG_LOGGER_INFO(
        MAIN_LOG,
        "Committing boundary condition (type:{}) for tag {}. Consists of "
        "{} border elements",
        static_cast<int>(condition.GetBCType()), n + 1, tdof_list.Size());

    switch (condition.GetBCType()) {
      case BoundaryConditionType::HOMOGENEOUS_DIRICHLET: {
        ess_tdof_list.Append(tdof_list);
        u.SetSubVector(tdof_list, 0.0);
        break;
      }
      case BoundaryConditionType::DIRICHLET: {
        ess_tdof_list.Append(tdof_list);
        if (condition.IsSpatiallyVarying()) {
          const auto& values = condition.GetValues();
          const auto& indices = condition.GetIndices();
          this->SetTrueDofsFromVertexData(indices.Size(), indices.Data(),
                                          values.Data(), boundary,
                                          temperature_gf);
          temperature_gf.GetTrueDofs(boundary_temp);
          for (int i = 0; i < tdof_list.Size(); ++i) {
            const int j = tdof_list[i];
            u(j) = boundary_temp(j);
          }
        } else {
          const auto& values = condition.GetValues();
          u.SetSubVector(tdof_list, *values.Data());
        }
        break;
      }
      case BoundaryConditionType::HOMOGENEOUS_NEUMANN: {
        // Nothing required, natural BC
        break;
      }
      case BoundaryConditionType::NEUMANN: {
        if (condition.IsSpatiallyVarying()) {
          const auto& values = condition.GetValues();
          const auto& indices = condition.GetIndices();
          this->SetTrueDofsFromVertexData(indices.Size(), indices.Data(),
                                          values.Data(), boundary,
                                          temperature_gf);
          this->AddNeumannCondition(n + 1, temperature_gf);
        } else {
          const auto& values = condition.GetValues();
          this->AddNeumannCondition(n + 1, *values.Data());
        }

        break;
      }
      default:
        DEBUG_ASSERT(
            false, global_assert{}, DebugLevel::ALWAYS{},
            "Unknown BoundaryConditionType. Set type is: " +
                std::to_string(static_cast<int>(condition.GetBCType())));
    }
  }
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Finished commiting boundary conditions");

  if (inhomogeneous_neumann_active_m) {
    this->FinalizeNeumannCondition();
  }

  SPDLOG_LOGGER_INFO(
      MAIN_LOG, "Forming system matrices with updated essential True DOF list.",
      number_of_boundary_conditions);
  M = new mfem::ParBilinearForm(&fespace);
  M->SetAssemblyLevel(mfem::AssemblyLevel::PARTIAL);
  M->AddDomainIntegrator(new mfem::MassIntegrator());
  M->Assemble();  // keep sparsity pattern of M and K the same
  M->FormSystemMatrix(ess_tdof_list, M_op);

  M_prec = new mfem::OperatorJacobiSmoother(*M, ess_tdof_list);
  M_solver.SetPreconditioner(*M_prec);
  M_solver.SetOperator(*M_op);

  SPDLOG_LOGGER_INFO(MAIN_LOG, "Finished constructing ConductionOperator");
}

void ConductionOperator::Mult(const mfem::Vector& u,
                              mfem::Vector& du_dt) const {
  // Compute:
  //    du_dt = M^{-1}*-K(u)
  // for du_dt
  mfem::ParGridFunction tmp_u(&fespace);
  mfem::ParLinearForm tmp_z(&fespace);
  tmp_u.SetFromTrueDofs(u);
  K->Mult(tmp_u, tmp_z);
  tmp_z.Neg();  // z = -z

  if (inhomogeneous_neumann_active_m) {
    tmp_z.Add(1.0, *neumann_m);
  }

  mfem::ParGridFunction tmp_du_dt(&fespace);
  tmp_du_dt = 0.0;
  mfem::Vector B;
  M->FormLinearSystem(ess_tdof_list, tmp_du_dt, tmp_z, A, du_dt, B, 0);
  M_solver.Mult(B, du_dt);
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

    dt_tensor_thermal_coeff_m->SetAConst(dt);
    T->AddDomainIntegrator(
        new mfem::DiffusionIntegrator(*dt_tensor_thermal_coeff_m));
    T->Assemble();
    T->FormSystemMatrix(ess_tdof_list, T_op);
    if (UsesTensorBasis(fespace)) {
      delete T_prec;
      T_prec = new mfem::OperatorJacobiSmoother(*T, ess_tdof_list);
      T_solver.SetPreconditioner(*T_prec);
    }
    T_solver.SetOperator(*T_op);
    T->Finalize();
    current_dt = dt;
  }
  MFEM_VERIFY(dt == current_dt, "");  // SDIRK methods use the same dt

  mfem::ParGridFunction tmp_u(&fespace);
  mfem::ParLinearForm tmp_z(&fespace);
  tmp_u.SetFromTrueDofs(u);
  K->Mult(tmp_u, tmp_z);
  tmp_z.Neg();  // z = -z

  // Add inhomogeneous Neumann to RHS
  if (inhomogeneous_neumann_active_m) {
    tmp_z.Add(1.0, *neumann_m);
  }

  mfem::ParGridFunction tmp_du_dt(&fespace);
  tmp_du_dt = 0.0;
  mfem::Vector B;
  T->FormLinearSystem(ess_tdof_list, tmp_du_dt, tmp_z, A, du_dt, B);
  T_solver.Mult(B, du_dt);
}

void ConductionOperator::SetParameters(const mfem::Vector& u) {
  if (K == nullptr) {
    mfem::DenseMatrix full_value(tensor_kappa_m.data(), mesh_m.GetDimension(),
                                 mesh_m.GetDimension());

    // Make assert that this (full_value) is positive definite
    DEBUG_ASSERT(full_value.Det() > 0.0, global_assert{}, DebugLevel::CHEAP{},
                 "Negative determinant value for tensor thermal coefficient. "
                 "Not Positive Definite. Determinant value: " +
                     std::to_string(full_value.Det()));

    delete tensor_thermal_coeff_m;
    tensor_thermal_coeff_m = new mfem::MatrixConstantCoefficient(full_value);

    // dt will be reset instead of 1.0 in below in implicit solve
    delete dt_tensor_thermal_coeff_m;
    dt_tensor_thermal_coeff_m =
        new mfem::ScalarMatrixProductCoefficient(1.0, *tensor_thermal_coeff_m);

    delete K;
    K = new mfem::ParBilinearForm(&fespace);
    K->SetAssemblyLevel(mfem::AssemblyLevel::PARTIAL);
    K->AddDomainIntegrator(
        new mfem::DiffusionIntegrator(*tensor_thermal_coeff_m));
    K->Assemble();  // keep sparsity pattern of M and K the same
    K->Finalize();
  }
  delete T;
  T = nullptr;  // Delete here to be reset on first entrance to SolveImplicit
}

void ConductionOperator::UpdateBoundaryConditions(mfem::Vector& u) {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Entered UpdateBoundaryConditions.");
  DEBUG_ASSERT(
      boundary_conditions_m.find("temperature") != boundary_conditions_m.end(),
      global_assert{}, DebugLevel::CHEAP{},
      "\"temperature\" not found in supplied boundary conditions.");
  const auto& bc_manager = boundary_conditions_m.at("temperature");
  if (bc_manager.GetNumberOfTimeVaryingNeumannConditions() == 0) {
    if (bc_manager.GetNumberOfTimeVaryingDirichletConditions() == 0) {
      SPDLOG_LOGGER_INFO(MAIN_LOG, "No time-varying conditions to be updated");
      return;
    }
  }

  if (inhomogeneous_neumann_active_m) {
    this->ResetNeumannCondition();
  }

  mfem::Array<int> tdof_list;

  const int number_of_boundary_conditions =
      bc_manager.GetNumberOfBoundaryConditions();

  int boundary_conditions_updated = 0;
  mfem::Vector boundary_temp;
  mfem::ParGridFunction temperature_gf(&fespace);
  for (int n = 0; n < number_of_boundary_conditions; ++n) {
    const auto& condition = bc_manager.GetBoundaryCondition(n + 1);

    mfem::Array<int>& boundary = boundary_marker_m[n];
    fespace.GetEssentialTrueDofs(boundary, tdof_list);

    if (condition.IsTimeVarying()) {
      SPDLOG_LOGGER_INFO(
          MAIN_LOG,
          "Updating boundary condition (type:{}) for tag {}. Consists of "
          "{} border elements",
          static_cast<int>(condition.GetBCType()), n + 1, tdof_list.Size());
    }

    // FIXME: For time-varying Dirichlet, should technically
    // be setting a Dirichlet condition on the acceleration too.
    // Will put off for now until discussing how to determine
    // the acceleration from discrete data points.

    switch (condition.GetBCType()) {
      case BoundaryConditionType::HOMOGENEOUS_DIRICHLET: {
        // Do nothing, cannot time-vary
        break;
      }
      case BoundaryConditionType::DIRICHLET: {
        if (condition.IsTimeVarying()) {
          ++boundary_conditions_updated;
          if (condition.IsSpatiallyVarying()) {
            const auto& values = condition.GetValues();
            const auto& indices = condition.GetIndices();
            this->SetTrueDofsFromVertexData(indices.Size(), indices.Data(),
                                            values.Data(), boundary,
                                            temperature_gf);
            temperature_gf.GetTrueDofs(boundary_temp);
            for (int i = 0; i < tdof_list.Size(); ++i) {
              const int j = tdof_list[i];
              u(j) = boundary_temp(j);
            }
          } else {
            const auto& values = condition.GetValues();
            u.SetSubVector(tdof_list, *values.Data());
          }
        }
        break;
      }
      case BoundaryConditionType::HOMOGENEOUS_NEUMANN: {
        // Nothing required, natural BC
        break;
      }
      case BoundaryConditionType::NEUMANN: {
        // Note: Need to add back whether time-varying or not
        // since reset entire ParLinearForm representing Neumann
        if (condition.IsTimeVarying()) {
          ++boundary_conditions_updated;
        }
        if (condition.IsSpatiallyVarying()) {
          const auto& values = condition.GetValues();
          const auto& indices = condition.GetIndices();
          this->SetTrueDofsFromVertexData(indices.Size(), indices.Data(),
                                          values.Data(), boundary,
                                          temperature_gf);
          this->AddNeumannCondition(n + 1, temperature_gf);
        } else {
          const auto& values = condition.GetValues();
          this->AddNeumannCondition(n + 1, *values.Data());
        }
        break;
      }
      default:
        DEBUG_ASSERT(
            false, global_assert{}, DebugLevel::ALWAYS{},
            "Unknown BoundaryConditionType. Set type is: " +
                std::to_string(static_cast<int>(condition.GetBCType())));
    }
  }

  SPDLOG_LOGGER_INFO(MAIN_LOG, "Updated {} boundary conditions.",
                     boundary_conditions_updated);

  if (inhomogeneous_neumann_active_m) {
    this->FinalizeNeumannCondition();
  }
}

ConductionOperator::~ConductionOperator(void) {
  delete M;
  delete K;
  delete T;
  delete T_prec;
  delete M_prec;
  for (auto& elem : neumann_coefficient_m) {
    delete elem;
    elem = nullptr;
  }
  delete neumann_m;
  delete tensor_thermal_coeff_m;
  delete dt_tensor_thermal_coeff_m;
}

void ConductionOperator::SetTrueDofsFromVertexData(
    const std::size_t a_size, const int* a_vertex_list,
    const double* a_vertex_data, mfem::Array<int>& a_boundary,
    mfem::ParGridFunction& a_temperature_gf) {
  mfem::ParGridFunction coarse_temperature_gf(&fespace_linear_m);
  coarse_temperature_gf = 0.0;
  // CASTING AWAY CONST TO WRAP IN MFEM VALUES. WILL NOT MODIFY
  mfem::Vector values(const_cast<double*>(a_vertex_data), a_size);
  mfem::Array<int> indices(const_cast<int*>(a_vertex_list), a_size);
  coarse_temperature_gf.SetSubVector(indices, values);

  mfem::GridFunctionCoefficient map_coeff(&coarse_temperature_gf);

  a_temperature_gf = 0.0;
  a_temperature_gf.ProjectBdrCoefficient(map_coeff, a_boundary);
}

void ConductionOperator::ResetNeumannCondition(void) {
  // FIXME: Make this more efficient. Maybe do not need full delete?
  delete neumann_m;
  neumann_m = new mfem::ParLinearForm(&fespace);
}

void ConductionOperator::AddNeumannCondition(
    const int a_tag, const mfem::ParGridFunction& a_grid_function) {
  DEBUG_ASSERT(neumann_m != nullptr, global_assert{}, DebugLevel::CHEAP{});
  DEBUG_ASSERT(a_tag > 0, global_assert{}, DebugLevel::CHEAP{},
               "Tag value must be strictly positive. Current tag value is: " +
                   std::to_string(a_tag));
  DEBUG_ASSERT(a_tag - 1 < static_cast<int>(boundary_marker_m.size()),
               global_assert{}, DebugLevel::CHEAP{},
               "Tag value must exist on mesh. Current tag value is: " +
                   std::to_string(a_tag));
  neumann_coefficient_m[a_tag - 1] =
      new mfem::GridFunctionCoefficient(&a_grid_function);
  neumann_m->AddBoundaryIntegrator(
      new mfem::BoundaryLFIntegrator(*(neumann_coefficient_m[a_tag - 1])),
      boundary_marker_m[a_tag - 1]);
}

void ConductionOperator::AddNeumannCondition(const int a_tag,
                                             const double a_value) {
  DEBUG_ASSERT(neumann_m != nullptr, global_assert{}, DebugLevel::CHEAP{});
  DEBUG_ASSERT(a_tag > 0, global_assert{}, DebugLevel::CHEAP{},
               "Tag value must be strictly positive. Current tag value is: " +
                   std::to_string(a_tag));
  DEBUG_ASSERT(a_tag - 1 < static_cast<int>(boundary_marker_m.size()),
               global_assert{}, DebugLevel::CHEAP{},
               "Tag value must exist on mesh. Current tag value is: " +
                   std::to_string(a_tag));

  neumann_coefficient_m[a_tag - 1] = new mfem::ConstantCoefficient(a_value);
  neumann_m->AddBoundaryIntegrator(
      new mfem::BoundaryLFIntegrator(*(neumann_coefficient_m[a_tag - 1])),
      boundary_marker_m[a_tag - 1]);
}

void ConductionOperator::FinalizeNeumannCondition(void) {
  neumann_m->Assemble();
}

}  // namespace chyps
