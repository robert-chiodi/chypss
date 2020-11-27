// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chyps/mfem_coefficients.hpp"

#include "chyps/debug_assert.hpp"

namespace chyps {

MaterialVaryingMatrixCoefficient::MaterialVaryingMatrixCoefficient(
    const mfem::ParFiniteElementSpace& a_fe_space)
    : MatrixCoefficient(a_fe_space.GetParMesh()->Dimension(), false),
      finite_element_space_m(&a_fe_space),
      material_matrices_m() {
  const auto& pmesh = *(finite_element_space_m->GetParMesh());
  material_matrices_m.resize(pmesh.attributes.Size(),
                             mfem::DenseMatrix(pmesh.Dimension()));
}

mfem::DenseMatrix& MaterialVaryingMatrixCoefficient::operator[](
    const int a_material_tag) {
  DEBUG_ASSERT(a_material_tag > 0, global_assert{}, DebugLevel::CHEAP{});
  DEBUG_ASSERT(
      static_cast<std::size_t>(a_material_tag - 1) < material_matrices_m.size(),
      global_assert{}, DebugLevel::CHEAP{});
  return material_matrices_m[a_material_tag - 1];
}

const mfem::DenseMatrix& MaterialVaryingMatrixCoefficient::operator[](
    const int a_material_tag) const {
  DEBUG_ASSERT(a_material_tag > 0, global_assert{}, DebugLevel::CHEAP{});
  DEBUG_ASSERT(
      static_cast<std::size_t>(a_material_tag - 1) < material_matrices_m.size(),
      global_assert{}, DebugLevel::CHEAP{});
  return material_matrices_m[a_material_tag - 1];
}

void MaterialVaryingMatrixCoefficient::Eval(
    mfem::DenseMatrix& a_K, mfem::ElementTransformation& a_T,
    const mfem::IntegrationPoint& a_ip) {
  DEBUG_ASSERT(
      static_cast<std::size_t>(a_T.Attribute - 1) < material_matrices_m.size(),
      global_assert{}, DebugLevel::CHEAP{});
  a_K = material_matrices_m[a_T.Attribute - 1];
}

ElementVaryingScalarCoefficient::ElementVaryingScalarCoefficient(
    const mfem::ParFiniteElementSpace& a_fe_space)
    : Coefficient(), finite_element_space_m(&a_fe_space), element_values_m() {
  const auto& pmesh = *(finite_element_space_m->GetParMesh());
  element_values_m.resize(pmesh.GetNE());
}

double& ElementVaryingScalarCoefficient::operator[](const int a_element_index) {
  DEBUG_ASSERT(a_element_index >= 0, global_assert{}, DebugLevel::CHEAP{});
  DEBUG_ASSERT(
      static_cast<std::size_t>(a_element_index) < element_values_m.size(),
      global_assert{}, DebugLevel::CHEAP{});
  return element_values_m[a_element_index];
}

double ElementVaryingScalarCoefficient::operator[](
    const int a_element_index) const {
  DEBUG_ASSERT(a_element_index >= 0, global_assert{}, DebugLevel::CHEAP{});
  DEBUG_ASSERT(
      static_cast<std::size_t>(a_element_index) < element_values_m.size(),
      global_assert{}, DebugLevel::CHEAP{});
  return element_values_m[a_element_index];
}

double ElementVaryingScalarCoefficient::Eval(
    mfem::ElementTransformation& a_T, const mfem::IntegrationPoint& a_ip) {
  DEBUG_ASSERT(
      static_cast<std::size_t>(a_T.ElementNo) < element_values_m.size(),
      global_assert{}, DebugLevel::CHEAP{});
  return element_values_m[a_T.ElementNo];
}

ElementVaryingMatrixCoefficient::ElementVaryingMatrixCoefficient(
    const mfem::ParFiniteElementSpace& a_fe_space)
    : MatrixCoefficient(a_fe_space.GetParMesh()->Dimension(), false),
      finite_element_space_m(&a_fe_space),
      element_matrices_m() {
  const auto& pmesh = *(finite_element_space_m->GetParMesh());
  element_matrices_m.resize(pmesh.GetNE(),
                            mfem::DenseMatrix(pmesh.Dimension()));
}

mfem::DenseMatrix& ElementVaryingMatrixCoefficient::operator[](
    const int a_element_index) {
  DEBUG_ASSERT(a_element_index >= 0, global_assert{}, DebugLevel::CHEAP{});
  DEBUG_ASSERT(
      static_cast<std::size_t>(a_element_index) < element_matrices_m.size(),
      global_assert{}, DebugLevel::CHEAP{});
  return element_matrices_m[a_element_index];
}

const mfem::DenseMatrix& ElementVaryingMatrixCoefficient::operator[](
    const int a_element_index) const {
  DEBUG_ASSERT(a_element_index >= 0, global_assert{}, DebugLevel::CHEAP{});
  DEBUG_ASSERT(
      static_cast<std::size_t>(a_element_index) < element_matrices_m.size(),
      global_assert{}, DebugLevel::CHEAP{});
  return element_matrices_m[a_element_index];
}

void ElementVaryingMatrixCoefficient::Eval(mfem::DenseMatrix& a_K,
                                           mfem::ElementTransformation& a_T,
                                           const mfem::IntegrationPoint& a_ip) {
  DEBUG_ASSERT(
      static_cast<std::size_t>(a_T.ElementNo) < element_matrices_m.size(),
      global_assert{}, DebugLevel::CHEAP{});
  a_K = element_matrices_m[a_T.ElementNo];
}

}  // namespace chyps
