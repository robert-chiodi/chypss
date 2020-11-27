// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef CHYPS_MFEM_COEFFICIENTS_HPP_
#define CHYPS_MFEM_COEFFICIENTS_HPP_

#include <vector>

#include <mfem/mfem.hpp>

namespace chyps {

class MaterialVaryingMatrixCoefficient : public mfem::MatrixCoefficient {
 public:
  MaterialVaryingMatrixCoefficient(
      const mfem::ParFiniteElementSpace& a_fe_space);

  mfem::DenseMatrix& operator[](const int a_material_tag);

  const mfem::DenseMatrix& operator[](const int a_material_tag) const;

  virtual void Eval(mfem::DenseMatrix& a_K, mfem::ElementTransformation& a_T,
                    const mfem::IntegrationPoint& a_ip) override final;

  virtual ~MaterialVaryingMatrixCoefficient(void) = default;

 private:
  const mfem::ParFiniteElementSpace* finite_element_space_m;
  std::vector<mfem::DenseMatrix> material_matrices_m;
};

class ElementVaryingScalarCoefficient : public mfem::Coefficient {
 public:
  ElementVaryingScalarCoefficient(
      const mfem::ParFiniteElementSpace& a_fe_space);

  double& operator[](const int a_element_index);

  double operator[](const int a_element_index) const;

  virtual double Eval(mfem::ElementTransformation& a_T,
                      const mfem::IntegrationPoint& a_ip) override final;

  virtual ~ElementVaryingScalarCoefficient(void) = default;

 private:
  const mfem::ParFiniteElementSpace* finite_element_space_m;
  std::vector<double> element_values_m;
};

class ElementVaryingMatrixCoefficient : public mfem::MatrixCoefficient {
 public:
  ElementVaryingMatrixCoefficient(
      const mfem::ParFiniteElementSpace& a_fe_space);

  mfem::DenseMatrix& operator[](const int a_element_index);

  const mfem::DenseMatrix& operator[](const int a_element_index) const;

  virtual void Eval(mfem::DenseMatrix& a_K, mfem::ElementTransformation& a_T,
                    const mfem::IntegrationPoint& a_ip) override final;

  virtual ~ElementVaryingMatrixCoefficient(void) = default;

 private:
  const mfem::ParFiniteElementSpace* finite_element_space_m;
  std::vector<mfem::DenseMatrix> element_matrices_m;
};

}  // namespace chyps

#endif  // CHYPS_MFEM_COEFFICIENTS_HPP_
