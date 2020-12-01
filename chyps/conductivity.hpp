// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef CHYPS_CONDUCTIVITY_HPP_
#define CHYPS_CONDUCTIVITY_HPP_

#include <vector>

#include <mfem/mfem.hpp>

#include "chyps/input_parser.hpp"
#include "chyps/mfem_coefficients.hpp"

namespace chyps {

// Forward declare Simulation
class Simulation;

// Scalar options are even, matrix options are positive. Keep this true, used in
// ConductivityHandler.
enum class ConductivityType {
  INVALID = -1,
  CONSTANT_SCALAR,
  CONSTANT_MATRIX,
  MATERIAL_VARYING_SCALAR,
  MATERIAL_VARYING_MATRIX,
  ELEMENT_VARYING_SCALAR,
  ELEMENT_VARYING_MATRIX,
  SIZE
};

// Define what each type means and how it must be provided in an InputFile.

// InputFile usage example for each type:
//
// For constant_scalar type.
//
// "HeatSolver" : {
// "Conductivity":{
//                 "type": "constant_scalar",
//                 "ConstantScalar": 1.0 // constant value
//                }
//                }
//
//
//
// For constant_matrix type.
//
// "HeatSolver" : {
// "Conductivity":{
//                 "type": "constant_matrix",
//                 // MESH_DIM*MESH_DIM
//                 // long array in column-major ordering for matrix
//                 "ConstantMatrix": [1.0, 0.0, 0.0, 1.0]
//                }
//                }
//
//
// For material_varying_scalar type.
//
// "HeatSolver" : {
// "Conductivity":{
//                 "type": "material_varying_scalar",
//                 // One vector entry per tag, tag value is entry index + 1
//                 // Below would be for tag values [1,4]
//                 "MaterialVaryingScalar": [0.25, 0.5, 0.75, 1.0]
//                }
//                }
//
//
// For material_varying_matrix type.
//
// "HeatSolver" : {
// "Conductivity":{
//                 "type": "material_varying_matrix",
//                 // Attribute (Material) tag used as key.
//                 // Each key gets MESH_DIM*MESH_DIM array in column-major
//                 "MaterialVaryingMatrix":{
//                                        "1": [0.25, 0.0, 0.0, 0.25],
//                                        "2": [0.5, 0.0, 0.0, 0.5],
//                                        "3": [0.75, 0.0, 0.0, 0.75],
//                                        "4": [1.0, 0.0, 0.0, 1.0],
//                                         }
//                }
//                }
//
//
// For element_varying_scalar type.
//
// "HeatSolver" : {
// "Conductivity":{
//                 "type": "element_varying_scalar",
//                 // Provide path to data providing one value per element
//                 "ElementVaryingScalar": "path_to_data_in_BP4_file"
//                }
//                }
//
//
// For element_varying_matrix type.
//
// "HeatSolver" : {
// "Conductivity":{
//                 "type": "element_varying_matrix",
//                 // Provide path to data providing MESH_DIM*MESH_DIM  values
//                 // per element in column-major order
//                 "ElementVaryingScalar":  "path_to_data_in_BP4_file"
//                }
//                }
//
//

/// \class Conductivity chyps/conductivity.hpp conductivity.hpp
/// \brief Class to control the conductivity used for the HeatSolver. Can
/// handle scalar or tensor values that are constant, vary based on element
/// attribute, or vary based on element index. It is assumed only those that
/// vary with element index change in time.
class Conductivity {
 public:
  Conductivity(void) = delete;

  Conductivity(Simulation& a_simulation,
               mfem::ParFiniteElementSpace& a_finite_element_space);

  bool IsScalarCoefficient(void) const;
  bool IsMatrixCoefficient(void) const;
  bool CanTimeVary(void) const;

  ConductivityType GetType(void) const;

  mfem::Coefficient& GetScalarCoefficient(void);
  const mfem::Coefficient& GetScalarCoefficient(void) const;
  mfem::MatrixCoefficient& GetMatrixCoefficient(void);
  const mfem::MatrixCoefficient& GetMatrixCoefficient(void) const;

  mfem::ConstantCoefficient& GetConstantScalarCoefficient(void);
  const mfem::ConstantCoefficient& GetConstantScalarCoefficient(void) const;

  mfem::MatrixConstantCoefficient& GetConstantMatrixCoefficient(void);
  const mfem::MatrixConstantCoefficient& GetConstantMatrixCoefficient(
      void) const;

  mfem::PWConstCoefficient& GetMaterialVaryingScalarCoefficient(void);
  const mfem::PWConstCoefficient& GetMaterialVaryingScalarCoefficient(
      void) const;

  MaterialVaryingMatrixCoefficient& GetMaterialVaryingMatrixCoefficient(void);
  const MaterialVaryingMatrixCoefficient& GetMaterialVaryingMatrixCoefficient(
      void) const;

  ElementVaryingScalarCoefficient& GetElementVaryingScalarCoefficient(void);
  const ElementVaryingScalarCoefficient& GetElementVaryingScalarCoefficient(
      void) const;

  ElementVaryingMatrixCoefficient& GetElementVaryingMatrixCoefficient(void);
  const ElementVaryingMatrixCoefficient& GetElementVaryingMatrixCoefficient(
      void) const;

  ~Conductivity(void);

 private:
  static ConductivityType StringToType(const std::string& a_name);

  Simulation& sim_m;
  // Could make this a tagged union
  ConductivityType type_m;
  mfem::Coefficient* scalar_coefficient_m;
  mfem::MatrixCoefficient* matrix_coefficient_m;
};

}  // namespace chyps

#endif  // CHYPS_CONDUCTIVITY_HPP_
