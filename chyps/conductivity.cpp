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
#include "chyps/simulation.hpp"

namespace chyps {

Conductivity::Conductivity(Simulation& a_simulation,
                           mfem::ParFiniteElementSpace& a_finite_element_space)
    : sim_m(a_simulation),
      type_m(ConductivityType::INVALID),
      scalar_coefficient_m(nullptr),
      matrix_coefficient_m(nullptr) {
  const InputParser& parser = sim_m.GetParser();
  const int number_of_attributes =
      sim_m.GetMesh().GetMfemMesh().attributes.Size();
  const int dimension = sim_m.GetMesh().GetDimension();

  type_m = this->StringToType(
      parser["HeatSolver/ConductionOperator/conductivity_type"]);
  switch (type_m) {
    case ConductivityType::CONSTANT_SCALAR: {
      auto value =
          parser["HeatSolver/ConductionOperator/ConstantScalar"].get<double>();
      scalar_coefficient_m = new mfem::ConstantCoefficient(value);

      break;
    }
    case ConductivityType::CONSTANT_MATRIX: {
      auto value = parser["HeatSolver/ConductionOperator/ConstantMatrix"]
                       .get<std::vector<double>>();
      DEBUG_ASSERT(
          value.size() == static_cast<std::size_t>(dimension * dimension),
          global_assert{}, DebugLevel::CHEAP{},
          "Thermal conductivity tensor of incorrect size. Provide as a "
          "column-major array of size MESH_DIM*MESH_DIM= " +
              std::to_string(dimension * dimension));
      mfem::DenseMatrix matrix(value.data(), dimension, dimension);
      matrix_coefficient_m = new mfem::MatrixConstantCoefficient(matrix);
      break;
    }

    case ConductivityType::MATERIAL_VARYING_SCALAR: {
      auto value = parser["HeatSolver/ConductionOperator/MaterialVaryingScalar"]
                       .get<std::vector<double>>();

      DEBUG_ASSERT(
          value.size() == static_cast<std::size_t>(number_of_attributes),
          global_assert{}, DebugLevel::CHEAP{},
          "Thermal conductivity tensor of incorrect size. "
          "MaterialVaryingScalar requires one entry per cell-material "
          "(attribute) in the domain. There are " +
              std::to_string(number_of_attributes) + " attributes");
      mfem::Vector conductivities(value.data(), number_of_attributes);
      scalar_coefficient_m = new mfem::PWConstCoefficient(conductivities);
      break;
    }

    case ConductivityType::MATERIAL_VARYING_MATRIX: {
      const auto value_list =
          parser["HeatSolver/ConductionOperator/MaterialVaryingMatrix"];
      matrix_coefficient_m =
          new MaterialVaryingMatrixCoefficient(a_finite_element_space);
      auto& material_varying_matrix =
          this->GetMaterialVaryingMatrixCoefficient();

      // Fill the matrix coefficient per material, as dictated in input file.
      for (int n = 1; n <= number_of_attributes; ++n) {
        DEBUG_ASSERT(
            value_list.contains(std::to_string(n)), global_assert{},
            DebugLevel::CHEAP{},
            "Thermal conductivity tensor not provided for mesh attribute tag " +
                std::to_string(n));

        auto value = value_list[std::to_string(n)].get<std::vector<double>>();
        DEBUG_ASSERT(
            value.size() == static_cast<std::size_t>(dimension * dimension),
            global_assert{}, DebugLevel::CHEAP{},
            "Thermal conductivity tensor of incorrect size for mesh "
            "attribute " +
                std::to_string(n) +
                ". Provide as a column-major array of size "
                "MESH_DIM*MESH_DIM= " +
                std::to_string(dimension * dimension));
        material_varying_matrix[n] =
            mfem::DenseMatrix(value.data(), dimension, dimension);
      }

      break;
    }

    case ConductivityType::ELEMENT_VARYING_SCALAR: {
      // Figure out way to write/read this
      DEBUG_ASSERT(false, global_assert{}, DebugLevel::ALWAYS{});
      const auto name_in_file =
          parser["HeatSolver/ConductionOperator/ElementVaryingScalar"]
              .get<std::string>();
      break;
    }

    case ConductivityType::ELEMENT_VARYING_MATRIX: {
      // Figure out way to write/read this
      DEBUG_ASSERT(false, global_assert{}, DebugLevel::ALWAYS{});
      const auto name_in_file =
          parser["HeatSolver/ConductionOperator/ElementVaryingMatrix"]
              .get<std::string>();
      break;
    }

    default:
      DEBUG_ASSERT(false, global_assert{}, DebugLevel::ALWAYS{},
                   "Unknown Conductivity type.");
  }
}

bool Conductivity::IsScalarCoefficient(void) const {
  // Convention, Scalar coefficients are even.
  return static_cast<int>(type_m) % 2 == 0;
}

bool Conductivity::IsMatrixCoefficient(void) const {
  // Convention, Matrix coefficients are odd.
  return static_cast<int>(type_m) % 2 == 1;
}

bool Conductivity::CanTimeVary(void) const {
  return this->GetType() == ConductivityType::ELEMENT_VARYING_SCALAR ||
         this->GetType() == ConductivityType::ELEMENT_VARYING_MATRIX;
}

ConductivityType Conductivity::GetType(void) const { return type_m; }

mfem::Coefficient& Conductivity::GetScalarCoefficient(void) {
  DEBUG_ASSERT(scalar_coefficient_m != nullptr, global_assert{},
               DebugLevel::CHEAP{});
  DEBUG_ASSERT(this->IsScalarCoefficient(), global_assert{},
               DebugLevel::CHEAP{});
  return *scalar_coefficient_m;
}

const mfem::Coefficient& Conductivity::GetScalarCoefficient(void) const {
  DEBUG_ASSERT(scalar_coefficient_m != nullptr, global_assert{},
               DebugLevel::CHEAP{});
  DEBUG_ASSERT(this->IsScalarCoefficient(), global_assert{},
               DebugLevel::CHEAP{});
  return *scalar_coefficient_m;
}

mfem::MatrixCoefficient& Conductivity::GetMatrixCoefficient(void) {
  DEBUG_ASSERT(matrix_coefficient_m != nullptr, global_assert{},
               DebugLevel::CHEAP{});
  DEBUG_ASSERT(this->IsMatrixCoefficient(), global_assert{},
               DebugLevel::CHEAP{});
  return *matrix_coefficient_m;
}

const mfem::MatrixCoefficient& Conductivity::GetMatrixCoefficient(void) const {
  DEBUG_ASSERT(matrix_coefficient_m != nullptr, global_assert{},
               DebugLevel::CHEAP{});
  DEBUG_ASSERT(this->IsMatrixCoefficient(), global_assert{},
               DebugLevel::CHEAP{});
  return *matrix_coefficient_m;
}

mfem::ConstantCoefficient& Conductivity::GetConstantScalarCoefficient(void) {
  DEBUG_ASSERT(scalar_coefficient_m != nullptr, global_assert{},
               DebugLevel::CHEAP{});
  DEBUG_ASSERT(type_m == ConductivityType::CONSTANT_SCALAR, global_assert{},
               DebugLevel::CHEAP{});
  return *dynamic_cast<mfem::ConstantCoefficient*>(scalar_coefficient_m);
}
const mfem::ConstantCoefficient& Conductivity::GetConstantScalarCoefficient(
    void) const {
  DEBUG_ASSERT(scalar_coefficient_m != nullptr, global_assert{},
               DebugLevel::CHEAP{});
  DEBUG_ASSERT(type_m == ConductivityType::CONSTANT_SCALAR, global_assert{},
               DebugLevel::CHEAP{});
  return *dynamic_cast<mfem::ConstantCoefficient*>(scalar_coefficient_m);
}

mfem::MatrixConstantCoefficient& Conductivity::GetConstantMatrixCoefficient(
    void) {
  DEBUG_ASSERT(matrix_coefficient_m != nullptr, global_assert{},
               DebugLevel::CHEAP{});
  DEBUG_ASSERT(type_m == ConductivityType::CONSTANT_MATRIX, global_assert{},
               DebugLevel::CHEAP{});
  return *dynamic_cast<mfem::MatrixConstantCoefficient*>(matrix_coefficient_m);
}
const mfem::MatrixConstantCoefficient&
Conductivity::GetConstantMatrixCoefficient(void) const {
  DEBUG_ASSERT(matrix_coefficient_m != nullptr, global_assert{},
               DebugLevel::CHEAP{});
  DEBUG_ASSERT(type_m == ConductivityType::CONSTANT_MATRIX, global_assert{},
               DebugLevel::CHEAP{});
  return *dynamic_cast<mfem::MatrixConstantCoefficient*>(matrix_coefficient_m);
}

mfem::PWConstCoefficient& Conductivity::GetMaterialVaryingScalarCoefficient(
    void) {
  DEBUG_ASSERT(scalar_coefficient_m != nullptr, global_assert{},
               DebugLevel::CHEAP{});
  DEBUG_ASSERT(type_m == ConductivityType::MATERIAL_VARYING_SCALAR,
               global_assert{}, DebugLevel::CHEAP{});
  return *dynamic_cast<mfem::PWConstCoefficient*>(scalar_coefficient_m);
}
const mfem::PWConstCoefficient&
Conductivity::GetMaterialVaryingScalarCoefficient(void) const {
  DEBUG_ASSERT(scalar_coefficient_m != nullptr, global_assert{},
               DebugLevel::CHEAP{});
  DEBUG_ASSERT(type_m == ConductivityType::MATERIAL_VARYING_SCALAR,
               global_assert{}, DebugLevel::CHEAP{});
  return *dynamic_cast<mfem::PWConstCoefficient*>(scalar_coefficient_m);
}

MaterialVaryingMatrixCoefficient&
Conductivity::GetMaterialVaryingMatrixCoefficient(void) {
  DEBUG_ASSERT(matrix_coefficient_m != nullptr, global_assert{},
               DebugLevel::CHEAP{});
  DEBUG_ASSERT(type_m == ConductivityType::MATERIAL_VARYING_MATRIX,
               global_assert{}, DebugLevel::CHEAP{});
  return *dynamic_cast<MaterialVaryingMatrixCoefficient*>(matrix_coefficient_m);
}
const MaterialVaryingMatrixCoefficient&
Conductivity::GetMaterialVaryingMatrixCoefficient(void) const {
  DEBUG_ASSERT(matrix_coefficient_m != nullptr, global_assert{},
               DebugLevel::CHEAP{});
  DEBUG_ASSERT(type_m == ConductivityType::MATERIAL_VARYING_MATRIX,
               global_assert{}, DebugLevel::CHEAP{});
  return *dynamic_cast<MaterialVaryingMatrixCoefficient*>(matrix_coefficient_m);
}

ElementVaryingScalarCoefficient&
Conductivity::GetElementVaryingScalarCoefficient(void) {
  DEBUG_ASSERT(scalar_coefficient_m != nullptr, global_assert{},
               DebugLevel::CHEAP{});
  DEBUG_ASSERT(type_m == ConductivityType::ELEMENT_VARYING_SCALAR,
               global_assert{}, DebugLevel::CHEAP{});
  return *dynamic_cast<ElementVaryingScalarCoefficient*>(scalar_coefficient_m);
}
const ElementVaryingScalarCoefficient&
Conductivity::GetElementVaryingScalarCoefficient(void) const {
  DEBUG_ASSERT(scalar_coefficient_m != nullptr, global_assert{},
               DebugLevel::CHEAP{});
  DEBUG_ASSERT(type_m == ConductivityType::ELEMENT_VARYING_SCALAR,
               global_assert{}, DebugLevel::CHEAP{});
  return *dynamic_cast<ElementVaryingScalarCoefficient*>(scalar_coefficient_m);
}

ElementVaryingMatrixCoefficient&
Conductivity::GetElementVaryingMatrixCoefficient(void) {
  DEBUG_ASSERT(matrix_coefficient_m != nullptr, global_assert{},
               DebugLevel::CHEAP{});
  DEBUG_ASSERT(type_m == ConductivityType::ELEMENT_VARYING_MATRIX,
               global_assert{}, DebugLevel::CHEAP{});
  return *dynamic_cast<ElementVaryingMatrixCoefficient*>(matrix_coefficient_m);
}
const ElementVaryingMatrixCoefficient&
Conductivity::GetElementVaryingMatrixCoefficient(void) const {
  DEBUG_ASSERT(matrix_coefficient_m != nullptr, global_assert{},
               DebugLevel::CHEAP{});
  DEBUG_ASSERT(type_m == ConductivityType::ELEMENT_VARYING_MATRIX,
               global_assert{}, DebugLevel::CHEAP{});
  return *dynamic_cast<ElementVaryingMatrixCoefficient*>(matrix_coefficient_m);
}

Conductivity::~Conductivity(void) {
  delete scalar_coefficient_m;
  delete matrix_coefficient_m;
  scalar_coefficient_m = nullptr;
  matrix_coefficient_m = nullptr;
}

ConductivityType Conductivity::StringToType(const std::string& a_name) {
  std::array<std::string, static_cast<std::size_t>(ConductivityType::SIZE)>
      name_to_enum({"constant_scalar", "constant_matrix",
                    "material_varying_scalar", "material_varying_matrix",
                    "material_varying_scalar", "material_varying_matrix"});
  const auto location =
      std::find(name_to_enum.begin(), name_to_enum.end(), a_name);
  DEBUG_ASSERT(
      location != name_to_enum.end(), global_assert{}, DebugLevel::CHEAP{},
      "ConductivityType with name of \"" + a_name + "\" is not known.");
  return static_cast<ConductivityType>(
      std::distance(name_to_enum.begin(), location));
}

}  // namespace chyps
