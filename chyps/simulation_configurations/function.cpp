// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chyps/simulation_configurations/function.hpp"

#include <array>
#include <functional>
#include <vector>

#include "chyps/debug_assert.hpp"

namespace chyps {

namespace function {

namespace scalar {

void AddParserOptions(InputParser& a_parser) {}

bool SupportedFieldType(const DataFieldType a_field_type) {
  static constexpr std::array<DataFieldType, 3> supported_types{
      {DataFieldType::GRID_FUNCTION, DataFieldType::TRUE_DOF,
       DataFieldType::ELEMENT}};
  return std::find(supported_types.begin(), supported_types.end(),
                   a_field_type) != supported_types.end();
}

void InitializeData(
    const DataFieldType a_field_type, const nlohmann::json& a_json_object,
    const InputParser& a_full_parser, const Mesh& a_mesh,
    mfem::ParFiniteElementSpace& a_finite_element_space, mfem::Vector& a_data,
    std::function<double(const mfem::Vector&)> a_setting_function) {
  switch (a_field_type) {
    case DataFieldType::GRID_FUNCTION: {
      mfem::ParGridFunction grid_function(&a_finite_element_space);

      // a_setting_function needs to be of form
      // double func(const mfem::Vector& x);
      mfem::FunctionCoefficient value_setter_function(a_setting_function);
      grid_function.ProjectCoefficient(value_setter_function);
      grid_function.GetTrueDofs(a_data);
      a_data = grid_function;
      break;
    }

    case DataFieldType::TRUE_DOF: {
      mfem::ParGridFunction grid_function(&a_finite_element_space);

      // a_setting_function needs to be of form
      // double func(const mfem::Vector& x);
      mfem::FunctionCoefficient value_setter_function(a_setting_function);
      grid_function.ProjectCoefficient(value_setter_function);
      grid_function.GetTrueDofs(a_data);
      break;
    }

    case DataFieldType::ELEMENT: {
      // a_setting_function needs to be of form
      // double func(const mfem::Vector& x);
      // Sets based on algebraic center of element
      const int nelem =
          static_cast<int>(a_mesh.GetLocalCount<MeshElement::ELEMENT>());
      a_data.SetSize(nelem);
      mfem::Vector center(a_mesh.GetDimension());
      const mfem::ParMesh& mfem_mesh = a_mesh.GetMfemMesh();
      const auto mfem_mesh_elements = mfem_mesh.GetElementsArray();

      for (int n = 0; n < nelem; ++n) {
        const mfem::Element& elem = *(mfem_mesh_elements[n]);

        // Compute algebraic center point for element
        auto nvert = elem.GetNVertices();
        const auto verts = elem.GetVertices();
        center = 0.0;
        for (int j = 0; j < nvert; ++j) {
          const double* pos = mfem_mesh.GetVertex(verts[j]);
          for (int d = 0; d < center.Size(); ++d) {
            center(d) += pos[d];
          }
        }
        center /= static_cast<double>(nvert);

        // Set value based on algebraic center
        a_data(n) = a_setting_function(center);
      }
      break;
    }

    default:
      DEBUG_ASSERT(false, global_assert{}, DebugLevel::ALWAYS{},
                   "Unknown field type.");
  }
}

}  // namespace scalar

namespace matrix {

void AddParserOptions(InputParser& a_parser) {}

bool SupportedFieldType(const DataFieldType a_field_type) {
  static constexpr std::array<DataFieldType, 1> supported_types{
      {DataFieldType::ELEMENT}};
  return std::find(supported_types.begin(), supported_types.end(),
                   a_field_type) != supported_types.end();
}

void InitializeData(const DataFieldType a_field_type,
                    const nlohmann::json& a_json_object,
                    const InputParser& a_full_parser, const Mesh& a_mesh,
                    mfem::ParFiniteElementSpace& a_finite_element_space,
                    const int a_number_of_rows, const int a_number_of_columns,
                    std::vector<mfem::DenseMatrix>& a_data,
                    std::function<void(const mfem::Vector&, mfem::DenseMatrix&)>
                        a_setting_function) {
  mfem::DenseMatrix base_matrix(a_number_of_rows, a_number_of_columns);
  switch (a_field_type) {
    case DataFieldType::ELEMENT: {
      // a_setting_function needs to be of form
      // void func(const mfem::Vector& x, mfem::DenseMatrix& K);
      // Sets based on algebraic center of element
      const int nelem =
          static_cast<int>(a_mesh.GetLocalCount<MeshElement::ELEMENT>());
      a_data.resize(nelem);
      mfem::Vector center(a_mesh.GetDimension());
      const mfem::ParMesh& mfem_mesh = a_mesh.GetMfemMesh();
      const auto mfem_mesh_elements = mfem_mesh.GetElementsArray();

      for (int n = 0; n < nelem; ++n) {
        const mfem::Element& elem = *(mfem_mesh_elements[n]);

        // Compute algebraic center point for element
        auto nvert = elem.GetNVertices();
        const auto verts = elem.GetVertices();
        center = 0.0;
        for (int j = 0; j < nvert; ++j) {
          const double* pos = mfem_mesh.GetVertex(verts[j]);
          for (int d = 0; d < center.Size(); ++d) {
            center(d) += pos[d];
          }
        }
        center /= static_cast<double>(nvert);
        a_setting_function(center, base_matrix);
        // Set value based on algebraic center
        a_data[n] = base_matrix;
      }
      break;
    }
    default:
      DEBUG_ASSERT(false, global_assert{}, DebugLevel::ALWAYS{},
                   "Unknown field type.");
  }
}

}  // namespace matrix

}  // namespace function

}  // namespace chyps
