// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chyps/simulation_configurations/attribute_varying.hpp"

#include <array>
#include <functional>
#include <vector>

#include "chyps/debug_assert.hpp"

namespace chyps {

namespace attribute_varying {

namespace scalar {

void AddParserOptions(InputParser& a_parser) {
  a_parser.AddOption(
      "SimulationInitializer/Initializers/attribute_varying/scalar/value",
      "Array of values with one entry per element attribute (material) on the "
      "mesh.");
}

bool SupportedFieldType(const DataFieldType a_field_type) {
  static constexpr std::array<DataFieldType, 3> supported_types{
      {DataFieldType::GRID_FUNCTION, DataFieldType::TRUE_DOF,
       DataFieldType::ELEMENT}};
  return std::find(supported_types.begin(), supported_types.end(),
                   a_field_type) != supported_types.end();
}

void InitializeData(const DataFieldType a_field_type,
                    const nlohmann::json& a_json_object,
                    const InputParser& a_full_parser, const Mesh& a_mesh,
                    mfem::ParFiniteElementSpace& a_finite_element_space,
                    mfem::Vector& a_data) {
  auto value = a_json_object["value"].get<std::vector<double>>();

  DEBUG_ASSERT(value.size() == static_cast<std::size_t>(
                                   a_mesh.GetMfemMesh().attributes.Size()),
               global_assert{}, DebugLevel::CHEAP{});

  switch (a_field_type) {
    case DataFieldType::GRID_FUNCTION: {
      mfem::ParGridFunction grid_function(&a_finite_element_space);
      mfem::Vector vec(value.data(), static_cast<int>(value.size()));
      mfem::PWConstCoefficient value_setter_function(vec);
      grid_function.ProjectCoefficient(value_setter_function);
      grid_function.GetTrueDofs(a_data);
      a_data = grid_function;
      break;
    }

    case DataFieldType::TRUE_DOF: {
      mfem::ParGridFunction grid_function(&a_finite_element_space);
      mfem::Vector vec(value.data(), static_cast<int>(value.size()));
      mfem::PWConstCoefficient value_setter_function(vec);
      grid_function.ProjectCoefficient(value_setter_function);
      grid_function.GetTrueDofs(a_data);
      break;
    }

    case DataFieldType::ELEMENT: {
      a_data.SetSize(a_mesh.GetLocalCount<MeshElement::ELEMENT>());
      for (int n = 0; n < a_data.Size(); ++n) {
        a_data(n) = value[a_mesh.GetMfemMesh().GetAttribute(n)];
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

void AddParserOptions(InputParser& a_parser) {
  a_parser.AddOption(
      "SimulationInitializer/Initializers/attribute_varying/matrix/value",
      "A JSON object with keys of attribute names with values "
      "that are a matrix of the appropriate size.");
}

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
                    std::vector<mfem::DenseMatrix>& a_data) {
  auto value_object = a_json_object["value"];

  const int number_of_attributes = a_mesh.GetMfemMesh().attributes.Size();

  std::vector<mfem::DenseMatrix> base_matrices(
      number_of_attributes,
      mfem::DenseMatrix(a_number_of_rows, a_number_of_columns));
  for (int n = 1; n <= number_of_attributes; ++n) {
    auto flat_matrix =
        value_object[std::to_string(n)].get<std::vector<double>>();
    DEBUG_ASSERT(
        flat_matrix.size() ==
            static_cast<std::size_t>(a_number_of_rows * a_number_of_columns),
        global_assert{}, DebugLevel::CHEAP{});
    for (int j = 0; j < a_number_of_columns; ++j) {
      for (int i = 0; i < a_number_of_rows; ++i) {
        base_matrices[n - 1](i, j) = flat_matrix[i + j * a_number_of_rows];
      }
    }
  }

  switch (a_field_type) {
    case DataFieldType::ELEMENT: {
      a_data.resize(a_mesh.GetLocalCount<MeshElement::ELEMENT>());
      for (int n = 0; n < static_cast<int>(a_data.size()); ++n) {
        a_data[n] = base_matrices[a_mesh.GetMfemMesh().GetAttribute(n)];
      }
      break;
    }
    default:
      DEBUG_ASSERT(false, global_assert{}, DebugLevel::ALWAYS{},
                   "Unknown field type.");
  }
}

}  // namespace matrix

}  // namespace attribute_varying

}  // namespace chyps
