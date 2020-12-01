// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chyps/simulation_configurations/constant.hpp"

#include <array>
#include <functional>
#include <vector>

#include "chyps/debug_assert.hpp"

namespace chyps {

namespace constant {

namespace scalar {

void AddParserOptions(InputParser& a_parser) {
  a_parser.AddOption("SimulationInitializer/Initializers/constant/scalar/value",
                     "Sets the value of all nodes to the supplied value.");
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
  const auto value = a_json_object["value"].get<double>();

  switch (a_field_type) {
    case DataFieldType::GRID_FUNCTION: {
      mfem::ParGridFunction grid_function(&a_finite_element_space);

      mfem::ConstantCoefficient value_setter_function(value);
      grid_function.ProjectCoefficient(value_setter_function);
      grid_function.GetTrueDofs(a_data);
      a_data = grid_function;
      break;
    }

    case DataFieldType::TRUE_DOF: {
      mfem::ParGridFunction grid_function(&a_finite_element_space);

      mfem::ConstantCoefficient value_setter_function(value);
      grid_function.ProjectCoefficient(value_setter_function);
      grid_function.GetTrueDofs(a_data);
      break;
    }

    case DataFieldType::ELEMENT: {
      a_data.SetSize(a_mesh.GetLocalCount<MeshElement::ELEMENT>());
      a_data = value;
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
  a_parser.AddOption("SimulationInitializer/Initializers/constant/matrix/value",
                     "Array of length row*column in column-major order.");
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
  auto values = a_json_object["value"].get<std::vector<double>>();
  DEBUG_ASSERT(values.size() == static_cast<std::size_t>(a_number_of_rows *
                                                         a_number_of_columns),
               global_assert{}, DebugLevel::CHEAP{},
               "Supplied array form of matrix is of incorrect size.");
  // Cast away const to use with mfem::DenseMatrix. Will not modify.
  mfem::DenseMatrix base_matrix(values.data(), a_number_of_rows,
                                a_number_of_columns);
  switch (a_field_type) {
    case DataFieldType::ELEMENT: {
      a_data.resize(a_mesh.GetLocalCount<MeshElement::ELEMENT>(), base_matrix);
      break;
    }
    default:
      DEBUG_ASSERT(false, global_assert{}, DebugLevel::ALWAYS{},
                   "Unknown field type.");
  }
}
}  // namespace matrix

}  // namespace constant

}  // namespace chyps
