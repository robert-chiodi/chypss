// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chyps/simulation_configurations/cooled_rod.hpp"

#include <functional>

#include "chyps/debug_assert.hpp"

namespace chyps {

namespace cooled_rod {

namespace scalar {

void AddParserOptions(InputParser& a_parser) {
  a_parser.AddOption(
      "SimulationInitializer/Initializers/cooled_rod/scalar/"
      "approximation_terms",
      "Number of terms to use in CooledRod solution approximation", 1);
}

bool SupportedFieldType(const DataFieldType a_field_type) {
  static constexpr std::array<DataFieldType, 2> supported_types{
      {DataFieldType::GRID_FUNCTION, DataFieldType::TRUE_DOF}};
  return std::find(supported_types.begin(), supported_types.end(),
                   a_field_type) != supported_types.end();
}

void InitializeData(const DataFieldType a_field_type,
                    const nlohmann::json& a_json_object,
                    const InputParser& a_full_parser, const Mesh& a_mesh,
                    mfem::ParFiniteElementSpace& a_finite_element_space,
                    mfem::Vector& a_data) {
  const auto approximation_terms =
      a_json_object["approximation_terms"].get<std::size_t>();
  auto value_setter = [=](const mfem::Vector& position) {
    // This is the analytical solution for an initial
    // condition of 25 and Neumann on left/Dirichelt on right.
    double result = 0.0;
    for (std::size_t n = 1; n <= approximation_terms; ++n) {
      const double dn = static_cast<double>(n);
      const double ic_factor = 8.0 / std::pow(M_PI * (1.0 - 2.0 * dn), 2);
      result += ic_factor * std::cos((dn - 0.5) * M_PI * position[0]);
    }
    return result + 24.0 + position[0];
  };

  mfem::FunctionCoefficient value_setter_function(value_setter);

  switch (a_field_type) {
    case DataFieldType::GRID_FUNCTION: {
      mfem::ParGridFunction grid_function(&a_finite_element_space);
      grid_function.ProjectCoefficient(value_setter_function);
      a_data = grid_function;
      break;
    }

    case DataFieldType::TRUE_DOF: {
      mfem::ParGridFunction grid_function(&a_finite_element_space);
      grid_function.ProjectCoefficient(value_setter_function);
      grid_function.GetTrueDofs(a_data);
      break;
    }

    default:
      DEBUG_ASSERT(false, global_assert{}, DebugLevel::ALWAYS{},
                   "Unknown field type.");
  }
}

}  // namespace scalar
}  // namespace cooled_rod

}  // namespace chyps
