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

#include <functional>

#include "chyps/debug_assert.hpp"

namespace chyps {

namespace constant {

void AddParserOptions(InputParser& a_parser) {
  a_parser.AddOption("SimulationInitializer/Initializers/constant/value",
                     "Sets the value of alll nodes to the supplied value.");
}

void InitializeData(const nlohmann::json& a_json_object,
                    const InputParser& a_full_parser,
                    mfem::ParFiniteElementSpace& a_finite_element_space,
                    mfem::Vector& a_data) {
  mfem::ParGridFunction grid_function(&a_finite_element_space);

  const auto value = a_json_object["value"].get<double>();
  auto value_setter = [=](const mfem::Vector& position) { return value; };

  mfem::FunctionCoefficient value_setter_function(value_setter);
  grid_function.ProjectCoefficient(value_setter_function);
  grid_function.GetTrueDofs(a_data);
}
}  // namespace constant

}  // namespace chyps
