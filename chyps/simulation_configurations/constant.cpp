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

Constant::Constant(InputParser& a_parser) : ConfigurationInitializer(a_parser) {
  parser_m.AddOptionNoDefault("SimulationInitializer/Constant/value",
                              "Value to set at all node locations", true);
}

void Constant::Initialize(void) {}

void Constant::FillRequiredData(RequiredData& a_data) {
  auto finite_element_collection = new mfem::H1_FECollection(
      parser_m["HeatSolver/order"].get<int>(), a_data.GetMesh().GetDimension());
  auto finite_element_space = new mfem::ParFiniteElementSpace(
      &(a_data.GetMesh().GetMfemMesh()), finite_element_collection);
  a_data.SetFiniteElementCollection(finite_element_collection);
  a_data.SetFiniteElementSpace(finite_element_space);

  auto temperature_field = new mfem::ParGridFunction(finite_element_space);

  const auto value =
      parser_m["SimulationInitializer/Constant/value"].get<double>();
  auto value_setter = [=](const mfem::Vector& position) { return value; };

  mfem::FunctionCoefficient temperature_setter(value_setter);
  temperature_field->ProjectCoefficient(temperature_setter);
  a_data.AddGridFunction("HeatSolver/temperature", temperature_field);
}

Constant::~Constant(void) {}

}  // namespace chyps
