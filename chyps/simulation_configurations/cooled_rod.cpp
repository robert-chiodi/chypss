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

CooledRod::CooledRod(InputParser& a_parser)
    : ConfigurationInitializer(a_parser) {
  parser_m.AddOption(
      "SimulationInitializer/CooledRod/approximation_terms",
      "Number of terms to use in CooledRod solution approximation", 1);
}

void CooledRod::Initialize(void) {}

void CooledRod::FillRequiredData(RequiredData& a_data) {
  auto finite_element_collection = new mfem::H1_FECollection(
      parser_m["HeatSolver/order"].get<int>(), a_data.GetMesh().GetDimension());
  auto finite_element_space = new mfem::ParFiniteElementSpace(
      &(a_data.GetMesh().GetMfemMesh()), finite_element_collection);
  a_data.SetFiniteElementCollection(finite_element_collection);
  a_data.SetFiniteElementSpace(finite_element_space);

  auto temperature_field = new mfem::ParGridFunction(finite_element_space);

  const auto approximation_terms =
      parser_m["SimulationInitializer/CooledRod/approximation_terms"]
          .get<double>();
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

  mfem::FunctionCoefficient temperature_setter(value_setter);
  temperature_field->ProjectCoefficient(temperature_setter);
  a_data.AddGridFunction("HeatSolver/temperature", temperature_field);
}

CooledRod::~CooledRod(void) {}

}  // namespace chyps
