// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chyps/simulation_configurations/quadratic_pulse.hpp"

#include <functional>

#include "chyps/debug_assert.hpp"

namespace chyps {

QuadraticPulse::QuadraticPulse(InputParser& a_parser)
    : ConfigurationInitializer(a_parser) {
  parser_m.AddOptionDefault(
      "SimulationInitializer/QuadraticPulse/pulse_amplitude",
      "Amplitude of quadratic pulse (value at center)", 1.0);
  parser_m.AddOptionDefault(
      "SimulationInitializer/QuadraticPulse/approximation_terms",
      "Number of terms to use in pulse approximation", 1);
  parser_m.AddOptionNoDefault(
      "SimulationInitializer/QuadraticPulse/boundary_configuration",
      "Configuration of homogeneous boundary conditions. Supplied by string of "
      "4 characters, D for Dirichlet and N for Neumann, correspoding to tag "
      "locations, (e.g.,\"HHNN\")",
      true);
}

void QuadraticPulse::Initialize(void) {}

void QuadraticPulse::FillRequiredData(RequiredData& a_data) {
  auto finite_element_collection = new mfem::H1_FECollection(
      parser_m["HeatSolver/order"].get<int>(), a_data.GetMesh().GetDimension());
  auto finite_element_space = new mfem::ParFiniteElementSpace(
      &(a_data.GetMesh().GetMfemMesh()), finite_element_collection);
  a_data.SetFiniteElementCollection(finite_element_collection);
  a_data.SetFiniteElementSpace(finite_element_space);

  auto temperature_field = new mfem::ParGridFunction(finite_element_space);
  std::function<double(const mfem::Vector&)> initializing_function =
      this->SelectInitialConditions();
  mfem::FunctionCoefficient temperature_setter(initializing_function);
  temperature_field->ProjectCoefficient(temperature_setter);
  a_data.AddGridFunction("HeatSolver/temperature", temperature_field);
}

QuadraticPulse::~QuadraticPulse(void) {}

std::function<double(const mfem::Vector&)>
QuadraticPulse::SelectInitialConditions(void) const {
  const std::string bc_config =
      parser_m["SimulationInitializer/QuadraticPulse/boundary_configuration"];
  const double domain_length = parser_m["Mesh/gen_bux"].get<double>() -
                               parser_m["Mesh/gen_blx"].get<double>();
  const double domain_height = parser_m["Mesh/gen_buy"].get<double>() -
                               parser_m["Mesh/gen_bly"].get<double>();
  const double amplitude =
      parser_m["SimulationInitializer/QuadraticPulse/pulse_amplitude"]
          .get<double>();
  const std::size_t approximation_terms =
      parser_m["SimulationInitializer/QuadraticPulse/approximation_terms"]
          .get<std::size_t>();
  const double rotation = parser_m["Mesh/rotation"].get<double>();
  const double mesh_rotation_rad = rotation * M_PI / 180.0;
  if (bc_config == "DDDN") {
    return [=](const mfem::Vector& position) {
      std::array<double, 2> rp;
      rp[0] = std::cos(-mesh_rotation_rad) * position[0] -
              std::sin(-mesh_rotation_rad) * position[1];
      rp[1] = std::sin(-mesh_rotation_rad) * position[0] +
              std::cos(-mesh_rotation_rad) * position[1];
      double result = 0.0;
      for (std::size_t m = 1; m <= approximation_terms; ++m) {
        for (std::size_t n = 1; n <= approximation_terms; ++n) {
          const double dm = static_cast<double>(m);
          const double dn = static_cast<double>(n);
          const double ic_factor =
              ((-1.0 + std::pow(-1.0, dn)) * amplitude *
               (-8.0 * std::pow(domain_length, 2) +
                (-4.0 + std::pow(domain_length, 2)) * std::pow(dn, 2) *
                    std::pow(M_PI, 2)) *
               (4.0 * std::pow(1.0 - 2.0 * dm, 2) * std::pow(M_PI, 2) +
                std::pow(domain_height, 2) *
                    (32.0 + (-1.0 + 2.0 * dm) * M_PI *
                                (8.0 * std::pow(-1.0, dm) + M_PI -
                                 2.0 * dm * M_PI)))) /
              (2.0 * std::pow(-1.0 + 2.0 * dm, 3) * std::pow(dn, 3) *
               std::pow(M_PI, 6));
          result += ic_factor * std::sin(dn * M_PI * rp[0] / domain_length) *
                    std::sin((dm - 0.5) * M_PI * rp[1] / domain_height);
        }
      }
      return result;
    };
  } else if (bc_config == "DDND") {
    return [=](const mfem::Vector& position) {
      std::array<double, 2> rp;
      rp[0] = std::cos(-mesh_rotation_rad) * position[0] -
              std::sin(-mesh_rotation_rad) * position[1];
      rp[1] = std::sin(-mesh_rotation_rad) * position[0] +
              std::cos(-mesh_rotation_rad) * position[1];
      double result = 0.0;
      for (std::size_t m = 1; m <= approximation_terms; ++m) {
        for (std::size_t n = 1; n <= approximation_terms; ++n) {
          const double dm = static_cast<double>(m);
          const double dn = static_cast<double>(n);
          const double ic_factor =
              ((-1.0 + std::pow(-1.0, dn)) * amplitude *
               (-8.0 * std::pow(domain_length, 2) +
                (-4.0 + std::pow(domain_length, 2)) * std::pow(dn, 2) *
                    std::pow(M_PI, 2)) *
               (8.0 * std::pow(domain_height, 2) * (1.0 - 2.0 * dm) * M_PI +
                std::pow(-1.0, dm) *
                    (-32.0 * std::pow(domain_height, 2) +
                     (-4.0 + std::pow(domain_height, 2)) *
                         std::pow(1.0 - 2.0 * dm, 2) * std::pow(M_PI, 2)))) /
              (2.0 * std::pow(-1.0 + 2.0 * dm, 3) * std::pow(dn, 3) *
               std::pow(M_PI, 6));
          result += ic_factor * std::sin(dn * M_PI * rp[0] / domain_length) *
                    std::cos((dm - 0.5) * M_PI * rp[1] / domain_height);
        }
      }
      return result;
    };
  } else if (bc_config == "DNDD") {
    return [=](const mfem::Vector& position) {
      std::array<double, 2> rp;
      rp[0] = std::cos(-mesh_rotation_rad) * position[0] -
              std::sin(-mesh_rotation_rad) * position[1];
      rp[1] = std::sin(-mesh_rotation_rad) * position[0] +
              std::cos(-mesh_rotation_rad) * position[1];
      double result = 0.0;
      for (std::size_t m = 1; m <= approximation_terms; ++m) {
        for (std::size_t n = 1; n <= approximation_terms; ++n) {
          const double dm = static_cast<double>(m);
          const double dn = static_cast<double>(n);
          const double ic_factor =
              ((-1.0 + std::pow(-1.0, dm)) * amplitude *
               (-8.0 * std::pow(domain_height, 2) +
                (-4.0 + std::pow(domain_height, 2)) * std::pow(dm, 2) *
                    std::pow(M_PI, 2)) *
               (4.0 * std::pow(1.0 - 2.0 * dn, 2) * std::pow(M_PI, 2) +
                std::pow(domain_length, 2) *
                    (32.0 + (-1.0 + 2.0 * dn) * M_PI *
                                (8.0 * std::pow(-1.0, dn) + M_PI -
                                 2.0 * dn * M_PI)))) /
              (2.0 * std::pow(-1.0 + 2.0 * dn, 3) * std::pow(dm, 3) *
               std::pow(M_PI, 6));
          result += ic_factor *
                    std::sin((dn - 0.5) * M_PI * rp[0] / domain_length) *
                    std::sin(dm * M_PI * rp[1] / domain_height);
        }
      }
      return result;
    };
  } else if (bc_config == "NDDD") {
    return [=](const mfem::Vector& position) {
      std::array<double, 2> rp;
      rp[0] = std::cos(-mesh_rotation_rad) * position[0] -
              std::sin(-mesh_rotation_rad) * position[1];
      rp[1] = std::sin(-mesh_rotation_rad) * position[0] +
              std::cos(-mesh_rotation_rad) * position[1];
      double result = 0.0;
      for (std::size_t m = 1; m <= approximation_terms; ++m) {
        for (std::size_t n = 1; n <= approximation_terms; ++n) {
          const double dm = static_cast<double>(m);
          const double dn = static_cast<double>(n);
          const double ic_factor =
              ((-1.0 + std::pow(-1.0, dm)) * amplitude *
               (-8.0 * std::pow(domain_height, 2) +
                (-4.0 + std::pow(domain_height, 2)) * std::pow(dm, 2) *
                    std::pow(M_PI, 2)) *
               (8.0 * std::pow(domain_length, 2) * (1.0 - 2.0 * dn) * M_PI +
                std::pow(-1.0, dn) *
                    (-32.0 * std::pow(domain_length, 2) +
                     (-4.0 + std::pow(domain_length, 2)) *
                         std::pow(1.0 - 2.0 * dn, 2) * std::pow(M_PI, 2)))) /
              (2.0 * std::pow(-1.0 + 2.0 * dn, 3) * std::pow(dm, 3) *
               std::pow(M_PI, 6));
          result += ic_factor *
                    std::cos((dn - 0.5) * M_PI * rp[0] / domain_length) *
                    std::sin(dm * M_PI * rp[1] / domain_height);
        }
      }
      return result;
    };
  } else if (bc_config == "DDNN") {
    return [=](const mfem::Vector& position) {
      std::array<double, 2> rp;
      rp[0] = std::cos(-mesh_rotation_rad) * position[0] -
              std::sin(-mesh_rotation_rad) * position[1];
      rp[1] = std::sin(-mesh_rotation_rad) * position[0] +
              std::cos(-mesh_rotation_rad) * position[1];
      double result = 0.0;
      for (std::size_t m = 1; m <= approximation_terms; ++m) {
        for (std::size_t n = 1; n <= approximation_terms; ++n) {
          const double dm = static_cast<double>(m);
          const double dn = static_cast<double>(n);
          const double ic_factor =
              -(((1.0 + std::pow(-1.0, m)) * (-1.0 + std::pow(-1.0, dn)) *
                 amplitude * std::pow(domain_height, 2) *
                 (-8.0 * std::pow(domain_length, 2) +
                  (-4.0 + std::pow(domain_length, 2)) * std::pow(dn, 2) *
                      std::pow(M_PI, 2))) /
                (std::pow(dm, 2) * std::pow(dn, 3) * std::pow(M_PI, 5)));
          result += ic_factor * std::sin(dn * M_PI * rp[0] / domain_length) *
                    std::cos(dm * M_PI * rp[1] / domain_height);
        }
      }
      return result;
    };
  } else {
    DEBUG_ASSERT(
        false, global_assert{}, DebugLevel::ALWAYS{},
        "Unknown boundary condition configuration of \"" + bc_config + "\"");
  }
  return [](const mfem::Vector& position) { return 0.0; };
}

}  // namespace chyps
