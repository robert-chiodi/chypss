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

namespace quadratic_pulse {

static std::function<double(const mfem::Vector&)> SelectInitialConditions(
    const nlohmann::json& a_json_object, const InputParser& a_full_parser);

void AddParserOptions(InputParser& a_parser) {
  a_parser.AddOption(
      "SimulationInitializer/Initializers/quadratic_pulse/"
      "pulse_amplitude",
      "Amplitude of quadratic pulse (value at center)", 1.0);
  a_parser.AddOption(
      "SimulationInitializer/Initializers/quadratic_pulse/"
      "approximation_terms",
      "Number of terms to use in pulse approximation", 1);
  a_parser.AddOption(
      "SimulationInitializer/Initializers/quadratic_pulse/"
      "boundary_configuration",
      "Configuration of homogeneous boundary conditions. Supplied by string of "
      "4 characters, D for Dirichlet and N for Neumann, correspoding to tag "
      "locations, (e.g.,\"HHNN\")");
}

void InitializeData(const nlohmann::json& a_json_object,
                    const InputParser& a_full_parser,
                    mfem::ParFiniteElementSpace& a_finite_element_space,
                    mfem::Vector& a_data) {
  mfem::ParGridFunction grid_function(&a_finite_element_space);

  std::function<double(const mfem::Vector&)> value_setter =
      SelectInitialConditions(a_json_object, a_full_parser);

  mfem::FunctionCoefficient value_setter_function(value_setter);
  grid_function.ProjectCoefficient(value_setter_function);
  grid_function.GetTrueDofs(a_data);
}

std::function<double(const mfem::Vector&)> SelectInitialConditions(
    const nlohmann::json& a_json_object, const InputParser& a_full_parser) {
  const auto bc_config =
      a_json_object["boundary_configuration"].get<std::string>();
  const double amplitude = a_json_object["pulse_amplitude"].get<double>();
  const std::size_t approximation_terms =
      a_json_object["approximation_terms"].get<std::size_t>();

  const double domain_length = a_full_parser["Mesh/gen_bux"].get<double>() -
                               a_full_parser["Mesh/gen_blx"].get<double>();
  const double domain_height = a_full_parser["Mesh/gen_buy"].get<double>() -
                               a_full_parser["Mesh/gen_bly"].get<double>();
  const double rotation = a_full_parser["Mesh/rotation"].get<double>();
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
              -(((1.0 + std::pow(-1.0, dm)) * (-1.0 + std::pow(-1.0, dn)) *
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

}  // namespace quadratic_pulse

}  // namespace chyps
