// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>
#include "tests/unit/parallel/mpi_session.hpp"

#include <array>
#include <cmath>

#include <mfem/mfem.hpp>
#include <nlohmann_json/json.hpp>

#include "tests/helper/command_line_input.hpp"
#include "tests/helper/compute_error.hpp"
#include "tests/helper/convergence_runner.hpp"

#include "chyps/debug_assert.hpp"
#include "chyps/io.hpp"
#include "chyps/logger.hpp"
#include "chyps/simulation.hpp"
#include "chyps/simulation_initializer.hpp"

namespace {

using namespace chyps;

auto TopLambda(const nlohmann::json& a_input_file, const double a_coefficient,
               const double a_pulse_amplitude,
               const std::size_t a_approximation_terms) {
  const double domain_length = a_input_file["Mesh"]["gen_bux"].get<double>() -
                               a_input_file["Mesh"]["gen_blx"].get<double>();
  const double domain_height = a_input_file["Mesh"]["gen_buy"].get<double>() -
                               a_input_file["Mesh"]["gen_bly"].get<double>();
  const double rotation = a_input_file["Mesh"]["rotation"].get<double>();
  const double mesh_rotation_rad = rotation * M_PI / 180.0;

  return [=](const double* a_position, const double a_time) {
    std::array<double, 2> rp;
    rp[0] = std::cos(-mesh_rotation_rad) * a_position[0] -
            std::sin(-mesh_rotation_rad) * a_position[1];
    rp[1] = std::sin(-mesh_rotation_rad) * a_position[0] +
            std::cos(-mesh_rotation_rad) * a_position[1];
    double result = 0.0;
    for (std::size_t m = 1; m <= a_approximation_terms; ++m) {
      for (std::size_t n = 1; n <= a_approximation_terms; ++n) {
        const double dm = static_cast<double>(m);
        const double dn = static_cast<double>(n);
        const double ic_factor =
            ((-1.0 + std::pow(-1.0, dn)) * a_pulse_amplitude *
             (-8.0 * std::pow(domain_length, 2) +
              (-4.0 + std::pow(domain_length, 2)) * std::pow(dn, 2) *
                  std::pow(M_PI, 2)) *
             (4.0 * std::pow(1.0 - 2.0 * dm, 2) * std::pow(M_PI, 2) +
              std::pow(domain_height, 2) *
                  (32.0 +
                   (-1.0 + 2.0 * dm) * M_PI *
                       (8.0 * std::pow(-1.0, dm) + M_PI - 2.0 * dm * M_PI)))) /
            (2.0 * std::pow(-1.0 + 2.0 * dm, 3) * std::pow(dn, 3) *
             std::pow(M_PI, 6));
        result += ic_factor * std::sin(dn * M_PI * rp[0] / domain_length) *
                  std::sin((dm - 0.5) * M_PI * rp[1] / domain_height) *
                  std::exp(-a_coefficient * M_PI * M_PI * a_time *
                           (std::pow(dm - 0.5, 2) / std::pow(domain_height, 2) +
                            std::pow(dn / domain_length, 2)));
      }
    }
    return result;
  };
}

TEST(ConstantCoefficientLinearHeatEquation, ScalarTop) {
  std::string file_name =
      "tests/verification/data/"
      "homogeneous_constant_coefficient_linear_heat_hn_top.json";
  static constexpr int number_of_refinements = 3;

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  const auto coefficient =
      input_file["HeatSolver"]["Conductivity"]["ConstantScalar"].get<double>();

  auto solution_lambda = TopLambda(input_file, coefficient, 1.0, 1);
  auto initial_condition_lambda = [=](const double* a_position) {
    return solution_lambda(a_position, 0.0);
  };

  FunctionInitializers initializers;
  initializers.AddScalarPositionFunction("HeatSolver/temperature",
                                         initial_condition_lambda);

  auto test_result = ConvergenceRunner(
      file_name, solution_lambda, number_of_refinements, true, &initializers);
}

TEST(ConstantCoefficientLinearHeatEquation, MatrixTop) {
  std::string file_name =
      "tests/verification/data/"
      "homogeneous_constant_coefficient_linear_heat_hn_top_matrix.json";
  static constexpr int number_of_refinements = 3;

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  const auto coefficient_vec =
      input_file["HeatSolver"]["Conductivity"]["ConstantMatrix"]
          .get<std::vector<double>>();
  const auto coefficient = coefficient_vec[0];

  auto solution_lambda = TopLambda(input_file, coefficient, 1.0, 1);
  auto initial_condition_lambda = [=](const double* a_position) {
    return solution_lambda(a_position, 0.0);
  };

  FunctionInitializers initializers;
  initializers.AddScalarPositionFunction("HeatSolver/temperature",
                                         initial_condition_lambda);

  auto test_result = ConvergenceRunner(
      file_name, solution_lambda, number_of_refinements, true, &initializers);
}

TEST(ConstantCoefficientLinearHeatEquation, AttributeScalarTop) {
  std::string file_name =
      "tests/verification/data/"
      "homogeneous_constant_coefficient_linear_heat_hn_top_attr_scalar.json";
  static constexpr int number_of_refinements = 3;

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  const auto coefficient_vec =
      input_file["HeatSolver"]["Conductivity"]["MaterialVaryingScalar"]
          .get<std::vector<double>>();
  const auto coefficient = coefficient_vec[0];

  auto solution_lambda = TopLambda(input_file, coefficient, 1.0, 1);
  auto initial_condition_lambda = [=](const double* a_position) {
    return solution_lambda(a_position, 0.0);
  };

  FunctionInitializers initializers;
  initializers.AddScalarPositionFunction("HeatSolver/temperature",
                                         initial_condition_lambda);

  auto test_result = ConvergenceRunner(
      file_name, solution_lambda, number_of_refinements, true, &initializers);
}

TEST(ConstantCoefficientLinearHeatEquation, AttributeMatrixTop) {
  std::string file_name =
      "tests/verification/data/"
      "homogeneous_constant_coefficient_linear_heat_hn_top_attr_matrix.json";
  static constexpr int number_of_refinements = 3;

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  const auto coefficient_vec =
      input_file["HeatSolver"]["Conductivity"]["MaterialVaryingMatrix"]["1"]
          .get<std::vector<double>>();
  const auto coefficient = coefficient_vec[0];

  auto solution_lambda = TopLambda(input_file, coefficient, 1.0, 1);
  auto initial_condition_lambda = [=](const double* a_position) {
    return solution_lambda(a_position, 0.0);
  };

  FunctionInitializers initializers;
  initializers.AddScalarPositionFunction("HeatSolver/temperature",
                                         initial_condition_lambda);

  auto test_result = ConvergenceRunner(
      file_name, solution_lambda, number_of_refinements, true, &initializers);
}

TEST(ConstantCoefficientLinearHeatEquation, ElementScalarTop) {
  std::string file_name =
      "tests/verification/data/"
      "homogeneous_constant_coefficient_linear_heat_hn_top_element_scalar.json";
  static constexpr int number_of_refinements = 3;

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  const auto coefficient =
      input_file["SimulationInitializer"]["FieldInitialization"]
                ["HeatSolver/Conductivity/kappa"]["Arguments"]["value"]
                    .get<double>();

  auto solution_lambda = TopLambda(input_file, coefficient, 1.0, 1);
  auto initial_condition_lambda = [=](const double* a_position) {
    return solution_lambda(a_position, 0.0);
  };

  FunctionInitializers initializers;
  initializers.AddScalarPositionFunction("HeatSolver/temperature",
                                         initial_condition_lambda);

  auto test_result = ConvergenceRunner(
      file_name, solution_lambda, number_of_refinements, true, &initializers);
}

TEST(ConstantCoefficientLinearHeatEquation, ElementMatrixTop) {
  std::string file_name =
      "tests/verification/data/"
      "homogeneous_constant_coefficient_linear_heat_hn_top_element_matrix.json";
  static constexpr int number_of_refinements = 3;

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  const auto coefficient_vec =
      input_file["SimulationInitializer"]["FieldInitialization"]
                ["HeatSolver/Conductivity/kappa"]["Arguments"]["value"]
                    .get<std::vector<double>>();
  const auto coefficient = coefficient_vec[0];

  auto solution_lambda = TopLambda(input_file, coefficient, 1.0, 1);
  auto initial_condition_lambda = [=](const double* a_position) {
    return solution_lambda(a_position, 0.0);
  };

  FunctionInitializers initializers;
  initializers.AddScalarPositionFunction("HeatSolver/temperature",
                                         initial_condition_lambda);

  auto test_result = ConvergenceRunner(
      file_name, solution_lambda, number_of_refinements, true, &initializers);
}

TEST(ConstantCoefficientLinearHeatEquation, ScalarBot) {
  std::string file_name =
      "tests/verification/data/"
      "homogeneous_constant_coefficient_linear_heat_hn_bot.json";
  static constexpr int number_of_refinements = 3;

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  const double domain_length = input_file["Mesh"]["gen_bux"].get<double>() -
                               input_file["Mesh"]["gen_blx"].get<double>();
  const double domain_height = input_file["Mesh"]["gen_buy"].get<double>() -
                               input_file["Mesh"]["gen_bly"].get<double>();
  const double rotation = input_file["Mesh"]["rotation"].get<double>();
  const double mesh_rotation_rad = rotation * M_PI / 180.0;

  static constexpr double amplitude = 1.0;
  static constexpr std::size_t approximation_terms = 1;

  const auto coefficient =
      input_file["HeatSolver"]["Conductivity"]["ConstantScalar"].get<double>();

  auto solution_lambda = [=](const double* a_position, const double a_time) {
    std::array<double, 2> rp;
    rp[0] = std::cos(-mesh_rotation_rad) * a_position[0] -
            std::sin(-mesh_rotation_rad) * a_position[1];
    rp[1] = std::sin(-mesh_rotation_rad) * a_position[0] +
            std::cos(-mesh_rotation_rad) * a_position[1];
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
                  std::cos((dm - 0.5) * M_PI * rp[1] / domain_height) *
                  std::exp(-coefficient * M_PI * M_PI * a_time *
                           (std::pow(dm - 0.5, 2) / std::pow(domain_height, 2) +
                            std::pow(dn / domain_length, 2)));
      }
    }
    return result;
  };

  auto initial_condition_lambda = [=](const double* a_position) {
    return solution_lambda(a_position, 0.0);
  };

  FunctionInitializers initializers;
  initializers.AddScalarPositionFunction("HeatSolver/temperature",
                                         initial_condition_lambda);

  auto test_result = ConvergenceRunner(
      file_name, solution_lambda, number_of_refinements, true, &initializers);
}

TEST(ConstantCoefficientLinearHeatEquation, ScalarRight) {
  std::string file_name =
      "tests/verification/data/"
      "homogeneous_constant_coefficient_linear_heat_hn_right.json";
  static constexpr int number_of_refinements = 3;

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  const double domain_length = input_file["Mesh"]["gen_bux"].get<double>() -
                               input_file["Mesh"]["gen_blx"].get<double>();
  const double domain_height = input_file["Mesh"]["gen_buy"].get<double>() -
                               input_file["Mesh"]["gen_bly"].get<double>();
  const double rotation = input_file["Mesh"]["rotation"].get<double>();
  const double mesh_rotation_rad = rotation * M_PI / 180.0;

  static constexpr double amplitude = 1.0;
  static constexpr std::size_t approximation_terms = 1;

  const auto coefficient =
      input_file["HeatSolver"]["Conductivity"]["ConstantScalar"].get<double>();

  auto solution_lambda = [=](const double* a_position, const double a_time) {
    std::array<double, 2> rp;
    rp[0] = std::cos(-mesh_rotation_rad) * a_position[0] -
            std::sin(-mesh_rotation_rad) * a_position[1];
    rp[1] = std::sin(-mesh_rotation_rad) * a_position[0] +
            std::cos(-mesh_rotation_rad) * a_position[1];
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
                  (32.0 +
                   (-1.0 + 2.0 * dn) * M_PI *
                       (8.0 * std::pow(-1.0, dn) + M_PI - 2.0 * dn * M_PI)))) /
            (2.0 * std::pow(-1.0 + 2.0 * dn, 3) * std::pow(dm, 3) *
             std::pow(M_PI, 6));
        result += ic_factor *
                  std::sin((dn - 0.5) * M_PI * rp[0] / domain_length) *
                  std::sin(dm * M_PI * rp[1] / domain_height) *
                  std::exp(-coefficient * M_PI * M_PI * a_time *
                           (std::pow(dn - 0.5, 2) / std::pow(domain_length, 2) +
                            std::pow(dm / domain_height, 2)));
      }
    }
    return result;
  };

  auto initial_condition_lambda = [=](const double* a_position) {
    return solution_lambda(a_position, 0.0);
  };

  FunctionInitializers initializers;
  initializers.AddScalarPositionFunction("HeatSolver/temperature",
                                         initial_condition_lambda);

  auto test_result = ConvergenceRunner(
      file_name, solution_lambda, number_of_refinements, true, &initializers);
}

TEST(ConstantCoefficientLinearHeatEquation, ScalarLeft) {
  std::string file_name =
      "tests/verification/data/"
      "homogeneous_constant_coefficient_linear_heat_hn_left.json";
  static constexpr int number_of_refinements = 3;

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  const double domain_length = input_file["Mesh"]["gen_bux"].get<double>() -
                               input_file["Mesh"]["gen_blx"].get<double>();
  const double domain_height = input_file["Mesh"]["gen_buy"].get<double>() -
                               input_file["Mesh"]["gen_bly"].get<double>();
  const double rotation = input_file["Mesh"]["rotation"].get<double>();
  const double mesh_rotation_rad = rotation * M_PI / 180.0;

  static constexpr double amplitude = 1.0;
  static constexpr std::size_t approximation_terms = 1;

  const auto coefficient =
      input_file["HeatSolver"]["Conductivity"]["ConstantScalar"].get<double>();

  auto solution_lambda = [=](const double* a_position, const double a_time) {
    std::array<double, 2> rp;
    rp[0] = std::cos(-mesh_rotation_rad) * a_position[0] -
            std::sin(-mesh_rotation_rad) * a_position[1];
    rp[1] = std::sin(-mesh_rotation_rad) * a_position[0] +
            std::cos(-mesh_rotation_rad) * a_position[1];
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
                  std::sin(dm * M_PI * rp[1] / domain_height) *
                  std::exp(-coefficient * M_PI * M_PI * a_time *
                           (std::pow(dm - 0.5, 2) / std::pow(domain_length, 2) +
                            std::pow(dn / domain_height, 2)));
      }
    }
    return result;
  };

  auto initial_condition_lambda = [=](const double* a_position) {
    return solution_lambda(a_position, 0.0);
  };

  FunctionInitializers initializers;
  initializers.AddScalarPositionFunction("HeatSolver/temperature",
                                         initial_condition_lambda);

  auto test_result = ConvergenceRunner(
      file_name, solution_lambda, number_of_refinements, true, &initializers);
}

TEST(ConstantCoefficientLinearHeatEquation, HomogeneousAndTensorCoefficient) {
  std::string file_name =
      "tests/verification/data/"
      "homogeneous_tensor_coefficient.json";
  static constexpr int number_of_refinements = 3;

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  const double domain_length = input_file["Mesh"]["gen_bux"].get<double>() -
                               input_file["Mesh"]["gen_blx"].get<double>();
  const double domain_height = input_file["Mesh"]["gen_buy"].get<double>() -
                               input_file["Mesh"]["gen_bly"].get<double>();
  const double rotation = input_file["Mesh"]["rotation"].get<double>();
  const double mesh_rotation_rad = rotation * M_PI / 180.0;

  static constexpr double amplitude = 1.0;
  static constexpr std::size_t approximation_terms = 2;

  const std::vector<double> tensor_kappa =
      input_file["HeatSolver"]["Conductivity"]["ConstantMatrix"]
          .get<std::vector<double>>();
  const double coefficient_x = tensor_kappa[0];
  const double coefficient_y = tensor_kappa[3];

  auto solution_lambda = [=](const double* a_position, const double a_time) {
    std::array<double, 2> rp;
    rp[0] = std::cos(-mesh_rotation_rad) * a_position[0] -
            std::sin(-mesh_rotation_rad) * a_position[1];
    rp[1] = std::sin(-mesh_rotation_rad) * a_position[0] +
            std::cos(-mesh_rotation_rad) * a_position[1];
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
                  std::cos(dm * M_PI * rp[1] / domain_height) *
                  std::exp(-M_PI * M_PI * a_time *
                           (coefficient_x * std::pow(dn / domain_length, 2) +
                            coefficient_y * std::pow(dm / domain_height, 2)));
      }
    }
    return result;
  };

  auto initial_condition_lambda = [=](const double* a_position) {
    return solution_lambda(a_position, 0.0);
  };

  FunctionInitializers initializers;
  initializers.AddScalarPositionFunction("HeatSolver/temperature",
                                         initial_condition_lambda);

  auto test_result = ConvergenceRunner(
      file_name, solution_lambda, number_of_refinements, true, &initializers);
}

TEST(ConstantCoefficientLinearHeatEquation, CooledRod) {
  std::string file_name =
      "tests/verification/data/"
      "neumann_cooled_rod.json";
  static constexpr int number_of_refinements = 2;

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  const double domain_length = input_file["Mesh"]["gen_bux"].get<double>() -
                               input_file["Mesh"]["gen_blx"].get<double>();
  DEBUG_ASSERT(std::fabs(domain_length - 1.0) < 1.0e-15, global_assert{},
               DebugLevel::ALWAYS{}, "Test requires domain length of 1.0");

  const auto coefficient =
      input_file["HeatSolver"]["Conductivity"]["ConstantScalar"].get<double>();
  DEBUG_ASSERT(std::fabs(coefficient - 1.0) < 1.0e-15, global_assert{},
               DebugLevel::ALWAYS{},
               "Test requires thermal coefficient of 1.0");

  static constexpr int approximation_terms = 1;

  auto solution_lambda = [=](const double* a_position, const double a_time) {
    double result = 0.0;
    for (int n = 1; n <= approximation_terms; ++n) {
      const double dn = static_cast<double>(n);
      const double ic_factor = 8.0 / std::pow(M_PI * (1.0 - 2.0 * dn), 2);
      result += ic_factor * std::cos((dn - 0.5) * M_PI * a_position[0]) *
                std::exp(-std::pow(dn - 0.5, 2) * M_PI * M_PI * a_time);
    }
    return result + 24.0 + a_position[0];
  };

  auto initial_condition_lambda = [=](const double* a_position) {
    return solution_lambda(a_position, 0.0);
  };

  FunctionInitializers initializers;
  initializers.AddScalarPositionFunction("HeatSolver/temperature",
                                         initial_condition_lambda);

  auto test_result = ConvergenceRunner(
      file_name, solution_lambda, number_of_refinements, true, &initializers);
}

}  // namespace
