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

#include "chyps/debug_assert.hpp"
#include "chyps/io.hpp"
#include "chyps/logger.hpp"
#include "chyps/simulation.hpp"
#include "chyps/simulation_initializer.hpp"

namespace {

using namespace chyps;

template <class SolutionFunctor>
int ConvergenceRunner(const std::string& a_input_name,
                      const std::string& a_configuration,
                      const SolutionFunctor& a_analytical_solution_lambda,
                      const int a_number_of_refines) {
  std::vector<std::string> init_input_string;
  init_input_string.push_back("Executable_name");
  init_input_string.push_back(a_input_name);
  init_input_string.push_back(a_configuration);
  init_input_string.push_back("-Mesh/parallel_refine");
  init_input_string.push_back("0");

  std::vector<std::string> chyps_input_string;
  chyps_input_string.push_back("Executable_name");
  chyps_input_string.push_back(a_input_name);
  chyps_input_string.push_back("-Mesh/parallel_refine");
  chyps_input_string.push_back("0");

  std::ifstream myfile(a_input_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  const double final_time = input_file["Simulation"]["end_time"].get<double>();
  const std::string& solution_file_name =
      input_file["Simulation"]["out_data"].get<std::string>();

  std::vector<std::array<double, 4>> error_metrics(a_number_of_refines + 1);
  for (int r = 0; r <= a_number_of_refines; ++r) {
    // Create file with initial conditions
    init_input_string[4] = std::to_string(r);
    auto init_input_char = FakeCommandLineInput(init_input_string);
    int argc = static_cast<int>(init_input_char.size());
    char** argv = init_input_char.data();
    chyps::SimulationInitializer(argc, argv, *mpi_session, SpdlogLevel::OFF);

    // Perform simulation
    chyps_input_string[3] = std::to_string(r);
    auto chyps_input_char = FakeCommandLineInput(chyps_input_string);
    argc = static_cast<int>(chyps_input_char.size());
    argv = chyps_input_char.data();
    chyps::main(argc, argv, *mpi_session, SpdlogLevel::OFF);

    IO run_file(*mpi_session, "Run");
    run_file.SetRead(solution_file_name);
    std::vector<double> temperature_field;
    run_file.GetImmediateBlock("HeatSolver/temperature", temperature_field);
    std::vector<double> point_field;
    run_file.GetImmediateBlock("vertices", point_field);
    DEBUG_ASSERT(point_field.size() / 2 == temperature_field.size(),
                 global_assert{}, DebugLevel::CHEAP{});
    std::vector<uint32_t> nvert(1);
    run_file.GetImmediateBlock("NumOfVertices", {0}, {1}, nvert.data());
    DEBUG_ASSERT(nvert[0] == point_field.size() / 2, global_assert{},
                 DebugLevel::CHEAP{});

    std::vector<double> correct_solution(temperature_field.size());
    for (std::size_t n = 0; n < nvert[0]; ++n) {
      const double* position_2d = &(point_field[2 * n]);
      correct_solution[n] =
          a_analytical_solution_lambda(position_2d, final_time);
    }

    const double global_l1 = GlobalL1Diff_Normalized(
        temperature_field, correct_solution, *mpi_session);
    const double global_l2 = GlobalL2Diff_Normalized(
        temperature_field, correct_solution, *mpi_session);
    const double global_linf =
        GlobalLinfDiff(temperature_field, correct_solution, *mpi_session);
    error_metrics[r][0] = static_cast<double>(nvert[0]);
    error_metrics[r][1] = global_l1;
    error_metrics[r][2] = global_l2;
    error_metrics[r][3] = global_linf;

    DeleteCommandLineInput(init_input_char);
    DeleteCommandLineInput(chyps_input_char);
  }

  if (mpi_session->IAmRoot()) {
    printf("%16s %16s %16s %16s %16s %16s %16s\n", "Resolution", "L1_err",
           "L1_conv", "L2_err", "L2_conv", "Linf_err", "Linf_conv");
    for (int r = 0; r <= a_number_of_refines; ++r) {
      if (r == 0) {
        printf("%16d %16.8E %16s %16.8E %16s %16.8E %16s\n",
               static_cast<int>(error_metrics[r][0]), error_metrics[r][1],
               "N/A", error_metrics[r][2], "N/A", error_metrics[r][3], "N/A");
      } else {
        const double l1_conv =
            log(error_metrics[r - 1][1] / error_metrics[r][1] /
                (error_metrics[r - 1][0] / error_metrics[r][0]));
        const double l2_conv =
            log(error_metrics[r - 1][2] / error_metrics[r][2] /
                (error_metrics[r - 1][0] / error_metrics[r][0]));
        const double linf_conv =
            log(error_metrics[r - 1][3] / error_metrics[r][3] /
                (error_metrics[r - 1][0] / error_metrics[r][0]));
        printf("%16d %16.8E %16.8f %16.8E %16.8f %16.8E %16.8f\n",
               static_cast<int>(error_metrics[r][0]), error_metrics[r][1],
               l1_conv, error_metrics[r][2], l2_conv, error_metrics[r][3],
               linf_conv);
      }
    }
  }
  return 1;  // Come up with actual enforceable error measure
}
TEST(LinearHeatEquation, HomogeneousAndConstantCoefficientTop) {
  std::string file_name =
      "tests/verification/data/"
      "homogeneous_constant_coefficient_linear_heat_hn_top.json";
  static constexpr int number_of_refinements = 5;

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  const double domain_length = input_file["Mesh"]["gen_bux"].get<double>() -
                               input_file["Mesh"]["gen_blx"].get<double>();
  const double domain_height = input_file["Mesh"]["gen_buy"].get<double>() -
                               input_file["Mesh"]["gen_bly"].get<double>();
  const double amplitude =
      input_file["SimulationInitializer"]["QuadraticPulse"]["pulse_amplitude"]
          .get<double>();
  const std::size_t approximation_terms =
      input_file["SimulationInitializer"]["QuadraticPulse"]
                ["approximation_terms"]
                    .get<std::size_t>();
  const double rotation = input_file["Mesh"]["rotation"].get<double>();
  const double mesh_rotation_rad = rotation * M_PI / 180.0;

  const std::vector<double> tensor_kappa =
      input_file["HeatSolver"]["kappa"].get<std::vector<double>>();
  const double coefficient = tensor_kappa[0];

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
             (4.0 * std::pow(1.0 - 2.0 * dm, 2) * std::pow(M_PI, 2) +
              std::pow(domain_height, 2) *
                  (32.0 +
                   (-1.0 + 2.0 * dm) * M_PI *
                       (8.0 * std::pow(-1.0, dm) + M_PI - 2.0 * dm * M_PI)))) /
            (2.0 * std::pow(-1.0 + 2.0 * dm, 3) * std::pow(dn, 3) *
             std::pow(M_PI, 6));
        result += ic_factor * std::sin(dn * M_PI * rp[0] / domain_length) *
                  std::sin((dm - 0.5) * M_PI * rp[1] / domain_height) *
                  std::exp(-coefficient * M_PI * M_PI * a_time *
                           (std::pow(dm - 0.5, 2) / std::pow(domain_height, 2) +
                            std::pow(dn / domain_length, 2)));
      }
    }
    return result;
  };

  int test_result = ConvergenceRunner(file_name, "quadratic_pulse",
                                      solution_lambda, number_of_refinements);
  EXPECT_EQ(test_result, 1);
}

TEST(LinearHeatEquation, HomogeneousAndConstantCoefficientBot) {
  std::string file_name =
      "tests/verification/data/"
      "homogeneous_constant_coefficient_linear_heat_hn_bot.json";
  static constexpr int number_of_refinements = 5;

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  const double domain_length = input_file["Mesh"]["gen_bux"].get<double>() -
                               input_file["Mesh"]["gen_blx"].get<double>();
  const double domain_height = input_file["Mesh"]["gen_buy"].get<double>() -
                               input_file["Mesh"]["gen_bly"].get<double>();
  const double amplitude =
      input_file["SimulationInitializer"]["QuadraticPulse"]["pulse_amplitude"]
          .get<double>();
  const std::size_t approximation_terms =
      input_file["SimulationInitializer"]["QuadraticPulse"]
                ["approximation_terms"]
                    .get<std::size_t>();
  const double rotation = input_file["Mesh"]["rotation"].get<double>();
  const double mesh_rotation_rad = rotation * M_PI / 180.0;

  const std::vector<double> tensor_kappa =
      input_file["HeatSolver"]["kappa"].get<std::vector<double>>();
  const double coefficient = tensor_kappa[0];

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

  int test_result = ConvergenceRunner(file_name, "quadratic_pulse",
                                      solution_lambda, number_of_refinements);
  EXPECT_EQ(test_result, 1);
}

TEST(LinearHeatEquation, HomogeneousAndConstantCoefficientRight) {
  std::string file_name =
      "tests/verification/data/"
      "homogeneous_constant_coefficient_linear_heat_hn_right.json";
  static constexpr int number_of_refinements = 5;

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  const double domain_length = input_file["Mesh"]["gen_bux"].get<double>() -
                               input_file["Mesh"]["gen_blx"].get<double>();
  const double domain_height = input_file["Mesh"]["gen_buy"].get<double>() -
                               input_file["Mesh"]["gen_bly"].get<double>();
  const double amplitude =
      input_file["SimulationInitializer"]["QuadraticPulse"]["pulse_amplitude"]
          .get<double>();
  const std::size_t approximation_terms =
      input_file["SimulationInitializer"]["QuadraticPulse"]
                ["approximation_terms"]
                    .get<std::size_t>();
  const double rotation = input_file["Mesh"]["rotation"].get<double>();
  const double mesh_rotation_rad = rotation * M_PI / 180.0;

  const std::vector<double> tensor_kappa =
      input_file["HeatSolver"]["kappa"].get<std::vector<double>>();
  const double coefficient = tensor_kappa[0];

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

  int test_result = ConvergenceRunner(file_name, "quadratic_pulse",
                                      solution_lambda, number_of_refinements);
  EXPECT_EQ(test_result, 1);
}

TEST(LinearHeatEquation, HomogeneousAndConstantCoefficientLeft) {
  std::string file_name =
      "tests/verification/data/"
      "homogeneous_constant_coefficient_linear_heat_hn_left.json";
  static constexpr int number_of_refinements = 5;

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  const double domain_length = input_file["Mesh"]["gen_bux"].get<double>() -
                               input_file["Mesh"]["gen_blx"].get<double>();
  const double domain_height = input_file["Mesh"]["gen_buy"].get<double>() -
                               input_file["Mesh"]["gen_bly"].get<double>();
  const double amplitude =
      input_file["SimulationInitializer"]["QuadraticPulse"]["pulse_amplitude"]
          .get<double>();
  const std::size_t approximation_terms =
      input_file["SimulationInitializer"]["QuadraticPulse"]
                ["approximation_terms"]
                    .get<std::size_t>();
  const double rotation = input_file["Mesh"]["rotation"].get<double>();
  const double mesh_rotation_rad = rotation * M_PI / 180.0;

  const std::vector<double> tensor_kappa =
      input_file["HeatSolver"]["kappa"].get<std::vector<double>>();
  const double coefficient = tensor_kappa[0];

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

  int test_result = ConvergenceRunner(file_name, "quadratic_pulse",
                                      solution_lambda, number_of_refinements);
  EXPECT_EQ(test_result, 1);
}

TEST(LinearHeatEquation, HomogeneousAndTensorCoefficient) {
  std::string file_name =
      "tests/verification/data/"
      "homogeneous_tensor_coefficient.json";
  static constexpr int number_of_refinements = 5;

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  const double domain_length = input_file["Mesh"]["gen_bux"].get<double>() -
                               input_file["Mesh"]["gen_blx"].get<double>();
  const double domain_height = input_file["Mesh"]["gen_buy"].get<double>() -
                               input_file["Mesh"]["gen_bly"].get<double>();
  const double amplitude =
      input_file["SimulationInitializer"]["QuadraticPulse"]["pulse_amplitude"]
          .get<double>();
  const std::size_t approximation_terms =
      input_file["SimulationInitializer"]["QuadraticPulse"]
                ["approximation_terms"]
                    .get<std::size_t>();
  const double rotation = input_file["Mesh"]["rotation"].get<double>();
  const double mesh_rotation_rad = rotation * M_PI / 180.0;

  const std::vector<double> tensor_kappa =
      input_file["HeatSolver"]["kappa"].get<std::vector<double>>();
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

  int test_result = ConvergenceRunner(file_name, "quadratic_pulse",
                                      solution_lambda, number_of_refinements);
  EXPECT_EQ(test_result, 1);
}

TEST(LinearHeatEquation, HomogeneousAndTensorCoefficientRot45) {
  std::string file_name =
      "tests/verification/data/"
      "homogeneous_tensor_coefficient_rotated.json";
  static constexpr int number_of_refinements = 5;

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  const double domain_length = input_file["Mesh"]["gen_bux"].get<double>() -
                               input_file["Mesh"]["gen_blx"].get<double>();
  const double domain_height = input_file["Mesh"]["gen_buy"].get<double>() -
                               input_file["Mesh"]["gen_bly"].get<double>();
  const double amplitude =
      input_file["SimulationInitializer"]["QuadraticPulse"]["pulse_amplitude"]
          .get<double>();
  const std::size_t approximation_terms =
      input_file["SimulationInitializer"]["QuadraticPulse"]
                ["approximation_terms"]
                    .get<std::size_t>();
  const double rotation = input_file["Mesh"]["rotation"].get<double>();
  const double mesh_rotation_rad = rotation * M_PI / 180.0;

  const std::vector<double> tensor_kappa =
      input_file["HeatSolver"]["kappa"].get<std::vector<double>>();
  const double coefficient_x = 0.5;
  const double coefficient_y = 0.5;

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

  int test_result = ConvergenceRunner(file_name, "quadratic_pulse",
                                      solution_lambda, number_of_refinements);
  EXPECT_EQ(test_result, 1);
}

TEST(LinearHeatEquation, CooledRod) {
  std::string file_name =
      "tests/verification/data/"
      "neumann_cooled_rod.json";
  static constexpr int number_of_refinements = 3;

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  const double domain_length = input_file["Mesh"]["gen_bux"].get<double>() -
                               input_file["Mesh"]["gen_blx"].get<double>();
  DEBUG_ASSERT(std::fabs(domain_length - 1.0) < 1.0e-15, global_assert{},
               DebugLevel::ALWAYS{}, "Test requires domain length of 1.0");

  const std::vector<double> tensor_kappa =
      input_file["HeatSolver"]["kappa"].get<std::vector<double>>();
  const double coefficient = tensor_kappa[0];
  DEBUG_ASSERT(std::fabs(coefficient - 1.0) < 1.0e-15, global_assert{},
               DebugLevel::ALWAYS{},
               "Test requires thermal coefficient of 1.0");
  const int approximation_terms =
      input_file["SimulationInitializer"]["CooledRod"]["approximation_terms"]
          .get<int>();

  auto solution_lambda = [=](const double* a_position, const double a_time) {
    double result = 0.0;
    for (std::size_t n = 1; n <= approximation_terms; ++n) {
      const double dn = static_cast<double>(n);
      const double ic_factor = 8.0 / std::pow(M_PI * (1.0 - 2.0 * dn), 2);
      result += ic_factor * std::cos((dn - 0.5) * M_PI * a_position[0]) *
                std::exp(-std::pow(dn - 0.5, 2) * M_PI * M_PI * a_time);
    }
    return result + 24.0 + a_position[0];
  };

  int test_result = ConvergenceRunner(file_name, "cooled_rod", solution_lambda,
                                      number_of_refinements);
  EXPECT_EQ(test_result, 1);
}

}  // namespace
