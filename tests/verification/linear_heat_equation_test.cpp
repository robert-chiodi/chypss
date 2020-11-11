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

#include "chyps/io.hpp"
#include "chyps/logger.hpp"
#include "chyps/simulation.hpp"

namespace {

using namespace chyps;

class DecayingPulse_CC {
 public:
  DecayingPulse_CC(const double a_coefficient, const double a_amplitude,
                   const double a_domain_length, const double a_domain_height)
      : coefficient_m(a_coefficient),
        amplitude_m(a_amplitude),
        domain_length_m(a_domain_length),
        domain_height_m(a_domain_height) {}

  double NeumannTop(const double* a_position, const double a_time) {
    // Use just first term so we are certain thereA are enough terms included
    // (i.e., 1 will be exact).
    const std::size_t m = 1;
    const std::size_t n = 1;
    const auto dn = static_cast<double>(n);
    const auto dm = static_cast<double>(m);
    const double domain_length = domain_length_m;
    const double domain_height = domain_height_m;
    const double amplitude = amplitude_m;
    const double ic_factor =
        (amplitude *
         (-4.0 * std::pow(M_PI, 2) +
          std::pow(domain_length, 2) * (-8.0 + std::pow(M_PI, 2))) *
         (-4.0 * std::pow(M_PI, 2) +
          std::pow(domain_height, 2) *
              (-32.0 + 8.0 * M_PI + std::pow(M_PI, 2)))) /
        std::pow(M_PI, 6);
    return ic_factor * std::sin(dn * M_PI * a_position[0] / domain_length) *
           std::sin((dm - 0.5) * M_PI * a_position[1] / domain_height) *
           std::exp(-coefficient_m * M_PI * M_PI * a_time *
                    (std::pow(dm - 0.5, 2) / std::pow(domain_height, 2) +
                     std::pow(dn / domain_length, 2)));
  }

  double NeumannBot(const double* a_position, const double a_time) {
    // Use just first term so we are certain thereA are enough terms included
    // (i.e., 1 will be exact).
    const std::size_t m = 1;
    const std::size_t n = 1;
    const auto dn = static_cast<double>(n);
    const auto dm = static_cast<double>(m);
    const double domain_length = domain_length_m;
    const double domain_height = domain_height_m;
    const double amplitude = amplitude_m;
    const double ic_factor =
        (amplitude *
         (-4.0 * std::pow(M_PI, 2) +
          std::pow(domain_length, 2) * (-8.0 + std::pow(M_PI, 2))) *
         (-4.0 * std::pow(M_PI, 2) +
          std::pow(domain_height, 2) *
              (-32.0 + 8.0 * M_PI + std::pow(M_PI, 2)))) /
        std::pow(M_PI, 6);
    return ic_factor * std::sin(dn * M_PI * a_position[0] / domain_length) *
           std::cos((dm - 0.5) * M_PI * a_position[1] / domain_height) *
           std::exp(-coefficient_m * M_PI * M_PI * a_time *
                    (std::pow(dm - 0.5, 2) / std::pow(domain_height, 2) +
                     std::pow(dn / domain_length, 2)));
  }

  double NeumannRight(const double* a_position, const double a_time) {
    // Use just first term so we are certain thereA are enough terms included
    // (i.e., 1 will be exact).
    const std::size_t m = 1;
    const std::size_t n = 1;
    const auto dn = static_cast<double>(n);
    const auto dm = static_cast<double>(m);
    const double domain_length = domain_length_m;
    const double domain_height = domain_height_m;
    const double amplitude = amplitude_m;
    // const double ic_factor =
    //     (-2.0 * amplitude *
    //      (-8.0 * std::pow(domain_height, 2) +
    //       (-4.0 + std::pow(domain_height, 2)) * std::pow(dn, 2) *
    //           std::pow(M_PI, 2)) *
    //      (4.0 * std::pow(1.0 - 2.0 * dm, 2) * std::pow(M_PI, 2) +
    //       std::pow(domain_length, 2) * (32.0 + (-1.0 + 2.0 * dm) * M_PI *
    //                                                (8.0 * std::pow(-1.0, m) +
    //                                                 M_PI - 2.0 * dm *
    //                                                 M_PI)))) /
    //     (2.0 * std::pow(-1.0 + 2.0 * dm, 3) * std::pow(dn, 3) *
    //      std::pow(M_PI, 6));
    // return ic_factor * std::sin(dn * M_PI * a_position[1] / domain_height) *
    //        std::sin((dm - 0.5) * M_PI * a_position[0] / domain_length) *
    //        std::exp(-coefficient_m * M_PI * M_PI * a_time *
    //                 (std::pow(dm - 0.5, 2) / std::pow(domain_length, 2) +
    //                  std::pow(dn / domain_height, 2)));

    const double ic_factor =
        (amplitude *
         (-4.0 * std::pow(M_PI, 2) +
          std::pow(domain_height, 2) * (-8.0 + std::pow(M_PI, 2))) *
         (-4.0 * std::pow(M_PI, 2) +
          std::pow(domain_length, 2) *
              (-32.0 + 8.0 * M_PI + std::pow(M_PI, 2)))) /
        std::pow(M_PI, 6);
    return ic_factor * std::sin(dm * M_PI * a_position[1] / domain_height) *
           std::sin((dn - 0.5) * M_PI * a_position[0] / domain_length) *
           std::exp(-coefficient_m * M_PI * M_PI * a_time *
                    (std::pow(dn - 0.5, 2) / std::pow(domain_length, 2) +
                     std::pow(dm / domain_height, 2)));
  }

  double NeumannLeft(const double* a_position, const double a_time) {
    // Use just first term so we are certain thereA are enough terms included
    // (i.e., 1 will be exact).
    const std::size_t m = 1;
    const std::size_t n = 1;
    const auto dn = static_cast<double>(n);
    const auto dm = static_cast<double>(m);
    const double domain_length = domain_length_m;
    const double domain_height = domain_height_m;
    const double amplitude = amplitude_m;
    const double ic_factor =
        (amplitude *
         (-4.0 * std::pow(M_PI, 2) +
          std::pow(domain_height, 2) * (-8.0 + std::pow(M_PI, 2))) *
         (-4.0 * std::pow(M_PI, 2) +
          std::pow(domain_length, 2) *
              (-32.0 + 8.0 * M_PI + std::pow(M_PI, 2)))) /
        std::pow(M_PI, 6);
    return ic_factor * std::sin(dn * M_PI * a_position[1] / domain_height) *
           std::cos((dm - 0.5) * M_PI * a_position[0] / domain_length) *
           std::exp(-coefficient_m * M_PI * M_PI * a_time *
                    (std::pow(dm - 0.5, 2) / std::pow(domain_length, 2) +
                     std::pow(dn / domain_height, 2)));
  }

 private:
  const double coefficient_m;
  const double amplitude_m;
  const double domain_length_m;
  const double domain_height_m;
};

class DecayingPulse_TC {
 public:
  DecayingPulse_TC(const double a_coefficient_x, const double a_coefficient_y,
                   const double a_amplitude, const double a_domain_length,
                   const double a_domain_height)
      : coefficient_x_m(a_coefficient_x),
        coefficient_y_m(a_coefficient_y),
        amplitude_m(a_amplitude),
        domain_length_m(a_domain_length),
        domain_height_m(a_domain_height) {}

  double NeumannTopBot(const double* a_position, const double a_time) {
    // Use just first term so we are certain thereA are enough terms included
    // (i.e., 1 will be exact).
    const std::size_t m = 2;
    const std::size_t n = 1;
    const auto dn = static_cast<double>(n);
    const auto dm = static_cast<double>(m);
    const double domain_length = domain_length_m;
    const double domain_height = domain_height_m;
    const double amplitude = amplitude_m;
    const double ic_factor =
        (amplitude * std::pow(domain_height, 2) *
         (-4.0 * std::pow(M_PI, 2) +
          std::pow(domain_length, 2) * (-8.0 + std::pow(M_PI, 2)))) /
        std::pow(M_PI, 5);
    return ic_factor * std::sin(dn * M_PI * a_position[0] / domain_length) *
           std::cos(dm * M_PI * a_position[1] / domain_height) *
           std::exp(-M_PI * M_PI * a_time *
                    (coefficient_x_m * std::pow(dn / domain_length, 2) +
                     coefficient_y_m * std::pow(dm / domain_height, 2)));
  }

 private:
  const double coefficient_x_m;
  const double coefficient_y_m;
  const double amplitude_m;
  const double domain_length_m;
  const double domain_height_m;
};

TEST(LinearHeatEquation, HomogeneousAndConstantCoefficientTop) {
  std::vector<std::string> input_string;
  input_string.push_back("Executable_name");
  std::string file_name =
      "tests/verification/data/"
      "homogeneous_constant_coefficient_linear_heat_hn_top.json";
  input_string.push_back(file_name);
  input_string.push_back("-Mesh/parallel_refine");
  input_string.push_back("0");

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  const double final_time = input_file["Simulation"]["end_time"].get<double>();

  static constexpr std::size_t res_levels = 5;
  std::array<std::array<double, 4>, res_levels> error_metrics;
  for (std::size_t r = 0; r < res_levels; ++r) {
    input_string[3] = std::to_string(r);
    auto input_char = FakeCommandLineInput(input_string);
    int argc = static_cast<int>(input_char.size());
    char** argv = input_char.data();
    chyps::main(argc, argv, *mpi_session, SpdlogLevel::OFF);

    IO run_file(*mpi_session, "Run");
    run_file.SetRead("tests/verification/data/HomogeneousCCSolution");
    std::vector<double> temperature_field;
    run_file.GetImmediateBlock("HeatSolver/Temperature", temperature_field);
    std::vector<double> point_field;
    run_file.GetImmediateBlock("vertices", point_field);
    ASSERT_EQ(point_field.size() / 2, temperature_field.size());
    std::vector<uint32_t> nvert(1);
    run_file.GetImmediateBlock("NumOfVertices", {0}, {1}, nvert.data());
    assert(nvert[0] == point_field.size() / 2);

    std::vector<double> correct_solution(temperature_field.size());
    DecayingPulse_CC analytical_solution(0.5, 1.0, 2.0, 2.0);
    for (std::size_t n = 0; n < nvert[0]; ++n) {
      const double* position_2d = &(point_field[2 * n]);
      correct_solution[n] =
          analytical_solution.NeumannTop(position_2d, final_time);
    }

    const double global_l1 = GlobalL1Diff_Normalized(
        temperature_field, correct_solution, *mpi_session);
    const double global_l2 = GlobalL2Diff_Normalized(
        temperature_field, correct_solution, *mpi_session);
    const double global_linf =
        GlobalLinfDiff(temperature_field, correct_solution, *mpi_session);
    error_metrics[r][0] = static_cast<double>(std::sqrt(nvert[0]) - 1);
    error_metrics[r][1] = global_l1;
    error_metrics[r][2] = global_l2;
    error_metrics[r][3] = global_linf;

    DeleteCommandLineInput(input_char);
  }

  if (mpi_session->IAmRoot()) {
    printf("%16s %16s %16s %16s %16s %16s %16s\n", "Resolution", "L1_err",
           "L1_conv", "L2_err", "L2_conv", "Linf_err", "Linf_conv");
    for (std::size_t r = 0; r < res_levels; ++r) {
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
}

TEST(LinearHeatEquation, HomogeneousAndConstantCoefficientBot) {
  std::cout << "ANALYTICAL SOLUTION WRONG? " << std::endl;

  std::vector<std::string> input_string;
  input_string.push_back("Executable_name");
  std::string file_name =
      "tests/verification/data/"
      "homogeneous_constant_coefficient_linear_heat_hn_bot.json";
  input_string.push_back(file_name);
  input_string.push_back("-Mesh/parallel_refine");
  input_string.push_back("0");

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  const double final_time = input_file["Simulation"]["end_time"].get<double>();

  static constexpr std::size_t res_levels = 5;
  std::array<std::array<double, 4>, res_levels> error_metrics;
  for (std::size_t r = 0; r < res_levels; ++r) {
    input_string[3] = std::to_string(r);
    auto input_char = FakeCommandLineInput(input_string);
    int argc = static_cast<int>(input_char.size());
    char** argv = input_char.data();
    chyps::main(argc, argv, *mpi_session, SpdlogLevel::OFF);

    IO run_file(*mpi_session, "Run");
    run_file.SetRead("tests/verification/data/HomogeneousCCSolution");
    std::vector<double> temperature_field;
    run_file.GetImmediateBlock("HeatSolver/Temperature", temperature_field);
    std::vector<double> point_field;
    run_file.GetImmediateBlock("vertices", point_field);
    ASSERT_EQ(point_field.size() / 2, temperature_field.size());
    std::vector<uint32_t> nvert(1);
    run_file.GetImmediateBlock("NumOfVertices", {0}, {1}, nvert.data());
    assert(nvert[0] == point_field.size() / 2);

    std::vector<double> correct_solution(temperature_field.size());
    DecayingPulse_CC analytical_solution(0.5, 1.0, 2.0, 2.0);
    for (std::size_t n = 0; n < nvert[0]; ++n) {
      const double* position_2d = &(point_field[2 * n]);
      correct_solution[n] =
          analytical_solution.NeumannBot(position_2d, final_time);
    }

    const double global_l1 = GlobalL1Diff_Normalized(
        temperature_field, correct_solution, *mpi_session);
    const double global_l2 = GlobalL2Diff_Normalized(
        temperature_field, correct_solution, *mpi_session);
    const double global_linf =
        GlobalLinfDiff(temperature_field, correct_solution, *mpi_session);
    error_metrics[r][0] = static_cast<double>(std::sqrt(nvert[0]) - 1);
    error_metrics[r][1] = global_l1;
    error_metrics[r][2] = global_l2;
    error_metrics[r][3] = global_linf;

    DeleteCommandLineInput(input_char);
  }

  if (mpi_session->IAmRoot()) {
    printf("%16s %16s %16s %16s %16s %16s %16s\n", "Resolution", "L1_err",
           "L1_conv", "L2_err", "L2_conv", "Linf_err", "Linf_conv");
    for (std::size_t r = 0; r < res_levels; ++r) {
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
}

TEST(LinearHeatEquation, HomogeneousAndConstantCoefficientRight) {
  std::vector<std::string> input_string;
  input_string.push_back("Executable_name");
  std::string file_name =
      "tests/verification/data/"
      "homogeneous_constant_coefficient_linear_heat_hn_right.json";
  input_string.push_back(file_name);
  input_string.push_back("-Mesh/parallel_refine");
  input_string.push_back("0");

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  const double final_time = input_file["Simulation"]["end_time"].get<double>();

  static constexpr std::size_t res_levels = 5;
  std::array<std::array<double, 4>, res_levels> error_metrics;
  for (std::size_t r = 0; r < res_levels; ++r) {
    input_string[3] = std::to_string(r);
    auto input_char = FakeCommandLineInput(input_string);
    int argc = static_cast<int>(input_char.size());
    char** argv = input_char.data();
    chyps::main(argc, argv, *mpi_session, SpdlogLevel::OFF);

    IO run_file(*mpi_session, "Run");
    run_file.SetRead("tests/verification/data/HomogeneousCCSolution");
    std::vector<double> temperature_field;
    run_file.GetImmediateBlock("HeatSolver/Temperature", temperature_field);
    std::vector<double> point_field;
    run_file.GetImmediateBlock("vertices", point_field);
    ASSERT_EQ(point_field.size() / 2, temperature_field.size());
    std::vector<uint32_t> nvert(1);
    run_file.GetImmediateBlock("NumOfVertices", {0}, {1}, nvert.data());
    assert(nvert[0] == point_field.size() / 2);

    std::vector<double> correct_solution(temperature_field.size());
    DecayingPulse_CC analytical_solution(0.5, 1.0, 2.0, 2.0);
    for (std::size_t n = 0; n < nvert[0]; ++n) {
      const double* position_2d = &(point_field[2 * n]);
      correct_solution[n] =
          analytical_solution.NeumannRight(position_2d, final_time);
    }

    const double global_l1 = GlobalL1Diff_Normalized(
        temperature_field, correct_solution, *mpi_session);
    const double global_l2 = GlobalL2Diff_Normalized(
        temperature_field, correct_solution, *mpi_session);
    const double global_linf =
        GlobalLinfDiff(temperature_field, correct_solution, *mpi_session);
    error_metrics[r][0] = static_cast<double>(std::sqrt(nvert[0]) - 1);
    error_metrics[r][1] = global_l1;
    error_metrics[r][2] = global_l2;
    error_metrics[r][3] = global_linf;

    DeleteCommandLineInput(input_char);
  }

  if (mpi_session->IAmRoot()) {
    printf("%16s %16s %16s %16s %16s %16s %16s\n", "Resolution", "L1_err",
           "L1_conv", "L2_err", "L2_conv", "Linf_err", "Linf_conv");
    for (std::size_t r = 0; r < res_levels; ++r) {
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
}

TEST(LinearHeatEquation, HomogeneousAndConstantCoefficientLeft) {
  std::cout << "ANALYTICAL SOLUTION WRONG? " << std::endl;

  std::vector<std::string> input_string;
  input_string.push_back("Executable_name");
  std::string file_name =
      "tests/verification/data/"
      "homogeneous_constant_coefficient_linear_heat_hn_left.json";
  input_string.push_back(file_name);
  input_string.push_back("-Mesh/parallel_refine");
  input_string.push_back("0");

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  const double final_time = input_file["Simulation"]["end_time"].get<double>();

  static constexpr std::size_t res_levels = 5;
  std::array<std::array<double, 4>, res_levels> error_metrics;
  for (std::size_t r = 0; r < res_levels; ++r) {
    input_string[3] = std::to_string(r);
    auto input_char = FakeCommandLineInput(input_string);
    int argc = static_cast<int>(input_char.size());
    char** argv = input_char.data();
    chyps::main(argc, argv, *mpi_session, SpdlogLevel::OFF);

    IO run_file(*mpi_session, "Run");
    run_file.SetRead("tests/verification/data/HomogeneousCCSolution");
    std::vector<double> temperature_field;
    run_file.GetImmediateBlock("HeatSolver/Temperature", temperature_field);
    std::vector<double> point_field;
    run_file.GetImmediateBlock("vertices", point_field);
    ASSERT_EQ(point_field.size() / 2, temperature_field.size());
    std::vector<uint32_t> nvert(1);
    run_file.GetImmediateBlock("NumOfVertices", {0}, {1}, nvert.data());
    assert(nvert[0] == point_field.size() / 2);

    std::vector<double> correct_solution(temperature_field.size());
    DecayingPulse_CC analytical_solution(0.5, 1.0, 2.0, 2.0);
    for (std::size_t n = 0; n < nvert[0]; ++n) {
      const double* position_2d = &(point_field[2 * n]);
      correct_solution[n] =
          analytical_solution.NeumannLeft(position_2d, final_time);
    }

    const double global_l1 = GlobalL1Diff_Normalized(
        temperature_field, correct_solution, *mpi_session);
    const double global_l2 = GlobalL2Diff_Normalized(
        temperature_field, correct_solution, *mpi_session);
    const double global_linf =
        GlobalLinfDiff(temperature_field, correct_solution, *mpi_session);
    error_metrics[r][0] = static_cast<double>(std::sqrt(nvert[0]) - 1);
    error_metrics[r][1] = global_l1;
    error_metrics[r][2] = global_l2;
    error_metrics[r][3] = global_linf;

    DeleteCommandLineInput(input_char);
  }

  if (mpi_session->IAmRoot()) {
    printf("%16s %16s %16s %16s %16s %16s %16s\n", "Resolution", "L1_err",
           "L1_conv", "L2_err", "L2_conv", "Linf_err", "Linf_conv");
    for (std::size_t r = 0; r < res_levels; ++r) {
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
}

TEST(LinearHeatEquation, HomogeneousAndTensorCoefficient) {
  std::vector<std::string> input_string;
  input_string.push_back("Executable_name");
  std::string file_name =
      "tests/verification/data/"
      "homogeneous_tensor_coefficient.json";
  input_string.push_back(file_name);
  input_string.push_back("-Mesh/parallel_refine");
  input_string.push_back("0");

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  const double final_time = input_file["Simulation"]["end_time"].get<double>();

  static constexpr std::size_t res_levels = 5;
  std::array<std::array<double, 4>, res_levels> error_metrics;
  for (std::size_t r = 0; r < res_levels; ++r) {
    input_string[3] = std::to_string(r);
    auto input_char = FakeCommandLineInput(input_string);
    int argc = static_cast<int>(input_char.size());
    char** argv = input_char.data();
    chyps::main(argc, argv, *mpi_session, SpdlogLevel::OFF);

    IO run_file(*mpi_session, "Run");
    run_file.SetRead("tests/verification/data/HomogeneousTCSolution");
    std::vector<double> temperature_field;
    run_file.GetImmediateBlock("HeatSolver/Temperature", temperature_field);
    std::vector<double> point_field;
    run_file.GetImmediateBlock("vertices", point_field);
    ASSERT_EQ(point_field.size() / 2, temperature_field.size());
    std::vector<uint32_t> nvert(1);
    run_file.GetImmediateBlock("NumOfVertices", {0}, {1}, nvert.data());
    assert(nvert[0] == point_field.size() / 2);

    std::vector<double> correct_solution(temperature_field.size());
    DecayingPulse_TC analytical_solution(0.5, 0.25, 1.0, 2.0, 2.0);
    for (std::size_t n = 0; n < nvert[0]; ++n) {
      const double* position_2d = &(point_field[2 * n]);
      correct_solution[n] =
          analytical_solution.NeumannTopBot(position_2d, final_time);
    }

    const double global_l1 = GlobalL1Diff_Normalized(
        temperature_field, correct_solution, *mpi_session);
    const double global_l2 = GlobalL2Diff_Normalized(
        temperature_field, correct_solution, *mpi_session);
    const double global_linf =
        GlobalLinfDiff(temperature_field, correct_solution, *mpi_session);
    error_metrics[r][0] = static_cast<double>(std::sqrt(nvert[0]) - 1);
    error_metrics[r][1] = global_l1;
    error_metrics[r][2] = global_l2;
    error_metrics[r][3] = global_linf;

    DeleteCommandLineInput(input_char);
  }

  if (mpi_session->IAmRoot()) {
    printf("%16s %16s %16s %16s %16s %16s %16s\n", "Resolution", "L1_err",
           "L1_conv", "L2_err", "L2_conv", "Linf_err", "Linf_conv");
    for (std::size_t r = 0; r < res_levels; ++r) {
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
}

}  // namespace
