// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef TESTS_HELPER_CONVERGENCE_RUNNER_HPP_
#define TESTS_HELPER_CONVERGENCE_RUNNER_HPP_

#include <array>
#include <vector>

#include "chyps/debug_assert.hpp"
#include "chyps/io.hpp"
#include "chyps/logger.hpp"
#include "chyps/simulation.hpp"
#include "chyps/simulation_initializer.hpp"

namespace chyps {

// Performs uniform refinement and reports errors.
// Error metrics returned are
// refinement_level : {sqrt(Nodes), L1, L2, Linf}
template <class SolutionFunctor>
std::vector<std::array<double, 4>> ConvergenceRunner(
    const std::string& a_input_name,
    const SolutionFunctor& a_analytical_solution_lambda,
    const int a_number_of_refines, const bool a_print = true) {
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
    chyps_input_string[3] = std::to_string(r);
    auto chyps_input_char = FakeCommandLineInput(chyps_input_string);
    int argc = static_cast<int>(chyps_input_char.size());
    char** argv = chyps_input_char.data();
    chyps::SimulationInitializer(argc, argv, *mpi_session, SpdlogLevel::OFF);

    // Perform simulation
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
    int send_temp = static_cast<int>(nvert[0]);
    int total_nodes = 0;
    MPI_Allreduce(&send_temp, &total_nodes, 1, MPI_INT, MPI_SUM,
                  mpi_session->GetComm());
    DEBUG_ASSERT(nvert[0] == point_field.size() / 2, global_assert{},
                 DebugLevel::CHEAP{});

    std::vector<double> correct_solution(temperature_field.size());
    for (std::size_t n = 0; n < nvert[0]; ++n) {
      const double* position_2d = &(point_field[2 * n]);
      correct_solution[n] =
          a_analytical_solution_lambda(position_2d, final_time);
    }

    const double global_l1 =
        GlobalL1Diff(temperature_field, correct_solution, *mpi_session);
    const double global_l2 =
        GlobalL2Diff(temperature_field, correct_solution, *mpi_session);
    const double global_linf =
        GlobalLinfDiff(temperature_field, correct_solution, *mpi_session);
    error_metrics[r][0] = static_cast<double>(total_nodes);
    error_metrics[r][1] = global_l1;
    error_metrics[r][2] = global_l2;
    error_metrics[r][3] = global_linf;

    DeleteCommandLineInput(chyps_input_char);
  }

  if (mpi_session->IAmRoot() && a_print) {
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
                (std::sqrt(error_metrics[r - 1][0] / error_metrics[r][0])));
        const double l2_conv =
            log(error_metrics[r - 1][2] / error_metrics[r][2] /
                (std::sqrt(error_metrics[r - 1][0] / error_metrics[r][0])));
        const double linf_conv =
            log(error_metrics[r - 1][3] / error_metrics[r][3] /
                (std::sqrt(error_metrics[r - 1][0] / error_metrics[r][0])));
        printf("%16d %16.8E %16.8f %16.8E %16.8f %16.8E %16.8f\n",
               static_cast<int>(error_metrics[r][0]), error_metrics[r][1],
               l1_conv, error_metrics[r][2], l2_conv, error_metrics[r][3],
               linf_conv);
      }
    }
  }
  return error_metrics;  // Come up with actual enforceable error measure
}

}  // namespace chyps

#endif  // TESTS_HELPER_CONVERGENCE_RUNNER_HPP_
