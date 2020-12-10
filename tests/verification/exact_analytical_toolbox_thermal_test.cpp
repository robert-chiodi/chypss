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

// Solution with 10 terms to the problem below.
// Analytical solution/problem from the
// Exact Analytical Conduction Toolbox from
// University of Nebraska Lincoln, problem X1C11B00T11.
//
//                  Perfect
//                  Contact
//
//             |      |     |
//             |            |
// T(0,t) = 0  |      |     | T(L1+L2, t) = 0
//             |            |
//             |      |     |
//
// Properties left:
// Initial temperature = 1.0
// k = 1.0
// cp = 1.0
// rho = 1.0
// L1 = 1.0
//
// Properties right:
// Initial temperature = 0.0
// k = 5.0
// cp = 7/5
// rho = 1.0
// L2 = 2.0
//
TEST(ExactAnalyticalConductionToolbox, X1C11B00T11) {
  std::string file_name =
      "tests/verification/data/"
      "exact_analytical_conduction_toolbox_X1C11B00T11.json";
  static constexpr int number_of_refinements = 3;

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  static constexpr double length_1 = 1.0;
  static constexpr double kappa_left_value = 1.0;
  static constexpr double kappa_right_value = 5.0;
  static constexpr double cp_left_value = 1.0;
  static constexpr double cp_right_value = 7.0 / 5.0;

  auto solution_lambda = [=](const double* a_position, const double a_time) {
    if (a_position[0] <= length_1) {
      return (0.3278290555730198650731943502 *
              std::sin(1.5070544363619441150401104054 * a_position[0])) /
                 std::exp(2.2712130741582170644415889513 * a_time) +
             (0.9196048788528327407288338611 *
              std::sin(3.0919721252949537119373024282 * a_position[0])) /
                 std::exp(9.560291623600992936037835404 * a_time) +
             (0.1416859170035374784952016363 *
              std::sin(4.521632744028237985738644879 * a_position[0])) /
                 std::exp(20.445162671868333137881720806 * a_time) +
             (0.00114123180146417721220158409 *
              std::sin(6.1829020225210679084026041094 * a_position[0])) /
                 std::exp(38.228277420095112133195055241 * a_time) +
             (0.0512342948788687584205187662 *
              std::sin(7.537632137297866136574637139 * a_position[0])) /
                 std::exp(56.81589823722559749622675623 * a_time) +
             (0.2946165934372050951158346242 *
              std::sin(9.2717260253001723036996203861 * a_position[0])) /
                 std::exp(85.96490348822853134548205164 * a_time) +
             (0.0802527801078554126066576483 *
              std::sin(10.556038977683732724587807407 * a_position[0])) /
                 std::exp(111.42995889837822511066609293 * a_time);
    } else {
      return (0.3278290555730198650731943502 *
              (0.99796917345058996294120293796 *
                   std::cos(0.7974582501700642497677609187 *
                            (-1.0 + a_position[0])) +
               0.0240758588108096525514621011 *
                   std::sin(0.7974582501700642497677609187 *
                            (-1.0 + a_position[0])))) /
                 std::exp(2.2712130741582170644415889513 * a_time) +
             (0.9196048788528327407288338611 *
              (0.0496001682169631132552739126 *
                   std::cos(1.6361178608548584700859542323 *
                            (-1.0 + a_position[0])) -
               0.37749925701210685811786374484 *
                   std::sin(1.6361178608548584700859542323 *
                            (-1.0 + a_position[0])))) /
                 std::exp(9.560291623600992936037835404 * a_time) +
             (0.1416859170035374784952016363 *
              (-0.9818611323673927160345014702 *
                   std::cos(2.3926231521330586001512488153 *
                            (-1.0 + a_position[0])) -
               0.0716626190126058874726484512 *
                   std::sin(2.3926231521330586001512488153 *
                            (-1.0 + a_position[0])))) /
                 std::exp(20.445162671868333137881720806 * a_time) +
             (0.00114123180146417721220158409 *
              (-0.1001152820524177911549541854 *
                   std::cos(3.2716842264538048472848975719 *
                            (-1.0 + a_position[0])) +
               0.37606552056401114622714484049 *
                   std::sin(3.2716842264538048472848975719 *
                            (-1.0 + a_position[0])))) /
                 std::exp(38.228277420095112133195055241 * a_time) +
             (0.0512342948788687584205187662 *
              (0.9503774169755365102104257122 *
                   std::cos(3.9885400219156842878153752928 *
                            (-1.0 + a_position[0])) +
               0.1175844774151194836392120073 *
                   std::sin(3.9885400219156842878153752928 *
                            (-1.0 + a_position[0])))) /
                 std::exp(56.81589823722559749622675623 * a_time) +
             (0.2946165934372050951158346242 *
              (0.1524550973566883791663907735 *
                   std::cos(4.9061362574539232599964995361 *
                            (-1.0 + a_position[0])) -
               0.3735462188442707416183106567 *
                   std::sin(4.9061362574539232599964995361 *
                            (-1.0 + a_position[0])))) /
                 std::exp(85.96490348822853134548205164 * a_time) +
             (0.0802527801078554126066576483 *
              (-0.9049494953891065577090267614 *
                   std::cos(5.5857307929711312780756915917 *
                            (-1.0 + a_position[0])) -
               0.1608310607860671603053276498 *
                   std::sin(5.5857307929711312780756915917 *
                            (-1.0 + a_position[0])))) /
                 std::exp(111.42995889837822511066609293 * a_time);
    }
  };
  auto initial_condition_lambda = [=](const double* a_position) {
    return solution_lambda(a_position, 0.0);
  };

  auto cp_initial_lambda = [=](const double* a_position) {
    if (a_position[0] <= length_1) {
      return cp_left_value;
    } else {
      return cp_right_value;
    }
  };

  auto rho_initial_lambda = [=](const double* a_position) { return 1.0; };
  auto kappa_initial_lambda = [=](const double* a_position) {
    if (a_position[0] <= length_1) {
      return kappa_left_value;
    } else {
      return kappa_right_value;
    }
  };

  FunctionInitializers initializers;
  initializers.AddScalarPositionFunction("HeatSolver/temperature",
                                         initial_condition_lambda);
  initializers.AddScalarPositionFunction("HeatSolver/cp", cp_initial_lambda);
  initializers.AddScalarPositionFunction("HeatSolver/rho", rho_initial_lambda);
  initializers.AddScalarPositionFunction("HeatSolver/Conductivity/kappa",
                                         kappa_initial_lambda);

  auto test_result =
      ConvergenceRunner(file_name, solution_lambda, number_of_refinements, true,
                        &initializers, 1);
}

// Solution with 10 terms to the problem below.
// Analytical solution/problem from the
// Exact Analytical Conduction Toolbox from
// University of Nebraska Lincoln, problem X1C12B00T11.
//
//                  Perfect
//                  Contact
//
//             |      |     |
//             |            |
// T(0,t) = 0  |      |     | dT/dn(L1+L2, t) = 0
//             |            |
//             |      |     |
//
// Properties left:
// Initial temperature = 1.0
// k = 1.0
// cp = 1.0
// rho = 1.0
// L1 = 1.0
//
// Properties right:
// Initial temperature = 0.0
// k = 5.0
// cp = 7/5
// rho = 1.0
// L2 = 2.0
//
TEST(ExactAnalyticalConductionToolbox, X1C12B00T11) {
  std::string file_name =
      "tests/verification/data/"
      "exact_analytical_conduction_toolbox_X1C12B00T11.json";
  static constexpr int number_of_refinements = 3;

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  static constexpr double length_1 = 1.0;
  static constexpr double kappa_left_value = 1.0;
  static constexpr double kappa_right_value = 5.0;
  static constexpr double cp_left_value = 1.0;
  static constexpr double cp_right_value = 7.0 / 5.0;

  auto solution_lambda = [=](const double* a_position, const double a_time) {
    if (a_position[0] <= length_1) {
      return (0.2582410291292657096936748776 *
              std::sin(0.53572895989671380442028885075 * a_position[0])) /
                 std::exp(0.28700551847201478773011753602 * a_time) +
             (0.646456047624290515064273895 *
              std::sin(2.5143668216738443521856073599 * a_position[0])) /
                 std::exp(6.32204051393422980476152617 * a_time) +
             (0.569064052365579255374072461 *
              std::sin(3.5935769764261007220206213528 * a_position[0])) /
                 std::exp(12.913795485499756064402124007 * a_time) +
             (0.0369436966392649960765965023 *
              std::sin(5.5569042555230729855216326159 * a_position[0])) /
                 std::exp(30.8791849050504380231150327 * a_time) +
             (0.0121609130227759446780826774 *
              std::sin(6.658949040203146081275418557 * a_position[0])) /
                 std::exp(44.34160232002240040581868908 * a_time) +
             (0.1435785715396410378698976265 *
              std::sin(8.59245865169831364821691615 * a_position[0])) /
                 std::exp(73.83034568114520209695326574 * a_time) +
             (0.248062826679199791486475849 *
              std::sin(9.731234023993870102612662702 * a_position[0])) /
                 std::exp(94.69691562973592964396146674 * a_time);
    } else {
      return (0.2582410291292657096936748776 *
              (0.51046800524399462180239831924 *
                   std::cos(0.2834811196044000045010643983 *
                            (-1.0 + a_position[0])) +
               0.3250104208488040551026880062 *
                   std::sin(0.2834811196044000045010643983 *
                            (-1.0 + a_position[0])))) /
                 std::exp(0.28700551847201478773011753602 * a_time) +
             (0.646456047624290515064273895 *
              (0.5869008896628118215718253424 *
                   std::cos(1.330477862988176269390920158 *
                            (-1.0 + a_position[0])) -
               0.3060222554262072862382146011 *
                   std::sin(1.330477862988176269390920158 *
                            (-1.0 + a_position[0])))) /
                 std::exp(6.32204051393422980476152617 * a_time) +
             (0.569064052365579255374072461 *
              (-0.4367514543405448432083784577 *
                   std::cos(1.9015421993581766677648956547 *
                            (-1.0 + a_position[0])) -
               0.3400101189946765966769250784 *
                   std::sin(1.9015421993581766677648956547 *
                            (-1.0 + a_position[0])))) /
                 std::exp(12.913795485499756064402124007 * a_time) +
             (0.0369436966392649960765965023 *
              (-0.6640937646919869222088810938 *
                   std::cos(2.9404373439021146138835440542 *
                            (-1.0 + a_position[0])) +
               0.2825848726355686536833688032 *
                   std::sin(2.9404373439021146138835440542 *
                            (-1.0 + a_position[0])))) /
                 std::exp(30.8791849050504380231150327 * a_time) +
             (0.0121609130227759446780826774 *
              (0.3669830815952850583680981 *
                   std::cos(3.5235846306859541802971679731 *
                            (-1.0 + a_position[0])) +
               0.3515929907116279629210623524 *
                   std::sin(3.5235846306859541802971679731 *
                            (-1.0 + a_position[0])))) /
                 std::exp(44.34160232002240040581868908 * a_time) +
             (0.1435785715396410378698976265 *
              (0.7394946304415241657424895748 *
                   std::cos(4.5467017485998195378866019837 *
                            (-1.0 + a_position[0])) -
               0.2544314927576362072976198936 *
                   std::sin(4.5467017485998195378866019837 *
                            (-1.0 + a_position[0])))) /
                 std::exp(73.83034568114520209695326574 * a_time) +
             (0.248062826679199791486475849 *
              (-0.301681717859184671165389737 *
                   std::cos(5.1492850354516266863740670769 *
                            (-1.0 + a_position[0])) -
               0.3603546152471659996837230831 *
                   std::sin(5.1492850354516266863740670769 *
                            (-1.0 + a_position[0])))) /
                 std::exp(94.69691562973592964396146674 * a_time);
    }
  };
  auto initial_condition_lambda = [=](const double* a_position) {
    return solution_lambda(a_position, 0.0);
  };

  auto cp_initial_lambda = [=](const double* a_position) {
    if (a_position[0] <= length_1) {
      return cp_left_value;
    } else {
      return cp_right_value;
    }
  };

  auto rho_initial_lambda = [=](const double* a_position) { return 1.0; };
  auto kappa_initial_lambda = [=](const double* a_position) {
    if (a_position[0] <= length_1) {
      return kappa_left_value;
    } else {
      return kappa_right_value;
    }
  };

  FunctionInitializers initializers;
  initializers.AddScalarPositionFunction("HeatSolver/temperature",
                                         initial_condition_lambda);
  initializers.AddScalarPositionFunction("HeatSolver/cp", cp_initial_lambda);
  initializers.AddScalarPositionFunction("HeatSolver/rho", rho_initial_lambda);
  initializers.AddScalarPositionFunction("HeatSolver/Conductivity/kappa",
                                         kappa_initial_lambda);

  auto test_result =
      ConvergenceRunner(file_name, solution_lambda, number_of_refinements, true,
                        &initializers, 1);
}

// Solution with 10 terms to the problem below.
// Analytical solution/problem from the
// Exact Analytical Conduction Toolbox from
// University of Nebraska Lincoln, problem X1C12B00T11.
//
//                  Perfect
//                  Contact
//
//             |      |     |
//             |            |
// T(0,t) = 1  |      |     | dT/dn(L1+L2, t) = 0
//             |            |
//             |      |     |
//
// Properties left:
// Initial temperature = 1.0
// k = 1.0
// cp = 1.0
// rho = 1.0
// L1 = 1.0
//
// Properties right:
// Initial temperature = 0.0
// k = 5.0
// cp = 7/5
// rho = 1.0
// L2 = 2.0
//
TEST(ExactAnalyticalConductionToolbox, X1C12B00T0) {
  std::string file_name =
      "tests/verification/data/"
      "exact_analytical_conduction_toolbox_X1C12B00T0.json";
  static constexpr int number_of_refinements = 3;

  std::ifstream myfile(file_name.c_str());
  nlohmann::json input_file =
      nlohmann::json::parse(myfile, nullptr, true, true);
  myfile.close();

  static constexpr double length_1 = 1.0;
  static constexpr double kappa_left_value = 1.0;
  static constexpr double kappa_right_value = 5.0;
  static constexpr double cp_left_value = 1.0;
  static constexpr double cp_right_value = 7.0 / 5.0;

  auto solution_lambda = [=](const double* a_position, const double a_time) {
    if (a_position[0] <= length_1) {
      return 1.0 -
             (1.8432193666408816 *
              std::sin(0.53572895989671380442028885075 * a_position[0])) /
                 std::exp(0.28700551847201478773011753602 * a_time) -
             (0.35722538056088504 *
              std::sin(2.5143668216738443521856073599 * a_position[0])) /
                 std::exp(6.32204051393422980476152617 * a_time) -
             (0.2995732677120724 *
              std::sin(3.5935769764261007220206213528 * a_position[0])) /
                 std::exp(12.913795485499756064402124007 * a_time) -
             (0.1463982316613783 *
              std::sin(5.5569042555230729855216326159 * a_position[0])) /
                 std::exp(30.8791849050504380231150327 * a_time) -
             (0.17429407380649936 *
              std::sin(6.658949040203146081275418557 * a_position[0])) /
                 std::exp(44.34160232002240040581868908 * a_time) -
             (0.08581269025268283 *
              std::sin(8.59245865169831364821691615 * a_position[0])) /
                 std::exp(73.83034568114520209695326574 * a_time) -
             (0.12698972172100167 *
              std::sin(9.731234023993870102612662702 * a_position[0])) /
                 std::exp(94.69691562973592964396146674 * a_time);
    } else {
      return 1.0 -
             (1.8432193666408816 *
              (0.51046800524399462180239831924 *
                   std::cos(0.2834811196044000045010643983 *
                            (-1.0 + a_position[0])) +
               0.3250104208488040551026880062 *
                   std::sin(0.2834811196044000045010643983 *
                            (-1.0 + a_position[0])))) /
                 std::exp(0.28700551847201478773011753602 * a_time) -
             (0.35722538056088504 *
              (0.5869008896628118215718253424 *
                   std::cos(1.330477862988176269390920158 *
                            (-1.0 + a_position[0])) -
               0.3060222554262072862382146011 *
                   std::sin(1.330477862988176269390920158 *
                            (-1.0 + a_position[0])))) /
                 std::exp(6.32204051393422980476152617 * a_time) -
             (0.2995732677120724 *
              (-0.4367514543405448432083784577 *
                   std::cos(1.9015421993581766677648956547 *
                            (-1.0 + a_position[0])) -
               0.3400101189946765966769250784 *
                   std::sin(1.9015421993581766677648956547 *
                            (-1.0 + a_position[0])))) /
                 std::exp(12.913795485499756064402124007 * a_time) -
             (0.1463982316613783 *
              (-0.6640937646919869222088810938 *
                   std::cos(2.9404373439021146138835440542 *
                            (-1.0 + a_position[0])) +
               0.2825848726355686536833688032 *
                   std::sin(2.9404373439021146138835440542 *
                            (-1.0 + a_position[0])))) /
                 std::exp(30.8791849050504380231150327 * a_time) -
             (0.17429407380649936 *
              (0.3669830815952850583680981 *
                   std::cos(3.5235846306859541802971679731 *
                            (-1.0 + a_position[0])) +
               0.3515929907116279629210623524 *
                   std::sin(3.5235846306859541802971679731 *
                            (-1.0 + a_position[0])))) /
                 std::exp(44.34160232002240040581868908 * a_time) -
             (0.08581269025268283 *
              (0.7394946304415241657424895748 *
                   std::cos(4.5467017485998195378866019837 *
                            (-1.0 + a_position[0])) -
               0.2544314927576362072976198936 *
                   std::sin(4.5467017485998195378866019837 *
                            (-1.0 + a_position[0])))) /
                 std::exp(73.83034568114520209695326574 * a_time) -
             (0.12698972172100167 *
              (-0.301681717859184671165389737 *
                   std::cos(5.1492850354516266863740670769 *
                            (-1.0 + a_position[0])) -
               0.3603546152471659996837230831 *
                   std::sin(5.1492850354516266863740670769 *
                            (-1.0 + a_position[0])))) /
                 std::exp(94.69691562973592964396146674 * a_time);
    }
  };
  auto initial_condition_lambda = [=](const double* a_position) {
    return solution_lambda(a_position, 0.0);
  };

  auto cp_initial_lambda = [=](const double* a_position) {
    if (a_position[0] <= length_1) {
      return cp_left_value;
    } else {
      return cp_right_value;
    }
  };

  auto rho_initial_lambda = [=](const double* a_position) { return 1.0; };
  auto kappa_initial_lambda = [=](const double* a_position) {
    if (a_position[0] <= length_1) {
      return kappa_left_value;
    } else {
      return kappa_right_value;
    }
  };

  FunctionInitializers initializers;
  initializers.AddScalarPositionFunction("HeatSolver/temperature",
                                         initial_condition_lambda);
  initializers.AddScalarPositionFunction("HeatSolver/cp", cp_initial_lambda);
  initializers.AddScalarPositionFunction("HeatSolver/rho", rho_initial_lambda);
  initializers.AddScalarPositionFunction("HeatSolver/Conductivity/kappa",
                                         kappa_initial_lambda);

  auto test_result =
      ConvergenceRunner(file_name, solution_lambda, number_of_refinements, true,
                        &initializers, 1);
}

}  // namespace
