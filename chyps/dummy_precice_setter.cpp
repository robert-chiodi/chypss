// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cmath>
#include <utility>

#include <mpi.h>
#include <mfem/mfem.hpp>

#include "chyps/debug_assert.hpp"
#include "chyps/git.hpp"
#include "chyps/input_parser.hpp"
#include "chyps/io.hpp"
#include "chyps/logger.hpp"
#include "chyps/mesh.hpp"
#include "chyps/mpi_parallel.hpp"
#include "chyps/precice_adapter.hpp"
#include "chyps/string_manipulation.hpp"

namespace chyps {

int main(int argc, char** argv, MPIParallel& a_mpi_session) {
  StartLogger(a_mpi_session, SpdlogLevel::INFO);
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Beginning simulation.");

  DEBUG_ASSERT(argc >= 4, global_assert{}, DebugLevel::CHEAP{},
               "The command line argument -bc_tag and preCICE boundary number "
               "on for the chyps_heat executable must be supplied.");

  InputParser input_parser;
  // Related to mesh generation and use with chyps_heat executable.
  // Should be given same command line options as chyps_heat to
  // behave correctly.
  IO file_io(a_mpi_session, "preCICE_IO");
  Mesh mesh(a_mpi_session, input_parser, file_io);

  input_parser.AddOption("Simulation/precice_config",
                         "XML File holding the configuration for preCICE",
                         std::string("../data/precice-config.xml"));
  input_parser.AddOption("Simulation/end_time", "Final time; start time is 0.",
                         0.5);
  input_parser.AddOption("Simulationtime_step", "Time step.", 1.0e-2);
  input_parser.AddOption(
      "Simulation/option_output_name",
      "Name of file to write the options used for the simulation. Provides way "
      "to rerun the same simulation. Should end in the extension .json",
      "simulation_configuration.json");
  input_parser.AddOption(
      "bc_tag",
      "The boundary attribute tag for the boundary in the CHyPS mesh. .Must "
      "be supplied on the command line.");
  input_parser.AddOption(
      "bc_val",
      "Value to assigne to the boundary in the CHyPS mesh. Must "
      "be supplied on the command line.",
      0.0);

  const std::string input_file_name = argv[1];
  if (input_file_name == "help") {
    if (a_mpi_session.IAmRoot()) {
      input_parser.PrintOptions();
    }
    DEBUG_ASSERT(false, global_assert{}, DebugLevel::ALWAYS{});
  }
  if (input_file_name != "ignore") {
    input_parser.ParseFromFile(input_file_name, a_mpi_session);
  }
  argc -= 2;  // Skip executable and input file name
  input_parser.ParseCL(argc, argv + 2);
  if (a_mpi_session.IAmRoot()) {
    std::ofstream options_used(
        input_parser["Simulation/option_output_name"].get<std::string>() +
        "_precice");
    options_used << "/*\n" << GetGitDescriptionString() << "\n*/\n";
    const auto now = std::chrono::system_clock::now();
    const auto time = std::chrono::system_clock::to_time_t(now);
    options_used << "/*\n"
                 << "Simulation start time:\n"
                 << std::ctime(&time) << "*/\n";
    options_used << input_parser.WriteToString();
    options_used.close();
  }
  input_parser.ClearOptions();

  DEBUG_ASSERT(input_parser["Mesh/mesh_file"].get<std::string>() == "generate",
               global_assert{}, DebugLevel::CHEAP{},
               "preCICE dummy solver requires a generated mesh.");

  // Modify variables in parser to be for this mesh
  // as opposed to the base mesh used in chyps_heat
  int my_tag = -1;
  double shift;
  std::array<std::array<double, 3>, 2> bounding_box{
      {{input_parser["Mesh/gen_blx"].get<double>(),
        input_parser["Mesh/gen_bly"].get<double>(),
        input_parser["Mesh/gen_blz"].get<double>()},
       {input_parser["Mesh/gen_bux"].get<double>(),
        input_parser["Mesh/gen_buy"].get<double>(),
        input_parser["Mesh/gen_buz"].get<double>()}}};

  const int main_tag = input_parser["bc_tag"].get<int>();
  switch (main_tag) {
    case 1: {
      my_tag = 2;
      shift = bounding_box[0][0] - bounding_box[1][0];
      input_parser.DirectSet("Mesh/gen_nx", 1);
      input_parser.DirectSet("Mesh/gen_blx", bounding_box[0][0] + shift);
      input_parser.DirectSet("Mesh/gen_bux", bounding_box[1][0] + shift);
      break;
    }
    case 2: {
      my_tag = 1;
      shift = bounding_box[1][0] - bounding_box[0][0];
      input_parser.DirectSet("Mesh/gen_nx", 1);
      input_parser.DirectSet("Mesh/gen_blx", bounding_box[0][0] + shift);
      input_parser.DirectSet("Mesh/gen_bux", bounding_box[1][0] + shift);
      break;
    }
    case 3: {
      my_tag = 4;
      shift = bounding_box[0][1] - bounding_box[1][1];
      input_parser.DirectSet("Mesh/gen_ny", 1);
      input_parser.DirectSet("Mesh/gen_bly", bounding_box[0][1] + shift);
      input_parser.DirectSet("Mesh/gen_buy", bounding_box[1][1] + shift);
      break;
    }
    case 4: {
      my_tag = 3;
      shift = bounding_box[1][1] - bounding_box[0][1];
      input_parser.DirectSet("Mesh/gen_ny", 1);
      input_parser.DirectSet("Mesh/gen_bly", bounding_box[0][1] + shift);
      input_parser.DirectSet("Mesh/gen_buy", bounding_box[1][1] + shift);
      break;
    }
    case 5: {
      my_tag = 6;
      shift = bounding_box[0][2] - bounding_box[1][2];
      input_parser.DirectSet("Mesh/gen_nz", 1);
      input_parser.DirectSet("Mesh/gen_blz", bounding_box[0][2] + shift);
      input_parser.DirectSet("Mesh/gen_buz", bounding_box[1][2] + shift);
      break;
    }
    case 6: {
      my_tag = 5;
      shift = bounding_box[1][2] - bounding_box[0][2];
      input_parser.DirectSet("Mesh/gen_nz", 1);
      input_parser.DirectSet("Mesh/gen_blz", bounding_box[0][2] + shift);
      input_parser.DirectSet("Mesh/gen_buz", bounding_box[1][2] + shift);
      break;
    }
  }
  mesh.Initialize();

  PreciceAdapter precice(
      "DummySolver",
      input_parser["Simulation/precice_config"].get<std::string>(),
      a_mpi_session.MyRank(), a_mpi_session.NumberOfRanks());
  precice.AddMesh("DummyMesh");
  const std::vector<double>* vertex_positions;
  const std::vector<int>* vertex_indices;
  std::tie(vertex_positions, vertex_indices) = mesh.GetBoundaryVertices(my_tag);
  precice.SetVertexPositions("DummyMesh", *vertex_positions);
  precice.AddData("temperature_" + ZeroFill(main_tag, 3), "DummyMesh",
                  DataOperation::WRITE);
  double* temperature_bc = new double[vertex_indices->size()];
  const double value = input_parser["bc_val"];
  std::fill(temperature_bc, temperature_bc + vertex_indices->size(), value);

  double time = 0.0;
  bool last_step = false;
  double dt = input_parser["Simulation/time_step"];
  double max_dt = precice.Initialize();
  const double final_time = input_parser["Simulation/end_time"];
  for (int ti = 1; !last_step; ++ti) {
    for (std::size_t n = 0; n < vertex_indices->size(); ++n) {
      const double vertex_x = (*vertex_positions)[2 * n];
      const double center_x_position =
          time / final_time * (bounding_box[1][0] - bounding_box[0][0]) +
          bounding_box[0][0];
      temperature_bc[n] =
          (1.0 - std::fabs(vertex_x - center_x_position) /
                     (bounding_box[1][0] - bounding_box[0][0])) *
          10.0;
    }
    precice.WriteBlockScalarData("temperature_" + ZeroFill(main_tag, 3),
                                 temperature_bc);
    double dt_adj = dt;
    dt_adj = std::min(max_dt, dt_adj);

    if (time + dt_adj > final_time) {
      dt_adj = final_time - time;
    }
    if (std::fabs(time + dt_adj - final_time) < 1.0e-14) {
      last_step = true;
    }
    dt = dt_adj;
    time += dt;
    max_dt = precice.Advance(dt);
  }

  precice.Finalize();

  return 0;
}
}  // namespace chyps

int main(int argc, char** argv) {
  chyps::MPIParallel mpi_session(&argc, &argv);
  return chyps::main(argc, argv, mpi_session);
}
