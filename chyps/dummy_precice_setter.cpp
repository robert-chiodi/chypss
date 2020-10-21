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

#include "chyps/input_parser.hpp"
#include "chyps/logger.hpp"
#include "chyps/mesh.hpp"
#include "chyps/precice_adapter.hpp"

namespace chyps {

int main(int argc, char** argv) {
  int num_procs, myid;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  StartLogger(myid, num_procs, SpdlogLevel::OFF);
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Beginning simulation.");

  InputParser input_parser;
  // Related to mesh generation and use with chyps_heat executable.
  // Should be given same command line options as chyps_heat to
  // behave correctly.
  Mesh mesh(MPI_COMM_WORLD, input_parser);

  input_parser.AddOption<int>(
      "bc_tag", "-bt", "--boundary-condition-tag",
      "Tag value for  boundary condition. Value for mesh on chyps_heat mesh.",
      OptionType::COMMAND_LINE);
  input_parser.AddOption<double>("bc_val", "-bv", "--boundary-condition-value",
                                 "Value for boundary condition.",
                                 OptionType::COMMAND_LINE);

  input_parser.AddOption("precice_config", "-pc", "--preCICE-config",
                         "XML File holding the configuration for preCICE",
                         std::string("../data/precice-config.xml"));
  input_parser.AddOption("end_time", "-tf", "--t-final",
                         "Final time; start time is 0.", 0.5);
  input_parser.AddOption("time_step", "-dt", "--time-step", "Time step.",
                         1.0e-2);

  input_parser.ParseCL(argc, argv);
  assert(input_parser.AllOptionsSet());
  assert(static_cast<std::string>(input_parser["mesh_file"]) == "generate");

  // Modify variables in parser to be for this mesh
  // as opposed to the base mesh used in chyps_heat
  int my_tag = -1;
  double shift;
  std::array<std::array<double, 3>, 2> bounding_box{
      {{input_parser["gen_blx"], input_parser["gen_bly"],
        input_parser["gen_blz"]},
       {input_parser["gen_bux"], input_parser["gen_buy"],
        input_parser["gen_buz"]}}};

  const int main_tag = input_parser["bc_tag"];
  switch (main_tag) {
    case 1: {
      my_tag = 2;
      shift = bounding_box[0][0] - bounding_box[1][0];
      input_parser.DirectSet("gen_nx", 1);
      input_parser.DirectSet("gen_blx", bounding_box[0][0] + shift);
      input_parser.DirectSet("gen_bux", bounding_box[1][0] + shift);
      break;
    }
    case 2: {
      my_tag = 1;
      shift = bounding_box[1][0] - bounding_box[0][0];
      input_parser.DirectSet("gen_nx", 1);
      input_parser.DirectSet("gen_blx", bounding_box[0][0] + shift);
      input_parser.DirectSet("gen_bux", bounding_box[1][0] + shift);
      break;
    }
    case 3: {
      my_tag = 4;
      shift = bounding_box[0][1] - bounding_box[1][1];
      input_parser.DirectSet("gen_ny", 1);
      input_parser.DirectSet("gen_bly", bounding_box[0][1] + shift);
      input_parser.DirectSet("gen_buy", bounding_box[1][1] + shift);
      break;
    }
    case 4: {
      my_tag = 3;
      shift = bounding_box[1][1] - bounding_box[0][1];
      input_parser.DirectSet("gen_ny", 1);
      input_parser.DirectSet("gen_bly", bounding_box[0][1] + shift);
      input_parser.DirectSet("gen_buy", bounding_box[1][1] + shift);
      break;
    }
    case 5: {
      my_tag = 6;
      shift = bounding_box[0][2] - bounding_box[1][2];
      input_parser.DirectSet("gen_nz", 1);
      input_parser.DirectSet("gen_blz", bounding_box[0][2] + shift);
      input_parser.DirectSet("gen_buz", bounding_box[1][2] + shift);
      break;
    }
    case 6: {
      my_tag = 5;
      shift = bounding_box[1][2] - bounding_box[0][2];
      input_parser.DirectSet("gen_nz", 1);
      input_parser.DirectSet("gen_blz", bounding_box[0][2] + shift);
      input_parser.DirectSet("gen_buz", bounding_box[1][2] + shift);
      break;
    }
  }
  mesh.Initialize();

  PreciceAdapter precice("DummySolver", "DummyMesh",
                         input_parser["precice_config"], myid, num_procs);
  std::vector<double> vertex_positions;
  std::vector<int> vertex_indices;
  std::tie(vertex_positions, vertex_indices) = mesh.GetBoundaryVertices(my_tag);
  precice.SetVertexPositions(vertex_positions);
  precice.AddData("Temperature", DataOperation::WRITE);
  double* temperature_bc = new double[precice.ScalarDataSize()];
  const double value = input_parser["bc_val"];
  for (std::size_t n = 0; n < precice.ScalarDataSize(); ++n) {
    temperature_bc[n] = value;
  }

  double time = 0.0;
  bool last_step = false;
  double dt = input_parser["time_step"];
  double max_dt = precice.Initialize();
  const double final_time = input_parser["end_time"];
  for (int ti = 1; !last_step; ++ti) {
    for (std::size_t n = 0; n < precice.ScalarDataSize(); ++n) {
      const double vertex_x = vertex_positions[2 * n];
      const double center_x_position =
          time / final_time * (bounding_box[1][0] - bounding_box[0][0]) +
          bounding_box[0][0];
      temperature_bc[n] =
          (1.0 - std::fabs(vertex_x - center_x_position) /
                     (bounding_box[1][0] - bounding_box[0][0])) *
          10.0;
    }
    precice.WriteBlockScalarData("Temperature", temperature_bc);
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

  MPI_Finalize();
  return 0;
}
}  // namespace chyps

int main(int argc, char** argv) { return chyps::main(argc, argv); }
