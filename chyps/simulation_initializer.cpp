// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chyps/simulation_initializer.hpp"

#include <vector>

#include "chyps/debug_assert.hpp"
#include "chyps/io.hpp"

// List of configuration header files
#include "chyps/simulation_configurations/constant.hpp"
#include "chyps/simulation_configurations/cooled_rod.hpp"
#include "chyps/simulation_configurations/quadratic_pulse.hpp"

namespace chyps {

namespace chyps_details {
ConfigurationInitializer* GetConfigurationInitializer(
    const std::string& a_configuration_name, InputParser& a_parser) {
  if (a_configuration_name == "constant") {
    return new Constant(a_parser);
  } else if (a_configuration_name == "cooled_rod") {
    return new CooledRod(a_parser);

  } else if (a_configuration_name == "quadratic_pulse") {
    return new QuadraticPulse(a_parser);
  } else {
    DEBUG_ASSERT(false, global_assert{}, DebugLevel::ALWAYS{},
                 "Unknown configuration type of : " + a_configuration_name);
  }
  return nullptr;  // Should never reach
}
}  // namespace chyps_details

RequiredData::RequiredData(const MPIParallel& a_mpi_session,
                           InputParser& a_parser, IO& a_file_io)
    : parser_m(a_parser),
      file_io_m(a_file_io),
      mesh_m(a_mpi_session, parser_m, &file_io_m),
      finite_element_collection_m(nullptr),
      finite_element_space_m(nullptr) {}

void RequiredData::Initialize(void) { mesh_m.Initialize(); }

void RequiredData::SetFiniteElementCollection(
    mfem::FiniteElementCollection* a_collection) {
  finite_element_collection_m = a_collection;
}

void RequiredData::SetFiniteElementSpace(
    mfem::ParFiniteElementSpace* a_finite_element_space) {
  finite_element_space_m = a_finite_element_space;
}

void RequiredData::AddGridFunction(const std::string& a_name,
                                   mfem::ParGridFunction* a_grid_function) {
  DEBUG_ASSERT(
      data_fields_m.find(a_name) == data_fields_m.end(), global_assert{},
      DebugLevel::CHEAP{},
      "Grid function under name \"" + a_name + "\" has already been added.");
  data_fields_m[a_name] = a_grid_function;
}

void RequiredData::WriteToFile(void) {
  file_io_m.BeginWriteStep(0, 0.0, -1.0);
  this->CheckAllRequiredFieldsAvailable();
  this->WriteDataFields();
  mesh_m.WriteMesh();  // Need to write mesh after adding all variables.
  file_io_m.EndWriteStep();
}

void RequiredData::CheckAllRequiredFieldsAvailable(void) const {
  static const std::vector<std::string> required_variables{
      {"HeatSolver/temperature"}};

  DEBUG_ASSERT(finite_element_collection_m != nullptr, global_assert{},
               DebugLevel::CHEAP{},
               "The finite element collection was not set through the "
               "SetFiniteElementCollection method.");

  DEBUG_ASSERT(finite_element_space_m != nullptr, global_assert{},
               DebugLevel::CHEAP{},
               "The finite element space was not set through the "
               "SetFiniteElementSpace method.");

  for (const auto elem : required_variables) {
    DEBUG_ASSERT(
        data_fields_m.find(elem) != data_fields_m.end(), global_assert{},
        DebugLevel::ALWAYS{},
        "Required variable \"" + elem + "\" missing in initialization.");

    DEBUG_ASSERT(data_fields_m.at(elem) != nullptr, global_assert{},
                 DebugLevel::CHEAP{},
                 "Required variable \"" + elem + "\" points to a nullptr.");
  }
}

void RequiredData::WriteDataFields(void) {
  for (const auto& elem : data_fields_m) {
    DEBUG_ASSERT(
        elem.second != nullptr, global_assert{}, DebugLevel::CHEAP{},
        "Variable \"" + elem.first + "\" is a nullptr in initializaton.");
    file_io_m.AddVariableForGridFunction(elem.first, *finite_element_space_m,
                                         true);
    file_io_m.PutDeferred(elem.first, *elem.second);
  }
}

Mesh& RequiredData::GetMesh(void) { return mesh_m; }

RequiredData::~RequiredData(void) {
  for (auto& elem : data_fields_m) {
    delete elem.second;
    elem.second = nullptr;
  }
  delete finite_element_space_m;
  delete finite_element_collection_m;
}

SimulationInitializer::SimulationInitializer(int argc, char** argv,
                                             MPIParallel& a_mpi_session,
                                             SpdlogLevel a_log_level)
    : mpi_session_m(a_mpi_session),
      parser_m(),
      configuration_initializer_m(nullptr) {
  StartLogger(mpi_session_m, a_log_level);

  DEBUG_ASSERT(argc >= 3, global_assert{}, DebugLevel::ALWAYS{},
               "Command line options should be provided as [input_file_name] "
               "[configuration_name] [...command line parser flags...]");

  parser_m.AddOptionNoDefault(
      "SimulationInitializer/simulation_configuration",
      "Name of the simulation configuration to intialize.", true);
  parser_m.AddOptionDefault(
      "SimulationInitializer/out_data",
      "Name of file (or BP4 directory) to write that holds "
      "data to start from. Do not include extension.",
      std::string("CHyPSInitData"));

  const std::string input_file_name = argv[1];
  const std::string configuration_name = argv[2];

  // Should develop and construct a RequiredData class here that holds a mesh
  // and whatever ParGridFunctions need to be filled in by the configuration
  // initializer. Can then write that to file after being set.

  IO file_io(mpi_session_m, "FILEIO");
  RequiredData required_data(mpi_session_m, parser_m, file_io);
  configuration_initializer_m =
      chyps_details::GetConfigurationInitializer(configuration_name, parser_m);

  if (input_file_name == "help") {
    if (a_mpi_session.IAmRoot()) {
      parser_m.PrintOptions();
    }
    DEBUG_ASSERT(false, global_assert{}, DebugLevel::ALWAYS{});
  }
  parser_m.ParseFromFile(input_file_name, a_mpi_session);
  argc += -3;  // Skip executable and input file name
  parser_m.ParseCL(argc, argv + 3);

  const std::string out_data_name =
      parser_m["SimulationInitializer/out_data"].get<std::string>();
  file_io.SetWrite(out_data_name);
  file_io.RootWriteAttribute("InputFile", parser_m.WriteToString());

  required_data.Initialize();
  configuration_initializer_m->Initialize();
  configuration_initializer_m->FillRequiredData(required_data);
  required_data.WriteToFile();
}

SimulationInitializer::~SimulationInitializer(void) {
  delete configuration_initializer_m;
}

}  // namespace chyps
