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

RequiredData::RequiredData(const MPIParallel& a_mpi_session,
                           InputParser& a_parser, IO& a_file_io)
    : mpi_session_m(a_mpi_session),
      parser_m(a_parser),
      file_io_m(a_file_io),
      mesh_m(a_mpi_session, parser_m, file_io_m),
      finite_element_collection_m(nullptr),
      finite_element_space_m(nullptr) {
  parser_m.AddOption(
      "SimulationInitializer/FieldInitialization",
      "A JSON object that lists the field names to be initialized as keys and "
      "values of individual initializers for how to initialize them.  As an "
      "example:\n"
      "\"HeatSolver/temperature\": {\n"
      " \"Initializer\": "
      "\"constant\",\n"
      "\" Arguments\":{"
      "\n \"value\": 1.0\n"
      "}\n"
      "See the options themselves in Initializers for the arguments that must "
      "be supplied to them");

  // Add parser options from all Initializers
  constant::AddParserOptions(parser_m);
  quadratic_pulse::AddParserOptions(parser_m);
  cooled_rod::AddParserOptions(parser_m);
}

void RequiredData::Initialize(void) { mesh_m.Initialize(); }

void RequiredData::InitializeFields(void) {
  // For now, require all added fields to be on H1 elements and
  // of the same order
  auto finite_element_collection = new mfem::H1_FECollection(
      parser_m["HeatSolver/order"].get<int>(), mesh_m.GetDimension());
  auto finite_element_space = new mfem::ParFiniteElementSpace(
      &(mesh_m.GetMfemMesh()), finite_element_collection);
  this->SetFiniteElementCollection(finite_element_collection);
  this->SetFiniteElementSpace(finite_element_space);

  const auto& field_initialization_arguments =
      parser_m["SimulationInitializer/FieldInitialization"];
  for (const auto& field : field_initialization_arguments.items()) {
    const std::string& field_name = field.key();
    const auto& field_initialization_params = field.value();
    const auto& initializer =
        field_initialization_params["Initializer"].get<std::string>();
    const auto& arguments_for_initializer =
        field_initialization_params["Arguments"];

    auto true_dof_field = new mfem::Vector;
    this->ApplyInitializer(initializer, arguments_for_initializer, parser_m,
                           *finite_element_space, *true_dof_field);

    // RequiredData destructor will
    // handle freeing the memory allocated
    this->AddTrueDofField(field_name, true_dof_field);
  }
}

void RequiredData::SetFiniteElementCollection(
    mfem::FiniteElementCollection* a_collection) {
  finite_element_collection_m = a_collection;
}

void RequiredData::SetFiniteElementSpace(
    mfem::ParFiniteElementSpace* a_finite_element_space) {
  finite_element_space_m = a_finite_element_space;
}

void RequiredData::AddTrueDofField(const std::string& a_name,
                                   mfem::Vector* a_grid_function) {
  DEBUG_ASSERT(
      data_fields_m.find(a_name) == data_fields_m.end(), global_assert{},
      DebugLevel::CHEAP{},
      "mfem::Vector under name \"" + a_name + "\" has already been added.");
  data_fields_m[a_name] = a_grid_function;
}

void RequiredData::WriteToFile(void) {
  file_io_m.BeginWriteStep(0, 0.0, -1.0);
  this->CheckAllRequiredFieldsAvailable();
  this->WriteDataFields();
  mesh_m.WriteMesh();  // Need to write mesh after adding all variables.
  file_io_m.EndWriteStep();
}

void RequiredData::ApplyInitializer(
    const std::string& a_initializer,
    const nlohmann::json& a_initializer_arguments,
    const InputParser& a_full_parser,
    mfem::ParFiniteElementSpace& a_finite_element_space,
    mfem::Vector& a_field) {
  if (a_initializer == "constant") {
    constant::InitializeData(a_initializer_arguments, a_full_parser,
                             a_finite_element_space, a_field);
  } else if (a_initializer == "cooled_rod") {
    cooled_rod::InitializeData(a_initializer_arguments, a_full_parser,
                               a_finite_element_space, a_field);

  } else if (a_initializer == "quadratic_pulse") {
    quadratic_pulse::InitializeData(a_initializer_arguments, a_full_parser,
                                    a_finite_element_space, a_field);
  } else {
    DEBUG_ASSERT(false, global_assert{}, DebugLevel::ALWAYS{},
                 "Unknown configuration type of : " + a_initializer);
  }
}

void RequiredData::CheckAllRequiredFieldsAvailable(void) const {
  static const std::vector<std::string> required_variables{
      {"HeatSolver/temperature", "HeatSolver/rho", "HeatSolver/cp"}};

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
    file_io_m.AddVariableForTrueDofs(elem.first, *finite_element_space_m, true);
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
    : mpi_session_m(a_mpi_session), parser_m() {
  StartLogger(mpi_session_m, a_log_level);

  DEBUG_ASSERT(argc >= 2, global_assert{}, DebugLevel::ALWAYS{},
               "Command line options should be provided as [input_file_name] "
               "[...command line parser flags...]");

  parser_m.AddOption("SimulationInitializer/initializer",
                     "Name of initializer to use.");
  parser_m.AddOption("SimulationInitializer/out_data",
                     "Name of file (or BP4 directory) to write that holds "
                     "data to start from. Do not include extension.");

  const std::string input_file_name = argv[1];

  // Should develop and construct a RequiredData class here that holds a mesh
  // and whatever ParGridFunctions need to be filled in by the configuration
  // initializer. Can then write that to file after being set.

  IO file_io(mpi_session_m, "FILEIO");
  RequiredData required_data(mpi_session_m, parser_m, file_io);

  if (input_file_name == "help") {
    if (a_mpi_session.IAmRoot()) {
      parser_m.PrintOptions();
    }
    DEBUG_ASSERT(false, global_assert{}, DebugLevel::ALWAYS{});
  }
  parser_m.ParseFromFile(input_file_name, a_mpi_session);
  argc += -2;  // Skip executable and input file name
  parser_m.ParseCL(argc, argv + 2);

  const std::string out_data_name =
      parser_m["SimulationInitializer/out_data"].get<std::string>();
  file_io.SetWrite(out_data_name);
  file_io.RootWriteAttribute("InputFile", parser_m.WriteToString());
  required_data.Initialize();
  required_data.InitializeFields();
  required_data.WriteToFile();
}

}  // namespace chyps
