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
#include "chyps/simulation_configurations/attribute_varying.hpp"
#include "chyps/simulation_configurations/constant.hpp"
#include "chyps/simulation_configurations/function.hpp"

namespace chyps {

InitializerType FunctionInitializers::Contains(
    const std::string& a_name) const {
  if (scalar_position_functions_m.find(a_name) !=
      scalar_position_functions_m.end()) {
    return InitializerType::SCALAR_POSITION;
  }

  if (matrix_position_functions_m.find(a_name) !=
      matrix_position_functions_m.end()) {
    return InitializerType::MATRIX_POSITION;
  }

  return InitializerType::INVALID;
}

void FunctionInitializers::AddScalarPositionFunction(
    const std::string& a_name, ScalarPositionFunction a_function) {
  DEBUG_ASSERT(this->Contains(a_name) == InitializerType::INVALID,
               global_assert{}, DebugLevel::CHEAP{},
               "Function for name \"" + a_name + "\" already added.");
  scalar_position_functions_m[a_name] = a_function;
}

void FunctionInitializers::AddMatrixPositionFunction(
    const std::string& a_name, MatrixPositionFunction a_function) {
  DEBUG_ASSERT(this->Contains(a_name) == InitializerType::INVALID,
               global_assert{}, DebugLevel::CHEAP{},
               "Function for name \"" + a_name + "\" already added.");
  matrix_position_functions_m[a_name] = a_function;
}

FunctionInitializers::ScalarPositionFunction
FunctionInitializers::GetScalarPositionFunction(
    const std::string& a_name) const {
  DEBUG_ASSERT(this->Contains(a_name) == InitializerType::SCALAR_POSITION,
               global_assert{}, DebugLevel::CHEAP{},
               "Cannot find ScalarPosition function with name" + a_name);
  return scalar_position_functions_m.at(a_name);
}

FunctionInitializers::MatrixPositionFunction
FunctionInitializers::GetMatrixPositionFunction(
    const std::string& a_name) const {
  DEBUG_ASSERT(this->Contains(a_name) == InitializerType::MATRIX_POSITION,
               global_assert{}, DebugLevel::CHEAP{},
               "Cannot find MatrixPosition function with name" + a_name);
  return matrix_position_functions_m.at(a_name);
}

RequiredData::RequiredData(const MPIParallel& a_mpi_session,
                           InputParser& a_parser, IO& a_file_io,
                           const FunctionInitializers* a_initializers)
    : mpi_session_m(a_mpi_session),
      parser_m(a_parser),
      file_io_m(a_file_io),
      initializers_m(a_initializers),
      mesh_m(a_mpi_session, parser_m, file_io_m),
      finite_element_collection_m(nullptr),
      finite_element_space_m(nullptr),
      scalar_fields_m(),
      matrix_fields_m() {
  parser_m.AddOption(
      "SimulationInitializer/FieldInitialization/DataName/Dimensionality",
      "The type of object to initialize. Value options are either \"scalar\" "
      "or \"matrix\".");

  parser_m.AddOption(
      "SimulationInitializer/FieldInitialization/DataName/NumberOfRows",
      "Number of rows in matrix if \"Dimensionality\" is of matrix type. "
      "Otherwise not required.");

  parser_m.AddOption(
      "SimulationInitializer/FieldInitialization/DataName/NumberOfColumns",
      "Number of columns in matrix if \"Dimensionality\" is of matrix type. "
      "Otherwise not required.");

  parser_m.AddOption(
      "SimulationInitializer/FieldInitialization/DataName/FieldType",
      "Location/type of field. Valid values are \"grid_function\", "
      "\"true_dof\", and "
      "\"element\".");

  parser_m.AddOption(
      "SimulationInitializer/FieldInitialization/DataName/Initializer",
      "Name of initializer to use for the field. See section "
      "SimulationInitializer/Initializers for list and options. Note: Not all "
      "initializers support all dimensionalities or field types.");

  parser_m.AddOption(
      "SimulationInitializer/FieldInitialization/DataName/Arguments",
      "Arguments passed to initializer. See description of individual "
      "initializers for expected arguments.");

  // Add parser options from all Initializers
  constant::scalar::AddParserOptions(parser_m);
  function::scalar::AddParserOptions(parser_m);
  attribute_varying::scalar::AddParserOptions(parser_m);

  constant::matrix::AddParserOptions(parser_m);
  function::matrix::AddParserOptions(parser_m);
  attribute_varying::matrix::AddParserOptions(parser_m);
}

void RequiredData::Initialize(void) { mesh_m.Initialize(); }

void RequiredData::InitializeFields(void) {
  // For now, require all added grid function or True DOF fields to be on
  // H1 elements and of the same order
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

    const auto field_type =
        field_initialization_params["FieldType"].get<std::string>();
    const auto field_enum = this->FieldStringToEnum(field_type);
    const auto initializer =
        field_initialization_params["Initializer"].get<std::string>();
    const auto& arguments_for_initializer =
        field_initialization_params["Arguments"];

    const std::string dimensionality =
        field_initialization_params["Dimensionality"].get<std::string>();

    if (dimensionality == "scalar") {
      auto scalar_field = new mfem::Vector;
      this->ApplyScalarInitializer(field_name, field_enum, initializer,
                                   arguments_for_initializer, parser_m, mesh_m,
                                   *finite_element_space, *scalar_field);
      this->AddScalarField(field_name, field_enum, scalar_field);
    } else if (dimensionality == "matrix") {
      auto matrix_list = new std::vector<mfem::DenseMatrix>;
      const auto number_of_rows =
          field_initialization_params["NumberOfRows"].get<int>();
      const auto number_of_columns =
          field_initialization_params["NumberOfColumns"].get<int>();
      this->ApplyMatrixInitializer(field_name, field_enum, initializer,
                                   arguments_for_initializer, parser_m, mesh_m,
                                   *finite_element_space, number_of_rows,
                                   number_of_columns, *matrix_list);
      this->AddMatrixField(field_name, field_enum, matrix_list);

    } else {
      DEBUG_ASSERT(false, global_assert{}, DebugLevel::ALWAYS{},
                   "Unknown dimensionality of \"" + dimensionality +
                       "\" for field named \"" + field_name + "\".");
    }
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

void RequiredData::AddScalarField(const std::string& a_name,
                                  const DataFieldType a_field_type,
                                  mfem::Vector* a_grid_function) {
  DEBUG_ASSERT(
      scalar_fields_m.find(a_name) == scalar_fields_m.end(), global_assert{},
      DebugLevel::CHEAP{},
      "Scalar field under name \"" + a_name + "\" has already been added.");
  scalar_fields_m[a_name] = std::make_pair(a_field_type, a_grid_function);
}

void RequiredData::AddMatrixField(
    const std::string& a_name, const DataFieldType a_field_type,
    std::vector<mfem::DenseMatrix>* a_matrix_list) {
  DEBUG_ASSERT(
      matrix_fields_m.find(a_name) == matrix_fields_m.end(), global_assert{},
      DebugLevel::CHEAP{},
      "Matrix field under name \"" + a_name + "\" has already been added.");
  matrix_fields_m[a_name] = std::make_pair(a_field_type, a_matrix_list);
}

void RequiredData::WriteToFile(void) {
  file_io_m.BeginWriteStep(0, 0.0, -1.0);
  this->CheckAllRequiredFieldsAvailable();
  this->WriteDataFields();
  mesh_m.WriteMesh();  // Need to write mesh after adding all variables.
  file_io_m.EndWriteStep();
}

void RequiredData::ApplyScalarInitializer(
    const std::string& a_field_name, const DataFieldType a_field_enum,
    const std::string& a_initializer,
    const nlohmann::json& a_initializer_arguments,
    const InputParser& a_full_parser, const Mesh& a_mesh,
    mfem::ParFiniteElementSpace& a_finite_element_space,
    mfem::Vector& a_field) {
  if (a_initializer == "constant") {
    DEBUG_ASSERT(constant::scalar::SupportedFieldType(a_field_enum),
                 global_assert{}, DebugLevel::CHEAP{});
    constant::scalar::InitializeData(a_field_enum, a_initializer_arguments,
                                     a_full_parser, a_mesh,
                                     a_finite_element_space, a_field);
  } else if (a_initializer == "function") {
    DEBUG_ASSERT(function::scalar::SupportedFieldType(a_field_enum),
                 global_assert{}, DebugLevel::CHEAP{});
    DEBUG_ASSERT(initializers_m != nullptr, global_assert{},
                 DebugLevel::CHEAP{});
    DEBUG_ASSERT(static_cast<int>(initializers_m->Contains(a_field_name)) > 0,
                 global_assert{}, DebugLevel::CHEAP{});
    auto setting_function =
        initializers_m->GetScalarPositionFunction(a_field_name);
    function::scalar::InitializeData(
        a_field_enum, a_initializer_arguments, a_full_parser, a_mesh,
        a_finite_element_space, a_field, setting_function);
  } else if (a_initializer == "attribute_varying") {
    DEBUG_ASSERT(attribute_varying::scalar::SupportedFieldType(a_field_enum),
                 global_assert{}, DebugLevel::CHEAP{});
    attribute_varying::scalar::InitializeData(
        a_field_enum, a_initializer_arguments, a_full_parser, a_mesh,
        a_finite_element_space, a_field);
  } else {
    DEBUG_ASSERT(false, global_assert{}, DebugLevel::ALWAYS{},
                 "Unknown scalar configuration type of : " + a_initializer);
  }
}

void RequiredData::ApplyMatrixInitializer(
    const std::string& a_field_name, const DataFieldType a_field_enum,
    const std::string& a_initializer,
    const nlohmann::json& a_initializer_arguments,
    const InputParser& a_full_parser, const Mesh& a_mesh,
    mfem::ParFiniteElementSpace& a_finite_element_space,
    const int a_number_of_rows, const int a_number_of_columns,
    std::vector<mfem::DenseMatrix>& a_field) {
  if (a_initializer == "constant") {
    DEBUG_ASSERT(constant::matrix::SupportedFieldType(a_field_enum),
                 global_assert{}, DebugLevel::CHEAP{});
    constant::matrix::InitializeData(
        a_field_enum, a_initializer_arguments, a_full_parser, a_mesh,
        a_finite_element_space, a_number_of_rows, a_number_of_columns, a_field);
  } else if (a_initializer == "function") {
    DEBUG_ASSERT(function::matrix::SupportedFieldType(a_field_enum),
                 global_assert{}, DebugLevel::CHEAP{});
    DEBUG_ASSERT(initializers_m != nullptr, global_assert{},
                 DebugLevel::CHEAP{});
    DEBUG_ASSERT(static_cast<int>(initializers_m->Contains(a_field_name)) > 0,
                 global_assert{}, DebugLevel::CHEAP{});
    auto setting_function =
        initializers_m->GetMatrixPositionFunction(a_field_name);
    function::matrix::InitializeData(
        a_field_enum, a_initializer_arguments, a_full_parser, a_mesh,
        a_finite_element_space, a_number_of_rows, a_number_of_columns, a_field,
        setting_function);
  } else if (a_initializer == "attribute_varying") {
    DEBUG_ASSERT(attribute_varying::matrix::SupportedFieldType(a_field_enum),
                 global_assert{}, DebugLevel::CHEAP{});
    attribute_varying::matrix::InitializeData(
        a_field_enum, a_initializer_arguments, a_full_parser, a_mesh,
        a_finite_element_space, a_number_of_rows, a_number_of_columns, a_field);
  } else {
    DEBUG_ASSERT(false, global_assert{}, DebugLevel::ALWAYS{},
                 "Unknown matrix configuration type of : " + a_initializer);
  }
}

DataFieldType RequiredData::FieldStringToEnum(const std::string& a_field_type) {
  if (a_field_type == "grid_function") {
    return DataFieldType::GRID_FUNCTION;
  } else if (a_field_type == "true_dof") {
    return DataFieldType::TRUE_DOF;
  } else if (a_field_type == "element") {
    return DataFieldType::ELEMENT;
  }
  DEBUG_ASSERT(false, global_assert{}, DebugLevel::ALWAYS{},
               "Unknown FieldType string of \"" + a_field_type + "\"");
  return DataFieldType::TRUE_DOF;  // Will never reach.
}

void RequiredData::CheckAllRequiredFieldsAvailable(void) const {
  static constexpr std::array<const char*, 3> required_variables{
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
    DEBUG_ASSERT(scalar_fields_m.find(elem) != scalar_fields_m.end(),
                 global_assert{}, DebugLevel::ALWAYS{},
                 "Required variable \"" + std::string(elem) +
                     "\" missing in initialization.");

    DEBUG_ASSERT(
        scalar_fields_m.at(elem).second != nullptr, global_assert{},
        DebugLevel::CHEAP{},
        "Required variable \"" + std::string(elem) + "\" points to a nullptr.");
  }

  if (parser_m.OptionSet("HeatSolver/Conductivity")) {
    if (parser_m["HeatSolver/Conductivity/type"] == "element_varying_scalar") {
      DEBUG_ASSERT(scalar_fields_m.find("HeatSolver/Conductivity/kappa") !=
                       scalar_fields_m.end(),
                   global_assert{}, DebugLevel::ALWAYS{},
                   "Element varying scalar kappa is being used but was not "
                   "initialized.");
    }
    if (parser_m["HeatSolver/Conductivity/type"] == "element_varying_matrix") {
      DEBUG_ASSERT(matrix_fields_m.find("HeatSolver/Conductivity/kappa") !=
                       matrix_fields_m.end(),
                   global_assert{}, DebugLevel::ALWAYS{},
                   "Element varying matrix kappa is being used but was not "
                   "initialized.");
    }
  }
}

void RequiredData::WriteDataFields(void) {
  // Write scalar fields
  for (const auto& elem : scalar_fields_m) {
    DEBUG_ASSERT(
        elem.second.second != nullptr, global_assert{}, DebugLevel::CHEAP{},
        "Variable \"" + elem.first + "\" is a nullptr in initializaton.");
    switch (elem.second.first) {
      case DataFieldType::GRID_FUNCTION: {
        file_io_m.AddVariableForGridFunction(elem.first,
                                             *finite_element_space_m, true);
        break;
      }

      case DataFieldType::TRUE_DOF: {
        file_io_m.AddVariableForTrueDofs(elem.first, *finite_element_space_m,
                                         true);
        break;
      }

      case DataFieldType::ELEMENT: {
        file_io_m.AddVariableForMesh<double>(elem.first, mesh_m,
                                             MeshElement::ELEMENT);
        break;
      }

      default:
        DEBUG_ASSERT(false, global_assert{}, DebugLevel::CHEAP{},
                     "Unhandled field type for writing scalar field to file.");
    }

    file_io_m.PutDeferred(elem.first, *(elem.second.second));
  }

  // Write matrix fields
  for (const auto& elem : matrix_fields_m) {
    DEBUG_ASSERT(
        elem.second.second != nullptr, global_assert{}, DebugLevel::CHEAP{},
        "Variable \"" + elem.first + "\" is a nullptr in initializaton.");
    switch (elem.second.first) {
      case DataFieldType::ELEMENT: {
        const mfem::DenseMatrix& matrix = (*elem.second.second)[0];
        const std::size_t number_of_rows =
            static_cast<std::size_t>(matrix.NumRows());
        const std::size_t number_of_columns =
            static_cast<std::size_t>(matrix.NumCols());
        file_io_m.AddMatrixForMesh(elem.first, mesh_m, MeshElement::ELEMENT,
                                   number_of_rows, number_of_columns);
        break;
      }

      default:
        DEBUG_ASSERT(false, global_assert{}, DebugLevel::CHEAP{},
                     "Unhandled field type for writing matrix field to file.");
    }
    file_io_m.PutDeferred(elem.first, *(elem.second.second));
  }
}

Mesh& RequiredData::GetMesh(void) { return mesh_m; }

RequiredData::~RequiredData(void) {
  for (auto& elem : scalar_fields_m) {
    delete elem.second.second;
    elem.second.second = nullptr;
  }
  for (auto& elem : matrix_fields_m) {
    delete elem.second.second;
    elem.second.second = nullptr;
  }
  delete finite_element_space_m;
  delete finite_element_collection_m;
}

SimulationInitializer::SimulationInitializer(
    int argc, char** argv, MPIParallel& a_mpi_session,
    const FunctionInitializers* a_initializers, SpdlogLevel a_log_level)
    : mpi_session_m(a_mpi_session), initializers_m(a_initializers), parser_m() {
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

  IO file_io(mpi_session_m, "FILEIO");
  RequiredData required_data(mpi_session_m, parser_m, file_io, a_initializers);

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
