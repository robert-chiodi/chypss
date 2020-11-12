// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef CHYPS_SIMULATION_INITIALIZER_HPP_
#define CHYPS_SIMULATION_INITIALIZER_HPP_

#include <string>
#include <unordered_map>

#include <mfem/mfem.hpp>

#include "chyps/input_parser.hpp"
#include "chyps/io.hpp"
#include "chyps/logger.hpp"
#include "chyps/mpi_parallel.hpp"

namespace chyps {

// Forward declaration to preveAnt circular referencing.
class ConfigurationInitializer;

class RequiredData {
 public:
  RequiredData(const MPIParallel& a_mpi_session, InputParser& a_parser,
               IO& a_file_io);

  /// \brief Initialize objects using parsed parameters (mostly constructs
  /// mesh).
  void Initialize(void);

  /// \brief Set the finite element collection that the finite element space
  /// will live on.
  ///
  /// NOTE: RequiredData object will take ownership of the finite element
  /// collection.
  void SetFiniteElementCollection(mfem::FiniteElementCollection* a_collection);

  /// \brief Set the finite element space that the data fields will live on.
  ///
  /// NOTE: RequiredData object will take ownership of the finite element space.
  void SetFiniteElementSpace(
      mfem::ParFiniteElementSpace* a_finite_element_space);

  /// \brief Add the grid function (with data correctly set) to the RequiredData
  /// object under the name `a_name`.
  ///
  /// NOTE: RequiredData object will take ownership of the GridFunction object.
  void AddGridFunction(const std::string& a_name,
                       mfem::ParGridFunction* a_grid_function);

  /// \brief Write out mesh and data fields to file to be used for initializing
  /// a simulation.
  void WriteToFile(void);

  /// \brief Return a reference to the contained mesh data is placed on.
  Mesh& GetMesh(void);

  ~RequiredData(void);

 private:
  void CheckAllRequiredFieldsAvailable(void) const;
  void WriteDataFields(void);

  InputParser& parser_m;
  IO& file_io_m;
  Mesh mesh_m;
  mfem::FiniteElementCollection* finite_element_collection_m;
  mfem::ParFiniteElementSpace* finite_element_space_m;
  std::unordered_map<std::string, mfem::ParGridFunction*> data_fields_m;
};

/// Simulation configuration types must contain no spaces and should be all
/// lowercase letters. The configuration name should be provided after the input
/// file on the command line.
class SimulationInitializer {
 public:
  SimulationInitializer(int argc, char** argv, MPIParallel& a_mpi_session,
                        SpdlogLevel a_log_level = SpdlogLevel::OFF);

  ~SimulationInitializer(void);

 private:
  MPIParallel& mpi_session_m;
  InputParser parser_m;
  ConfigurationInitializer* configuration_initializer_m;
};

}  // namespace chyps

#endif  // CHYPS_SIMULATION_INITIALIZER_HPP_
