// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// \file precice_adapter.hpp
/// \brief File handling all use of Precice for coupling between codes.

#ifndef CHYPS_PRECICE_ADAPTER_H_
#define CHYPS_PRECICE_ADAPTER_H_

#include <string>
#include <unordered_map>
#include <vector>

#include <precice/SolverInterface.hpp>

namespace chyps {

/// enum DataOperation precice_adapter.hpp chyps/precice_adapter.hpp
/// \brief Simple enum to denote whether Precice data should be written to or
/// read into the supplied buffer.
enum class DataOperation { READ = 0, WRITE };

/// \class PreciceData precice_adapter.hpp chyps/precice_adapter.hpp
/// \brief Wrapper for data information used in Precice for exchaning data
/// between solvers.
class PreciceData {
 public:
  /// \brief Default constructor.
  PreciceData(void) = default;

  /// \brief Constructor that stores reference ID used in Precice and whether
  /// data is for reading or writing.
  PreciceData(const int a_data_id, const DataOperation a_operation);

  /// \brief Return data ID used by Precice.
  int GetDataID(void) const;

  /// \brief Return true if data is intended to be read.
  bool ForReading(void) const;

  /// \brief Return true if data is intended to be written.
  bool ForWriting(void) const;

 private:
  int data_id_m;
  DataOperation data_operation_m;
};

/// \class PreciceAdapter precice_adapter.hpp chyps/precice_adapter.hpp
/// \brief Adapter to abstract away use of Precice for coupling of two separate
/// solvers.
class PreciceAdapter {
 public:
  /// \brief Do not allow default construction since precice::SolverInterface
  /// does not allow it.
  PreciceAdapter(void) = delete;

  /// \brief Constructor that sets up initial parallel interface in Precice.
  ///
  /// Assumed that the mesh used will be "mfem_mesh", which must match in the
  /// Precice XML config file.
  PreciceAdapter(const std::string& a_solver_name,
                 const std::string& a_config_file_name, const int a_proc_rank,
                 const int a_proc_size);

  /// \brief Initializes Precice SolverInterface and returns initial timestep.
  double Initialize(void);

  /// \brief Set the vertices for the surface mesh coupling solvers in Precice.
  ///
  /// The number of vertices will be assumed to be a_positions.size().
  /// Each individual vertex should be supplied as dimension_m doubles.
  /// For example, if the mesh consistents of 10 vertices for a 2D surface,
  /// a_positions.size() = 20, ordered as V0_x V0_y V1_x V1_y ...
  void SetVertexPositions(const std::vector<double>& a_positions);

  /// \brief Return the number of vertices that exist in the coupling mesh.
  int NumberOfVertices(void) const;

  /// \brief Add Data to be tracked by Precice. This does not do any update
  /// itself, but registers it inside Precice.
  void AddData(const std::string& a_name, const DataOperation a_operation);

  /// \brief Write the data to be read by the other solver(s).
  ///
  /// The supplied a_name must match one previously added for writing via
  /// AddData with DataOperation::WRITE. The values written will be taken from
  /// a_data. The size of a_data must be atleast this->NumberOfVertices() long.
  void WriteBlockVectorData(const std::string& a_name, const double* a_data);

  /// \brief Read the data written by the other solver(s).
  ///
  /// The supplied a_name must match one previously added for reading via
  /// AddData with DataOperation::READ. The values will be placed in
  /// a_data. The size of a_data must be atleast this->NumberOfVertices() long.
  void ReadBlockVectorData(const std::string& a_name, double* a_data);

  /// \brief Custom destructor that finalizes the precice::SolverInterface,
  /// freeing any used memory by Precice.
  ~PreciceAdapter(void);

 private:
  precice::SolverInterface interface_m;
  std::vector<int> vertex_ids_m;
  std::unordered_map<std::string, PreciceData> data_m;
  int dimension_m;
  int mesh_id_m;
};

}  // namespace chyps

#endif  // CHYPS_PRECICE_ADAPTER_H_
