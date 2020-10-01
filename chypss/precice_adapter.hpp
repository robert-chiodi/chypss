// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef CHYPSS_PRECICE_ADAPTER_H_
#define CHYPSS_PRECICE_ADAPTER_H_

#include <string>
#include <unordered_map>
#include <vector>

#include <precice/SolverInterface.hpp>

namespace chyps {

enum class DataOperation { READ = 0, WRITE };

class PreciceData {
 public:
  PreciceData(void) = default;
  PreciceData(const int a_data_id, const DataOperation a_operation);

  int GetDataID(void) const;
  bool ForReading(void) const;
  bool ForWriting(void) const;

 private:
  int data_id_m;
  DataOperation data_operation_m;
};

class PreciceAdapter {
 public:
  PreciceAdapter(void) = delete;
  PreciceAdapter(const std::string& a_solver_name,
                 const std::string& a_config_file_name, const int a_proc_rank,
                 const int a_proc_size);

  double Initialize(void);

  void SetVertexPositions(const std::vector<double>& a_positions);

  void AddData(const std::string& a_name, const DataOperation a_operation);

  void WriteBlockVectorData(const std::string& a_name, const double* a_data);

  void ReadBlockVectorData(const std::string& a_name, double* a_data);

  int NumberOfVertices(void) const;

  ~PreciceAdapter(void);

 private:
  precice::SolverInterface interface_m;
  std::vector<int> vertex_ids_m;
  std::unordered_map<std::string, PreciceData> data_m;
  int dimension_m;
  int mesh_id_m;
};

}  // namespace chyps

#endif  // CHYPSS_PRECICE_ADAPTER_H_
