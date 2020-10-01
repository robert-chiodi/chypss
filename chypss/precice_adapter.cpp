// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chypss/precice_adapter.hpp"

#include <cassert>

namespace chyps {

PreciceData::PreciceData(const int a_data_id, const DataOperation a_operation)
    : data_id_m(a_data_id), data_operation_m(a_operation) {}

int PreciceData::GetDataID(void) const { return data_id_m; }

bool PreciceData::ForReading(void) const {
  return data_operation_m == DataOperation::READ;
}

bool PreciceData::ForWriting(void) const {
  return data_operation_m == DataOperation::WRITE;
}

PreciceAdapter::PreciceAdapter(const std::string& a_solver_name,
                               const std::string& a_config_file_name,
                               const int a_proc_rank, const int a_proc_size)
    : interface_m(a_solver_name, a_config_file_name, a_proc_rank, a_proc_size) {
  dimension_m = interface_m.getDimensions();
  const std::string mesh_name("mfem_mesh");
  mesh_id_m = interface_m.getMeshID(mesh_name);
}

double PreciceAdapter::Initialize(void) { return interface_m.initialize(); }

void PreciceAdapter::SetVertexPositions(
    const std::vector<double>& a_positions) {
  assert(a_positions.size() % dimension_m == 0);
  vertex_ids_m.resize(a_positions.size() / dimension_m);
  interface_m.setMeshVertices(mesh_id_m, a_positions.size() / dimension_m,
                              a_positions.data(), vertex_ids_m.data());
}

void PreciceAdapter::AddData(const std::string& a_name,
                             const DataOperation a_operation) {
  assert(data_m.find(a_name) ==
         data_m.end());  // Data not added with name already
  const int data_id = interface_m.getDataID(a_name, mesh_id_m);
  data_m.emplace(a_name, PreciceData(data_id, a_operation));
}

void PreciceAdapter::WriteBlockVectorData(const std::string& a_name,
                                          const double* a_data) {
  assert(data_m.find(a_name) != data_m.end());
  const auto& data_details = data_m[a_name];
  assert(data_details.ForWriting());
  interface_m.writeBlockVectorData(data_details.GetDataID(),
                                   this->NumberOfVertices(),
                                   vertex_ids_m.data(), a_data);
}

void PreciceAdapter::ReadBlockVectorData(const std::string& a_name,
                                         double* a_data) {
  assert(data_m.find(a_name) != data_m.end());
  const auto data_details = data_m[a_name];
  assert(data_details.ForReading());
  interface_m.readBlockVectorData(data_details.GetDataID(),
                                  this->NumberOfVertices(), vertex_ids_m.data(),
                                  a_data);
}

int PreciceAdapter::NumberOfVertices(void) const {
  return static_cast<int>(vertex_ids_m.size());
}

PreciceAdapter::~PreciceAdapter(void) { interface_m.finalize(); }

}  // namespace chyps
