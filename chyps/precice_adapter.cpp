// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chyps/precice_adapter.hpp"

#include <cassert>

#include "chyps/logger.hpp"

namespace chyps {

PreciceData::PreciceData(const int a_data_id, const DataOperation a_operation)
    : data_id_m(a_data_id), data_operation_m(a_operation) {
  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Constructed preCICE data with data_id {} for operator {}",
                     a_data_id, a_operation);
}

int PreciceData::GetDataID(void) const { return data_id_m; }

bool PreciceData::ForReading(void) const {
  return data_operation_m == DataOperation::READ;
}

bool PreciceData::ForWriting(void) const {
  return data_operation_m == DataOperation::WRITE;
}

PreciceAdapter::PreciceAdapter(const std::string& a_solver_name,
                               const std::string& a_mesh_name,
                               const std::string& a_config_file_name,
                               const int a_proc_rank, const int a_proc_size)
    : interface_m(a_solver_name, a_config_file_name, a_proc_rank, a_proc_size) {
  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Constructed preCICE adapter for {} and read config file "
                     "{}. This is rank {} of world size {}",
                     a_solver_name, a_config_file_name, a_proc_rank,
                     a_proc_size);
  dimension_m = interface_m.getDimensions();
  mesh_id_m = interface_m.getMeshID(a_mesh_name);
  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Registered mesh {} and given ID {}. This mesh has {} "
                     "dimensions.",
                     a_mesh_name, mesh_id_m, dimension_m);
}

double PreciceAdapter::Initialize(void) {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Initializing PreciceAdapter...");
  return interface_m.initialize();
  SPDLOG_LOGGER_INFO(MAIN_LOG, "PreciceAdapter initialized.");
}

void PreciceAdapter::SetVertexPositions(
    const std::vector<double>& a_positions) {
  assert(a_positions.size() % dimension_m == 0);
  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Setting vertex positions for {} dimensional mesh.",
                     dimension_m);
  vertex_ids_m.resize(a_positions.size() / dimension_m);
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Setting the position of {} vertices.",
                     vertex_ids_m.size());
  interface_m.setMeshVertices(mesh_id_m, a_positions.size() / dimension_m,
                              a_positions.data(), vertex_ids_m.data());
  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Vertex position successfully copied to preCICE.");
}

int PreciceAdapter::NumberOfVertices(void) const {
  return static_cast<int>(vertex_ids_m.size());
}

bool PreciceAdapter::IsCouplingOngoing(void) const {
  return interface_m.isCouplingOngoing();
}

void PreciceAdapter::AddData(const std::string& a_name,
                             const DataOperation a_operation) {
  assert(data_m.find(a_name) ==
         data_m.end());  // Data not added with name already
  const int data_id = interface_m.getDataID(a_name, mesh_id_m);
  data_m.emplace(a_name, PreciceData(data_id, a_operation));
  SPDLOG_LOGGER_INFO(
      MAIN_LOG,
      "Registered data for {} on mesh with id {}. Data given id of {}", a_name,
      mesh_id_m, data_id);
}

std::size_t PreciceAdapter::ScalarDataSize(void) const {
  return this->NumberOfVertices();
}

std::size_t PreciceAdapter::VectorDataSize(void) const {
  return this->NumberOfVertices() * dimension_m;
}

void PreciceAdapter::WriteBlockScalarData(const std::string& a_name,
                                          const double* a_data) {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Writing data for {} to preCICE", a_name);
  assert(data_m.find(a_name) != data_m.end());
  const auto& data_details = data_m[a_name];
  assert(data_details.ForWriting());
  interface_m.writeBlockScalarData(data_details.GetDataID(),
                                   this->ScalarDataSize(), vertex_ids_m.data(),
                                   a_data);
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Wrote data {} to preCICE", a_name);
}

void PreciceAdapter::ReadBlockScalarData(const std::string& a_name,
                                         double* a_data) {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Reading data for {} from preCICE", a_name);
  assert(data_m.find(a_name) != data_m.end());
  const auto data_details = data_m[a_name];
  assert(data_details.ForReading());
  interface_m.readBlockScalarData(data_details.GetDataID(),
                                  this->ScalarDataSize(), vertex_ids_m.data(),
                                  a_data);
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Read data {} from preCICE", a_name);
}

void PreciceAdapter::WriteBlockVectorData(const std::string& a_name,
                                          const double* a_data) {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Writing data for {} to preCICE", a_name);
  assert(data_m.find(a_name) != data_m.end());
  const auto& data_details = data_m[a_name];
  assert(data_details.ForWriting());
  interface_m.writeBlockVectorData(data_details.GetDataID(),
                                   this->VectorDataSize(), vertex_ids_m.data(),
                                   a_data);
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Wrote data {} to preCICE", a_name);
}

void PreciceAdapter::ReadBlockVectorData(const std::string& a_name,
                                         double* a_data) {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Reading data for {} from preCICE", a_name);
  assert(data_m.find(a_name) != data_m.end());
  const auto data_details = data_m[a_name];
  assert(data_details.ForReading());
  interface_m.readBlockVectorData(data_details.GetDataID(),
                                  this->VectorDataSize(), vertex_ids_m.data(),
                                  a_data);
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Read data {} from preCICE", a_name);
}

double PreciceAdapter::Advance(const double a_dt) {
  return interface_m.advance(a_dt);
}

void PreciceAdapter::Finalize(void) {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Finalizing PreciceAdapter");
  interface_m.finalize();
}

}  // namespace chyps
