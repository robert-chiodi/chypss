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

#include "chyps/debug_assert.hpp"
#include "chyps/logger.hpp"

namespace chyps {

PreciceMesh::PreciceMesh(const int a_mesh_id)
    : mesh_id_m(a_mesh_id), vertex_ids_m() {}

void PreciceMesh::SetVertexPositions(const std::vector<double>& a_positions,
                                     const int a_dimension,
                                     precice::SolverInterface& a_interface) {
  DEBUG_ASSERT(a_positions.size() % a_dimension == 0, global_assert{},
               DebugLevel::CHEAP{},
               "Vector of vertex positions must be of incorrect size.");
  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Setting vertex positions for {} dimensional mesh.",
                     a_dimension);
  vertex_ids_m.resize(a_positions.size() / a_dimension);
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Setting the position of {} vertices.",
                     vertex_ids_m.size());
  a_interface.setMeshVertices(mesh_id_m, static_cast<int>(vertex_ids_m.size()),
                              a_positions.data(), vertex_ids_m.data());
  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Vertex position successfully copied to preCICE.");
}

const std::vector<int>& PreciceMesh::GetVertexIds(void) const {
  DEBUG_ASSERT(
      vertex_ids_m.size() > 0, global_assert{}, DebugLevel::CHEAP{},
      "Vertices not set for mesh with id " + std::to_string(mesh_id_m));
  return vertex_ids_m;
}

std::size_t PreciceMesh::NumberOfVertices(void) const {
  return vertex_ids_m.size();
}

int PreciceMesh::GetId(void) const { return mesh_id_m; }

PreciceData::PreciceData(const int a_data_id, const DataOperation a_operation,
                         const PreciceMesh& a_mesh)
    : associated_mesh_m(a_mesh),
      data_id_m(a_data_id),
      data_operation_m(a_operation) {
  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Constructed preCICE data with data_id {} for operator {}",
                     a_data_id, a_operation);
}

int PreciceData::GetId(void) const { return data_id_m; }

const PreciceMesh& PreciceData::GetMesh(void) const {
  return associated_mesh_m;
}

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
  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Constructed preCICE adapter for {} and read config file "
                     "{}. This is rank {} of world size {}",
                     a_solver_name, a_config_file_name, a_proc_rank,
                     a_proc_size);
  dimension_m = interface_m.getDimensions();
  SPDLOG_LOGGER_INFO(
      MAIN_LOG, "preCICE adapter setup for mesh of dimension {}.", dimension_m);
}

void PreciceAdapter::AddMesh(const std::string& a_mesh_name) {
  DEBUG_ASSERT(meshes_m.find(a_mesh_name) == meshes_m.end(), global_assert{},
               DebugLevel::CHEAP{},
               "Mesh with name \"" + a_mesh_name + "\" already added.");
  const int mesh_id = interface_m.getMeshID(a_mesh_name);
  meshes_m.emplace(a_mesh_name, PreciceMesh(mesh_id));
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Registered mesh with name {}. Given Id of {}",
                     a_mesh_name, mesh_id);
}

void PreciceAdapter::SetVertexPositions(
    const std::string& a_mesh_name, const std::vector<double>& a_positions) {
  DEBUG_ASSERT(meshes_m.find(a_mesh_name) != meshes_m.end(), global_assert{},
               DebugLevel::CHEAP{},
               "Mesh with name \"" + a_mesh_name + "\" does not exist.");

  meshes_m[a_mesh_name].SetVertexPositions(a_positions, dimension_m,
                                           interface_m);
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Set {} vertices for mesh {}",
                     a_positions.size() / dimension_m, a_mesh_name);
}

double PreciceAdapter::Initialize(void) {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Initializing PreciceAdapter...");
  return interface_m.initialize();
}

bool PreciceAdapter::IsCouplingOngoing(void) const {
  return interface_m.isCouplingOngoing();
}

void PreciceAdapter::AddData(const std::string& a_name,
                             const std::string& a_mesh_name,
                             const DataOperation a_operation) {
  DEBUG_ASSERT(
      data_m.find(a_name) == data_m.end(), global_assert{}, DebugLevel::CHEAP{},
      "Data with name \"" + a_name + "\"already added to preCICE adapter.");
  DEBUG_ASSERT(meshes_m.find(a_mesh_name) != meshes_m.end(), global_assert{},
               DebugLevel::CHEAP{},
               "Mesh with name \"" + a_mesh_name + "\" does not exist.");
  const auto& mesh = meshes_m[a_mesh_name];
  const int data_id = interface_m.getDataID(a_name, mesh.GetId());
  data_m.emplace(a_name, PreciceData(data_id, a_operation, mesh));
  SPDLOG_LOGGER_INFO(
      MAIN_LOG, "Registered data for {} on mesh named {}. Data given id of {}",
      a_name, a_mesh_name, data_id);
}

void PreciceAdapter::WriteBlockScalarData(const std::string& a_name,
                                          const double* a_data) {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Writing data for {} to preCICE", a_name);
  DEBUG_ASSERT(data_m.find(a_name) != data_m.end(), global_assert{},
               DebugLevel::CHEAP{},
               "Data with name \"" + a_name + "\" has not been added.");
  const auto& data_details = data_m.at(a_name);
  DEBUG_ASSERT(data_details.ForWriting(), global_assert{}, DebugLevel::CHEAP{},
               "Data with name \"" + a_name + "\" not marked for writing.");
  const auto& associated_mesh = data_details.GetMesh();
  interface_m.writeBlockScalarData(
      data_details.GetId(), associated_mesh.NumberOfVertices(),
      associated_mesh.GetVertexIds().data(), a_data);
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Wrote data {} to preCICE", a_name);
}

void PreciceAdapter::ReadBlockScalarData(const std::string& a_name,
                                         double* a_data) {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Reading data for {} from preCICE", a_name);
  DEBUG_ASSERT(data_m.find(a_name) != data_m.end(), global_assert{},
               DebugLevel::CHEAP{},
               "Data with name \"" + a_name + "\" has not been added.");
  const auto& data_details = data_m.at(a_name);
  DEBUG_ASSERT(data_details.ForReading(), global_assert{}, DebugLevel::CHEAP{},
               "Data with name \"" + a_name + "\" not marked for reading.");
  const auto& associated_mesh = data_details.GetMesh();
  interface_m.readBlockScalarData(
      data_details.GetId(), associated_mesh.NumberOfVertices(),
      associated_mesh.GetVertexIds().data(), a_data);
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Read data {} from preCICE", a_name);
}

void PreciceAdapter::WriteBlockVectorData(const std::string& a_name,
                                          const double* a_data) {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Writing data for {} to preCICE", a_name);
  DEBUG_ASSERT(data_m.find(a_name) != data_m.end(), global_assert{},
               DebugLevel::CHEAP{},
               "Data with name \"" + a_name + "\" has not been added.");
  const auto& data_details = data_m.at(a_name);
  DEBUG_ASSERT(data_details.ForWriting(), global_assert{}, DebugLevel::CHEAP{},
               "Data with name \"" + a_name + "\" not marked for writing.");
  const auto& associated_mesh = data_details.GetMesh();
  interface_m.writeBlockVectorData(
      data_details.GetId(), dimension_m * associated_mesh.NumberOfVertices(),
      associated_mesh.GetVertexIds().data(), a_data);
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Wrote data {} to preCICE", a_name);
}

void PreciceAdapter::ReadBlockVectorData(const std::string& a_name,
                                         double* a_data) {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Reading data for {} from preCICE", a_name);
  DEBUG_ASSERT(data_m.find(a_name) != data_m.end(), global_assert{},
               DebugLevel::CHEAP{},
               "Data with name \"" + a_name + "\" has not been added.");
  const auto& data_details = data_m.at(a_name);
  DEBUG_ASSERT(data_details.ForReading(), global_assert{}, DebugLevel::CHEAP{},
               "Data with name \"" + a_name + "\" not marked for reading.");
  const auto& associated_mesh = data_details.GetMesh();
  interface_m.readBlockVectorData(
      data_details.GetId(), dimension_m * associated_mesh.NumberOfVertices(),
      associated_mesh.GetVertexIds().data(), a_data);
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
