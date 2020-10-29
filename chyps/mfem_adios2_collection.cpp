// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chyps/mfem_adios2_collection.hpp"
#include "chyps/logger.hpp"

namespace chyps {

MfemAdios2Collection::MfemAdios2Collection(const MPI_Comm& a_mpi_comm,
                                           const std::string& a_collection_name,
                                           mfem::ParMesh& a_mesh,
                                           const std::string& a_engine_type)
    : adios2_collection_m(a_mpi_comm, a_collection_name, &a_mesh,
                          a_engine_type),
      name_data_map_m() {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Made new ADIOS2 collection: {}",
                     a_collection_name);
  // adios2_collection_m.SetLevelsOfDetail(1);
  //  adios2_collection_m.SetParameter("RefinedData", "false");
  // adios2_collection_m.SetParameter("FullData", "true");
}

void MfemAdios2Collection::RegisterField(
    const std::string& a_name, mfem::ParFiniteElementSpace* a_element_space) {
  auto grid_function = new mfem::ParGridFunction(a_element_space);
  // FIXME : Make this an exception
  assert(name_data_map_m.find(a_name) ==
         name_data_map_m.end());  // Make sure field with name not already added
  name_data_map_m[a_name] = grid_function;
  name_update_map_m[a_name] = false;
  adios2_collection_m.RegisterField(a_name, grid_function);
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Registered field {} in ADIOS2 collection {}",
                     a_name, adios2_collection_m.GetCollectionName());
}

void MfemAdios2Collection::UpdateField(const std::string& a_name,
                                       const mfem::Vector& a_data) {
  assert(name_data_map_m.find(a_name) !=
         name_data_map_m.end());  // Make sure field with name exists
  name_data_map_m[a_name]->SetFromTrueDofs(a_data);
  name_update_map_m[a_name] = true;
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Updated field {} in ADIOS2 collection {}",
                     a_name, adios2_collection_m.GetCollectionName());
}

void MfemAdios2Collection::WriteOutFields(const int a_cycle,
                                          const double a_time) {
  assert(this->AllFieldsUpdatedSinceLastWrite());
  adios2_collection_m.SetCycle(a_cycle);
  adios2_collection_m.SetTime(a_time);
  adios2_collection_m.Save();
  for (auto& elem : name_update_map_m) {
    elem.second = false;
  }
  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Wrote ADIOS2 collection {} at cycle {} and time {:8.6E}",
                     adios2_collection_m.GetCollectionName(), a_cycle, a_time);
}

MfemAdios2Collection::~MfemAdios2Collection(void) {
  for (auto& elem : name_data_map_m) {
    delete elem.second;
    elem.second = nullptr;
  }
}

bool MfemAdios2Collection::AllFieldsUpdatedSinceLastWrite(void) const {
  for (const auto& elem : name_update_map_m) {
    if (!elem.second) {
      return false;
    }
  }
  return true;
}

}  // namespace chyps
