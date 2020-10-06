// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chyps/mfem_visit_collection.hpp"

namespace chyps {

MfemVisItCollection::MfemVisItCollection(const MPI_Comm& a_mpi_comm,
                                         const std::string& a_collection_name,
                                         mfem::ParMesh& a_mesh)
    : visit_collection_m(a_mpi_comm, a_collection_name, &a_mesh),
      name_data_map_m() {}

void MfemVisItCollection::RegisterField(
    const std::string& a_name, mfem::ParFiniteElementSpace* a_element_space) {
  auto grid_function = new mfem::ParGridFunction(a_element_space);
  // FIXME : Make this an exception
  assert(name_data_map_m.find(a_name) ==
         name_data_map_m.end());  // Make sure field with name not already added
  name_data_map_m[a_name] = grid_function;
  name_update_map_m[a_name] = false;
  visit_collection_m.RegisterField(a_name, grid_function);
}

void MfemVisItCollection::UpdateField(const std::string& a_name,
                                      const mfem::Vector& a_data) {
  assert(name_data_map_m.find(a_name) !=
         name_data_map_m.end());  // Make sure field with name exists
  name_data_map_m[a_name]->SetFromTrueDofs(a_data);
  name_update_map_m[a_name] = true;
}

void MfemVisItCollection::WriteOutFields(const int a_cycle,
                                         const double a_time) {
  assert(this->AllFieldsUpdatedSinceLastWrite());
  visit_collection_m.SetCycle(a_cycle);
  visit_collection_m.SetTime(a_time);
  visit_collection_m.Save();
  for (auto& elem : name_update_map_m) {
    elem.second = false;
  }
}

MfemVisItCollection::~MfemVisItCollection(void) {
  for (auto& elem : name_data_map_m) {
    delete elem.second;
    elem.second = nullptr;
  }
}

bool MfemVisItCollection::AllFieldsUpdatedSinceLastWrite(void) const {
  for (const auto& elem : name_update_map_m) {
    if (!elem.second) {
      return false;
    }
  }
  return true;
}

}  // namespace chyps