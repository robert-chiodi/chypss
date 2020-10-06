// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// \file mfem_visit_collection.hpp
/// \brief File handling export of variables in the MFEM framework for viewing
/// in VisIt.

#ifndef CHYPS_MFEM_VISIT_COLLECTION_HPP_
#define CHYPS_MFEM_VISIT_COLLECTION_HPP_

#include <string>
#include <unordered_map>

#include <mfem/mfem.hpp>

namespace chyps {

/// \class MfemVisItCollection mfem_visit_collection.hpp
/// chyps/mfem_visit_collection.hpp
/// \brief Class that manages organization of variables and writing to disk for
/// viewing in VisIt. Leans heavily on MFEM's VisItDataCollection.
///
///
class MfemVisItCollection {
 public:
  MfemVisItCollection(const MPI_Comm& a_mpi_comm,
                      const std::string& a_collection_name,
                      mfem::ParMesh& a_mesh);

  void RegisterField(const std::string& a_name,
                     mfem::ParFiniteElementSpace* a_element_space);
  void UpdateField(const std::string& a_name, const mfem::Vector& a_data);

  void WriteOutFields(const int a_cycle, const double a_time);

  ~MfemVisItCollection(void);

 private:
  bool AllFieldsUpdatedSinceLastWrite(void) const;

  mfem::VisItDataCollection visit_collection_m;
  std::unordered_map<std::string, mfem::ParGridFunction*> name_data_map_m;
  std::unordered_map<std::string, bool> name_update_map_m;
};

}  // namespace chyps

#endif  // CHYPS_MFEM_VISIT_COLLECTION_HPP_