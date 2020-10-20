// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef CHYPS_MESH_HPP_
#define CHYPS_MESH_HPP_

#include <vector>

#include <mfem/mfem.hpp>

#include <mpi.h>

#include "chyps/boundary_condition.hpp"
#include "chyps/input_parser.hpp"

namespace chyps {

class Mesh {
 public:
  Mesh(void) = delete;

  /// \brief Initialize Mesh for parallel use and collect options in a_parser.
  Mesh(const MPI_Comm& a_mpi_comm, InputParser& a_parser);

  /// \brief Construct mesh for use by flow solvers.
  void Initialize(void);

  /// \brief Return MPI Communicator used for this mesh.
  const MPI_Comm& GetMPIComm(void) const;

  /// \brief Return the underlying MFEM mesh.
  mfem::ParMesh& GetMfemMesh(void);

  /// \brief Return spatial dimension of the mesh.
  int GetDimension(void) const;

  /// \brief Return number of boundary conditions.
  int GetNumberOfBoundaryConditions(void) const;

  /// \brief Return boundary condition for boundary element tag a_tag.
  const BoundaryCondition& GetBoundaryCondition(const int a_tag) const;

  /// \brief Return number of homogeneous Dirichlet conditions set.
  ///
  /// Note: Must be called after Mesh::CommitBoundaryConditions.
  int GetNumberOfHomogeneousDirichletConditions(void) const;

  /// \brief Return number of inhomogeneous Dirichlet conditions set.
  ///
  /// Note: Must be called after Mesh::CommitBoundaryConditions.
  int GetNumberOfDirichletConditions(void) const;

  /// \brief Return number of time-varying Dirichlet conditions.
  ///
  /// Note: This must be <= the number of inhomogeneous Dirichlet conditions.
  int GetNumberOfTimeVaryingDirichletConditions(void) const;

  /// \brief Return number of homogeneous Neumann conditions set.
  ///
  /// Note: Must be called after Mesh::CommitBoundaryConditions.
  int GetNumberOfHomogeneousNeumannConditions(void) const;

  /// \brief Return number of inhomogeneous Neumann conditions set.
  ///
  /// Note: Must be called after Mesh::CommitBoundaryConditions.
  int GetNumberOfNeumannConditions(void) const;

  /// \brief Return number of time-varying Neumann conditions.
  ///
  /// Note: This must be <= the number of inhomogeneous Neumann conditions.
  int GetNumberOfTimeVaryingNeumannConditions(void) const;

  /// \brief Function computes and returns a std::vector of vertex positions for
  /// border vertices in a_mesh tagged with a_tag. Number of vertices is
  /// vector.size()/Mesh::GetDimension();
  std::pair<std::vector<double>, std::vector<int>> GetBoundaryVertices(
      const int a_tag);

  /// \brief Sets boundary condition for boundary elements tagged with a_tag.
  ///
  /// These boundary conditions are stored for calling again later and passed to
  /// a solver to be used.
  void SetBoundaryCondition(const int a_tag,
                            const BoundaryCondition& a_condition);

  /// \brief Commits the boundary conditions used for this mesh.
  ///
  /// If changes are made to the boundary conditions (beyond the values
  /// at the values) they must be recommitted by calling
  /// CommitBoundaryConditions.
  void CommitBoundaryConditions(void);

  /// \brief Destructor to free all heap-allocated objects.
  ~Mesh(void);

 private:
  bool AllOptionsSupplied(void);
  void GatherOptions(void);
  void ReadAndRefineMesh(void);
  void AllocateVariables(void);

  /// \brief Generate a line mesh of a 1D domain for use by MFEM.
  ///
  /// This function returns a pointer to an MFEM mesh object and an array of
  /// vertices. Neither are stored and both must be deleted explicitly by the
  /// user.
  std::pair<mfem::Mesh*, double*> GenerateLineMesh(
      const std::array<std::array<double, 1>, 2>& a_bounding_box,
      const int a_nx);

  /// \brief Generate quad mesh of a rectangular domain for use by MFEM.
  ///
  /// This function returns a pointer to an MFEM mesh object and an array of
  /// vertices. Neither are stored and both must be deleted explicitly by the
  /// user.
  std::pair<mfem::Mesh*, double*> GenerateQuadMesh(
      const std::array<std::array<double, 2>, 2>& a_bounding_box,
      const int a_nx, const int a_ny);

  /// \brief Generate hex mesh of a cuboid domain for use by MFEM.
  ///
  /// This function returns a pointer to an MFEM mesh object and an array of
  /// vertices. Neither are stored and both must be deleted explicitly by the
  /// user.
  std::pair<mfem::Mesh*, double*> GenerateHexMesh(
      const std::array<std::array<double, 3>, 2>& a_bounding_box,
      const int a_nx, const int a_ny, const int a_nz);

  InputParser& parser_m;
  const MPI_Comm& mpi_comm_m;
  mfem::ParMesh* parallel_mesh_m;
  std::vector<BoundaryCondition> boundary_conditions_m;
  std::vector<int> boundary_condition_counts_m;
};

}  // namespace chyps

#endif  // CHYPS_MESH_HPP_
