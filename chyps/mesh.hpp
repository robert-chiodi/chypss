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

#include <cstddef>
#include <vector>

#include <mfem/mfem.hpp>

#include <mpi.h>

#include "chyps/boundary_condition.hpp"
#include "chyps/boundary_condition_manager.hpp"
#include "chyps/input_parser.hpp"
#include "chyps/mpi_parallel.hpp"

namespace chyps {

// Forward declare IO to prevent cyclical reference
class IO;

enum class MeshElement { ELEMENT = 0, VERTEX };

class Mesh {
 public:
  Mesh(void) = delete;

  /// \brief Initialize Mesh for parallel use and collect options in a_parser.
  Mesh(const MPIParallel& a_mpi_session, InputParser& a_parser, IO* a_file_io);

  /// \brief Construct mesh for use by flow solvers.
  void Initialize(void);

  /// \brief Return MPI Communicator used for this mesh.
  const MPI_Comm& GetMPIComm(void) const;

  /// \brief Return the underlying MFEM mesh.
  mfem::ParMesh& GetMfemMesh(void);

  /// \brief Return spatial dimension of the mesh.
  int GetDimension(void) const;

  /// \brief Whether this mesh uses AMR or not. Currently AMR is not
  /// implemented, so this is always false.
  static constexpr bool UsesAMR(void) { return false; }

  template <enum MeshElement>
  std::size_t GetGlobalCount(void) const;
  template <enum MeshElement>
  std::size_t GetOffsetStart(void) const;
  template <enum MeshElement>
  std::size_t GetLocalCount(void) const;

  /// \brief Return the number of boundary tag values that exist on the mesh.
  ///
  /// It is assumed that all tags from [1,ntags] are valid and exist on the
  /// mesh.
  int GetNumberOfBoundaryTagValues(void) const;

  /// \brief Function computes and returns a std::vector of vertex positions for
  /// border vertices in a_mesh tagged with a_tag. Number of vertices is
  /// vector.size()/Mesh::GetDimension();
  std::pair<std::vector<double>, std::vector<int>> GetBoundaryVertices(
      const int a_tag) const;

  /// \brief Writes the mesh to the file opened in file_io_m.
  ///
  /// Note: Mesh must be initialized before calling this and
  /// file_io_m must be in Write mode.
  void WriteMesh(void);

  /// \brief Destructor to free all heap-allocated objects.
  ~Mesh(void);

 private:
  void GatherOptions(void);
  void AddIOVariables(void);
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

  bool FileWritingEnabled(void) const;

  uint32_t GLVISToVTKType(const int glvisType) const noexcept;

  InputParser& parser_m;
  const MPIParallel& mpi_session_m;
  IO* file_io_m;
  mfem::ParMesh* parallel_mesh_m;
  std::size_t element_offset_m;
};

}  // namespace chyps

#endif  // CHYPS_MESH_HPP_
