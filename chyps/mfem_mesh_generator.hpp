// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <utility>

#include <mfem/mfem.hpp>

#ifndef CHYPS_MFEM_MESH_GENERATOR_HPP_
#define CHYPS_MFEM_MESH_GENERATOR_HPP_

namespace chyps {

/// \function GenerateLineMesh mfem_mesh_generator.hpp
/// \brief Generate a line mesh of a 1D domain for use by MFEM.
///
/// This function returns a pointer to an MFEM mesh object and an array of
/// vertices. Neither are stored and both must be deleted explicitly by the
/// user.
std::pair<mfem::Mesh*, double*> GenerateLineMesh(
    const std::array<std::array<double, 1>, 2>& a_bounding_box, const int a_nx);

/// \function GenerateQuadMesh mfem_mesh_generator.hpp
/// \brief Generate quad mesh of a rectangular domain for use by MFEM.
///
/// This function returns a pointer to an MFEM mesh object and an array of
/// vertices. Neither are stored and both must be deleted explicitly by the
/// user.
std::pair<mfem::Mesh*, double*> GenerateQuadMesh(
    const std::array<std::array<double, 2>, 2>& a_bounding_box, const int a_nx,
    const int a_ny);

/// \function GenerateHexMesh mfem_mesh_generator.hpp
/// \brief Generate hex mesh of a cuboid domain for use by MFEM.
///
/// This function returns a pointer to an MFEM mesh object and an array of
/// vertices. Neither are stored and both must be deleted explicitly by the
/// user.
std::pair<mfem::Mesh*, double*> GenerateHexMesh(
    const std::array<std::array<double, 3>, 2>& a_bounding_box, const int a_nx,
    const int a_ny, const int a_nz);

}  // namespace chyps

#endif  // CHYPS_MFEM_MESH_GENERATOR_HPP_
