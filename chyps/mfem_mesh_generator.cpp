// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "mfem_mesh_generator.hpp"

#include <vector>

std::pair<mfem::Mesh*, double*> GenerateQuadMesh(
    const std::array<std::array<double, 2>, 2>& a_bounding_box, const int a_nx,
    const int a_ny) {
  const int nvert = (a_nx + 1) * (a_ny + 1);
  const int nelem = a_nx * a_ny;
  double* vertices = new double[nvert * 3];

  std::cout << nelem << std::endl;

  std::vector<int> elem_indices(4 * nelem);

  std::vector<int> elem_attrib(nelem);
  std::fill(elem_attrib.begin(), elem_attrib.end(), 1);

  const int nboundary = 2 * a_nx + 2 * a_ny;
  std::vector<int> boundary_indices(2 * nboundary);
  std::vector<int> boundary_attrib(nboundary);
  int nbi = 0;
  for (int j = 0; j < a_ny; ++j) {
    for (int i = 0; i < a_nx; ++i) {
      const int elem_index = j * a_nx + i;
      elem_indices[4 * elem_index] = i + j * (a_nx + 1);
      elem_indices[4 * elem_index + 1] = elem_indices[4 * elem_index] + 1;
      elem_indices[4 * elem_index + 2] =
          elem_indices[4 * elem_index] + (a_nx + 1) + 1;
      elem_indices[4 * elem_index + 3] =
          elem_indices[4 * elem_index] + (a_nx + 1);
      if (i == 0) {
        boundary_attrib[nbi] = 1;
        boundary_indices[2 * nbi] = i + j * (a_nx + 1);
        boundary_indices[2 * nbi + 1] = boundary_indices[2 * nbi] + (a_nx + 1);
        ++nbi;
      }
      if (i == a_nx - 1) {
        boundary_attrib[nbi] = 2;
        boundary_indices[2 * nbi] = i + j * (a_nx + 1) + 1;
        boundary_indices[2 * nbi + 1] = boundary_indices[2 * nbi] + (a_nx + 1);
        ++nbi;
      }
      if (j == 0) {
        boundary_attrib[nbi] = 3;
        boundary_indices[2 * nbi] = i + j * (a_nx + 1);
        boundary_indices[2 * nbi + 1] = boundary_indices[2 * nbi] + 1;
        ++nbi;
      }
      if (j == a_ny - 1) {
        boundary_attrib[nbi] = 4;
        boundary_indices[2 * nbi] = i + j * (a_nx + 1) + (a_nx + 1);
        boundary_indices[2 * nbi + 1] = boundary_indices[2 * nbi] + 1;
        ++nbi;
      }
    }
  }

  const double dx =
      (a_bounding_box[1][0] - a_bounding_box[0][0]) / static_cast<double>(a_nx);
  const double dy =
      (a_bounding_box[1][1] - a_bounding_box[0][1]) / static_cast<double>(a_ny);
  // Apparently a vertex is always 3D in MFEM but they ignore the third
  // dimension of each vertex in 2D
  for (int j = 0; j < a_ny + 1; ++j) {
    for (int i = 0; i < a_nx + 1; ++i) {
      const int node_index = j * (a_nx + 1) + i;
      const double x_loc = a_bounding_box[0][0] + static_cast<double>(i) * dx;
      const double y_loc = a_bounding_box[0][1] + static_cast<double>(j) * dy;
      vertices[3 * node_index] = x_loc;
      vertices[3 * node_index + 1] = y_loc;
      vertices[3 * node_index + 2] = 0.0;
    }
  }
  return std::make_pair(
      new mfem::Mesh(vertices, nvert, elem_indices.data(),
                     mfem::Geometry::Type::SQUARE, elem_attrib.data(), nelem,
                     boundary_indices.data(), mfem::Geometry::Type::SEGMENT,
                     boundary_attrib.data(), nboundary, 2, -1),
      vertices);
}

std::pair<mfem::Mesh*, double*> GenerateHexMesh(
    const std::array<std::array<double, 3>, 2>& a_bounding_box, const int a_nx,
    const int a_ny, const int a_nz) {
  const int nvert = (a_nx + 1) * (a_ny + 1) * (a_nz + 1);
  const int nelem = a_nx * a_ny * a_nz;
  double* vertices = new double[nvert * 3];

  std::vector<int> elem_indices(8 * nelem);

  std::vector<int> elem_attrib(nelem);
  std::fill(elem_attrib.begin(), elem_attrib.end(), 1);

  const int nboundary =
      2 * (a_ny * a_nz) + 2 * (a_nx * a_nz) + 2 * (a_nx * a_ny);
  std::vector<int> boundary_indices(4 * nboundary);
  std::vector<int> boundary_attrib(nboundary);
  int nbi = 0;
  const int x_jump = 1;
  const int y_jump = a_nx + 1;
  const int z_jump = (a_ny + 1) * (a_nx + 1);
  for (int k = 0; k < a_nz; ++k) {
    for (int j = 0; j < a_ny; ++j) {
      for (int i = 0; i < a_nx; ++i) {
        const int elem_index = k * a_ny * a_nx + j * a_nx + i;
        elem_indices[8 * elem_index] =
            k * (a_ny + 1) * (a_nx + 1) + j * (a_nx + 1) + i;
        elem_indices[8 * elem_index + 1] =
            elem_indices[8 * elem_index] + y_jump;
        elem_indices[8 * elem_index + 2] =
            elem_indices[8 * elem_index + 1] + z_jump;
        elem_indices[8 * elem_index + 3] =
            elem_indices[8 * elem_index + 2] - y_jump;
        elem_indices[8 * elem_index + 4] =
            elem_indices[8 * elem_index] + x_jump;
        elem_indices[8 * elem_index + 5] =
            elem_indices[8 * elem_index + 1] + x_jump;
        elem_indices[8 * elem_index + 6] =
            elem_indices[8 * elem_index + 2] + x_jump;
        elem_indices[8 * elem_index + 7] =
            elem_indices[8 * elem_index + 3] + x_jump;
        if (i == 0) {
          boundary_attrib[nbi] = 1;
          boundary_indices[4 * nbi] =
              k * (a_ny + 1) * (a_nx + 1) + j * (a_nx + 1) + i;
          boundary_indices[4 * nbi + 1] = boundary_indices[4 * nbi] + y_jump;
          boundary_indices[4 * nbi + 2] =
              boundary_indices[4 * nbi + 1] + z_jump;
          boundary_indices[4 * nbi + 3] =
              boundary_indices[4 * nbi + 2] - y_jump;
          ++nbi;
        }
        if (i == a_nx - 1) {
          boundary_attrib[nbi] = 2;
          boundary_indices[4 * nbi] =
              k * (a_ny + 1) * (a_nx + 1) + j * (a_nx + 1) + i;
          boundary_indices[4 * nbi + 1] = boundary_indices[4 * nbi] + z_jump;
          boundary_indices[4 * nbi + 2] =
              boundary_indices[4 * nbi + 1] + y_jump;
          boundary_indices[4 * nbi + 3] =
              boundary_indices[4 * nbi + 2] - z_jump;
          ++nbi;
        }
        if (j == 0) {
          boundary_attrib[nbi] = 3;
          boundary_indices[4 * nbi] =
              k * (a_ny + 1) * (a_nx + 1) + j * (a_nx + 1) + i;
          boundary_indices[4 * nbi + 1] = boundary_indices[4 * nbi] + z_jump;
          boundary_indices[4 * nbi + 2] =
              boundary_indices[4 * nbi + 1] + x_jump;
          boundary_indices[4 * nbi + 3] =
              boundary_indices[4 * nbi + 2] - z_jump;
          ++nbi;
        }
        if (j == a_ny - 1) {
          boundary_attrib[nbi] = 4;
          boundary_indices[4 * nbi] =
              k * (a_ny + 1) * (a_nx + 1) + j * (a_nx + 1) + i;
          boundary_indices[4 * nbi + 1] = boundary_indices[4 * nbi] + x_jump;
          boundary_indices[4 * nbi + 2] =
              boundary_indices[4 * nbi + 1] + z_jump;
          boundary_indices[4 * nbi + 3] =
              boundary_indices[4 * nbi + 2] - x_jump;
          ++nbi;
        }
        if (k == 0) {
          boundary_attrib[nbi] = 5;
          boundary_indices[4 * nbi] =
              k * (a_ny + 1) * (a_nx + 1) + j * (a_nx + 1) + i;
          boundary_indices[4 * nbi + 1] = boundary_indices[4 * nbi] + x_jump;
          boundary_indices[4 * nbi + 2] =
              boundary_indices[4 * nbi + 1] + y_jump;
          boundary_indices[4 * nbi + 3] =
              boundary_indices[4 * nbi + 2] - x_jump;
          ++nbi;
        }
        if (k == a_nz - 1) {
          boundary_attrib[nbi] = 6;
          boundary_indices[4 * nbi] =
              k * (a_ny + 1) * (a_nx + 1) + j * (a_nx + 1) + i;
          boundary_indices[4 * nbi + 1] = boundary_indices[4 * nbi] + y_jump;
          boundary_indices[4 * nbi + 2] =
              boundary_indices[4 * nbi + 1] + x_jump;
          boundary_indices[4 * nbi + 3] =
              boundary_indices[4 * nbi + 2] - y_jump;
          ++nbi;
        }
      }
    }
  }

  const double dx =
      (a_bounding_box[1][0] - a_bounding_box[0][0]) / static_cast<double>(a_nx);
  const double dy =
      (a_bounding_box[1][1] - a_bounding_box[0][1]) / static_cast<double>(a_ny);
  const double dz =
      (a_bounding_box[1][2] - a_bounding_box[0][2]) / static_cast<double>(a_nz);
  for (int k = 0; k < a_nz + 1; ++k) {
    for (int j = 0; j < a_ny + 1; ++j) {
      for (int i = 0; i < a_nx + 1; ++i) {
        const int node_index = k * (a_ny + 1) * (a_nx + 1) + j * (a_nx + 1) + i;
        const double x_loc = a_bounding_box[0][0] + static_cast<double>(i) * dx;
        const double y_loc = a_bounding_box[0][1] + static_cast<double>(j) * dy;
        const double z_loc = a_bounding_box[0][2] + static_cast<double>(k) * dz;
        vertices[3 * node_index] = x_loc;
        vertices[3 * node_index + 1] = y_loc;
        vertices[3 * node_index + 2] = z_loc;
      }
    }
  }
  return std::make_pair(
      new mfem::Mesh(vertices, nvert, elem_indices.data(),
                     mfem::Geometry::Type::CUBE, elem_attrib.data(), nelem,
                     boundary_indices.data(), mfem::Geometry::Type::SQUARE,
                     boundary_attrib.data(), nboundary, 3, -1),
      vertices);
}
