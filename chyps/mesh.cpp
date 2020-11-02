// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chyps/mesh.hpp"

#include <algorithm>
#include <set>

#include "chyps/io.hpp"
#include "chyps/logger.hpp"

namespace chyps {

Mesh::Mesh(const MPIParallel& a_mpi_session, InputParser& a_parser,
           IO* a_file_io)
    : parser_m(a_parser),
      mpi_session_m(a_mpi_session),
      file_io_m(a_file_io),
      parallel_mesh_m(nullptr),
      boundary_condition_manager_m(nullptr),
      element_offset_m(0) {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Constructing mesh object");
  this->GatherOptions();
}

void Mesh::Initialize(BoundaryConditionManager& a_boundary_condition_manager) {
  if (!this->AllOptionsSupplied()) {
    std::cout << "Not all options needed in parser are supplied" << std::endl;
    std::cout << "Make sure that the InputParser has been parsed before "
                 "calling Initialize and that all required options are "
                 "specified or have a valid default value."
              << std::endl;
    std::exit(-1);
  }
  boundary_condition_manager_m = &a_boundary_condition_manager;
  this->ReadAndRefineMesh();
  this->AllocateVariables();
}

const MPI_Comm& Mesh::GetMPIComm(void) const { return mpi_session_m.GetComm(); }

mfem::ParMesh& Mesh::GetMfemMesh(void) {
  assert(parallel_mesh_m != nullptr);
  return *parallel_mesh_m;
}

int Mesh::GetDimension(void) const {
  assert(parallel_mesh_m != nullptr);
  return parallel_mesh_m->Dimension();
}

template <>
std::size_t Mesh::GetGlobalCount<MeshElement::ELEMENT>(void) const {
  assert(parallel_mesh_m != nullptr);
  return static_cast<std::size_t>(parallel_mesh_m->GetGlobalNE());
}

template <>
std::size_t Mesh::GetOffsetStart<MeshElement::ELEMENT>(void) const {
  assert(parallel_mesh_m != nullptr);
  return element_offset_m;
}

template <>
std::size_t Mesh::GetLocalCount<MeshElement::ELEMENT>(void) const {
  assert(parallel_mesh_m != nullptr);
  return static_cast<std::size_t>(parallel_mesh_m->GetNE());
}

template <>
std::size_t Mesh::GetLocalCount<MeshElement::VERTEX>(void) const {
  assert(parallel_mesh_m != nullptr);
  return static_cast<std::size_t>(parallel_mesh_m->GetNV());
}

int Mesh::GetNumberOfBoundaryTagValues(void) const {
  assert(parallel_mesh_m != nullptr);
  return parallel_mesh_m->bdr_attributes.Max();
}

const BoundaryConditionManager& Mesh::GetBoundaryConditionManager(void) const {
  assert(boundary_condition_manager_m != nullptr);
  return *boundary_condition_manager_m;
}

std::pair<std::vector<double>, std::vector<int>> Mesh::GetBoundaryVertices(
    const int a_tag) const {
  SPDLOG_LOGGER_INFO(
      MAIN_LOG,
      "Building boundary vertex position list for boundaries tagged with {}",
      a_tag);
  const int dimension = this->GetDimension();

  mfem::Array<int> element_vertices;
  std::set<int> unique_border_vertices;
  const int number_of_border_elements = parallel_mesh_m->GetNBE();
  for (int n = 0; n < number_of_border_elements; ++n) {
    const mfem::Element* element = parallel_mesh_m->GetBdrElement(n);
    if (element->GetAttribute() == a_tag) {
      element->GetVertices(element_vertices);
      for (int v = 0; v < element_vertices.Size(); ++v) {
        unique_border_vertices.insert(element_vertices[v]);
      }
    }
  }

  const std::size_t number_of_vertices = unique_border_vertices.size();
  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "{} unique boundary vertices found in mesh with tag {}",
                     number_of_vertices, a_tag);

  std::vector<double> boundary_vertices(number_of_vertices *
                                        static_cast<std::size_t>(dimension));
  std::vector<int> boundary_vertex_indices(number_of_vertices);
  std::size_t nvertex = 0;
  for (const auto& elem : unique_border_vertices) {
    boundary_vertex_indices[nvertex] = elem;
    const double* vertex = parallel_mesh_m->GetVertex(elem);
    for (int d = 0; d < dimension; ++d) {
      boundary_vertices[dimension * nvertex + d] = vertex[d];
    }
    ++nvertex;
  }
  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Built boundary vertex position list. Has length {} "
                     "consisting of {} vertices.",
                     boundary_vertices.size(),
                     boundary_vertices.size() / dimension);

  return std::make_pair(boundary_vertices, boundary_vertex_indices);
}

Mesh::~Mesh(void) {
  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Destructing Mesh and freeing associated memory");
  delete parallel_mesh_m;
  parallel_mesh_m = nullptr;
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Mesh successfully destructed.");
}

bool Mesh::AllOptionsSupplied(void) { return parser_m.AllOptionsSet(); }
void Mesh::GatherOptions(void) {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Adding options to look for in parser");
  parser_m.AddOption("mesh_file", "-m", "--mesh", "Mesh file to use.",
                     std::string("../data/star.mesh"));
  parser_m.AddOption("serial_refine", "-rs", "--refine-serial",
                     "Number of times to refine the mesh uniformly in serial.",
                     2);
  parser_m.AddOption(
      "parallel_refine", "-rp", "--refine-parallel",
      "Number of times to refine the mesh uniformly in parallel.", 1);
  parser_m.AddOption(
      "gen_nx", "-nx", "--nx",
      "If using generated mesh, number of elements in x direction.", -1);
  parser_m.AddOption(
      "gen_ny", "-ny", "--ny",
      "If using generated mesh, number of elements in y direction.", -1);
  parser_m.AddOption(
      "gen_nz", "-nz", "--nz",
      "If using generated mesh, number of elements in z direction.", -1);
  parser_m.AddOption(
      "gen_blx", "-blx", "--blx",
      "If using generated mesh, lower x dimension of the cuboid domain", -1.0);
  parser_m.AddOption(
      "gen_bly", "-bly", "--bly",
      "If using generated mesh, lower y dimension of the cuboid domain", -1.0);
  parser_m.AddOption(
      "gen_blz", "-blz", "--blz",
      "If using generated mesh, lower z dimension of the cuboid domain", -1.0);
  parser_m.AddOption(
      "gen_bux", "-bux", "--bux",
      "If using generated mesh, upper x dimension of the cuboid domain", 1.0);
  parser_m.AddOption(
      "gen_buy", "-buy", "--buy",
      "If using generated mesh, upper y dimension of the cuboid domain", 1.0);
  parser_m.AddOption(
      "gen_buz", "-buz", "--buz",
      "If using generated mesh, upper z dimension of the cuboid domain", 1.0);
  SPDLOG_LOGGER_INFO(MAIN_LOG, "All options added to parser for Mesh class");
}

void Mesh::ReadAndRefineMesh(void) {
  mfem::Mesh* serial_mesh = nullptr;
  double* vertices = nullptr;
  const auto file_name = parser_m["mesh_file"].get<std::string>();

  SPDLOG_LOGGER_INFO(MAIN_LOG, "Mesh reading mesh named {}", file_name);

  if (file_name == "generate") {
    const auto nx = parser_m["gen_nx"].get<int>();
    const auto ny = parser_m["gen_ny"].get<int>();
    const auto nz = parser_m["gen_nz"].get<int>();
    if (ny <= 0) {
      assert(nx > 0);
      assert(nz <= 0);
      std::array<std::array<double, 1>, 2> bounding_box{
          {parser_m["gen_blx"].get<double>(),
           parser_m["gen_bly"].get<double>()}};
      auto mesh_and_vertices = this->GenerateLineMesh(bounding_box, nx);
      serial_mesh = mesh_and_vertices.first;
      vertices = mesh_and_vertices.second;
    } else if (nz <= 0) {
      assert(nx > 0);
      assert(ny > 0);
      // 2D Mesh with Quads
      std::array<std::array<double, 2>, 2> bounding_box{
          {{parser_m["gen_blx"].get<double>(),
            parser_m["gen_bly"].get<double>()},
           {parser_m["gen_bux"].get<double>(),
            parser_m["gen_buy"].get<double>()}}};
      auto mesh_and_vertices = this->GenerateQuadMesh(bounding_box, nx, ny);
      serial_mesh = mesh_and_vertices.first;
      vertices = mesh_and_vertices.second;
    } else {
      assert(nx > 0);
      assert(ny > 0);
      assert(nz > 0);
      // 3D Mesh with Hexs
      std::array<std::array<double, 3>, 2> bounding_box{
          {{parser_m["gen_blx"].get<double>(),
            parser_m["gen_bly"].get<double>(),
            parser_m["gen_blz"].get<double>()},
           {parser_m["gen_bux"].get<double>(),
            parser_m["gen_buy"].get<double>(),
            parser_m["gen_buz"].get<double>()}}};
      auto mesh_and_vertices = this->GenerateHexMesh(bounding_box, nx, ny, nz);
      serial_mesh = mesh_and_vertices.first;
      vertices = mesh_and_vertices.second;
    }
  } else {
    SPDLOG_LOGGER_INFO(MAIN_LOG, "Reading mesh from file.");
    serial_mesh = new mfem::Mesh(file_name.c_str(), 1, 1);
    SPDLOG_LOGGER_INFO(MAIN_LOG, "Mesh successfully read from file");
  }
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Mesh consists of {} dimensions",
                     serial_mesh->Dimension());

  const auto serial_refinement = parser_m["serial_refine"].get<int>();
  const auto parallel_refinement = parser_m["parallel_refine"].get<int>();

  SPDLOG_LOGGER_INFO(MAIN_LOG, "Performing {} levels of serial refinement",
                     serial_refinement);
  for (int lev = 0; lev < serial_refinement; ++lev) {
    serial_mesh->UniformRefinement();
  }

  // FIXME:  1D mesh fails in conversion to parallel_mesh_m
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Converting serial mesh to parallel mesh.");
  parallel_mesh_m = new mfem::ParMesh(this->GetMPIComm(), *serial_mesh);
  delete serial_mesh;
  delete[] vertices;
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Parallel mesh successfully built");

  SPDLOG_LOGGER_INFO(MAIN_LOG, "Performing {} levels of parallel refinement",
                     parallel_refinement);
  for (int lev = 0; lev < parallel_refinement; ++lev) {
    parallel_mesh_m->UniformRefinement();
  }

  SPDLOG_LOGGER_INFO(MAIN_LOG, "Final parallel mesh successfully created");
}

void Mesh::AllocateVariables(void) {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Allocating space for {} boundary conditions",
                     parallel_mesh_m->bdr_attributes.Max());

  SPDLOG_LOGGER_INFO(MAIN_LOG, "Computing offsets for mesh elements.");
  std::size_t local_elems = this->GetLocalCount<MeshElement::ELEMENT>();
  MPI_Scan(&local_elems, &element_offset_m, 1, MPI_LONG, MPI_SUM,
           this->GetMPIComm());
  element_offset_m -= local_elems;
}

std::pair<mfem::Mesh*, double*> Mesh::GenerateLineMesh(
    const std::array<std::array<double, 1>, 2>& a_bounding_box,
    const int a_nx) {
  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Creating 1D Line Mesh from {} to {} with {} elements.",
                     a_bounding_box[0][0], a_bounding_box[1][0], a_nx);
  const int nvert = a_nx + 1;
  const int nelem = a_nx;
  double* vertices = new double[nvert * 3];

  std::vector<int> elem_indices(2 * nelem);

  std::vector<int> elem_attrib(nelem);
  std::fill(elem_attrib.begin(), elem_attrib.end(), 1);

  const int nboundary = 2;
  std::vector<int> boundary_indices(nboundary);
  std::vector<int> boundary_attrib(nboundary);
  for (int i = 0; i < a_nx; ++i) {
    elem_indices[2 * i] = i;
    elem_indices[2 * i + 1] = i + 1;
  }
  boundary_indices[0] = 0;
  boundary_attrib[0] = 1;
  boundary_indices[1] = a_nx;
  boundary_attrib[1] = 2;

  const double dx =
      (a_bounding_box[1][0] - a_bounding_box[0][0]) / static_cast<double>(a_nx);
  // Apparently a vertex is always 3D in MFEM
  for (int i = 0; i < a_nx + 1; ++i) {
    const int node_index = i;
    const double x_loc = a_bounding_box[0][0] + static_cast<double>(i) * dx;
    vertices[3 * node_index] = x_loc;
    vertices[3 * node_index + 1] = 0.0;
    vertices[3 * node_index + 2] = 0.0;
  }
  return std::make_pair(
      new mfem::Mesh(vertices, nvert, elem_indices.data(),
                     mfem::Geometry::Type::SEGMENT, elem_attrib.data(), nelem,
                     boundary_indices.data(), mfem::Geometry::Type::POINT,
                     boundary_attrib.data(), nboundary, 1, -1),
      vertices);
}

std::pair<mfem::Mesh*, double*> Mesh::GenerateQuadMesh(
    const std::array<std::array<double, 2>, 2>& a_bounding_box, const int a_nx,
    const int a_ny) {
  SPDLOG_LOGGER_INFO(
      MAIN_LOG,
      "Creating 2D Quad Mesh from ({},{}) to ({},{}) with ({},{}) elements.",
      a_bounding_box[0][0], a_bounding_box[0][1], a_bounding_box[1][0],
      a_bounding_box[1][1], a_nx, a_ny);

  const int nvert = (a_nx + 1) * (a_ny + 1);
  const int nelem = a_nx * a_ny;
  double* vertices = new double[nvert * 3];

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

std::pair<mfem::Mesh*, double*> Mesh::GenerateHexMesh(
    const std::array<std::array<double, 3>, 2>& a_bounding_box, const int a_nx,
    const int a_ny, const int a_nz) {
  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Creating 3D Hex Mesh from ({},{},{}) to ({},{},{}) with "
                     "({},{},{}) elements.",
                     a_bounding_box[0][0], a_bounding_box[0][1],
                     a_bounding_box[0][2], a_bounding_box[1][0],
                     a_bounding_box[1][1], a_bounding_box[1][2], a_nx, a_ny,
                     a_nz);

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

bool Mesh::FileWritingEnabled(void) const {
  return file_io_m != nullptr && file_io_m->IsWriteModeActive();
}

// Write now we will handle IO with independent (rank local) writes.
// In future, would like to figure out how to handle global writes into
// one-indexed array/space. Only hold-up is how to handle vertices (get a global
// vertex number, which I don't believe MFEM stores).
void Mesh::AddIOVariables(void) {
  if (!this->FileWritingEnabled()) {
    return;
  }
  file_io_m->WriteAttribute("format", std::string("MFEM ADIOS2 BP v0.2"));
  file_io_m->WriteAttribute("format/version", std::string("0.2"));
  file_io_m->WriteAttribute("format/mfem_mesh", std::string("MFEM mesh v1.0"));
  std::vector<std::string> export_types{"Paraview: ADIOS2VTXReader",
                                        "VTK: vtkADIOS2VTXReader.h"};
  file_io_m->WriteAttribute("format/viz_tools", export_types);

  file_io_m->RootAddVariable<uint32_t>("dimension");
  file_io_m->AddVariable<uint32_t>("NumOfElements");
  file_io_m->AddVariable<uint32_t>("NumOfVertices");
  file_io_m->RootAddVariable<uint32_t>("types");
  // Need to compute and include sizes for these.
  // Assumes same element type across all of mesh.
  std::size_t element_nvertices =
      static_cast<std::size_t>(parallel_mesh_m->GetElement(0)->GetNVertices());
  file_io_m->AddVariable<uint64_t>(
      "connectivity",
      {this->GetLocalCount<MeshElement::ELEMENT>(), 1 + element_nvertices},
      true);
  file_io_m->AddVariable<int32_t>(
      "attribute", {this->GetLocalCount<MeshElement::ELEMENT>()}, true);
  file_io_m->MarkAsElementVariable("attribute");

  // Will handle vertices in WriteMesh function.
}

void Mesh::WriteMesh(void) {
  if (!this->FileWritingEnabled()) {
    return;
  }
  assert(file_io_m->IsWriteModeActive());
  assert(file_io_m->OngoingWriteStep());
  auto min_max_order = file_io_m->MinMaxVariableOrder();
  // Handle only constant order mesh/simulation for now
  assert(min_max_order[0] == min_max_order[1]);
  mfem::H1_FECollection fec(min_max_order[0], this->GetDimension());
  mfem::ParFiniteElementSpace fes(parallel_mesh_m, &fec, this->GetDimension(),
                                  mfem::Ordering::byVDIM);
  mfem::ParGridFunction nodes(&fes);
  parallel_mesh_m->GetNodes(nodes);

  //*****************
  // This part inherently assumes mesh will only be written once.
  this->AddIOVariables();
  // Handle vertices
  const auto ndofs = static_cast<std::size_t>(fes.GetNDofs());
  const auto dim = static_cast<std::size_t>(this->GetDimension());
  // True here will be invalid for AMR and H/P refinement.
  file_io_m->AddVariable<double>("vertices", {ndofs, dim}, true);
  //*****************

  uint32_t dimension = static_cast<uint32_t>(this->GetDimension());
  uint32_t number_of_elements =
      static_cast<uint32_t>(this->GetLocalCount<MeshElement::ELEMENT>());
  uint32_t number_of_vertices = static_cast<uint32_t>(fes.GetNDofs());
  // static_cast<uint32_t>(this->GetLocalCount<MeshElement::VERTEX>());
  uint32_t vtk_type = static_cast<uint32_t>(this->GLVISToVTKType(
      static_cast<int>(parallel_mesh_m->GetElement(0)->GetGeometryType())));
  file_io_m->RootPutDeferred("dimension", &dimension);
  file_io_m->PutDeferred("NumOfElements", &number_of_elements);
  file_io_m->PutDeferred("NumOfVertices", &number_of_vertices);
  file_io_m->RootPutDeferred("types", &vtk_type);
  auto connectivity_span = file_io_m->PutSpan<uint64_t>("connectivity");
  auto attribute_span = file_io_m->PutSpan<int32_t>("attribute");

  std::size_t span_position = 0;
  for (uint32_t n = 0; n < number_of_elements; ++n) {
    attribute_span[n] = static_cast<uint32_t>(parallel_mesh_m->GetAttribute(n));
    auto element = parallel_mesh_m->GetElement(n);
    const auto nvertices = element->GetNVertices();
    const int* vertex_buffer = element->GetVertices();
    connectivity_span[span_position] = nvertices;
    std::copy(vertex_buffer, vertex_buffer + nvertices,
              &connectivity_span[span_position + 1]);
    span_position += (1 + nvertices);
  }
  assert(span_position ==
         static_cast<std::size_t>(
             (parallel_mesh_m->GetElement(0)->GetNVertices() + 1) *
             number_of_elements));

  file_io_m->PutDeferred("vertices", nodes.GetData());

  file_io_m->PerformPuts();
}

// Note: Taken from adios2stream.cpp in MFEM, added
// by William F Godoy godoywf@ornl.gov on Jan 22, 2020
uint32_t Mesh::GLVISToVTKType(const int glvisType) const noexcept {
  uint32_t vtkType = 0;
  switch (glvisType) {
    case mfem::Geometry::Type::POINT:
      vtkType = 1;
      break;
    case mfem::Geometry::Type::SEGMENT:
      vtkType = 3;
      break;
    case mfem::Geometry::Type::TRIANGLE:
      vtkType = 5;
      break;
    case mfem::Geometry::Type::SQUARE:
      // vtkType = 8;
      vtkType = 9;
      break;
    case mfem::Geometry::Type::TETRAHEDRON:
      vtkType = 10;
      break;
    case mfem::Geometry::Type::CUBE:
      // vtkType = 11;
      vtkType = 12;
      break;
    case mfem::Geometry::Type::PRISM:
      vtkType = 13;
      break;
    default:
      vtkType = 0;
      break;
  }
  return vtkType;
}

}  // namespace chyps
