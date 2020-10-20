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

#include "chyps/logger.hpp"

namespace chyps {

namespace details {
// H is homogeneous, TV is time-varying
enum BCType {
  H_DIRICHLET = 0,
  DIRICHLET,
  TV_DIRICHLET,
  H_NEUMANN,
  NEUMANN,
  TV_NEUMANN,
  COUNT
};
}  // namespace details

Mesh::Mesh(const MPI_Comm& a_mpi_comm, InputParser& a_parser)
    : parser_m(a_parser),
      mpi_comm_m(a_mpi_comm),
      parallel_mesh_m(nullptr),
      boundary_conditions_m() {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Constructing mesh object");
  this->GatherOptions();
}

void Mesh::Initialize(void) {
  if (!this->AllOptionsSupplied()) {
    std::cout << "Not all options needed in parser are supplied" << std::endl;
    std::cout << "Make sure that the InputParser has been parsed before "
                 "calling Initialize and that all required options are "
                 "specified or have a valid default value."
              << std::endl;
    std::exit(-1);
  }
  this->ReadAndRefineMesh();
  this->AllocateVariables();
}

const MPI_Comm& Mesh::GetMPIComm(void) const { return mpi_comm_m; }

mfem::ParMesh& Mesh::GetMfemMesh(void) {
  assert(parallel_mesh_m != nullptr);
  return *parallel_mesh_m;
}

int Mesh::GetDimension(void) const {
  assert(parallel_mesh_m != nullptr);
  return parallel_mesh_m->Dimension();
}

int Mesh::GetNumberOfBoundaryConditions(void) const {
  return static_cast<int>(boundary_conditions_m.size());
}

const BoundaryCondition& Mesh::GetBoundaryCondition(const int a_tag) const {
  assert(a_tag >= 1);
  assert(a_tag - 1 < this->GetNumberOfBoundaryConditions());
  return boundary_conditions_m[a_tag - 1];
}

int Mesh::GetNumberOfHomogeneousDirichletConditions(void) const {
  return boundary_condition_counts_m[details::BCType::H_DIRICHLET];
}

int Mesh::GetNumberOfDirichletConditions(void) const {
  return boundary_condition_counts_m[details::BCType::DIRICHLET];
}

int Mesh::GetNumberOfTimeVaryingDirichletConditions(void) const {
  assert(boundary_condition_counts_m[details::BCType::TV_DIRICHLET] <=
         boundary_condition_counts_m[details::BCType::DIRICHLET]);
  return boundary_condition_counts_m[details::BCType::TV_DIRICHLET];
}

int Mesh::GetNumberOfHomogeneousNeumannConditions(void) const {
  return boundary_condition_counts_m[details::BCType::H_NEUMANN];
}

int Mesh::GetNumberOfNeumannConditions(void) const {
  return boundary_condition_counts_m[details::BCType::NEUMANN];
}

int Mesh::GetNumberOfTimeVaryingNeumannConditions(void) const {
  assert(boundary_condition_counts_m[details::BCType::TV_NEUMANN] <=
         boundary_condition_counts_m[details::BCType::NEUMANN]);
  return boundary_condition_counts_m[details::BCType::TV_NEUMANN];
}

void Mesh::SetBoundaryCondition(const int a_tag,
                                const BoundaryCondition& a_condition) {
  assert(a_tag >= 1 && a_tag <= static_cast<int>(boundary_conditions_m.size()));
  boundary_conditions_m[a_tag - 1] = a_condition;
  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Added BoundaryCondition of type {} for boundary tag {}",
                     static_cast<int>(a_condition.GetBCType()), a_tag);
}

void Mesh::CommitBoundaryConditions(void) {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Commiting {} boundary conditions",
                     boundary_conditions_m.size());
  boundary_condition_counts_m.resize(details::BCType::COUNT);
  std::fill(boundary_condition_counts_m.begin(),
            boundary_condition_counts_m.end(), 0);
  for (const auto& condition : boundary_conditions_m) {
    if (condition.GetBCType() == BoundaryConditionType::HOMOGENEOUS_DIRICHLET) {
      ++boundary_condition_counts_m[details::BCType::H_DIRICHLET];
    } else if (condition.GetBCType() == BoundaryConditionType::DIRICHLET) {
      ++boundary_condition_counts_m[details::BCType::DIRICHLET];
      if (condition.IsTimeVarying()) {
        ++boundary_condition_counts_m[details::BCType::TV_DIRICHLET];
      }
    } else if (condition.GetBCType() ==
               BoundaryConditionType::HOMOGENEOUS_NEUMANN) {
      ++boundary_condition_counts_m[details::BCType::H_NEUMANN];
    } else if (condition.GetBCType() == BoundaryConditionType::NEUMANN) {
      ++boundary_condition_counts_m[details::BCType::NEUMANN];
      if (condition.IsTimeVarying()) {
        ++boundary_condition_counts_m[details::BCType::TV_NEUMANN];
      }

    } else {
      std::cout << "Unknown boundary condition type of : "
                << static_cast<int>(condition.GetBCType()) << std::endl;
      std::exit(-1);
    }
  }
}

std::pair<std::vector<double>, std::vector<int>> Mesh::GetBoundaryVertices(
    const int a_tag) {
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
  const std::string file_name = parser_m["mesh_file"];

  SPDLOG_LOGGER_INFO(MAIN_LOG, "Mesh reading mesh named {}", file_name);

  if (file_name == "generate") {
    const int nx = parser_m["gen_nx"];
    const int ny = parser_m["gen_ny"];
    const int nz = parser_m["gen_nz"];
    if (ny <= 0) {
      assert(nx > 0);
      assert(nz <= 0);
      std::array<std::array<double, 1>, 2> bounding_box{
          {parser_m["gen_blx"], parser_m["gen_bly"]}};
      auto mesh_and_vertices = this->GenerateLineMesh(bounding_box, nx);
      serial_mesh = mesh_and_vertices.first;
      vertices = mesh_and_vertices.second;
    } else if (nz <= 0) {
      assert(nx > 0);
      assert(ny > 0);
      // 2D Mesh with Quads
      std::array<std::array<double, 2>, 2> bounding_box{
          {{parser_m["gen_blx"], parser_m["gen_bly"]},
           {parser_m["gen_bux"], parser_m["gen_buy"]}}};
      auto mesh_and_vertices = this->GenerateQuadMesh(bounding_box, nx, ny);
      serial_mesh = mesh_and_vertices.first;
      vertices = mesh_and_vertices.second;
    } else {
      assert(nx > 0);
      assert(ny > 0);
      assert(nz > 0);
      // 3D Mesh with Hexs
      std::array<std::array<double, 3>, 2> bounding_box{
          {{parser_m["gen_blx"], parser_m["gen_bly"], parser_m["gen_blz"]},
           {parser_m["gen_bux"], parser_m["gen_buy"], parser_m["gen_buz"]}}};
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

  const auto serial_refinement = static_cast<int>(parser_m["serial_refine"]);
  const auto parallel_refinement =
      static_cast<int>(parser_m["parallel_refine"]);

  SPDLOG_LOGGER_INFO(MAIN_LOG, "Performing {} levels of serial refinement",
                     serial_refinement);
  for (int lev = 0; lev < serial_refinement; ++lev) {
    serial_mesh->UniformRefinement();
  }

  // FIXME:  1D mesh fails in conversion to parallel_mesh_m
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Converting serial mesh to parallel mesh.");
  parallel_mesh_m = new mfem::ParMesh(mpi_comm_m, *serial_mesh);
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
  boundary_conditions_m.resize(parallel_mesh_m->bdr_attributes.Max());
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
  boundary_indices[1] = a_nx + 1;
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

}  // namespace chyps
