// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chyps/io.hpp"

#include "chyps/debug_assert.hpp"
#include "chyps/logger.hpp"

namespace chyps {
IO::IO(const MPIParallel& a_mpi_session, const std::string a_io_name)
    : io_driver_m(a_mpi_session.GetComm()),
      mpi_session_m(a_mpi_session),
      read_m(io_driver_m.DeclareIO(a_io_name + "_READ")),
      write_m(io_driver_m.DeclareIO(a_io_name + "_WRITE")),
      read_engine_m(nullptr),
      write_engine_m(nullptr),
      element_data_variables_m(),
      point_data_variables_m(),
      point_data_orders_m(),
      active_write_step_m(false),
      active_read_step_m(false) {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Constructed new IO object");

  this->SetEngineType();
  this->SetDefaultParameters();
  this->RootAddVariable<std::string>("REAL_WORLD_TIME");
  this->RootAddVariable<uint64_t>("CYCLE");
  this->RootAddVariable<double>("TIME");
  this->RootAddVariable<double>("DT");
}

void IO::SetParameter(const std::string a_key, const std::string a_value) {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Setting parameter {} with value {}", a_key,
                     a_value);
  read_m.SetParameter(a_key, a_value);
  write_m.SetParameter(a_key, a_value);
}

void IO::SetWrite(std::string a_file_name) {
  if (write_engine_m != nullptr) {
    write_engine_m->Close();
    delete write_engine_m;
  }
  a_file_name.append(".bp");
  write_engine_m =
      new adios2::Engine(write_m.Open(a_file_name, adios2::Mode::Write));
}

void IO::SetRead(std::string a_file_name) {
  if (read_engine_m != nullptr) {
    read_engine_m->Close();
    delete read_engine_m;
  }
  a_file_name.append(".bp");
  read_engine_m =
      new adios2::Engine(read_m.Open(a_file_name, adios2::Mode::Read));
}

bool IO::IsWriteModeActive(void) const { return write_engine_m != nullptr; }

bool IO::IsReadModeActive(void) const { return read_engine_m != nullptr; }

void IO::MarkAsElementVariable(const std::string a_name) {
  auto result = std::find(element_data_variables_m.begin(),
                          element_data_variables_m.end(), a_name);
  if (result == element_data_variables_m.end()) {
    element_data_variables_m.push_back(a_name);
  }
}

void IO::MarkAsPointVariable(const std::string a_name, const int a_order) {
  auto result = std::find(point_data_variables_m.begin(),
                          point_data_variables_m.end(), a_name);
  if (result == point_data_variables_m.end()) {
    point_data_variables_m.push_back(a_name);
    point_data_orders_m.push_back(a_order);
  }
}

// Note: This routine is for all ranks writing their own section. Need to work
// on how to handle the all ranks to one section approach.
void IO::AddVariableForGridFunction(
    const std::string a_variable_name,
    const mfem::ParFiniteElementSpace& a_element_space,
    const bool a_dimensions_static) {
  const std::size_t ndofs =
      static_cast<std::size_t>(a_element_space.GetNDofs());
  std::string node_ordering;
  if (a_element_space.GetVDim() == 1) {
    write_m.DefineVariable<double>(a_variable_name, {}, {}, {ndofs},
                                   a_dimensions_static);
  } else {
    const std::size_t components =
        static_cast<std::size_t>(a_element_space.GetVDim());
    if (a_element_space.GetOrdering() == mfem::Ordering::byNODES) {
      node_ordering = "nodes";
      write_m.DefineVariable<double>(a_variable_name, {}, {},
                                     {components, ndofs}, a_dimensions_static);
    } else if (a_element_space.GetOrdering() == mfem::Ordering::byVDIM) {
      node_ordering = "vdim";
      write_m.DefineVariable<double>(a_variable_name, {}, {},
                                     {ndofs, components}, a_dimensions_static);
    } else {
      DEBUG_ASSERT(false, global_assert{}, DebugLevel::ALWAYS{},
                   "Unkown node ordering for Element Space when adding to IO");
    }
  }
  this->WriteAttributeForVariable(
      a_variable_name, "FiniteElementSpace",
      std::string(a_element_space.FEColl()->Name()));
  this->WriteAttributeForVariable(a_variable_name, "DataOrdering",
                                  node_ordering);
  // Assumes constant order elements below
  this->WriteAttributeForVariable(a_variable_name, "FiniteElementOrder",
                                  a_element_space.GetOrder(0));
  this->MarkAsPointVariable(a_variable_name, a_element_space.GetOrder(0));
}

// Note: This routine is for all ranks writing their own section. Need to work
// on how to handle the all ranks to one section approach.
void IO::AddVariableForTrueDofs(
    const std::string a_variable_name,
    const mfem::ParFiniteElementSpace& a_element_space,
    const bool a_dimensions_static) {
  const std::size_t components =
      static_cast<std::size_t>(a_element_space.GetVDim());
  const std::size_t true_dofs =
      static_cast<std::size_t>(a_element_space.GetTrueVSize()) / components;

  std::string node_ordering;
  if (components == 1) {
    write_m.DefineVariable<double>(a_variable_name, {}, {}, {true_dofs},
                                   a_dimensions_static);
  } else {
    if (a_element_space.GetOrdering() == mfem::Ordering::byNODES) {
      node_ordering = "nodes";
      write_m.DefineVariable<double>(a_variable_name, {}, {},
                                     {components, true_dofs},
                                     a_dimensions_static);
    } else if (a_element_space.GetOrdering() == mfem::Ordering::byVDIM) {
      node_ordering = "vdim";
      write_m.DefineVariable<double>(a_variable_name, {}, {},
                                     {true_dofs, components},
                                     a_dimensions_static);
    } else {
      DEBUG_ASSERT(false, global_assert{}, DebugLevel::ALWAYS{},
                   "Unkown node ordering for Element Space when adding to IO");
    }
  }
  this->WriteAttributeForVariable(
      a_variable_name, "FiniteElementSpace",
      std::string(a_element_space.FEColl()->Name()));
  this->WriteAttributeForVariable(a_variable_name, "DataOrdering",
                                  node_ordering);
  // Assumes constant order elements below
  this->WriteAttributeForVariable(a_variable_name, "FiniteElementOrder",
                                  a_element_space.GetOrder(0));
  this->MarkAsPointVariable(a_variable_name, a_element_space.GetOrder(0));
}

void IO::BeginWriteStep(const uint64_t a_cycle, const double a_time,
                        const double a_dt) {
  DEBUG_ASSERT(!this->OngoingWriteStep(), global_assert{}, DebugLevel::CHEAP{},
               "Already in middle of write step.");
  DEBUG_ASSERT(write_engine_m != nullptr, global_assert{}, DebugLevel::CHEAP{});
  active_write_step_m = true;
  write_engine_m->BeginStep();
  const auto now = std::chrono::system_clock::now();
  const auto time = std::chrono::system_clock::to_time_t(now);
  std::string real_world_time(std::ctime(&time));
  this->RootPutDeferred("REAL_WORLD_TIME", &real_world_time);
  this->RootPutDeferred("CYCLE", &a_cycle);
  this->RootPutDeferred("TIME", &a_time);
  this->RootPutDeferred("DT", &a_dt);
  this->PerformPuts();
}

bool IO::OngoingWriteStep(void) const {
  DEBUG_ASSERT(write_engine_m != nullptr, global_assert{}, DebugLevel::CHEAP{});
  return active_write_step_m;
}

void IO::EndWriteStep(void) {
  DEBUG_ASSERT(
      this->OngoingWriteStep(), global_assert{}, DebugLevel::CHEAP{},
      "No write step has been started. Has BeginWriteStep been called?");
  DEBUG_ASSERT(write_engine_m != nullptr, global_assert{}, DebugLevel::CHEAP{});
  active_write_step_m = false;
  write_engine_m->EndStep();
}

void IO::CloseWriteEngine(void) {
  if (write_engine_m != nullptr) {
    write_engine_m->Close();
    delete write_engine_m;
    write_engine_m = nullptr;
    active_write_step_m = false;
  }
}

void IO::BeginReadStep(void) {
  DEBUG_ASSERT(!this->OngoingReadStep(), global_assert{}, DebugLevel::CHEAP{},
               "Already in middle of read step.");
  DEBUG_ASSERT(read_engine_m != nullptr, global_assert{}, DebugLevel::CHEAP{});
  active_read_step_m = true;
  read_engine_m->BeginStep();
}

bool IO::OngoingReadStep(void) const {
  DEBUG_ASSERT(read_engine_m != nullptr, global_assert{}, DebugLevel::CHEAP{});
  return active_read_step_m;
}

void IO::EndReadStep(void) {
  DEBUG_ASSERT(this->OngoingReadStep(), global_assert{}, DebugLevel::CHEAP{},
               "No read step has been started. Has BeginReadStep been called?");
  DEBUG_ASSERT(read_engine_m != nullptr, global_assert{}, DebugLevel::CHEAP{});
  active_read_step_m = false;
  read_engine_m->EndStep();
}

void IO::CloseReadEngine(void) {
  if (read_engine_m != nullptr) {
    read_engine_m->Close();
    delete read_engine_m;
    read_engine_m = nullptr;
    active_read_step_m = false;
  }
}

void IO::PutDeferred(const std::string& a_variable_name,
                     const mfem::Vector& a_vector) {
  DEBUG_ASSERT(this->IsWriteModeActive(), global_assert{}, DebugLevel::CHEAP{},
               "IO must be in write mode for writing to fields.");
  DEBUG_ASSERT(this->OngoingWriteStep(), global_assert{}, DebugLevel::CHEAP{},
               "An ongoing IO step is required for writing.");
  this->Put(a_variable_name, a_vector.GetData(), adios2::Mode::Deferred);
}

void IO::PerformPuts(void) {
  DEBUG_ASSERT(this->IsWriteModeActive(), global_assert{}, DebugLevel::CHEAP{},
               "IO must be in write mode for writing to fields.");
  DEBUG_ASSERT(this->OngoingWriteStep(), global_assert{}, DebugLevel::CHEAP{},
               "An ongoing IO step is required for writing.");
  write_engine_m->PerformPuts();
}

void IO::WriteXMLSchema(void) {
  this->WriteAttribute("vtk.xml", this->VTKSchema());
}

IO::~IO(void) {
  // All objects from ADIOS are freed in ADIOS besides io_driver_m.
  this->CloseWriteEngine();
  this->CloseReadEngine();
}

void IO::SetEngineType(void) {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Setting engine in ADIOS2");
  read_m.SetEngine("BP4");
  write_m.SetEngine("BP4");
}

void IO::SetDefaultParameters(void) {
  SPDLOG_LOGGER_INFO(MAIN_LOG, "Accepting all default parameters");
  // Do nothing for now, accepting all defaults for engine from ADIOS2
}

ExistenceLocation IO::FindExistence(const bool a_read, const bool a_write) {
  int val = 0;
  val += a_read ? 1 : 0;
  val += a_write ? 2 : 0;
  return static_cast<ExistenceLocation>(val);
}

std::array<std::size_t, 3> IO::GetMeshSizes(const Mesh& a_mesh,
                                            const MeshElement a_type) const {
  std::array<std::size_t, 3> mesh_sizes{{0, 0, 0}};

  switch (a_type) {
    case MeshElement::ELEMENT: {
      mesh_sizes[0] = a_mesh.GetGlobalCount<MeshElement::ELEMENT>();
      mesh_sizes[1] = a_mesh.GetOffsetStart<MeshElement::ELEMENT>();
      mesh_sizes[2] = a_mesh.GetLocalCount<MeshElement::ELEMENT>();
      break;
    }
    case MeshElement::VERTEX: {
      DEBUG_ASSERT(false, global_assert{}, DebugLevel::ALWAYS{},
                   "Output of vertex data is not yet handled.");
      // global_length = a_mesh.GetGlobalCount<MeshElement::VERTEX>();
      // local_offset = a_mesh.GetOffsetStart<MeshElement::VERTEX>();
      // local_count = a_mesh.GetLocalCount<MeshElement::VERTEX>();
      break;
    }
    default:
      DEBUG_ASSERT(false, global_assert{}, DebugLevel::ALWAYS{},
                   "Unknown mesh element type of : " +
                       std::to_string(static_cast<int>(a_type)));
  }
  DEBUG_ASSERT(mesh_sizes[1] < mesh_sizes[0], global_assert{},
               DebugLevel::CHEAP{},
               "Local offset must be less than the global size of the data.\n "
               "Local offset: " +
                   std::to_string(mesh_sizes[1]) +
                   "\nGlobal size: " + std::to_string(mesh_sizes[0]));
  DEBUG_ASSERT(
      mesh_sizes[2] <= mesh_sizes[0], global_assert{}, DebugLevel::CHEAP{},
      "Local size must be less than or equal to the global size of the data.\n"
      "Local size: " +
          std::to_string(mesh_sizes[2]) +
          "\nGlobal size: " + std::to_string(mesh_sizes[0]));

  return mesh_sizes;
}

// Note: Taken from adios2stream.cpp in MFEM, added
// by William F Godoy godoywf@ornl.gov on Jan 22, 2020
std::string IO::VTKSchema(void) const {
  std::string VTKSchema = R"(
<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.2" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="NumOfVertices" NumberOfCells="NumOfElements">
      <Points>
        <DataArray Name="vertices" />
      </Points>
      <Cells>
        <DataArray Name="connectivity" />
        <DataArray Name="types" />
      </Cells>)";
  VTKSchema += '\n';

  // Add Cell Data
  VTKSchema += R"(<CellData>)";
  VTKSchema += '\n';

  for (const auto& name : element_data_variables_m) {
    VTKSchema += R"(<DataArray Name=")";
    VTKSchema += name;
    VTKSchema += R"(" />)";
    VTKSchema += '\n';
  }
  VTKSchema += R"(</CellData>)";
  VTKSchema += '\n';

  // Add Point Data
  VTKSchema += R"(<PointData>)";
  VTKSchema += '\n';
  for (const auto& name : point_data_variables_m) {
    VTKSchema += R"(<DataArray Name=")";
    VTKSchema += name;
    VTKSchema += R"(" />)";
    VTKSchema += '\n';
  }

  VTKSchema += "        <DataArray Name=\"TIME\">\n";
  VTKSchema += "          TIME\n";
  VTKSchema += "        </DataArray>\n";

  VTKSchema += R"(</PointData>)";

  // Wrap up rest of XML
  VTKSchema += R"(
       </Piece>
     </UnstructuredGrid>
   </VTKFile>)";

  return VTKSchema;
}

std::array<int, 2> IO::MinMaxVariableOrder(void) const {
  if (point_data_orders_m.size() == 0) {
    return std::array<int, 2>{1, 1};
  } else {
    std::array<int, 2> min_max{{INT_MAX, -INT_MAX}};
    for (const auto element : point_data_orders_m) {
      min_max[0] = std::min(min_max[0], element);
      min_max[1] = std::max(min_max[1], element);
    }
    DEBUG_ASSERT(min_max[0] <= min_max[1], global_assert{}, DebugLevel::CHEAP{},
                 "minimum value= " + std::to_string(min_max[0]) +
                     " maximum value= " + std::to_string(min_max[1]));
    return min_max;
  }
}

}  // namespace chyps
