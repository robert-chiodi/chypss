// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef CHYPS_IO_TPP_
#define CHYPS_IO_TPP_

#include <cassert>
#include <utility>

namespace chyps {

template <class T>
void IO::AddVariableForMesh(const std::string a_variable_name,
                            const Mesh& a_mesh, const MeshElement a_type) {
  auto io_mesh_sizes = this->GetMeshSizes(a_mesh, a_type);
  this->AddVariable<T>(a_variable_name, {io_mesh_sizes[0]}, {io_mesh_sizes[1]},
                       {io_mesh_sizes[2]}, !a_mesh.UsesAMR());
}

template <class T>
void IO::AddVariable(const std::string a_variable_name,
                     const adios2::Dims& a_global_shape,
                     const adios2::Dims& a_local_data_start,
                     const adios2::Dims& a_local_data_count,
                     const bool a_dimensions_static) {
  write_m.DefineVariable<T>(a_variable_name, a_global_shape, a_local_data_start,
                            a_local_data_count, a_dimensions_static);
}

template <class T>
void IO::AddVariable(const std::string a_variable_name) {
  write_m.DefineVariable<T>(a_variable_name, {adios2::LocalValueDim});
}

template <class T>
void IO::AddVariable(const std::string a_variable_name,
                     const adios2::Dims& a_local_data_count,
                     const bool a_dimensions_static) {
  write_m.DefineVariable<T>(a_variable_name, {}, {}, a_local_data_count,
                            a_dimensions_static);
}

template <class T>
void IO::RootAddVariable(const std::string a_variable_name,
                         const adios2::Dims& a_local_data_count,
                         const bool a_dimensions_static) {
  if (mpi_session_m.IAmRoot()) {
    write_m.DefineVariable<T>(a_variable_name, {a_local_data_count}, {0},
                              {a_local_data_count}, a_dimensions_static);
  }
}

template <class T>
void IO::RootAddVariable(const std::string a_variable_name) {
  if (mpi_session_m.IAmRoot()) {
    write_m.DefineVariable<T>(a_variable_name);
  }
}

template <class T>
void IO::WriteAttribute(const std::string& a_name, const T& a_data) {
  auto adios_attribute = write_m.DefineAttribute(a_name, a_data);
}

template <class T>
void IO::WriteAttribute(const std::string& a_name,
                        const std::vector<T>& a_data) {
  auto adios_attribute =
      write_m.DefineAttribute(a_name, a_data.data(), a_data.size());
}

template <class T>
void IO::WriteAttributeForVariable(const std::string& a_variable_name,
                                   const std::string& a_name,
                                   const std::vector<T>& a_data) {
  auto adios_attribute = write_m.DefineAttribute(
      a_name, a_data.data(), a_data.size(), a_variable_name);
}

template <class T>
void IO::RootWriteAttribute(const std::string& a_name, const T& a_data) {
  if (mpi_session_m.IAmRoot()) {
    auto adios_attribute = write_m.DefineAttribute(a_name, a_data);
  }
}

template <class T>
void IO::RootWriteAttribute(const std::string& a_name,
                            const std::vector<T>& a_data) {
  if (mpi_session_m.IAmRoot()) {
    auto adios_attribute =
        write_m.DefineAttribute(a_name, a_data.data(), a_data.size());
  }
}

template <class T>
void IO::RootWriteAttributeForVariable(const std::string& a_variable_name,
                                       const std::string& a_name,
                                       const std::vector<T>& a_data) {
  if (mpi_session_m.IAmRoot()) {
    auto adios_attribute = write_m.DefineAttribute(
        a_name, a_data.data(), a_data.size(), a_variable_name);
  }
}

template <class T>
bool IO::ReadAttribute(const std::string& a_name, std::vector<T>& a_data) {
  assert(this->IsReadModeActive());
  auto adios2_attribute = read_m.InquireAttribute<T>(a_name);
  if (!adios2_attribute) {
    return false;
  }
  a_data = adios2_attribute.Data();
  return true;
}

template <class T>
bool IO::ReadAttributeForVariable(const std::string& a_var_name,
                                  const std::string& a_name,
                                  std::vector<T>& a_data) {
  assert(this->IsReadModeActive());
  auto adios2_attribute = read_m.InquireAttribute<T>(a_name, a_var_name);
  if (!adios2_attribute) {
    return false;
  }
  a_data = adios2_attribute.Data();
  return true;
}

template <class T>
ExistenceLocation IO::DoesVariableExist(const std::string& a_name) {
  const auto in_read = static_cast<bool>(read_m.InquireVariable<T>(a_name));
  const auto in_write = static_cast<bool>(write_m.InquireVariable<T>(a_name));
  return this->FindExistence(in_read, in_write);
}

template <class T>
ExistenceLocation IO::DoesAttributeExist(const std::string& a_name) {
  const auto in_read = static_cast<bool>(read_m.InquireAttribute<T>(a_name));
  const auto in_write = static_cast<bool>(write_m.InquireAttribute<T>(a_name));
  return this->FindExistence(in_read, in_write);
}

template <class T>
ExistenceLocation IO::DoesAttributeExist(const std::string& a_variable_name,
                                         const std::string& a_name) {
  const auto in_read =
      static_cast<bool>(read_m.InquireAttribute<T>(a_name, a_variable_name));
  const auto in_write =
      static_cast<bool>(write_m.InquireAttribute<T>(a_name, a_variable_name));
  return this->FindExistence(in_read, in_write);
}

template <class T>
void IO::PutDeferred(const std::string& a_variable_name, const T* a_data) {
  assert(this->IsWriteModeActive());
  assert(this->OngoingWriteStep());
  this->Put(a_variable_name, a_data, adios2::Mode::Deferred);
}

template <class T>
typename adios2::Variable<T>::Span IO::PutSpan(
    const std::string& a_variable_name) {
  assert(this->IsWriteModeActive());
  assert(this->OngoingWriteStep());
  auto adios2_variable = write_m.InquireVariable<T>(a_variable_name);
  assert(adios2_variable);
  return write_engine_m->Put(adios2_variable);
}

template <class T>
void IO::RootPutDeferred(const std::string& a_variable_name, const T* a_data) {
  if (mpi_session_m.IAmRoot()) {
    assert(this->IsWriteModeActive());
    assert(this->OngoingWriteStep());
    this->Put(a_variable_name, a_data, adios2::Mode::Deferred);
  }
}

template <class T>
void IO::GetDeferred(const std::string& a_variable_name,
                     const adios2::Dims& a_local_start,
                     const adios2::Dims& a_local_count, T* a_data) {
  assert(this->IsReadModeActive());
  assert(this->OngoingReadStep());
  this->Get(a_variable_name, a_local_start, a_local_count, a_data,
            adios2::Mode::Deferred);
}

template <class T>
void IO::GetDeferred(const std::string& a_variable_name, const Mesh& a_mesh,
                     const MeshElement a_type, T* a_data) {
  assert(this->IsReadModeActive());
  assert(this->OngoingReadStep());
  auto io_mesh_sizes = this->GetMeshSizes(a_mesh, a_type);
  assert(io_mesh_sizes[1] < io_mesh_sizes[0]);
  assert(io_mesh_sizes[2] <= io_mesh_sizes[0]);
  this->Get(a_variable_name, {io_mesh_sizes[1]}, {io_mesh_sizes[2]}, a_data,
            adios2::Mode::Deferred);
}

template <class T>
void IO::GetDeferredBlock(const std::string& a_variable_name,
                          const adios2::Dims& a_local_start,
                          const adios2::Dims& a_local_count, T* a_data) {
  assert(this->IsReadModeActive());
  this->GetBlock(a_variable_name, a_local_start, a_local_count, a_data,
                 adios2::Mode::Deferred);
}

template <class T>
void IO::PutImmediate(const std::string& a_variable_name, const T* a_data) {
  assert(this->IsWriteModeActive());
  this->Put(a_variable_name, a_data, adios2::Mode::Sync);
}

template <class T>
void IO::RootPutImmediate(const std::string& a_variable_name, const T* a_data) {
  if (mpi_session_m.IAmRoot()) {
    assert(this->IsWriteModeActive());
    this->Put(a_variable_name, a_data, adios2::Mode::Sync);
  }
}

template <class T>
void IO::GetImmediate(const std::string& a_variable_name,
                      const adios2::Dims& a_local_start,
                      const adios2::Dims& a_local_count, T* a_data) {
  assert(this->IsReadModeActive());
  this->Get(a_variable_name, a_local_start, a_local_count, a_data,
            adios2::Mode::Sync);
}

template <class T>
void IO::GetImmediate(const std::string& a_variable_name, const Mesh& a_mesh,
                      const MeshElement a_type, T* a_data) {
  assert(this->IsReadModeActive());
  auto io_mesh_sizes = this->GetMeshSizes(a_mesh, a_type);
  assert(io_mesh_sizes[1] < io_mesh_sizes[0]);
  assert(io_mesh_sizes[2] <= io_mesh_sizes[0]);
  this->Get(a_variable_name, {io_mesh_sizes[1]}, {io_mesh_sizes[2]}, a_data,
            adios2::Mode::Sync);
}

template <class T>
void IO::GetImmediateBlock(const std::string& a_variable_name,
                           const adios2::Dims& a_local_start,
                           const adios2::Dims& a_local_count, T* a_data) {
  assert(this->IsReadModeActive());
  this->GetBlock(a_variable_name, {a_local_start}, {a_local_count}, a_data,
                 adios2::Mode::Sync);
}

template <class T>
void IO::GetImmediateBlock(const std::string& a_variable_name,
                           std::vector<T>& a_data) {
  assert(this->IsReadModeActive());
  this->GetBlock(a_variable_name, a_data, adios2::Mode::Sync);
}

template <class T>
void IO::GetImmediateValue(const std::string& a_variable_name, T* a_data) {
  assert(this->IsReadModeActive());
  this->GetSingleValue(a_variable_name, a_data, adios2::Mode::Sync);
}

template <class T>
void IO::Put(const std::string& a_variable_name, const T* a_data,
             const adios2::Mode a_mode) {
  auto adios_variable = write_m.InquireVariable<T>(a_variable_name);
  assert(adios_variable);  // Was found
  assert(this->IsWriteModeActive());
  write_engine_m->Put(adios_variable, a_data, a_mode);
}

template <class T>
void IO::Get(const std::string& a_variable_name,
             const adios2::Dims& a_local_start,
             const adios2::Dims& a_local_count, T* a_data,
             const adios2::Mode a_mode) {
  auto adios_variable = read_m.InquireVariable<T>(a_variable_name);
  assert(adios_variable);  // Was found
  assert(this->IsReadModeActive());
  adios_variable.SetSelection({a_local_start, a_local_count});
  read_engine_m->Get(adios_variable, a_data, a_mode);
}

template <class T>
void IO::GetBlock(const std::string& a_variable_name,
                  const adios2::Dims& a_local_start,
                  const adios2::Dims& a_local_count, T* a_data,
                  const adios2::Mode a_mode) {
  auto adios_variable = read_m.InquireVariable<T>(a_variable_name);
  assert(adios_variable);  // Was found
  assert(this->IsReadModeActive());
  const auto block_id = static_cast<std::size_t>(mpi_session_m.MyRank());
  adios_variable.SetBlockSelection(block_id);
  // FIXME: Might be better to allow users to find the step to load by being
  // closest to a time?
  const auto step_to_load = adios_variable.Steps();
  adios_variable.SetStepSelection({step_to_load - 1, 1});
  // FIXME: Below only correct for 1D arrays
  assert(adios_variable.SelectionSize() == a_local_count[0]);
  read_engine_m->Get(adios_variable, a_data, a_mode);
}

template <class T>
void IO::GetBlock(const std::string& a_variable_name, std::vector<T>& a_data,
                  const adios2::Mode a_mode) {
  auto adios_variable = read_m.InquireVariable<T>(a_variable_name);
  assert(adios_variable);  // Was found
  assert(this->IsReadModeActive());
  const auto block_id = static_cast<std::size_t>(mpi_session_m.MyRank());
  adios_variable.SetBlockSelection(block_id);
  // FIXME: Might be better to allow users to find the step to load by being
  // closest to a time?
  const auto step_to_load = adios_variable.Steps();
  adios_variable.SetStepSelection({step_to_load - 1, 1});
  read_engine_m->Get(adios_variable, a_data, a_mode);
}

template <class T>
void IO::GetSingleValue(const std::string& a_variable_name, T* a_data,
                        const adios2::Mode a_mode) {
  auto adios_variable = read_m.InquireVariable<T>(a_variable_name);
  assert(adios_variable);  // Was found
  assert(this->IsReadModeActive());
  // FIXME: Might be better to allow users to find the step to load by being
  // closest to a time?
  const auto step_to_load = adios_variable.Steps();
  adios_variable.SetStepSelection({step_to_load - 1, 1});
  read_engine_m->Get(adios_variable, a_data, a_mode);
}

}  // namespace chyps

#endif  // CHYPS_IO_TPP_
