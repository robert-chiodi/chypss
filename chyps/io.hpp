// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// \file io.hpp
/// \brief I/O using ADIOS2. https://github.com/ornladios/ADIOS2.

#ifndef CHYPS_IO_HPP_
#define CHYPS_IO_HPP_

#include <array>
#include <string>
#include <vector>

#include <adios2.h>
#include <mfem/mfem.hpp>

#include "chyps/mesh.hpp"
#include "chyps/mpi_parallel.hpp"

namespace chyps {

enum class ExistenceLocation { NONE = 0, READ, WRITE, BOTH };

/// \brief IO class is currently written to write out chunks for each MPI Rank.
/// This means the same number of ranks must be used to read and write data.
class IO {
 public:
  IO(void) = delete;

  /// \brief Construct IO for processes on the provided communicator.
  IO(const MPIParallel& a_mpi_comm, const std::string a_io_name);

  /// \brief Set Parameter for the ADIOS2 engine selected. See
  /// ADIOS2 documentation
  /// (https://adios2.readthedocs.io/en/latest/engines/engines.html)
  /// for available parameters.
  void SetParameter(const std::string a_key, const std::string a_value);

  /// \brief Setup the IO for writing. The provided a_file_name
  /// will have the correct extension appended to it for the
  /// given engine.
  void SetWrite(std::string a_file_name);

  /// \brief Setup the IO for reading. The provided a_file_name
  /// will have the correct extension appended to it for the
  /// given engine. (i.e., to read "test.bp", a_file_name should be
  /// "test".
  void SetRead(std::string a_file_name);

  /// \brief Return whether the IO object is acrively in write mode.
  bool IsWriteModeActive(void) const;

  /// \brief Return whether the IO object is acrively in read mode.
  bool IsReadModeActive(void) const;

  /// \brief Add mesh variable for later export (through PutDeferred or
  /// PutImmediate).
  template <class T>
  void AddVariableForMesh(const std::string a_variable_name, const Mesh& a_mesh,
                          const MeshElement a_type);

  /// \brief Add a mfem::DenseMatrix variable for later export (through
  /// PutDeferred or PutImmediate). Column-major ordering will be used.
  void AddMatrixForMesh(const std::string a_variable_name, const Mesh& a_mesh,
                        const MeshElement a_type,
                        const std::size_t a_number_of_rows,
                        const std::size_t a_number_of_columns);

  /// \brief Add ParGridFunction variable for later export (through PutDeferred
  /// or PutImmediate). Underlying data assumed to be of type double. The data
  /// is also marked as a point variable using the MarkAsPointVariableMethod.
  /// It is assumed the GridFunction is on a constant element-order mesh.
  ///
  /// NOTE: Currently only implemented as all ranks write to their own section.
  void AddVariableForGridFunction(
      const std::string a_variable_name,
      const mfem::ParFiniteElementSpace& a_element_space,
      const bool a_dimensions_static);

  /// \brief Add mfem::Vector of True Dofs variable for later export (through
  /// PutDeferred or PutImmediate). Underlying data assumed to be of type
  /// double. The data is also marked as a point variable using the
  /// MarkAsPointVariableMethod. It is assumed the GridFunction is on a constant
  /// element-order mesh.
  ///
  /// NOTE: Currently only implemented as all ranks write to their own section.
  void AddVariableForTrueDofs(
      const std::string a_variable_name,
      const mfem::ParFiniteElementSpace& a_element_space,
      const bool a_dimensions_static);

  /// \brief Add variable for later export (through PutDeferred or
  /// PutImmediate). The global shape, local start of the data, and local count
  /// of the data must be specified. The boolean a_dimensions_static indicate
  /// whether the dimensions might change (false) or will always remain the same
  /// (true).
  template <class T>
  void AddVariable(const std::string a_variable_name,
                   const adios2::Dims& a_global_shape,
                   const adios2::Dims& a_local_data_start,
                   const adios2::Dims& a_local_data_count,
                   const bool a_dimensions_static);

  /// \brief Add variable for later export (through PutDeferred or
  /// PutImmediate). The data is assumed to be a scalar.
  template <class T>
  void AddVariable(const std::string a_variable_name);

  /// \brief Add local variable for each rank to later export (through
  /// PutDeferred or PutImmediate).
  template <class T>
  void AddVariable(const std::string a_variable_name,
                   const adios2::Dims& a_local_data_count,
                   const bool a_dimensions_static);

  /// \brief Add variable for later export (through PutDeferred or
  /// PutImmediate). This will only be performed on the Root MPI rank. The
  /// boolean a_dimensions_static indicate whether the dimensions might change
  /// (false) or will always remain the same (true).
  template <class T>
  void RootAddVariable(const std::string a_variable_name,
                       const adios2::Dims& a_local_data_count,
                       const bool a_dimensions_static);

  /// \brief Add a variable for later export (through PutDeferred or
  /// PutImmediate). This will only be performed on the Root MPI rank. It is
  /// also assumed that the data is a scalar.
  template <class T>
  void RootAddVariable(const std::string a_variable_name);

  /// \brief Mark the variable with a_name (and separately added with an
  /// AddVariable call, as an element variable. This will be used when exporting
  /// to VTK XML for reading in Paraview.
  void MarkAsElementVariable(const std::string a_name);

  /// \brief Mark the variable with a_name (and separately added with an
  /// AddVariable call, as a point variable. This will be used when exporting
  /// to VTK XML for reading in Paraview. The order is also given so
  /// the appropriate underlying mesh can also be created for visualization.
  ///
  /// NOTE: MFEM ParGridFunction variables are considered point variables.
  void MarkAsPointVariable(const std::string a_name, const int a_order);

  /// \brief Writes attribute (with a_name) and the associated a_data.
  ///
  /// NOTE: a_name must be a unique attribute name, not used before.
  /// NOTE: Attribute will be written once per MPI process.
  template <class T>
  void WriteAttribute(const std::string& a_name, const T& a_data);

  /// \brief Writes attribute (with a_name) and the associated a_data.
  ///
  /// NOTE: a_name must be a unique attribute name, not used before.
  /// NOTE: Attribute will be written once per MPI process.
  template <class T>
  void WriteAttribute(const std::string& a_name, const std::vector<T>& a_data);

  /// \brief Writes attribute (with a_name) for the given variable (with name
  /// a_variable_name). This variable must have been previously added with the
  /// AddVariable method. The data a_data is written with the attribute.
  ///
  /// NOTE: a_name must be a unique attribute name, not used before.
  /// NOTE: Attribute will be written once per MPI process.
  template <class T>
  void WriteAttributeForVariable(const std::string& a_variable_name,
                                 const std::string& a_name,
                                 const std::vector<T>& a_data);

  /// \brief Writes attribute (with a_name) for the given variable (with name
  /// a_variable_name). This variable must have been previously added with the
  /// AddVariable method. The data a_data is written with the attribute.
  ///
  /// NOTE: a_name must be a unique attribute name, not used before.
  /// NOTE: Attribute will be written once per MPI process.
  template <class T>
  void WriteAttributeForVariable(const std::string& a_variable_name,
                                 const std::string& a_name, const T& a_data);

  /// \brief Writes attribute (with a_name) and the associated a_data from the
  /// Root MPI rank. Does nothing for all other ranks.
  ///
  /// NOTE: a_name must be a unique attribute name, not used before.
  template <class T>
  void RootWriteAttribute(const std::string& a_name, const T& a_data);

  /// \brief Writes attribute (with a_name) and the associated a_data from the
  /// Root MPI rank. Does nothing for all other ranks.
  ///
  /// NOTE: a_name must be a unique attribute name, not used before.
  template <class T>
  void RootWriteAttribute(const std::string& a_name,
                          const std::vector<T>& a_data);

  /// \brief Writes attribute (with a_name) for the given variable (with name
  /// a_variable_name). This variable must have been previously added with the
  /// AddVariable method. The data a_data is written with the attribute. Write
  /// is performed only from the Root MPI rank. Nothing is done for other ranks.
  ///
  /// NOTE: a_name must be a unique attribute name, not used before.
  template <class T>
  void RootWriteAttributeForVariable(const std::string& a_variable_name,
                                     const std::string& a_name,
                                     const std::vector<T>& a_data);

  /// \brief Reads attribute (with a_name) from the file currently open for
  /// reading. Requires IO to be in read mode and for attribute to exist in the
  /// file.  Returns true if the attribute was found, and a_data is filled with
  /// the attribute data. Otherwise returns false and a_data is left unchanged.
  template <class T>
  bool ReadAttribute(const std::string& a_name, std::vector<T>& a_data);

  /// \brief Reads attribute (with a_name) from the file currently open for
  /// reading. Requires IO to be in read mode and for attribute to exist in the
  /// file. Returns true if the attribute was found, and a_data is filled with
  /// the attribute data. Otherwise returns false and a_data is left unchanged.
  template <class T>
  bool ReadAttributeForVariable(const std::string& a_variable_name,
                                const std::string& a_name,
                                std::vector<T>& a_data);

  /// \brief Begin step inside the ADIOS2 file. This, along with EndStep, can be
  /// seen as Gates between which all calls to PutDeferred and GetDeferred are
  /// buffered. All calls will then be fulfilled during call to EndStep().
  void BeginWriteStep(const uint64_t a_cycle, const double a_time,
                      const double a_dt);

  /// \brief Returns if an IO step is currently ongoing (i.e., BeginStep has
  /// been called for this IO object but EndStep has not been called yet).
  bool OngoingWriteStep(void) const;

  /// \brief End step for the ADIOS2 file. See documentation for BeginStep
  /// method for usage.
  void EndWriteStep(void);

  /// \brief Close current engine contained by the IO object (either a read
  /// engine or write engine. If the IO object has yet to be set to read or
  /// write, this method does nothing.
  void CloseWriteEngine(void);

  /// \brief See BeginWriteStep. Does same but for engine open to read files.
  void BeginReadStep(void);

  /// \brief See OngoingWriteStep. Does same but for engine open to read files.
  bool OngoingReadStep(void) const;

  /// \brief See EndWriteStep. Does same but for engine open to read files.
  void EndReadStep(void);

  /// \brief See CloseWriteEngine. Does same but for engine open to read files.
  void CloseReadEngine(void);

  /// \brief Checks if a variable exists in either writing or reading.
  ///
  /// NOTE: If file is not open for reading/writing, it is considered as not
  /// existing there.
  template <class T>
  ExistenceLocation DoesVariableExist(const std::string& a_name);

  /// \brief Checks if an attribute exists in either writing or reading.
  ///
  /// NOTE: If file is not open for reading/writing, it is considered as not
  /// existing there.
  template <class T>
  ExistenceLocation DoesAttributeExist(const std::string& a_name);

  /// \brief Checks if an attribute for the variable exists in either writing or
  /// reading.
  ///
  /// NOTE: If file is not open for reading/writing, it is considered as not
  /// existing there.
  template <class T>
  ExistenceLocation DoesAttributeExist(const std::string& a_variable_name,
                                       const std::string& a_name);

  /// \brief Taking an mfem::Vector, provide the underlying pointer
  /// for writing. The underlying type is assumed to be double.
  ///
  /// The variable with name a_variable_name must have been previously
  /// added through a call to AddVariable. Furthermore, the mfem::Vector
  /// must be of the correct length, depending on whether it was
  /// added as a TrueDof field or a ParGridFunction.
  ///
  /// Note: DO NOT modify data contained in
  /// a_data or the validity of the pointer until after the call to EndStep.
  void PutDeferred(const std::string& a_variable_name,
                   const mfem::Vector& a_vector);

  /// \brief Taking a std::vector of mfem::DenseMatrix objects, flatten
  /// the matrix into an ADIOS2::Span for writing. The matrix will be written
  /// in column major ordering to match that used by mfem::DenseMatrix.
  /// The underlying type is assumed to be double and all matrices
  /// in `a_matrix_list` are assumed to be the same shape.
  ///
  /// The variable with name a_variable_name must have been previously
  /// added through a call to AddMatrixForMesh.
  ///
  /// Note: DO NOT modify data contained in
  /// a_data or the validity of the pointer until after the call to EndStep.
  // FIXME : Check if filling a span would make this a PutImmediate or
  // PutDeferred
  void PutDeferred(const std::string& a_variable_name,
                   const std::vector<mfem::DenseMatrix>& a_matrix_list);

  /// \brief Deferred Put of data for I/O. Will hold reference to data
  /// given to it between BeginStep and EndStep and write during EndStep.
  /// This should be used by default for better performance.
  ///
  /// The variable of type T and name a_variable_name must have been previously
  /// added through a call to AddVariable. Note: DO NOT modify data contained in
  /// a_data or the validity of the pointer until after the call to EndStep.
  template <class T>
  void PutDeferred(const std::string& a_variable_name, const T* a_data);

  /// \brief Deferred put of data that returns a span to be filled. This allows
  /// copying the data to be Put directly into the buffer being written without
  /// an auxiliary storage array being defined first by user.
  ///
  /// The variable of type T and name a_variable_name must have been previously
  /// added through a call to AddVariable.
  template <class T>
  typename adios2::Variable<T>::Span PutSpan(
      const std::string& a_variable_name);

  /// \brief Deferred Put of data for I/O. Will hold reference to data
  /// given to it between BeginStep and EndStep and write during EndStep. Only
  /// performed on the Root MPI rank. This should be used by default for better
  /// performance.
  ///
  /// The variable of type T and name a_variable_name must have been previously
  /// added through a call to AddVariable. Note: DO NOT modify data contained in
  /// a_data or the validity of the pointer until after the call to EndStep.
  template <class T>
  void RootPutDeferred(const std::string& a_variable_name, const T* a_data);

  /// \brief Deferred Get of data for I/O. Will hold reference to data
  /// given to it between BeginStep and EndStep and populate it during EndStep.
  /// The a_local_start and a_local_count indicate the (potentially
  /// multidimensional) slice of data that should be read from the larger global
  /// size. This method should be used over the immediate version  for better
  /// performance.
  ///
  // Note: DO NOT modify validity of the
  /// pointer until after the call to EndStep. Do not try to use data until
  /// after EndStep.
  ///
  /// Note: When perforAming a Get, the variable SHOULD NOT be added first
  /// via a AddVariable method. The variables will be populated from the file
  /// open for reading.
  template <class T>
  void GetDeferred(const std::string& a_variable_name,
                   const adios2::Dims& a_local_start,
                   const adios2::Dims& a_local_count, T* a_data);

  /// \brief Deferred Get of data for I/O. Will hold reference to data
  /// given to it between BeginStep and EndStep and populate it during EndStep.
  /// Offset and local size determined via the passed mesh object.  This method
  /// should be used over the immediate version for better performance.

  /// Note: DO NOT modify validity of the
  /// pointer until after the call to EndStep. Do not try to use data until
  /// after EndStep.
  ///
  /// Note: When perforAming a Get, the variable SHOULD NOT be added first
  /// via a AddVariable method. The variables will be populated from the file
  /// open for reading.
  template <class T>
  void GetDeferred(const std::string& a_variable_name, const Mesh& a_mesh,
                   const MeshElement a_type, T* a_data);

  /// \brief Deferred Get of data for I/O. Will hold reference to data
  /// given to it between BeginStep and EndStep and populate it during EndStep.
  /// This should be used to read blocks of data from locally written files
  /// (where one rank writes to one block). This method should be used over the
  /// immediate version  for better performance.
  ///
  /// Note: DO NOT modify validity of the
  /// pointer until after the call to EndStep. Do not try to use data until
  /// after EndStep.
  ///
  /// Note: When perforAming a Get, the variable SHOULD NOT be added first
  /// via a AddVariable method. The variables will be populated from the file
  /// open for reading.
  template <class T>
  void GetDeferredBlock(const std::string& a_variable_name,
                        const adios2::Dims& a_local_start,
                        const adios2::Dims& a_local_count, T* a_data);

  /// \brief Immediate Put of data for I/O. Will write data from a_data before
  /// returning.
  ///
  /// The variable of type T and name a_variable_name must have been previously
  /// added through a call to AddVariable.
  ///
  // Note: Prefer PutDeferred for better performance.
  template <class T>
  void PutImmediate(const std::string& a_variable_name, const T* a_data);

  /// \brief Immediate Put of data for I/O. Will write data from a_data before
  /// returning. Only performed on the Root MPI rank.
  ///
  /// The variable of type T and name a_variable_name must have been previously
  /// added through a call to AddVariable.
  ///
  // Note: Prefer PutDeferred for better performance.
  template <class T>
  void RootPutImmediate(const std::string& a_variable_name, const T* a_data);

  /// \brief Immediate Get of data for I/O. Will read data into a_data before
  /// returning. The a_local_start and a_local_count indicate the
  /// (potentially multidimensional) slice of data that should be read from the
  /// larger global size.
  ///
  ///
  /// Note: Prefer GetDeferred for better performance.
  ///
  /// Note: When performing a Get, the variable SHOULD NOT be added first
  /// via a AddVariable method. The variables will be populated from the file
  /// open for reading.
  template <class T>
  void GetImmediate(const std::string& a_variable_name,
                    const adios2::Dims& a_local_start,
                    const adios2::Dims& a_local_count, T* a_data);

  /// \brief Immediate Get of data for I/O. Will read data into a_data before
  /// returning. The local start and local count are determined from the passed
  /// Mesh.
  ///
  ///
  /// Note: Prefer GetDeferred for better performance.
  ///
  /// Note: When performing a Get, the variable SHOULD NOT be added first
  /// via a AddVariable method. The variables will be populated from the file
  /// open for reading.
  template <class T>
  void GetImmediate(const std::string& a_variable_name, const Mesh& a_mesh,
                    const MeshElement a_type, T* a_data);

  /// \brief Immediate Get of data for I/O. Will read data into a_data before
  /// returning. Used to read in values written with a RootAddVariable call.
  ///
  ///
  /// Note: Prefer GetDeferred for better performance.
  ///
  /// Note: When performing a Get, the variable SHOULD NOT be added first
  /// via a AddVariable method. The variables will be populated from the file
  /// open for reading.
  template <class T>
  void GetImmediateValue(const std::string& a_variable_name, T* a_data);

  /// \brief Immediate Get of data for I/O. Will read data into a_data before
  /// returning.
  ///
  ///
  /// Note: Prefer GetDeferred for better performance.
  ///
  /// Note: When perforAming a Get, the variable SHOULD NOT be added first
  /// via a AddVariable method. The variables will be populated from the file
  /// open for reading.
  template <class T>
  void GetImmediateBlock(const std::string& a_variable_name,
                         const adios2::Dims& a_local_start,
                         const adios2::Dims& a_local_count, T* a_data);

  /// \brief Immediate Get of data for I/O. Will resize data automatically and
  /// return it.
  ///
  ///
  /// Note: Prefer GetDeferred for better performance.
  ///
  /// Note: When perforAming a Get, the variable SHOULD NOT be added first
  /// via a AddVariable method. The variables will be populated from the file
  /// open for reading.
  template <class T>
  void GetImmediateBlock(const std::string& a_variable_name,
                         std::vector<T>& a_data);

  /// \brief Immediate Get of data for I/O. Will resize a_matrix_list and return
  /// filled. Only reads a rank-local block of the data.
  ///
  ///
  /// Note: Prefer GetDeferred for better performance.
  ///
  /// Note: When performing a Get, the variable SHOULD NOT be added first
  /// via a AddVariable method. The variables will be populated from the file
  /// open for reading.
  void GetImmediateBlock(const std::string& a_variable_name,
                         std::vector<mfem::DenseMatrix>& a_matrix_list);

  /// \brief Effectively performs a flush on all PutDeferred variables
  /// that have yet to be written.
  void PerformPuts(void);

  /// \brief Writes the XML Schema to the current file open for writing
  /// as an attribute. Used to allow viewing of data in Paraview.
  void WriteXMLSchema(void);

  /// \brief Return the minimum and maximum orders used for the
  /// point_data_variables addd.
  std::array<int, 2> MinMaxVariableOrder(void) const;

  /// \brief Destructor that closes all open engines and frees memory.
  ~IO(void);

 private:
  void SetEngineType(void);
  void SetDefaultParameters(void);
  static ExistenceLocation FindExistence(const bool a_read, const bool a_write);
  std::array<std::size_t, 3> GetMeshSizes(const Mesh& a_mesh,
                                          const MeshElement a_type) const;

  template <class T>
  void Put(const std::string& a_variable_name, const T* a_data,
           const adios2::Mode a_mode);
  template <class T>
  void Get(const std::string& a_variable_name,
           const adios2::Dims& a_local_start, const adios2::Dims& a_local_count,
           T* a_data, const adios2::Mode a_mode);
  template <class T>
  void GetBlock(const std::string& a_variable_name,
                const adios2::Dims& a_local_start,
                const adios2::Dims& a_local_count, T* a_data,
                const adios2::Mode a_mode);
  template <class T>
  void GetBlock(const std::string& a_variable_name, std::vector<T>& a_data,
                const adios2::Mode a_mode);

  template <class T>
  void GetSingleValue(const std::string& a_variable_name, T* a_data,
                      const adios2::Mode a_mode);

  std::string VTKSchema(void) const;

  adios2::ADIOS io_driver_m;
  const MPIParallel& mpi_session_m;
  adios2::IO read_m;
  adios2::IO write_m;
  adios2::Engine* read_engine_m;
  adios2::Engine* write_engine_m;
  std::vector<std::string> element_data_variables_m;
  std::vector<std::string> point_data_variables_m;
  std::vector<int> point_data_orders_m;
  bool active_write_step_m;
  bool active_read_step_m;
};

}  // namespace chyps

#include "chyps/io.tpp"

#endif  // CHYPS_IO_HPP_
