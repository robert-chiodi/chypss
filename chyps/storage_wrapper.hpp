// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// \file storage_wrapper.hpp
/// \brief File that contains implementation of StorageWrapper.

#ifndef CHYPS_STORAGE_WRAPPER_HPP_
#define CHYPS_STORAGE_WRAPPER_HPP_

#include <cstddef>

namespace chyps {
/// \class StorageWrapper storage_wrapper.hpp chyps/storage_wrapper.hpp
/// \brief Class that stores either a reference to storage somewhere else
/// or internally stores the variable. Useful to respect constness of
/// given variables.
template <class Type>
class StorageWrapper {
 public:
  /// \brief Default constructor that holds no data.
  StorageWrapper(void);

  /// \brief Copy constructor. If SetValues was used with deep_copy = false in
  /// a_other, a shallow copy will be used, and referenced_values_m will store
  /// the same address as in a_other. If deep_copy = true was used for SetValues
  /// in a_other, a deep_copy of the data will be used in this copy constructor
  /// as well.
  StorageWrapper(const StorageWrapper& a_other);

  /// \brief Move constructor. If SetValues was used with deep_copy = false in
  /// a_other, a shallow copy will be used, and referenced_values_m will store
  /// the same address as in a_other. If deep_copy = true was used for SetValues
  /// in a_other, a deep_copy of the data will be used in this copy constructor
  /// as well.
  StorageWrapper(StorageWrapper&& a_other);

  /// \brief Copy assignment. If SetValues was used with deep_copy = false in
  /// a_other, a shallow copy will be used, and referenced_values_m will store
  /// the same address as in a_other. If deep_copy = true was used for SetValues
  /// in a_other, a deep_copy of the data will be used in this copy constructor
  /// as well.
  StorageWrapper& operator=(const StorageWrapper& a_other);

  /// \brief Move assignment. If SetValues was used with deep_copy = false in
  /// a_other, a shallow copy will be used, and referenced_values_m will store
  /// the same address as in a_other. If deep_copy = true was used for SetValues
  /// in a_other, a deep_copy of the data will be used in this copy constructor
  /// as well.
  StorageWrapper& operator=(StorageWrapper&& a_other);

  /// \brief Sets values wrapped by this class. If a_deep_copy = true, internal
  /// storage is allocated and the values in a_values are copied over.
  /// Otherwise, the pointer is directed to a_values.
  void SetValues(const std::size_t a_size, const Type* a_values,
                 const bool a_deep_copy = false);

  /// \brief Returns size of current wrapped data.
  std::size_t Size(void) const;

  /// \brief Returns pointer to current wrapped data.
  const Type* Data(void) const;

  /// \brief If any data is owned, it is deleted.
  ~StorageWrapper(void);

 private:
  bool OwnsData(void) const;
  void ClearData(void);
  void AllocateAndCopyData(const std::size_t a_size, const Type* a_values);

  Type* stored_values_m;
  const Type* referenced_values_m;
  std::size_t size_m;
};

}  // namespace chyps

#include "chyps/storage_wrapper.tpp"

#endif  // CHYPS_STORAGE_WRAPPER_HPP_
