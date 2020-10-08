// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef CHYPS_STORAGE_WRAPPER_TPP_
#define CHYPS_STORAGE_WRAPPER_TPP_

#include <algorithm>
#include <cassert>

namespace chyps {

template <class Type>
StorageWrapper<Type>::StorageWrapper(void)
    : stored_values_m(nullptr), referenced_values_m(nullptr), size_m(0) {}

template <class Type>
StorageWrapper<Type>::StorageWrapper(const StorageWrapper& a_other) {
  if (a_other.OwnsData()) {
    stored_values_m = nullptr;
    referenced_values_m = nullptr;
    size_m = 0;
    this->AllocateAndCopyData(a_other.Size(), a_other.Data());
  } else {
    stored_values_m = nullptr;
    referenced_values_m = a_other.Data();
    size_m = a_other.Size();
  }
}

template <class Type>
StorageWrapper<Type>::StorageWrapper(StorageWrapper&& a_other)
    : stored_values_m(a_other.stored_values_m),
      referenced_values_m(a_other.referenced_values_m),
      size_m(a_other.size_m) {
  a_other.stored_values_m = nullptr;
  a_other.referenced_values_m = nullptr;
  a_other.size_m = 0;
}

// FIXME: Make more efficient. Do not delete and reallocate if already have
// allocated amount same size as the other.
template <class Type>
StorageWrapper<Type>& StorageWrapper<Type>::operator=(
    const StorageWrapper& a_other) {
  if (this != &a_other) {
    this->ClearData();
    if (a_other.OwnsData()) {
      this->AllocateAndCopyData(a_other.Size(), a_other.Data());
    } else {
      referenced_values_m = a_other.Data();
      size_m = a_other.Size();
    }
  }
  return *this;
}

template <class Type>
StorageWrapper<Type>& StorageWrapper<Type>::operator=(
    StorageWrapper&& a_other) {
  if (this != &a_other) {
    this->ClearData();
    stored_values_m = a_other.stored_values_m;
    referenced_values_m = a_other.referenced_values_m;
    size_m = a_other.size_m;

    a_other.stored_values_m = nullptr;
    a_other.referenced_values_m = nullptr;
    a_other.size_m = 0;
  }
  return *this;
}

// FIXME: Make it so that if allocated storage of correct size already, just use
// that storage
template <class Type>
void StorageWrapper<Type>::SetValues(const std::size_t a_size,
                                     const Type* a_values,
                                     const bool a_deep_copy) {
  if (this->OwnsData()) {
    this->ClearData();
  }
  if (a_deep_copy) {
    this->AllocateAndCopyData(a_size, a_values);
  } else {
    size_m = a_size;
    referenced_values_m = a_values;
  }
}

template <class Type>
std::size_t StorageWrapper<Type>::Size(void) const {
  return size_m;
}

template <class Type>
const Type* StorageWrapper<Type>::Data(void) const {
  return this->OwnsData() ? stored_values_m : referenced_values_m;
}

template <class Type>
bool StorageWrapper<Type>::OwnsData(void) const {
  return stored_values_m != nullptr;
}

template <class Type>
StorageWrapper<Type>::~StorageWrapper(void) {
  this->ClearData();
}

template <class Type>
void StorageWrapper<Type>::ClearData(void) {
  delete[] stored_values_m;
  stored_values_m = nullptr;
  referenced_values_m = nullptr;
  size_m = 0;
}

template <class Type>
void StorageWrapper<Type>::AllocateAndCopyData(const std::size_t a_size,
                                               const Type* a_values) {
  assert(stored_values_m == nullptr);
  assert(referenced_values_m == nullptr);
  size_m = a_size;
  stored_values_m = new Type[size_m];
  std::copy(a_values, a_values + a_size, stored_values_m);
}

}  // namespace chyps
#endif  // CHYPS_STORAGE_WRAPPER_TPP_
