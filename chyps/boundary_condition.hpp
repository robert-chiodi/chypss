// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// \file boundary_condition.hpp
/// \brief File to handle boundary condition details used to configure boundary
/// conditions for solvers.

#ifndef CHYPS_BOUNDARY_CONDITION_HPP_
#define CHYPS_BOUNDARY_CONDITION_HPP_

#include <cstddef>

#include "chyps/storage_wrapper.hpp"

namespace chyps {

/// \enum BoundaryConditionType boundary_condition.hpp
/// chyps/boundary_condition.hpp
/// \brief Enumerator that provides options for
/// boundary condition configurations.
enum class BoundaryConditionType {
  HOMOGENEOUS_DIRICHLET = 0,
  HOMOGENEOUS_NEUMANN,
  DIRICHLET,
  NEUMANN
};

/// \class BoundaryCondition boundary_condition.hpp chyps/boundary_condition.tpp
class BoundaryCondition {
 public:
  /// \brief Default constructor. Default boundary condition type is
  /// BoundaryConditionType::HOMOGENEOUS_DIRICHLET.
  BoundaryCondition(void);

  /// \brief Construct while setting BoundaryConditionType to a_type. Data must
  /// be separately allocated through SetValues(...).
  BoundaryCondition(const BoundaryConditionType a_type,
                    const bool a_spatial_varying = false,
                    const bool a_time_varying = false);

  /// \brief Copy constructor. Please note, if underlying data
  /// supplied via SetValues was set with deep_copy = true, the copy
  /// will involve a deep copy. However, if deep_copy = false was
  /// used, a shallow copy will be performed of the underlying values.
  BoundaryCondition(const BoundaryCondition& a_other) = default;

  /// \brief Move constructor. Please note, if underlying data
  /// supplied via SetValues was set with deep_copy = true, the copy
  /// will involve a deep copy. However, if deep_copy = false was
  /// used, a shallow copy will be performed of the underlying values.
  BoundaryCondition(BoundaryCondition&& a_other) = default;

  /// \brief Copy assignment. Please note, if underlying data
  /// supplied via SetValues was set with deep_copy = true, the copy
  /// will involve a deep copy. However, if deep_copy = false was
  /// used, a shallow copy will be performed of the underlying values.
  BoundaryCondition& operator=(const BoundaryCondition& a_other) = default;

  /// \brief Move assignment. Please note, if underlying data
  /// supplied via SetValues was set with deep_copy = true, the copy
  /// will involve a deep copy. However, if deep_copy = false was
  /// used, a shallow copy will be performed of the underlying values.
  BoundaryCondition& operator=(BoundaryCondition&& a_other) = default;

  /// \brief Sets a single value for the boundary condition.
  void SetValues(const double a_value);

  /// \brief Single value for boundary condition. Can point to different
  /// (external) storage if a_deep_copy = false;
  void SetValues(const double* a_value, const bool a_deep_copy = false);

  /// \brief Set a_size values for boundary condition. Can point to different
  /// (external) storage if a_deep_copy = false; The vertices and values
  /// will be assumed to use the same value of a_deep_copy.
  void SetValues(const std::size_t a_size, const double* a_values,
                 const int* a_vertices, const bool a_deep_copy = false);

  /// \brief Return type of this boundary condition, which is one of the options
  /// in the BoundaryConditionType enum. Individual solver will use this to
  /// properly handle their boundary conditions.
  BoundaryConditionType GetBCType(void) const;

  /// \brief Returns whether the BC value is constant or spatially varying,
  /// as set during construction.
  bool IsSpatiallyVarying(void) const;

  /// \brief Return if boundary condition might vary in time between iterations.
  ///
  /// This essentially implies the BC should be updated in the solver before
  /// each time advancement.
  bool IsTimeVarying(void) const;

  /// \brief Return values in form of a StorageWrapper.
  const StorageWrapper<double>& GetValues(void) const;

  /// \brief Return the index corresponding to each entry in GetValues.
  const StorageWrapper<int>& GetIndices(void) const;

  /// \brief Default destructor.
  ~BoundaryCondition(void) = default;

 private:
  bool LogicalConditionsSet(void) const;

  BoundaryConditionType bc_type_m;
  StorageWrapper<double> values_m;
  StorageWrapper<int> indices_m;
  bool spatial_varying_m;
  bool time_varying_m;
};

}  // namespace chyps

#endif  // CHYPS_BOUNDARY_CONDITION_HPP_
