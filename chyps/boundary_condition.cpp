// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chyps/boundary_condition.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>

#include "chyps/logger.hpp"

namespace chyps {

BoundaryCondition::BoundaryCondition(void)
    : bc_type_m(BoundaryConditionType::HOMOGENEOUS_DIRICHLET),
      values_m(),
      spatial_varying_m(false),
      time_varying_m(false) {
  SPDLOG_LOGGER_INFO(
      MAIN_LOG, "Constructed boundary condition: type({}), time-varying({})",
      static_cast<int>(BoundaryConditionType::HOMOGENEOUS_DIRICHLET), false);
}

BoundaryCondition::BoundaryCondition(const BoundaryConditionType a_type,
                                     const bool a_spatial_varying,
                                     const bool a_time_varying)
    : bc_type_m(a_type),
      values_m(),
      spatial_varying_m(a_spatial_varying),
      time_varying_m(a_time_varying) {
  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Constructed boundary condition: type({}), "
                     "spatial-varying({}), time-varying({})",
                     static_cast<int>(a_type), a_spatial_varying,
                     a_time_varying);
  assert(this->LogicalConditionsSet());
}

void BoundaryCondition::SetValues(const double a_value) {
  SPDLOG_LOGGER_INFO(
      MAIN_LOG,
      "Setting constant value of {:8.6E} for boundary condition of type {}",
      a_value, static_cast<int>(this->GetBCType()));
  assert(!this->IsSpatiallyVarying());
  values_m.SetValues(1, &a_value, true);
}

void BoundaryCondition::SetValues(const double* a_value,
                                  const bool a_deep_copy) {
  SPDLOG_LOGGER_INFO(
      MAIN_LOG, "Setting reference to double for boundary condition of type {}",
      static_cast<int>(this->GetBCType()));
  assert(!this->IsSpatiallyVarying());
  values_m.SetValues(1, a_value, a_deep_copy);
}

void BoundaryCondition::SetValues(const std::size_t a_size,
                                  const double* a_values, const int* a_indices,
                                  const bool a_deep_copy) {
  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Setting boundary condition of type {} with {} entries. "
                     "Are values deep copied? ({})",
                     static_cast<int>(this->GetBCType()), a_size, a_deep_copy);
  values_m.SetValues(a_size, a_values, a_deep_copy);
  indices_m.SetValues(a_size, a_indices, a_deep_copy);
}

BoundaryConditionType BoundaryCondition::GetBCType(void) const {
  return bc_type_m;
}

const StorageWrapper<double>& BoundaryCondition::GetValues(void) const {
  return values_m;
}

const StorageWrapper<int>& BoundaryCondition::GetIndices(void) const {
  return indices_m;
}

bool BoundaryCondition::IsSpatiallyVarying(void) const {
  return spatial_varying_m;
}

bool BoundaryCondition::IsTimeVarying(void) const { return time_varying_m; }

bool BoundaryCondition::LogicalConditionsSet(void) const {
  if (this->IsSpatiallyVarying() || this->IsTimeVarying()) {
    if (this->GetBCType() != BoundaryConditionType::DIRICHLET &&
        this->GetBCType() != BoundaryConditionType::NEUMANN) {
      std::cout << "Homogeneous condition used but boundary condition set to "
                   "spatially-varying and/or time-varying. "
                << std::endl;
      return false;
    }
  }
  return true;
}

}  // namespace chyps
