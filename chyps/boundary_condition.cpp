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

namespace chyps {

BoundaryCondition::BoundaryCondition(void)
    : bc_type_m(BoundaryConditionType::HOMOGENEOUS_DIRICHLET),
      values_m(),
      time_varying_m(false) {}

BoundaryCondition::BoundaryCondition(const BoundaryConditionType a_type,
                                     const bool a_time_varying)
    : bc_type_m(a_type), values_m(), time_varying_m(a_time_varying) {}

void BoundaryCondition::SetValues(const double a_value) {
  values_m.SetValues(1, &a_value, true);
}

void BoundaryCondition::SetValues(const double* a_value,
                                  const bool a_deep_copy) {
  values_m.SetValues(1, a_value, a_deep_copy);
}

void BoundaryCondition::SetValues(const std::size_t a_size,
                                  const double* a_values,
                                  const bool a_deep_copy) {
  values_m.SetValues(a_size, a_values, a_deep_copy);
}

BoundaryConditionType BoundaryCondition::GetBCType(void) const {
  return bc_type_m;
}

const StorageWrapper<double>& BoundaryCondition::GetValues(void) const {
  return values_m;
}

bool BoundaryCondition::IsTimeVarying(void) const { return time_varying_m; }

}  // namespace chyps
