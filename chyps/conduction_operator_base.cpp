// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chyps/conduction_operator_base.hpp"

namespace chyps {

ConductionOperatorBase::ConductionOperatorBase(mfem::ParFiniteElementSpace& f)
    : mfem::TimeDependentOperator(f.GetTrueVSize(), 0.0), fespace(f) {}

}  // namespace chyps
