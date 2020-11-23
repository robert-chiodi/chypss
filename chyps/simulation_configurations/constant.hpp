// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef CHYPS_SIMULATION_CONFIGURATIONS_CONSTANT_HPP_
#define CHYPS_SIMULATION_CONFIGURATIONS_CONSTANT_HPP_

#include <functional>

#include <mfem/mfem.hpp>

#include "chyps/input_parser.hpp"

namespace chyps {

namespace constant {
void AddParserOptions(InputParser& a_parser);
void InitializeData(const nlohmann::json& a_json_object,
                    const InputParser& a_full_parser,
                    mfem::ParFiniteElementSpace& a_finite_element_space,
                    mfem::Vector& a_data);
}  // namespace constant

}  // namespace chyps

#endif  // CHYPS_SIMULATION_CONFIGURATIONS_CONSTANT_HPP_
