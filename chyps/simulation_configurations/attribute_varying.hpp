// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef CHYPS_SIMULATION_CONFIGURATIONS_ATTRIBUTE_VARYING_HPP_
#define CHYPS_SIMULATION_CONFIGURATIONS_ATTRIBUTE_VARYING_HPP_

#include <functional>

#include <mfem/mfem.hpp>

#include "chyps/input_parser.hpp"
#include "chyps/mesh.hpp"
#include "chyps/simulation_initializer.hpp"

namespace chyps {

namespace attribute_varying {

namespace scalar {
void AddParserOptions(InputParser& a_parser);
bool SupportedFieldType(const DataFieldType a_field_type);
void InitializeData(const DataFieldType a_field_type,
                    const nlohmann::json& a_json_object,
                    const InputParser& a_full_parser, const Mesh& a_mesh,
                    mfem::ParFiniteElementSpace& a_finite_element_space,
                    mfem::Vector& a_data);
}  // namespace scalar

namespace matrix {
void AddParserOptions(InputParser& a_parser);
bool SupportedFieldType(const DataFieldType a_field_type);
void InitializeData(const DataFieldType a_field_type,
                    const nlohmann::json& a_json_object,
                    const InputParser& a_full_parser, const Mesh& a_mesh,
                    mfem::ParFiniteElementSpace& a_finite_element_space,
                    const int a_number_of_rows, const int a_number_of_columns,
                    std::vector<mfem::DenseMatrix>& a_data);
}  // namespace matrix

}  // namespace attribute_varying
}  // namespace chyps

#endif  // CHYPS_SIMULATION_CONFIGURATIONS_ATTRIBUTE_VARYING_HPP_
