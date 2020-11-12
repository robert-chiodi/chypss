// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef CHYPS_SIMULATION_CONFIGURATIONS_CONFIGURATION_INITIALIZER_INTERFACE_HPP_
#define CHYPS_SIMULATION_CONFIGURATIONS_CONFIGURATION_INITIALIZER_INTERFACE_HPP_

#include "chyps/input_parser.hpp"
#include "chyps/simulation_initializer.hpp"

namespace chyps {

class ConfigurationInitializer {
 public:
  ConfigurationInitializer(InputParser& a_parser) : parser_m(a_parser) {}

  virtual void Initialize(void) = 0;
  virtual void FillRequiredData(RequiredData& a_data) = 0;

  virtual ~ConfigurationInitializer(void) = default;

 protected:
  InputParser& parser_m;
};

}  // namespace chyps

#endif  // CHYPS_SIMULATION_CONFIGURATIONS_CONFIGURATION_INITIALIZER_INTERFACE_HPP_
