// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef CHYPS_SIMULATION_CONFIGURATIONS_COOLED_ROD_HPP_
#define CHYPS_SIMULATION_CONFIGURATIONS_COOLED_ROD_HPP_

#include <functional>

#include "chyps/simulation_configurations/configuration_initializer_interface.hpp"

namespace chyps {

class CooledRod : public ConfigurationInitializer {
 public:
  CooledRod(InputParser& a_parser);

  virtual void Initialize(void) override final;
  virtual void FillRequiredData(RequiredData& a_data) override final;

  virtual ~CooledRod(void) override final;

 private:
};

}  // namespace chyps

#endif  // CHYPS_SIMULATION_CONFIGURATIONS_COOLED_ROD_HPP_
