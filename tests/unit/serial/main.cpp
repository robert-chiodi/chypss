// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>

#define CHYPS_LOGGER_OFF
#include "chyps/logger.hpp"

using namespace chyps;

// Custom main function needed to turn off Spdlog.
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  StartLogger(0, 1, SpdlogLevel::OFF);
  return RUN_ALL_TESTS();
}
