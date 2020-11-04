// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef CHYPS_DEBUG_ASSERT_HPP_
#define CHYPS_DEBUG_ASSERT_HPP_

#include <mpi.h>
#include <iostream>

// This header contains the foonathan debug_assert code that we will be using.
#include <debug_assert/debug_assert.hpp>

#ifdef NDEBUG
#define CHYPS_DEBUG_ASSERT_LEVEL 1
#endif

#ifndef CHYPS_DEBUG_ASSERT_LEVEL
#define CHYPS_DEBUG_ASSERT_LEVEL 9
#endif

// NOTE: A level of 0 will turn off all assertions, including those used for
// fatal errors.

namespace chyps {

struct global_assert : debug_assert::set_level<CHYPS_DEBUG_ASSERT_LEVEL> {
  static void handle(const debug_assert::source_location& loc,
                     const char* expression,
                     std::string message = "") noexcept {
    std::cerr << "Assertion failure '" << loc.file_name << ':'
              << loc.line_number << ":\n";
    std::cerr << "Failing expression: " << expression;
    if (message != "") {
      std::cerr << "\nMessage: " << message;
    }
    std::cerr << '\n';
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
};

namespace DebugLevel {
// Largest amount of debugging. Used for full debug. Can include expensive
// checks.
using FULL = debug_assert::level<9>;

// Cheap debugging checks, such as checking bounds on an object or index. A
// check that is not expected to greatly increase computational cost. Allows
// debug runs to be done more efficiently if full option not needed.
using CHEAP = debug_assert::level<5>;

// To never be turned off. Used as critical or fatal failures.
using ALWAYS = debug_assert::level<1>;
}  // namespace DebugLevel

}  // namespace chyps

#endif  // CHYPS_DEBUG_ASSERT_HPP_
