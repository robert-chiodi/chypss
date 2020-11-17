// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// \file solver_interface.hpp
/// \brief An abstract interface for a solver class specifying classes
/// for driving the solution.
///
/// For an example of a solver that implements this interface, see the class
/// HeatSolver.

#ifndef CHYPS_SOLVER_INTERFACE_HPP_
#define CHYPS_SOLVER_INTERFACE_HPP_

#include "chyps/input_parser.hpp"
#include "chyps/mesh.hpp"

namespace chyps {

/// \class SolverInterface solver_interface.hpp chyps/solver_interface.hpp
/// \brief Abstract interface class for solvers.
class SolverInterface {
 public:
  /// \brief Default constructor. Does nothing.
  SolverInterface(void) = default;

  /// \brief Perform all initialization needed by the solver.
  ///
  /// This method should initialize all necessary operators and
  /// variables that will be used by the solver. It will be called
  /// separately from advancement of the solution, and will only
  /// be called once.
  virtual void Initialize(void) = 0;

  /// \brief Return the a timestep <= a_proposed_dt that
  /// meets any requirements of the solver (e.g., for stability).
  virtual double AdjustTimeStep(const double a_proposed_dt) const = 0;

  /// \brief Advance the solution from time ``a_time`` by the time step
  /// ``a_dt``.
  virtual double Advance(const double a_time, const double a_dt) = 0;

  /// \brief Write out data to disk for visualization in VisIt.
  virtual void WriteFields(const int a_cycle, const double a_time) = 0;

  /// \brief Empty virtual destructor to be overriden by inherited class.
  virtual ~SolverInterface(void) = default;

 private:
};
}  // namespace chyps

#endif  // CHYPS_SOLVER_INTERFACE_HPP_
