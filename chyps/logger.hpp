// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef CHYPS_LOGGER_HPP_
#define CHYPS_LOGGER_HPP_

#ifdef CHYPS_LOGGER_OFF
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_OFF
#else
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG
#endif

#include <memory>

#include <spdlog/spdlog.h>
#include "spdlog/sinks/basic_file_sink.h"  // support for basic file logging

/// \file logger.hpp
/// \brief File to include to use spdlog logger. Can be turned off via
/// compile time option, which adjusts the SPDLOG_ACTIVE_LEVEL
///
/// Note: Always used this when including the logger since
/// it will ensure that the macro SPDLOG_ACTIVE_LEVEL is
/// define prior to inclusion of spdlog (and therefore
/// make sure we respect the compile-time chosen level of logging.

namespace chyps {

/// \enum SpdlogLevel logger.hpp chyps/logger.hpp
/// \briefEnumerator for choosing level of spdlog output at runtime.
///
/// Note: If CHYPS_LOGGER_OFF is defined, the the debugger will always
/// be set to OFF.
enum class SpdlogLevel { TRACE = 0, DEBUG, INFO, WARN, ERROR, CRITICAL, OFF };

/// \brief Main logger to be used for all logging.
extern std::shared_ptr<spdlog::logger> MAIN_LOG;

/// \function StartLogger logger.hpp chyps/logger.hpp
/// \brief Start global logging object, with one file per
/// mpi_rank.
void StartLogger(const int a_rank_id, const int a_world_size,
                 const SpdlogLevel a_level = SpdlogLevel::INFO);

}  // namespace chyps

#endif  // CHYPS_LOGGER_HPP_
