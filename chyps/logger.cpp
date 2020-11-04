// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chyps/logger.hpp"

#include <iostream>
#include <string>

#include "chyps/mpi_parallel.hpp"

namespace chyps {

std::shared_ptr<spdlog::logger> MAIN_LOG;

void StartLogger(const MPIParallel& a_mpi_session, const SpdlogLevel a_level) {
  StartLogger(a_mpi_session.MyRank(), a_mpi_session.NumberOfRanks(), a_level);
}

void StartLogger(const int a_rank, const int a_world_size,
                 const SpdlogLevel a_level) {
#if CHYPS_LOGGER_OFF
  MAIN_LOG = nullptr;
  spdlog::set_level(spdlog::level::off);
  return;
#endif
  if (spdlog::get("MAIN_LOG") == nullptr) {
    // Log has yet to be created.
    if (a_world_size > 32) {
      if (a_rank == 0) {
        std::cout << "Too many processes " << a_world_size
                  << " in use for logger" << std::endl;
        std::cout << "Logger being DISABLED." << std::endl;
        std::cout << "For optimal performance, recompile with CHYPS_LOGGER_OFF "
                     "defined."
                  << std::endl;
        std::cout
            << "NOTE: This is not a steadfast requirement, however using "
               "the logger with many processes will lead to poor "
               "performance and excessive I/O use. Reconsider if necessary "
               "to log for this simulation."
            << std::endl;
      }
      spdlog::set_level(spdlog::level::off);
    }

    int world_digits = 0;
    int ws = a_world_size;
    do {
      ++world_digits;
      ws /= 10;
    } while (ws / 10 > 0);

    std::string number = std::to_string(a_rank);
    std::string id_suffix =
        std::string(world_digits - number.length(), '0') + number;
    std::string log_name =
        std::string("log_") + id_suffix + std::string(".txt");
    try {
      MAIN_LOG = spdlog::basic_logger_st("MAIN_LOG", "logs/" + log_name, true);
    } catch (const spdlog::spdlog_ex& ex) {
      std::cout << "Log initialization failed: " << ex.what() << std::endl;
    }
  }

  switch (a_level) {
    case SpdlogLevel::TRACE: {
      spdlog::set_level(spdlog::level::trace);
      break;
    }
    case SpdlogLevel::DEBUG: {
      spdlog::set_level(spdlog::level::debug);
      break;
    }
    case SpdlogLevel::INFO: {
      spdlog::set_level(spdlog::level::info);
      break;
    }
    case SpdlogLevel::WARN: {
      spdlog::set_level(spdlog::level::warn);
      break;
    }
    case SpdlogLevel::ERROR: {
      spdlog::set_level(spdlog::level::err);
      break;
    }
    case SpdlogLevel::CRITICAL: {
      spdlog::set_level(spdlog::level::critical);
      break;
    }
    case SpdlogLevel::OFF: {
      spdlog::set_level(spdlog::level::off);
      break;
    }
  }
}

}  // namespace chyps
