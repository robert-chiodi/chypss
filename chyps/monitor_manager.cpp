// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chyps/monitor_manager.hpp"

#include <filesystem>
#include <utility>

#include "chyps/debug_assert.hpp"
#include "chyps/string_manipulation.hpp"

namespace chyps {

MonitorManager::MonitorManager(void) {
  static constexpr uint32_t max_directory_attempt = 1000;

  std::filesystem::path dir_path;
  uint32_t n = 0;
  for (n = 0; n < max_directory_attempt; ++n) {
    directory_base_m = "./monitor_" + ZeroFill(n, 3);
    dir_path = std::filesystem::path(directory_base_m);
    if (!std::filesystem::is_directory(dir_path)) {
      break;
    }
  }
  DEBUG_ASSERT(n < max_directory_attempt, global_assert{}, DebugLevel::ALWAYS{},
               "Cannot create monitor directory. More than " +
                   std::to_string(max_directory_attempt) +
                   " directories in existence.");

  std::filesystem::create_directory(dir_path);
}

MonitorFile* MonitorManager::CreateMonitorFile(
    const std::string& a_name, const std::vector<std::string>& a_header,
    const std::vector<FieldType>& a_format) {
  auto format_copy = a_format;
  return CreateMonitorFile(a_name, a_header, std::move(format_copy));
}

MonitorFile* MonitorManager::CreateMonitorFile(
    const std::string& a_name, const std::vector<std::string>& a_header,
    std::vector<FieldType>&& a_format) {
  auto insert_location =
      file_list_m.insert({a_name, MonitorFile(directory_base_m + '/' + a_name,
                                              a_header, std::move(a_format))});
  DEBUG_ASSERT(insert_location.second, global_assert{}, DebugLevel::CHEAP{},
               "Monitor file with name \"" + a_name + "\" already added.");
  return &(insert_location.first->second);
}

void MonitorManager::WriteStepToFiles(const int a_iter, const double a_time,
                                      const double a_dt) {
  for (auto& elem : file_list_m) {
    elem.second.WriteLineToFile(a_iter, a_time, a_dt);
    elem.second.ClearLine();
  }
  if (a_iter % IO_FLUSH_FREQUENCY == 0) {
    for (auto& elem : file_list_m) {
      elem.second.Flush();
    }
  }
}

}  // namespace chyps
