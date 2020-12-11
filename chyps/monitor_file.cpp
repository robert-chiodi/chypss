// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chyps/monitor_file.hpp"

#include <utility>

#include "chyps/debug_assert.hpp"

namespace chyps {

MonitorFile::MonitorFile(const std::string& a_name,
                         std::initializer_list<std::string> a_header,
                         std::initializer_list<FieldType> a_format)
    : monitor_file_m(nullptr),
      updated_m(false),
      format_m(a_format),
      column_width_m(a_header.size()),
      data_m(a_header.size()) {
  DEBUG_ASSERT(
      a_header.size() == format_m.size(), global_assert{}, DebugLevel::CHEAP{},
      "Supplied header vector and format vector are of different length.");
  this->InitializeMonitorFile(a_name, a_header);
}

MonitorFile::MonitorFile(const std::string& a_name,
                         const std::vector<std::string>& a_header,
                         const std::vector<FieldType>& a_format)
    : monitor_file_m(nullptr),
      updated_m(false),
      format_m(a_format),
      column_width_m(a_header.size()),
      data_m(a_header.size()) {
  DEBUG_ASSERT(
      a_header.size() == a_format.size(), global_assert{}, DebugLevel::CHEAP{},
      "Supplied header vector and format vector are of different length.");
  this->InitializeMonitorFile(a_name, a_header);
}

MonitorFile::MonitorFile(MonitorFile&& a_other)
    : monitor_file_m(a_other.monitor_file_m),
      updated_m(a_other.updated_m),
      format_m(std::move(a_other.format_m)),
      column_width_m(std::move(a_other.column_width_m)),
      data_m(std::move(a_other.data_m)) {
  a_other.monitor_file_m = nullptr;
  a_other.updated_m = false;
}

MonitorFile& MonitorFile::operator=(MonitorFile&& a_other) {
  if (this != &a_other) {
    monitor_file_m = a_other.monitor_file_m;
    updated_m = a_other.updated_m;
    format_m = std::move(a_other.format_m);
    column_width_m = std::move(a_other.column_width_m);
    data_m = std::move(a_other.data_m);

    a_other.monitor_file_m = nullptr;
    a_other.updated_m = false;
  }
  return *this;
}

void MonitorFile::SetEntries(std::initializer_list<double> a_data) {
  DEBUG_ASSERT(a_data.size() == data_m.size(), global_assert{},
               DebugLevel::CHEAP{},
               "Supplied data and storage of different size.");
  std::copy(a_data.begin(), a_data.end(), data_m.begin());
  updated_m = true;
}

void MonitorFile::SetEntries(const std::vector<double>& a_data) {
  DEBUG_ASSERT(a_data.size() == data_m.size(), global_assert{},
               DebugLevel::CHEAP{},
               "Supplied data and storage of different size.");
  std::copy(a_data.begin(), a_data.end(), data_m.begin());
  updated_m = true;
}

void MonitorFile::SetEntry(const std::size_t a_index, const int a_value) {
  DEBUG_ASSERT(a_index < data_m.size(), global_assert{}, DebugLevel::CHEAP{});
  DEBUG_ASSERT(format_m[a_index] == FieldType::INT, global_assert{},
               DebugLevel::CHEAP{},
               "The column for a_index \"" + std::to_string(a_index) +
                   "\" was marked to be an integer.");
  data_m[a_index] = static_cast<double>(a_value);
  updated_m = true;
}

void MonitorFile::SetEntry(const std::size_t a_index, const double a_value) {
  DEBUG_ASSERT(a_index < data_m.size(), global_assert{}, DebugLevel::CHEAP{});
  DEBUG_ASSERT(format_m[a_index] == FieldType::DOUBLE, global_assert{},
               DebugLevel::CHEAP{},
               "The column for a_index \"" + std::to_string(a_index) +
                   "\" was marked to be a double.");
  data_m[a_index] = a_value;
  updated_m = true;
}

void MonitorFile::WriteLineToFile(const int a_iter, const double a_time,
                                  const double a_dt) {
  DEBUG_ASSERT(updated_m, global_assert{}, DebugLevel::CHEAP{},
               "No new data has been set since the last writing to file.");
  DEBUG_ASSERT(monitor_file_m != nullptr, global_assert{}, DebugLevel::CHEAP{});
  fprintf(monitor_file_m, "%-*d", MIN_WIDTH, a_iter);
  fprintf(monitor_file_m, "%-*.14E", MIN_WIDTH, a_time);
  fprintf(monitor_file_m, "%-*.14E", MIN_WIDTH, a_dt);
  for (std::size_t n = 0; n < format_m.size(); ++n) {
    switch (format_m[n]) {
      case FieldType::INT: {
        fprintf(monitor_file_m, "%-*d", column_width_m[n],
                static_cast<int>(data_m[n]));
        break;
      }
      case FieldType::DOUBLE: {
        fprintf(monitor_file_m, "%-*.14E", column_width_m[n], data_m[n]);
        break;
      }
      default:
        DEBUG_ASSERT(false, global_assert{}, DebugLevel::ALWAYS{},
                     "Unknown field type.");
        break;
    }
  }
  fprintf(monitor_file_m, "\n");
}

void MonitorFile::ClearLine(void) { updated_m = false; }

void MonitorFile::Flush(void) { fflush(monitor_file_m); }

MonitorFile::~MonitorFile(void) {
  if (monitor_file_m != nullptr) {
    fclose(monitor_file_m);
  }
  monitor_file_m = nullptr;
}

void MonitorFile::InitializeMonitorFile(
    const std::string& a_name, const std::vector<std::string>& a_header) {
  monitor_file_m = fopen(a_name.c_str(), "w");
  DEBUG_ASSERT(monitor_file_m != nullptr, global_assert{}, DebugLevel::CHEAP{},
               "Trouble opening file \"" + a_name + "\"");

  fprintf(monitor_file_m, "%-*s", MIN_WIDTH, "Iteration");
  fprintf(monitor_file_m, "%-*s", MIN_WIDTH, "Time [s]");
  fprintf(monitor_file_m, "%-*s", MIN_WIDTH, "dt [s]");
  auto cw = column_width_m.begin();
  for (const auto& string : a_header) {
    const int width = std::max(MIN_WIDTH, static_cast<int>(string.size() + 1));
    *cw++ = width;
    fprintf(monitor_file_m, "%-*s", width, string.c_str());
  }
  fprintf(monitor_file_m, "\n");
}

void MonitorFile::InitializeMonitorFile(
    const std::string& a_name, std::initializer_list<std::string> a_header) {
  monitor_file_m = fopen(a_name.c_str(), "w");
  DEBUG_ASSERT(monitor_file_m != nullptr, global_assert{}, DebugLevel::CHEAP{},
               "Trouble opening file \"" + a_name + "\"");

  fprintf(monitor_file_m, "%-*s", MIN_WIDTH, "Iteration");
  fprintf(monitor_file_m, "%-*s", MIN_WIDTH, "Time [s]");
  fprintf(monitor_file_m, "%-*s", MIN_WIDTH, "dt [s]");
  auto cw = column_width_m.begin();
  for (const auto& string : a_header) {
    const int width = std::max(MIN_WIDTH, static_cast<int>(string.size() + 1));
    *cw++ = width;
    fprintf(monitor_file_m, "%-*s", width, string.c_str());
  }
  fprintf(monitor_file_m, "\n");
}

}  // namespace chyps
