// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef CHYPS_MONITOR_FILE_HPP_
#define CHYPS_MONITOR_FILE_HPP_

#include <stdio.h>
#include <initializer_list>
#include <string>
#include <vector>

namespace chyps {

/// \enum FieldType chyps/monitor_file.hpp monitor_file.hpp
/// \brief Enum used to indicate the type for a column in a MonitorFile
enum class FieldType { INT = 0, DOUBLE };

/// \class MonitorFile chyps/monitor_file.hpp monitor_file.hpp
/// \breif Class to handle logic controlling a single monitor file.
class MonitorFile {
  static constexpr int MIN_WIDTH = 25;

 public:
  MonitorFile(void) = default;

  /// \brief Open the file with a_name and setup the header. Store the format to
  /// ensure proper writing of data later.
  MonitorFile(const std::string& a_name,
              const std::vector<std::string>& a_header,
              const std::vector<FieldType>& a_format);

  /// \brief Open the file with a_name and setup the header. Store the format to
  /// ensure proper writing of data later.
  MonitorFile(const std::string& a_name,
              std::initializer_list<std::string> a_header,
              std::initializer_list<FieldType> a_format);

  /// \brief Remove opy operator.
  MonitorFile(const MonitorFile& a_other) = delete;

  /// \brief Move operator
  MonitorFile(MonitorFile&& a_other);

  /// \brief Remove copy assignment.
  MonitorFile& operator=(const MonitorFile& a_other) = delete;

  /// \brief Move assignment.
  MonitorFile& operator=(MonitorFile&& a_other);

  /// \brief Set all entries that will be written to file. Data corresponding
  /// to columns that are formatted as an integer will be cast to int.
  void SetEntries(std::initializer_list<double> a_data);

  /// \brief Set all entries that will be written to file. Data corresponding
  /// to columns that are formatted as an integer will be cast to int.
  void SetEntries(const std::vector<double>& a_data);

  /// \brief Set entry at `a_index` to the int `a_value`. The column
  /// corresponding to `a_index` must have been formatted to be an int
  /// during construction.
  void SetEntry(const std::size_t a_index, const int a_value);

  /// \brief Set entry at `a_index` to the double `a_value`. The column
  /// corresponding to `a_index` must have been formatted to be a double
  /// during construction.
  void SetEntry(const std::size_t a_index, const double a_value);

  /// \brief Write the stored data for this iteration to disk.
  void WriteLineToFile(const int a_iter, const double a_time,
                       const double a_dt);

  /// \brief Clear the stored data in data_m to be updated again.
  void ClearLine(void);

  /// \brief Flush the file associated with this MonitorFile object.
  void Flush(void);

  /// \brief Close the file that is open.
  ~MonitorFile(void);

 private:
  void InitializeMonitorFile(const std::string& a_name,
                             const std::vector<std::string>& a_header);
  void InitializeMonitorFile(const std::string& a_name,
                             std::initializer_list<std::string> a_header);

  FILE* monitor_file_m;
  bool updated_m;
  std::vector<FieldType> format_m;
  std::vector<int> column_width_m;
  std::vector<double> data_m;
};

}  // namespace chyps

#endif  // CHYPS_MONITOR_FILE_HPP_
