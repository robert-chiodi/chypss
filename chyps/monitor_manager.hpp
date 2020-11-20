// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef CHYPS_MONITOR_MANAGER_HPP_
#define CHYPS_MONITOR_MANAGER_HPP_

#include <string>
#include <unordered_map>
#include <vector>

#include "chyps/monitor_file.hpp"

namespace chyps {

/// \class MonitorManager chyps/monitor_manager.hpp monitor_manager.hpp
/// \brief Class that manages a collection of files that monitor simulation
/// progress and results, written to a monitor directory.
///
/// This class stores a collection of MonitorFiles and will
/// drive their writing to disk and cleaning between iterations.
/// The storage of data itself should be handled by individual classes.
class MonitorManager {
  static constexpr int IO_FLUSH_FREQUENCY = 20;

 public:
  /// \brief Constructor that will determine the directory to use.
  ///
  /// The directory name will be monitor, HOWEVER, if a monitor directory
  /// already exists, an integer will be appended to keep this directory
  /// unique and avoid overwriting the previous directory and its files.
  MonitorManager(void);

  /// \brief Create a MonitorFile with the given name, header, and format field.
  ///
  /// `a_name` should be without the prefix "monitor", which will be
  /// automatically appended.
  ///
  /// `a_header` is a list of strings describing what  each column is. This is
  /// in addition to assumed (first three) entries of iteration#, time, and dt.
  ///
  /// `a_format` is a format specifier, marking a column to either use integer
  /// or double numbers.
  MonitorFile* CreateMonitorFile(const std::string& a_name,
                                 const std::vector<std::string>& a_header,
                                 const std::vector<FieldType>& a_format);

  /// \brief Create a MonitorFile with the given name, header, and format field.
  ///
  /// `a_name` should be without the prefix "monitor", which will be
  /// automatically appended.
  ///
  /// `a_header` is a list of strings describing what  each column is. This is
  /// in addition to assumed (first three) entries of iteration#, time, and dt.
  ///
  /// `a_format` is a format specifier, marking a column to either use integer
  /// or double numbers.
  MonitorFile* CreateMonitorFile(const std::string& a_name,
                                 const std::vector<std::string>& a_header,
                                 std::vector<FieldType>&& a_format);

  /// \brief Write the current MonitorFile::line_m stored in each
  /// MonitorFile to disk.
  void WriteStepToFiles(const int a_iter, const double a_time,
                        const double a_dt);

  ///\brief Default destructor. Will call destructor for all MonitorFile
  /// objects.
  ~MonitorManager(void) = default;

 private:
  std::string directory_base_m;
  std::unordered_map<std::string, MonitorFile> file_list_m;
};

}  // namespace chyps

#endif  // CHYPS_MONITOR_MANAGER_HPP_
