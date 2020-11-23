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
///
/// NOTE: If left uninitialized (by not calling the mehtod
/// MonitorManager::Initialize()), the monitor manager can be considered in an
/// off state. In this state,  `nullptr` will be handed back from the method
/// CreateMonitorFile(...). Ensure that anything that attempts to write to the
/// file checks for this by either checking for nullptr or checking if the
/// monitor manager is off.
class MonitorManager {
  static constexpr int IO_FLUSH_FREQUENCY = 20;

 public:
  /// \brief Constructor.
  MonitorManager(void);

  /// \brief Creates the directory for monitor files and turns the
  /// MonitorManager to the `on` state. If not called, MonitorManager will be
  /// considered off and return `nullptr` from calls to CreateMonitorFile().
  ///
  /// The directory name will be `a_directory_name`_000, HOWEVER, if this
  /// directory already exists, an integer will be appended to keep this
  /// directory unique and avoid overwriting the previous directory and its
  /// files.
  void Initialize(const std::string& a_directory_name);

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

  /// \brief Returns whether MonitorFileManager is on (has been initialized).
  bool IsOn(void) const;

  ///\brief Default destructor. Will call destructor for all MonitorFile
  /// objects.
  ~MonitorManager(void) = default;

 private:
  std::string directory_base_m;
  std::unordered_map<std::string, MonitorFile> file_list_m;
  bool on_m;
};

}  // namespace chyps

#endif  // CHYPS_MONITOR_MANAGER_HPP_
