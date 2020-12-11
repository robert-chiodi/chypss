// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef CHYPS_TIMER_MANAGER_HPP_
#define CHYPS_TIMER_MANAGER_HPP_

#include <chrono>
#include <unordered_map>

#include <mfem/mfem.hpp>

#include "chyps/monitor_file.hpp"
#include "chyps/monitor_manager.hpp"
#include "chyps/mpi_parallel.hpp"

namespace chyps {

/// \class TimerManager chyps/timer_manager.hpp timer_manager.hpp
/// \brief A class to create and manage a collection of timers. These
/// timers can be tied to a MonitorFile object to be written to a file.
///
/// The individual timers are referenced through the TimerManager
/// with a unique name. Timers can only be added before a
/// call to this->CreateMonitorFile(), after which no
/// additional timers can be created for this TimerManagerObject.
/// If a timer with name "total" or "Total" exists, it will
/// be moved to the first timer in the monitor.
///
/// NOTE: If using TimerManager in parallel (with a supplied MPIParallel object
/// during construction), all ranks must add the same timers. The reported time
/// will be the maximum time for that timer across all ranks.
class TimerManager {
  class TimerTracker {
   public:
    TimerTracker(void);
    TimerTracker(const std::size_t a_location);

    void Start(void);
    double Stop(void);  // Returns time in seconds since start
    bool Active(void) const;
    std::size_t GetTimeLocation(void) const;
    void SetTimeLocation(const std::size_t a_location);

   private:
    bool active_m;
    std::chrono::time_point<std::chrono::system_clock> time_m;
    std::size_t cumulative_time_location_m;
  };

 public:
  /// \brief Constructor that initializes an empty TimerManager. If using in
  /// parallel, a valid MPIParallel object must be supplied to perform an
  /// MPI_Reduce prior to pushing times to the monitor file.
  TimerManager(const MPIParallel* a_mpi_session);

  /// \brief Add a timer to be tracked by this manager. The timer will be
  /// interacted with through a_name. The name `a_name` must be unique to
  /// this TimerManager object. All timers must be added prior to creating
  /// a monitor file.
  void AddTimer(const std::string& a_name);

  /// \brief Start the timer with a_name.
  void StartTimer(const std::string& a_name);

  /// \brief Stop the timer with a_name. The timer must be active (been started
  /// with StartTime).
  void StopTimer(const std::string& a_name);

  /// \brief Reset any stored time for a_name to 0.
  void ResetTimer(const std::string& a_name);

  /// \brief Reset times for all timers to 0.
  void ResetAllTimers(void);

  /// \brief Returns the number of timers active.
  ///
  /// NOTE: This requires iterating over all timers and therefore
  /// could be expensive.
  std::size_t TimerActive(void) const;

  /// \brief Returns true if the timer with a_name is active.
  bool TimerActive(const std::string& a_name) const;

  /// \brief Returns the amount of time accumulated while this timer
  /// was active.
  ///
  /// NOTE: If this timer is currently active, the time elapsed from
  /// StartTimer() is NOT accounted for. Only time between matching StartTimer()
  /// and StopTimer() calls since the last reset are counted.
  double GetTime(const std::string& a_name) const;

  /// \brief Create a monitor file to write the timer values to. After creating
  /// the monitor file, no additional timers can be added to this TimerManager.
  void CreateMonitorFile(const std::string& a_file_name,
                         MonitorManager& a_monitor_manager);

  /// \brief Write the current stored times to the associated monitor file
  /// obtained by this->CreateMonitorFile(...). This call must come after
  /// this->CreateMonitorFile().
  ///
  /// Up until this, all ranks will store their own times for the respective
  /// timers. Before pushing times to the monitor (by the root rank), we will
  /// take the max values of the time across the ranks to report in the monitor.
  ///
  /// NOTE: This function features a collective call to all ranks in
  /// mpi_session_m, therefore all MPI ranks must call this function.
  void PushTimesToMonitor(void);

  /// \brief Returns whether there is a MPIParallel object associated
  /// with the timer.
  bool IsParallelTimer(void) const;

  /// \brief Default destructor.
  ~TimerManager(void) = default;

 private:
  bool TimerExists(const std::string& a_name) const;
  bool CanStillModify(void) const;
  bool MonitorFileOff(void) const;

  const MPIParallel* mpi_session_m;
  std::unordered_map<std::string, TimerTracker> timer_collection_m;
  std::vector<double> cumulative_time_m;
  MonitorFile* monitor_file_m;
  bool monitor_file_off_m;
};

}  // namespace chyps

#endif  // CHYPS_TIMER_MANAGER_HPP_
