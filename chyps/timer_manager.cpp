// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chyps/timer_manager.hpp"

#include "chyps/debug_assert.hpp"

namespace chyps {

TimerManager::TimerTracker::TimerTracker(void)
    : active_m(false), time_m(), cummulative_time_m(0.0) {}

void TimerManager::TimerTracker::Start(void) {
  DEBUG_ASSERT(!this->Active(), global_assert{}, DebugLevel::CHEAP{},
               "Timer already active.");
  active_m = true;
  time_m = std::chrono::system_clock::now();
}
void TimerManager::TimerTracker::Stop(void) {
  DEBUG_ASSERT(this->Active(), global_assert{}, DebugLevel::CHEAP{},
               "Trying to stop timer that has not been started.");
  const std::chrono::duration<double> duration =
      std::chrono::system_clock::now() - time_m;
  cummulative_time_m += duration.count();
  active_m = false;
}
void TimerManager::TimerTracker::Reset(void) { cummulative_time_m = 0.0; }

bool TimerManager::TimerTracker::Active(void) const { return active_m; }

double TimerManager::TimerTracker::GetTime(void) const {
  return cummulative_time_m;
}

TimerManager::TimerManager(const MPIParallel* a_mpi_session)
    : mpi_session_m(a_mpi_session),
      timer_collection_m(),
      monitor_file_m(nullptr),
      total_exists_m(-1) {}

void TimerManager::AddTimer(const std::string& a_name) {
  DEBUG_ASSERT(!this->TimerExists(a_name), global_assert{}, DebugLevel::CHEAP{},
               "Timer with name \"" + a_name + "\" already added.");
  DEBUG_ASSERT(
      this->CanStillModify(), global_assert{}, DebugLevel::CHEAP{},
      "Timers cannot be added after a call to this->CreateMonitorFile(...)");
  timer_collection_m.emplace(a_name, TimerTracker());
}
void TimerManager::StartTimer(const std::string& a_name) {
  DEBUG_ASSERT(this->TimerExists(a_name), global_assert{}, DebugLevel::CHEAP{},
               "Timer with name \"" + a_name +
                   "\" has not been added to this TimerManager.");
  DEBUG_ASSERT(!this->TimerActive(a_name), global_assert{}, DebugLevel::CHEAP{},
               "Timer with name \"" + a_name + "\" already started.");
  timer_collection_m[a_name].Start();
}
void TimerManager::StopTimer(const std::string& a_name) {
  DEBUG_ASSERT(
      this->TimerActive(a_name), global_assert{}, DebugLevel::CHEAP{},
      "Tried to stop timer \"" + a_name + "\" that has not been started.");
  timer_collection_m[a_name].Stop();
}
void TimerManager::ResetTimer(const std::string& a_name) {
  DEBUG_ASSERT(this->TimerExists(a_name), global_assert{}, DebugLevel::CHEAP{},
               "Tried to reset timer \"" + a_name + "\" that does not exist.");
  timer_collection_m[a_name].Reset();
}
void TimerManager::ResetAllTimers(void) {
  for (auto& elem : timer_collection_m) {
    elem.second.Reset();
  }
}

std::size_t TimerManager::TimerActive(void) const {
  std::size_t active_timers = 0;
  for (const auto& elem : timer_collection_m) {
    active_timers += elem.second.Active() ? 1 : 0;
  }
  return active_timers;
}
bool TimerManager::TimerActive(const std::string& a_name) const {
  DEBUG_ASSERT(this->TimerExists(a_name), global_assert{}, DebugLevel::CHEAP{},
               "Timer with \"" + a_name + "\" does not exist.");
  return timer_collection_m.at(a_name).Active();
}

double TimerManager::GetTime(const std::string& a_name) const {
  DEBUG_ASSERT(this->TimerExists(a_name), global_assert{}, DebugLevel::CHEAP{},
               "Attempted to get time spent in timer \"" + a_name +
                   "\" but timer does not exist.");
  return timer_collection_m.at(a_name).GetTime();
}

void TimerManager::CreateMonitorFile(const std::string& a_file_name,
                                     MonitorManager& a_monitor_manager) {
  DEBUG_ASSERT(monitor_file_m == nullptr, global_assert{}, DebugLevel::CHEAP{},
               "TimerManager object already associated with a monitor file "
               "from a monitor manager.");

  if (!this->IsParallelTimer() || mpi_session_m->IAmRoot()) {
    std::vector<std::string> header(timer_collection_m.size());
    std::size_t ind = 0;
    for (const auto& elem : timer_collection_m) {
      if (elem.first == "Total" || elem.first == "total") {
        total_exists_m = static_cast<int>(ind);
      }
      header[ind++] = elem.first;
    }

    // Moving Total to the front
    if (total_exists_m > 0) {
      const std::string total_name = header[total_exists_m];
      for (int n = 1; n <= total_exists_m; ++n) {
        header[n] = header[n - 1];
      }
      header[0] = total_name;
    }

    monitor_file_m = a_monitor_manager.CreateMonitorFile(
        a_file_name, header,
        std::vector<FieldType>(header.size(), FieldType::DOUBLE));
  }
}
void TimerManager::PushTimesToMonitor(void) {
  std::vector<double> times(timer_collection_m.size());
  std::size_t ind = 0;
  for (const auto& elem : timer_collection_m) {
    times[ind++] = elem.second.GetTime();
  }

  if (this->IsParallelTimer()) {
    if (mpi_session_m->IAmNotRoot()) {
      MPI_Reduce(times.data(), nullptr, static_cast<int>(times.size()),
                 MPI_DOUBLE, MPI_MAX, 0, mpi_session_m->GetComm());
    } else {
      DEBUG_ASSERT(monitor_file_m != nullptr, global_assert{},
                   DebugLevel::CHEAP{},
                   "TimerManager object is not associated with a monitor file. "
                   "Need to call this->CreateMonitorFile(...) first.");
      std::vector<double> global_times(times.size());
      MPI_Reduce(times.data(), global_times.data(),
                 static_cast<int>(times.size()), MPI_DOUBLE, MPI_MAX, 0,
                 mpi_session_m->GetComm());

      if (total_exists_m > 0) {
        const double total_time = global_times[total_exists_m];
        for (int n = 1; n <= total_exists_m; ++n) {
          global_times[n] = global_times[n - 1];
        }
        global_times[0] = total_time;
      }
      monitor_file_m->SetEntries(global_times);
    }
  } else {
    DEBUG_ASSERT(monitor_file_m != nullptr, global_assert{},
                 DebugLevel::CHEAP{},
                 "TimerManager object is not associated with a monitor file. "
                 "Need to call this->CreateMonitorFile(...) first.");
    if (total_exists_m >= 0) {
      const double total_time = times[total_exists_m];
      for (int n = 1; n >= total_exists_m; ++n) {
        times[n] = times[n - 1];
      }
      times[0] = total_time;
    }
    monitor_file_m->SetEntries(times);
  }
}

bool TimerManager::TimerExists(const std::string& a_name) const {
  return timer_collection_m.find(a_name) != timer_collection_m.end();
}

bool TimerManager::CanStillModify(void) const {
  return monitor_file_m == nullptr;
}

bool TimerManager::IsParallelTimer(void) const {
  return mpi_session_m != nullptr;
}

}  // namespace chyps
