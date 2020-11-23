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
    : active_m(false),
      time_m(),
      cumulative_time_location_m(static_cast<std::size_t>(-1)) {}

TimerManager::TimerTracker::TimerTracker(const std::size_t a_location)
    : active_m(false), time_m(), cumulative_time_location_m(a_location) {}

void TimerManager::TimerTracker::Start(void) {
  DEBUG_ASSERT(!this->Active(), global_assert{}, DebugLevel::CHEAP{},
               "Timer already active.");
  active_m = true;
  time_m = std::chrono::system_clock::now();
}
double TimerManager::TimerTracker::Stop(void) {
  DEBUG_ASSERT(this->Active(), global_assert{}, DebugLevel::CHEAP{},
               "Trying to stop timer that has not been started.");
  const std::chrono::duration<double> duration =
      std::chrono::system_clock::now() - time_m;
  active_m = false;
  return duration.count();
}

bool TimerManager::TimerTracker::Active(void) const { return active_m; }

std::size_t TimerManager::TimerTracker::GetTimeLocation(void) const {
  return cumulative_time_location_m;
}

void TimerManager::TimerTracker::SetTimeLocation(const std::size_t a_location) {
  cumulative_time_location_m = a_location;
}

TimerManager::TimerManager(const MPIParallel* a_mpi_session)
    : mpi_session_m(a_mpi_session),
      timer_collection_m(),
      cumulative_time_m(),
      monitor_file_m(nullptr),
      monitor_file_off_m(true) {}

void TimerManager::AddTimer(const std::string& a_name) {
  DEBUG_ASSERT(!this->TimerExists(a_name), global_assert{}, DebugLevel::CHEAP{},
               "Timer with name \"" + a_name + "\" already added.");
  DEBUG_ASSERT(
      this->CanStillModify(), global_assert{}, DebugLevel::CHEAP{},
      "Timers cannot be added after a call to this->CreateMonitorFile(...)");

  cumulative_time_m.push_back(0.0);
  if (a_name == "Total" || a_name == "total") {
    std::string alternative_name = a_name == "Total" ? "total" : "Total";
    DEBUG_ASSERT(
        !this->TimerExists(alternative_name), global_assert{},
        DebugLevel::CHEAP{},
        "A timer named \"Total\" and \"total\" cannot exist at the same time.");

    // Make total the first entry, so increment all others.
    for (auto& elem : timer_collection_m) {
      elem.second.SetTimeLocation(elem.second.GetTimeLocation() + 1);
    }
    timer_collection_m.emplace(a_name, TimerTracker(0));
  } else {
    timer_collection_m.emplace(
        a_name,
        TimerTracker(static_cast<std::size_t>(cumulative_time_m.size() - 1)));
  }
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
  auto& timer = timer_collection_m[a_name];
  DEBUG_ASSERT(timer.GetTimeLocation() < cumulative_time_m.size(),
               global_assert{}, DebugLevel::CHEAP{});
  cumulative_time_m[timer.GetTimeLocation()] += timer.Stop();
}
void TimerManager::ResetTimer(const std::string& a_name) {
  DEBUG_ASSERT(this->TimerExists(a_name), global_assert{}, DebugLevel::CHEAP{},
               "Tried to reset timer \"" + a_name + "\" that does not exist.");
  DEBUG_ASSERT(
      timer_collection_m[a_name].GetTimeLocation() < cumulative_time_m.size(),
      global_assert{}, DebugLevel::CHEAP{});
  cumulative_time_m[timer_collection_m[a_name].GetTimeLocation()] = 0.0;
}
void TimerManager::ResetAllTimers(void) {
  std::fill(cumulative_time_m.begin(), cumulative_time_m.end(), 0.0);
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
  DEBUG_ASSERT(timer_collection_m.at(a_name).GetTimeLocation() <
                   cumulative_time_m.size(),
               global_assert{}, DebugLevel::CHEAP{});
  return cumulative_time_m[timer_collection_m.at(a_name).GetTimeLocation()];
}

void TimerManager::CreateMonitorFile(const std::string& a_file_name,
                                     MonitorManager& a_monitor_manager) {
  DEBUG_ASSERT(monitor_file_m == nullptr, global_assert{}, DebugLevel::CHEAP{},
               "TimerManager object already associated with a monitor file "
               "from a monitor manager.");

  if (this->IsParallelTimer()) {
    if (mpi_session_m->IAmRoot()) {
      monitor_file_off_m = a_monitor_manager.IsOn() ? false : true;
    }
    int bool_as_int = monitor_file_off_m ? 1 : 0;
    MPI_Bcast(&bool_as_int, 1, MPI_INT, mpi_session_m->GetRootRank(),
              mpi_session_m->GetComm());
    monitor_file_off_m = bool_as_int == 1;
  } else {
    monitor_file_off_m = a_monitor_manager.IsOn() ? false : true;
  }

  if (this->MonitorFileOff()) {
    return;
  }

  if (!this->IsParallelTimer() || mpi_session_m->IAmRoot()) {
    std::vector<std::string> header(timer_collection_m.size());
    for (const auto& elem : timer_collection_m) {
      DEBUG_ASSERT(elem.second.GetTimeLocation() < header.size(),
                   global_assert{}, DebugLevel::CHEAP{});
      header[elem.second.GetTimeLocation()] = elem.first;
    }
    monitor_file_m = a_monitor_manager.CreateMonitorFile(
        a_file_name, header,
        std::vector<FieldType>(header.size(), FieldType::DOUBLE));
  }
}
void TimerManager::PushTimesToMonitor(void) {
  if (this->MonitorFileOff()) {
    // MonitorManager off.
    return;
  }

  if (this->IsParallelTimer()) {
    if (mpi_session_m->IAmNotRoot()) {
      MPI_Reduce(cumulative_time_m.data(), nullptr,
                 static_cast<int>(cumulative_time_m.size()), MPI_DOUBLE,
                 MPI_MAX, mpi_session_m->GetRootRank(),
                 mpi_session_m->GetComm());
    } else {
      DEBUG_ASSERT(monitor_file_m != nullptr, global_assert{},
                   DebugLevel::CHEAP{},
                   "TimerManager object is not associated with a monitor file. "
                   "Need to call this->CreateMonitorFile(...) first.");
      std::vector<double> global_times(cumulative_time_m.size());
      MPI_Reduce(cumulative_time_m.data(), global_times.data(),
                 static_cast<int>(cumulative_time_m.size()), MPI_DOUBLE,
                 MPI_MAX, mpi_session_m->GetRootRank(),
                 mpi_session_m->GetComm());
      monitor_file_m->SetEntries(global_times);
    }
  } else {
    DEBUG_ASSERT(monitor_file_m != nullptr, global_assert{},
                 DebugLevel::CHEAP{},
                 "TimerManager object is not associated with a monitor file. "
                 "Need to call this->CreateMonitorFile(...) first.");
    monitor_file_m->SetEntries(cumulative_time_m);
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

bool TimerManager::MonitorFileOff(void) const { return monitor_file_off_m; }

}  // namespace chyps
