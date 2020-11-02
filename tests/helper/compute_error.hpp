// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef TESTS_HELPER_COMPUTE_ERROR_HPP_
#define TESTS_HELPER_COMPUTE_ERROR_HPP_

#include <cmath>
#include <vector>

#include "chyps/mpi_parallel.hpp"

namespace chyps {

inline double GlobalMaxValue(const std::vector<double>& a_data,
                             const MPIParallel& a_mpi_session) {
  double value = -DBL_MAX;
  for (const auto& elem : a_data) {
    value = std::max(value, elem);
  }
  double max_value;
  MPI_Allreduce(&value, &max_value, 1, MPI_DOUBLE, MPI_MAX,
                a_mpi_session.GetComm());
  return max_value;
}

inline double GlobalL1Diff(const std::vector<double>& a_data_1,
                           const std::vector<double>& a_data_2,
                           const MPIParallel& a_mpi_session) {
  assert(a_data_1.size() == a_data_2.size());
  double diff = 0.0;
  for (std::size_t n = 0; n < a_data_1.size(); ++n) {
    diff += std::fabs(a_data_1[n] - a_data_2[n]);
  }
  double summed_diff;
  MPI_Allreduce(&diff, &summed_diff, 1, MPI_DOUBLE, MPI_SUM,
                a_mpi_session.GetComm());
  return summed_diff;
}

inline double GlobalL2Diff(const std::vector<double>& a_data_1,
                           const std::vector<double>& a_data_2,
                           const MPIParallel& a_mpi_session) {
  assert(a_data_1.size() == a_data_2.size());
  double diff = 0.0;
  for (std::size_t n = 0; n < a_data_1.size(); ++n) {
    diff += (a_data_1[n] - a_data_2[n]) * (a_data_1[n] - a_data_2[n]);
  }
  double summed_diff;
  MPI_Allreduce(&diff, &summed_diff, 1, MPI_DOUBLE, MPI_SUM,
                a_mpi_session.GetComm());
  summed_diff = std::sqrt(summed_diff);
  return summed_diff;
}

inline double GlobalLinfDiff(const std::vector<double>& a_data_1,
                             const std::vector<double>& a_data_2,
                             const MPIParallel& a_mpi_session) {
  assert(a_data_1.size() == a_data_2.size());
  double diff = -DBL_MAX;
  for (std::size_t n = 0; n < a_data_1.size(); ++n) {
    diff = std::max(diff, std::fabs(a_data_1[n] - a_data_2[n]));
  }
  double max_diff;
  MPI_Allreduce(&diff, &max_diff, 1, MPI_DOUBLE, MPI_MAX,
                a_mpi_session.GetComm());
  return max_diff;
}

inline double GlobalL2Diff_Normalized(const std::vector<double>& a_data_1,
                                      const std::vector<double>& a_data_2,
                                      const MPIParallel& a_mpi_session) {
  assert(a_data_1.size() == a_data_2.size());
  auto summed_diff = GlobalL2Diff(a_data_1, a_data_2, a_mpi_session);
  std::size_t local_count = a_data_1.size();
  std::size_t global_count;
  MPI_Allreduce(&local_count, &global_count, 1, my_MPI_SIZE_T, MPI_SUM,
                a_mpi_session.GetComm());
  return summed_diff / static_cast<double>(global_count);
}

inline double GlobalL1Diff_Normalized(const std::vector<double>& a_data_1,
                                      const std::vector<double>& a_data_2,
                                      const MPIParallel& a_mpi_session) {
  assert(a_data_1.size() == a_data_2.size());
  auto summed_diff = GlobalL1Diff(a_data_1, a_data_2, a_mpi_session);
  std::size_t local_count = a_data_1.size();
  std::size_t global_count;
  MPI_Allreduce(&local_count, &global_count, 1, my_MPI_SIZE_T, MPI_SUM,
                a_mpi_session.GetComm());
  return summed_diff / static_cast<double>(global_count);
}

}  // namespace chyps

#endif  // TESTS_HELPER_COMPUTE_ERROR_HPP_
