// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef TESTS_HELPER_COMMAND_LINE_INPUT_HPP_
#define TESTS_HELPER_COMMAND_LINE_INPUT_HPP_

#include <string>
#include <vector>

namespace chyps {

inline std::vector<char*> FakeCommandLineInput(
    const std::vector<std::string>& a_input) {
  std::vector<char*> vec;
  for (const auto& s : a_input) {
    char* c_str = new char[s.size() + 1];
    std::copy(s.begin(), s.end(), c_str);
    c_str[s.size()] = '\0';
    vec.push_back(c_str);
  }
  return vec;
}

inline void AddFakeFlag(const std::string a_flag, const std::string a_value,
                        std::vector<char*>& a_fake_cl) {
  char* c_str = new char[a_flag.size() + 1];
  std::copy(a_flag.begin(), a_flag.end(), c_str);
  c_str[a_flag.size()] = '\0';
  a_fake_cl.push_back(c_str);

  c_str = new char[a_value.size() + 1];
  std::copy(a_value.begin(), a_value.end(), c_str);
  c_str[a_value.size()] = '\0';
  a_fake_cl.push_back(c_str);
}

inline void DeleteCommandLineInput(std::vector<char*>& a_input) {
  for (auto& elem : a_input) {
    delete[] elem;
    elem = nullptr;
  }
}

inline bool StringEqual(const std::string& a_s1, const std::string& a_s2) {
  if (a_s1.size() != a_s2.size()) {
    return false;
  } else {
    for (std::size_t n = 0; n < a_s1.size(); ++n) {
      if (a_s1[n] != a_s2[n]) {
        return false;
      }
    }
  }
  return true;
}
}  // namespace chyps

#endif  // TESTS_HELPER_COMMAND_LINE_INPUT_HPP_
