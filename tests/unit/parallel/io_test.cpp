// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chyps/io.hpp"

#include <array>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include "tests/unit/parallel/mpi_session.hpp"

#include <fstream>

#include "chyps/boundary_condition_manager.hpp"
#include "chyps/input_parser.hpp"
#include "chyps/mesh.hpp"
#include "chyps/mpi_parallel.hpp"

namespace {

using namespace chyps;

TEST(IO, OpenAndCloseWrite) {
  std::string write_file_name = "tests/unit/parallel/data/io_test";

  IO file_io(*mpi_session, "HandleName");
  file_io.SetWrite(write_file_name);
}

TEST(IO, IsWriteMode) {
  std::string write_file_name = "tests/unit/parallel/data/io_test";
  IO file_io(*mpi_session, "HandleName");
  file_io.SetWrite(write_file_name);
  EXPECT_TRUE(file_io.IsWriteModeActive());
  EXPECT_FALSE(file_io.IsReadModeActive());
}

TEST(IO, IsReadMode) {
  std::string write_file_name = "tests/unit/parallel/data/io_test";
  IO file_io(*mpi_session, "HandleName");
  file_io.SetWrite(write_file_name);
  file_io.CloseWriteEngine();
  file_io.SetRead(write_file_name);
  EXPECT_FALSE(file_io.IsWriteModeActive());
  EXPECT_TRUE(file_io.IsReadModeActive());
  file_io.CloseReadEngine();
  EXPECT_FALSE(file_io.IsWriteModeActive());
  EXPECT_FALSE(file_io.IsReadModeActive());
}

TEST(IO, AttributeInt) {
  std::string write_file_name = "tests/unit/parallel/data/io_test";
  IO file_io(*mpi_session, "HandleName");
  file_io.SetWrite(write_file_name);
  std::vector<int> attribute{{-1, 0, 1, 2, 18}};
  file_io.WriteAttribute("Test int", attribute);
  file_io.CloseWriteEngine();

  file_io.SetRead(write_file_name);
  std::vector<int> read_in;
  bool found = file_io.ReadAttribute("Test int", read_in);
  EXPECT_TRUE(found);
  for (std::size_t n = 0; n < attribute.size(); ++n) {
    EXPECT_EQ(attribute[n], read_in[n]);
  }
}

TEST(IO, AttributeDouble) {
  std::string write_file_name = "tests/unit/parallel/data/io_test";
  IO file_io(*mpi_session, "HandleName");
  file_io.SetWrite(write_file_name);
  std::vector<double> attribute{{-1.32, 0.968, 3.94, -15.96, 18.42}};
  file_io.WriteAttribute("Test double", attribute);
  file_io.CloseWriteEngine();

  file_io.SetRead(write_file_name);
  std::vector<double> read_in;
  bool found = file_io.ReadAttribute("Test double", read_in);
  EXPECT_TRUE(found);
  for (std::size_t n = 0; n < attribute.size(); ++n) {
    EXPECT_DOUBLE_EQ(attribute[n], read_in[n]);
  }
}

TEST(IO, AttributeString) {
  std::string write_file_name = "tests/unit/parallel/data/io_test";
  IO file_io(*mpi_session, "HandleName");
  file_io.SetWrite(write_file_name);
  std::vector<std::string> attribute{
      {"The", "Best", "Series", "Of", "Strings"}};
  file_io.WriteAttribute("Test string", attribute);
  file_io.CloseWriteEngine();

  file_io.SetRead(write_file_name);
  std::vector<std::string> read_in;
  bool found = file_io.ReadAttribute("Test string", read_in);
  EXPECT_TRUE(found);
  for (std::size_t n = 0; n < attribute.size(); ++n) {
    EXPECT_EQ(attribute[n], read_in[n]);
  }
}

TEST(IO, AttributeForVariable) {
  std::string write_file_name = "tests/unit/parallel/data/io_test";
  IO file_io(*mpi_session, "HandleName");
  file_io.SetWrite(write_file_name);
  std::vector<std::string> attribute{
      {"The", "Best", "Series", "Of", "Strings"}};
  file_io.AddVariable<int>("VarName", {}, {}, {1}, true);
  file_io.WriteAttributeForVariable("VarName", "Test string", attribute);
  file_io.CloseWriteEngine();

  file_io.SetRead(write_file_name);
  std::vector<std::string> read_in;
  bool found =
      file_io.ReadAttributeForVariable("VarName", "Test string", read_in);
  EXPECT_TRUE(found);
  for (std::size_t n = 0; n < attribute.size(); ++n) {
    EXPECT_EQ(attribute[n], read_in[n]);
  }
}

TEST(IO, AttributeWithSeparateDrivers) {
  std::string write_file_name = "tests/unit/parallel/data/io_test";
  IO file_io(*mpi_session, "HandleName");
  file_io.SetWrite(write_file_name);
  std::vector<std::string> attribute{
      {"The", "Best", "Series", "Of", "Strings"}};
  file_io.AddVariable<int>("VarName", {}, {}, {1}, true);
  file_io.WriteAttributeForVariable("VarName", "Test string", attribute);
  file_io.CloseWriteEngine();

  IO file_in(*mpi_session, "ReadName");
  file_in.SetRead(write_file_name);
  std::vector<std::string> read_in;
  bool found =
      file_in.ReadAttributeForVariable("VarName", "Test string", read_in);
  EXPECT_TRUE(found);
  for (std::size_t n = 0; n < attribute.size(); ++n) {
    EXPECT_EQ(attribute[n], read_in[n]);
  }
}

TEST(IO, RootAttributeString) {
  std::string write_file_name = "tests/unit/parallel/data/io_test";
  IO file_io(*mpi_session, "HandleName");
  file_io.SetWrite(write_file_name);
  std::vector<std::string> attribute{
      {"The", "Best", "Series", "Of", "Strings"}};
  file_io.RootWriteAttribute("Test string", attribute);
  file_io.CloseWriteEngine();

  file_io.SetRead(write_file_name);
  std::vector<std::string> read_in;
  bool found = file_io.ReadAttribute("Test string", read_in);
  EXPECT_TRUE(found);
  for (std::size_t n = 0; n < attribute.size(); ++n) {
    EXPECT_EQ(attribute[n], read_in[n]);
  }
}

TEST(IO, RootAttributeForVariable) {
  std::string write_file_name = "tests/unit/parallel/data/io_test";
  IO file_io(*mpi_session, "HandleName");
  file_io.SetWrite(write_file_name);
  std::vector<std::string> attribute{
      {"The", "Best", "Series", "Of", "Strings"}};
  file_io.AddVariable<int>("VarName", {}, {}, {1}, true);
  file_io.RootWriteAttributeForVariable("VarName", "Test string", attribute);
  file_io.CloseWriteEngine();

  file_io.SetRead(write_file_name);
  std::vector<std::string> read_in;
  bool found =
      file_io.ReadAttributeForVariable("VarName", "Test string", read_in);
  EXPECT_TRUE(found);
  for (std::size_t n = 0; n < attribute.size(); ++n) {
    EXPECT_EQ(attribute[n], read_in[n]);
  }
}

TEST(IO, OngoingStep) {
  std::string write_file_name = "tests/unit/parallel/data/io_test";
  IO file_io(*mpi_session, "HandleName");
  file_io.SetWrite(write_file_name);
  EXPECT_FALSE(file_io.OngoingWriteStep());
  file_io.BeginWriteStep(0, 0.0, 0.0);
  EXPECT_TRUE(file_io.OngoingWriteStep());
  file_io.EndWriteStep();
  EXPECT_FALSE(file_io.OngoingWriteStep());
}

TEST(IO, DoesVariableExistWrite) {
  std::string write_file_name = "tests/unit/parallel/data/io_test";
  IO file_io(*mpi_session, "HandleName");
  file_io.SetWrite(write_file_name);
  file_io.BeginWriteStep(0, 0.0, 0.0);
  static constexpr std::size_t size_per_proc = 10;
  file_io.AddVariable<int>("VarName", {}, {}, {size_per_proc}, true);
  EXPECT_EQ(file_io.DoesVariableExist<int>("VarName"),
            ExistenceLocation::WRITE);
  EXPECT_EQ(file_io.DoesVariableExist<int>("NonExistent"),
            ExistenceLocation::NONE);
}

TEST(IO, DoesVariableExistRead) {
  std::string write_file_name = "tests/unit/parallel/data/io_test";
  IO file_io(*mpi_session, "HandleName");
  file_io.SetWrite(write_file_name);
  file_io.BeginWriteStep(0, 0.0, 0.0);
  static constexpr std::size_t size_per_proc = 10;
  file_io.AddVariable<int>(
      "VarName", {mpi_session->NumberOfRanks() * size_per_proc},
      {mpi_session->MyRank() * size_per_proc}, {size_per_proc}, true);
  std::vector<int> data(size_per_proc);
  std::fill(data.begin(), data.end(), mpi_session->MyRank());
  file_io.PutDeferred("VarName", data.data());
  file_io.EndWriteStep();
  file_io.CloseWriteEngine();

  IO file_in(*mpi_session, "ReadName");
  file_in.SetRead(write_file_name);
  EXPECT_EQ(file_in.DoesVariableExist<int>("VarName"), ExistenceLocation::READ);
  EXPECT_EQ(file_in.DoesVariableExist<int>("NonExistent"),
            ExistenceLocation::NONE);
}

TEST(IO, DoesAttributeExistWrite) {
  std::string write_file_name = "tests/unit/parallel/data/io_test";
  IO file_io(*mpi_session, "HandleName");
  file_io.SetWrite(write_file_name);
  static constexpr std::size_t size_per_proc = 10;
  file_io.AddVariable<int>(
      "VarName", {mpi_session->NumberOfRanks() * size_per_proc},
      {mpi_session->MyRank() * size_per_proc}, {size_per_proc}, true);
  std::vector<std::string> string_attr{"Test vec string"};
  file_io.WriteAttribute("Test string", string_attr);
  file_io.WriteAttributeForVariable("VarName", "TString", string_attr);

  EXPECT_EQ(file_io.DoesAttributeExist<std::string>("Test string"),
            ExistenceLocation::WRITE);
  EXPECT_EQ(file_io.DoesAttributeExist<int>("Test string"),
            ExistenceLocation::NONE);
  EXPECT_EQ(file_io.DoesAttributeExist<std::string>("VarName/TString"),
            ExistenceLocation::WRITE);
  EXPECT_EQ(file_io.DoesAttributeExist<std::string>("VarName", "TString"),
            ExistenceLocation::WRITE);

  EXPECT_EQ(file_io.DoesAttributeExist<std::string>("NonExistent", "TString"),
            ExistenceLocation::NONE);
  EXPECT_EQ(file_io.DoesAttributeExist<std::string>("VarName", "NonExistent"),
            ExistenceLocation::NONE);
}

TEST(IO, DoesAttributeExistRead) {
  std::string write_file_name = "tests/unit/parallel/data/io_test";
  IO file_io(*mpi_session, "HandleName");
  file_io.SetWrite(write_file_name);
  static constexpr std::size_t size_per_proc = 10;
  file_io.AddVariable<int>(
      "VarName", {mpi_session->NumberOfRanks() * size_per_proc},
      {mpi_session->MyRank() * size_per_proc}, {size_per_proc}, true);
  std::vector<std::string> string_attr{"Test vec string"};
  file_io.WriteAttribute("Test string", string_attr);
  file_io.WriteAttributeForVariable("VarName", "TString", string_attr);
  file_io.CloseWriteEngine();

  IO file_in(*mpi_session, "ReadName");
  file_in.SetRead(write_file_name);

  EXPECT_EQ(file_in.DoesAttributeExist<std::string>("Test string"),
            ExistenceLocation::READ);
  EXPECT_EQ(file_in.DoesAttributeExist<int>("Test string"),
            ExistenceLocation::NONE);
  EXPECT_EQ(file_in.DoesAttributeExist<std::string>("VarName/TString"),
            ExistenceLocation::READ);
  EXPECT_EQ(file_in.DoesAttributeExist<std::string>("VarName", "TString"),
            ExistenceLocation::READ);

  EXPECT_EQ(file_in.DoesAttributeExist<std::string>("NonExistent", "TString"),
            ExistenceLocation::NONE);
  EXPECT_EQ(file_in.DoesAttributeExist<std::string>("VarName", "NonExistent"),
            ExistenceLocation::NONE);
}

TEST(IO, DeferredPutGet) {
  std::string write_file_name = "tests/unit/parallel/data/io_test";
  IO file_io(*mpi_session, "HandleName");
  file_io.SetWrite(write_file_name);
  file_io.BeginWriteStep(0, 0.0, 0.0);
  static constexpr std::size_t size_per_proc = 10;
  file_io.AddVariable<int>(
      "VarName", {mpi_session->NumberOfRanks() * size_per_proc},
      {mpi_session->MyRank() * size_per_proc}, {size_per_proc}, true);
  std::vector<int> data(size_per_proc);
  std::fill(data.begin(), data.end(), mpi_session->MyRank());
  file_io.PutDeferred("VarName", data.data());
  file_io.EndWriteStep();

  IO file_in(*mpi_session, "ReadName");
  file_in.SetRead(write_file_name);
  std::vector<int> read_data(size_per_proc);
  file_in.BeginReadStep();
  file_in.GetDeferred("VarName", {mpi_session->MyRank() * size_per_proc},
                      {size_per_proc}, read_data.data());
  file_in.EndReadStep();

  for (std::size_t n = 0; n < data.size(); ++n) {
    EXPECT_EQ(read_data[n], data[n]);
  }
}

TEST(IO, ImmediatePutGet) {
  std::string write_file_name = "tests/unit/parallel/data/io_test";
  IO file_io(*mpi_session, "HandleName");
  file_io.SetWrite(write_file_name);
  file_io.BeginWriteStep(0, 0.0, 0.0);
  static constexpr std::size_t size_per_proc = 10;
  file_io.AddVariable<int>(
      "VarName", {mpi_session->NumberOfRanks() * size_per_proc},
      {mpi_session->MyRank() * size_per_proc}, {size_per_proc}, true);
  std::vector<int> data(size_per_proc);
  std::fill(data.begin(), data.end(), mpi_session->MyRank());
  file_io.PutImmediate("VarName", data.data());
  std::fill(data.begin(), data.end(),
            mpi_session->MyRank() + 20);  // Modify after put
  file_io.EndWriteStep();

  IO file_in(*mpi_session, "ReadName");
  file_in.SetRead(write_file_name);
  std::vector<int> read_data(size_per_proc);
  file_in.BeginReadStep();
  file_in.GetImmediate("VarName", {mpi_session->MyRank() * size_per_proc},
                       {size_per_proc}, read_data.data());

  // Check before ending step.
  for (std::size_t n = 0; n < data.size(); ++n) {
    EXPECT_EQ(read_data[n], data[n] - 20);
  }
  file_in.EndReadStep();
}

TEST(IO, PutGetMesh) {
  InputParser parser;
  IO file_in(*mpi_session, "ReadName");
  Mesh mesh(*mpi_session, parser, file_in);

  parser.DirectSet("mesh_file", std::string("generate"));
  parser.DirectSet("gen_nx", 10);
  parser.DirectSet("gen_ny", 10);
  parser.DirectSet("gen_nz", 10);
  parser.DirectSet("serial_refine", 0);
  parser.DirectSet("parallel_refine", 0);

  mesh.Initialize();

  IO file_io(*mpi_session, "HandleName");
  std::string write_file_name = "tests/unit/parallel/data/io_test";
  file_io.SetWrite(write_file_name);
  file_io.BeginWriteStep(0, 0.0, 0.0);
  file_io.AddVariableForMesh<int>("VarName", mesh, MeshElement::ELEMENT);
  std::vector<int> data(mesh.GetLocalCount<MeshElement::ELEMENT>());
  std::fill(data.begin(), data.end(), mpi_session->MyRank());
  file_io.PutDeferred("VarName", data.data());
  file_io.EndWriteStep();
  file_io.CloseWriteEngine();

  file_in.SetRead(write_file_name);
  std::vector<int> read_data(mesh.GetLocalCount<MeshElement::ELEMENT>());
  file_in.BeginReadStep();
  file_in.GetDeferred("VarName", mesh, MeshElement::ELEMENT, read_data.data());
  file_in.EndReadStep();

  // Check before ending step.
  for (std::size_t n = 0; n < data.size(); ++n) {
    EXPECT_EQ(read_data[n], data[n]);
  }
}

}  // namespace
