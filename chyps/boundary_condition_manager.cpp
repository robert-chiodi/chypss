// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "chyps/boundary_condition_manager.hpp"

#include <string>
#include <utility>

#include "chyps/logger.hpp"
#include "chyps/mesh.hpp"

namespace chyps {

namespace details {
// H is homogeneous, TV is time-varying
enum BCType {
  H_DIRICHLET = 0,
  DIRICHLET,
  TV_DIRICHLET,
  H_NEUMANN,
  NEUMANN,
  TV_NEUMANN,
  COUNT
};
}  // namespace details

BoundaryConditionManager::BoundaryConditionManager(InputParser& a_parser)
    : parser_m(a_parser),
      mesh_m(nullptr),
      boundary_condition_counts_m(),
      boundary_conditions_m(),
      precice_condition_m(),
      vertex_positions_m(),
      indices_m(),
      values_m() {
  this->GatherOptions();
}

void BoundaryConditionManager::Initialize(const Mesh& a_mesh) {
  // FIXME : make this an exception.
  if (!this->AllOptionsSupplied()) {
    std::cout << "Not all options needed for BoundaryConditionManager supplied"
              << std::endl;
    std::cout << "Make sure that the InputParser has been parsed before "
                 "calling Initialize and that all required options are "
                 "specified or have a valid default value."
              << std::endl;
    std::exit(-1);
  }
  mesh_m = &a_mesh;
  std::size_t number_of_bcs =
      static_cast<std::size_t>(mesh_m->GetNumberOfBoundaryTagValues());
  boundary_condition_counts_m.resize(details::BCType::COUNT, 0);
  boundary_conditions_m.resize(number_of_bcs);
  precice_condition_m.resize(number_of_bcs, false);
  values_m.resize(number_of_bcs);
  indices_m.resize(number_of_bcs);

  this->SetBoundaryConditionsFromInput();
}

int BoundaryConditionManager::GetNumberOfBoundaryConditions(void) const {
  return static_cast<int>(boundary_conditions_m.size());
}

void BoundaryConditionManager::SetBoundaryConditionType(
    const int a_tag, const BoundaryConditionType& a_bc_type,
    const bool a_spatially_varying, const bool a_time_varying) {
  assert(this->HasBeenInitialized());
  assert(a_tag > 0);
  assert(a_tag <= this->GetNumberOfBoundaryConditions());

  SPDLOG_LOGGER_INFO(MAIN_LOG,
                     "Setting boundary condition for tag {}: type({}), "
                     "spatially_varying({}), temporally_varying({})",
                     a_tag, static_cast<int>(a_bc_type), a_spatially_varying,
                     a_time_varying);

  // Remove current occupier of a_tag from boundary_condition_counts_m;
  const auto previous_condition_count_type =
      this->GetConditionCountType(boundary_conditions_m[a_tag - 1]);
  --boundary_condition_counts_m[previous_condition_count_type];

  boundary_conditions_m[a_tag - 1] =
      BoundaryCondition(a_bc_type, a_spatially_varying, a_time_varying);
  if (a_bc_type == BoundaryConditionType::DIRICHLET ||
      a_bc_type == BoundaryConditionType::NEUMANN) {
    if (a_spatially_varying) {
      std::tie(vertex_positions_m[a_tag - 1], indices_m[a_tag - 1]) =
          mesh_m->GetBoundaryVertices(a_tag);
      values_m[a_tag - 1].resize(indices_m[a_tag - 1].size());
      boundary_conditions_m[a_tag - 1].SetValues(
          values_m[a_tag - 1].size(), values_m[a_tag - 1].data(),
          indices_m[a_tag - 1].data(), false);
    } else {
      values_m[a_tag - 1].resize(1);
      boundary_conditions_m[a_tag - 1].SetValues(values_m[a_tag - 1].data(),
                                                 false);
    }
  }
  const auto current_condition_count_type =
      this->GetConditionCountType(boundary_conditions_m[a_tag - 1]);
  ++boundary_condition_counts_m[current_condition_count_type];
}

std::pair<const std::vector<double>*, const std::vector<int>*>
BoundaryConditionManager::GetBoundaryVertices(const int a_tag) const {
  assert(this->HasBeenInitialized());
  assert(a_tag > 0);
  assert(a_tag <= this->GetNumberOfBoundaryConditions());
  assert(boundary_conditions_m[a_tag - 1].GetBCType() ==
             BoundaryConditionType::DIRICHLET ||
         boundary_conditions_m[a_tag - 1].GetBCType() ==
             BoundaryConditionType::NEUMANN);
  assert(values_m[a_tag - 1].size() > 0);
  return std::make_pair(&(vertex_positions_m[a_tag - 1]),
                        &(indices_m[a_tag - 1]));
}

void BoundaryConditionManager::SetBoundaryConditionValues(
    const int a_tag, const std::vector<double>& a_values) {
  assert(this->HasBeenInitialized());
  assert(a_tag > 0);
  assert(a_tag <= this->GetNumberOfBoundaryConditions());
  assert(boundary_conditions_m[a_tag - 1].GetBCType() ==
             BoundaryConditionType::DIRICHLET ||
         boundary_conditions_m[a_tag - 1].GetBCType() ==
             BoundaryConditionType::NEUMANN);
  assert(boundary_conditions_m[a_tag - 1].IsSpatiallyVarying() ||
         a_values.size() == 1);
  std::copy(a_values.begin(), a_values.end(), values_m[a_tag - 1].begin());
}

void BoundaryConditionManager::SetBoundaryConditionValues(
    const int a_tag, const double a_value) {
  assert(this->HasBeenInitialized());
  assert(a_tag > 0);
  assert(a_tag <= this->GetNumberOfBoundaryConditions());
  assert(boundary_conditions_m[a_tag - 1].GetBCType() ==
             BoundaryConditionType::DIRICHLET ||
         boundary_conditions_m[a_tag - 1].GetBCType() ==
             BoundaryConditionType::NEUMANN);
  assert(!boundary_conditions_m[a_tag - 1].IsSpatiallyVarying());
  values_m[a_tag - 1][0] = a_value;
}

void BoundaryConditionManager::SetBoundaryConditionAsPrecice(const int a_tag) {
  assert(this->HasBeenInitialized());
  assert(a_tag > 0);
  assert(a_tag <= this->GetNumberOfBoundaryConditions());
  assert(boundary_conditions_m[a_tag - 1].GetBCType() ==
             BoundaryConditionType::DIRICHLET ||
         boundary_conditions_m[a_tag - 1].GetBCType() ==
             BoundaryConditionType::NEUMANN);
  precice_condition_m[a_tag - 1] = true;
}

double* BoundaryConditionManager::GetDataBuffer(const int a_tag) {
  assert(this->HasBeenInitialized());
  assert(a_tag > 0);
  assert(a_tag <= this->GetNumberOfBoundaryConditions());
  assert(boundary_conditions_m[a_tag - 1].GetBCType() ==
             BoundaryConditionType::DIRICHLET ||
         boundary_conditions_m[a_tag - 1].GetBCType() ==
             BoundaryConditionType::NEUMANN);
  assert(precice_condition_m[a_tag - 1]);
  return values_m[a_tag - 1].data();
}

const BoundaryCondition& BoundaryConditionManager::GetBoundaryCondition(
    const int a_tag) const {
  assert(this->HasBeenInitialized());
  assert(a_tag > 0);
  assert(a_tag <= this->GetNumberOfBoundaryConditions());
  return boundary_conditions_m[a_tag - 1];
}

int BoundaryConditionManager::GetNumberOfHomogeneousDirichletConditions(
    void) const {
  return boundary_condition_counts_m[details::BCType::H_DIRICHLET];
}

int BoundaryConditionManager::GetNumberOfDirichletConditions(void) const {
  return boundary_condition_counts_m[details::BCType::DIRICHLET];
}

int BoundaryConditionManager::GetNumberOfTimeVaryingDirichletConditions(
    void) const {
  assert(boundary_condition_counts_m[details::BCType::TV_DIRICHLET] <=
         boundary_condition_counts_m[details::BCType::DIRICHLET]);
  return boundary_condition_counts_m[details::BCType::TV_DIRICHLET];
}

int BoundaryConditionManager::GetNumberOfHomogeneousNeumannConditions(
    void) const {
  return boundary_condition_counts_m[details::BCType::H_NEUMANN];
}

int BoundaryConditionManager::GetNumberOfNeumannConditions(void) const {
  return boundary_condition_counts_m[details::BCType::NEUMANN];
}

int BoundaryConditionManager::GetNumberOfTimeVaryingNeumannConditions(
    void) const {
  assert(boundary_condition_counts_m[details::BCType::TV_NEUMANN] <=
         boundary_condition_counts_m[details::BCType::NEUMANN]);
  return boundary_condition_counts_m[details::BCType::TV_NEUMANN];
}

bool BoundaryConditionManager::AllOptionsSupplied(void) const {
  return parser_m.AllOptionsSet();
}

void BoundaryConditionManager::GatherOptions(void) {
  parser_m.AddOptionDefault(
      "BCManager/bc_default",
      "Default boundary condition type to use for the boundary "
      "conditions not "
      "explicitly supplied in the bc_list option. Options are "
      "HOMOGENEOUS_DIRICHLET or HOMOGENEOUS_NEUMANN",
      std::string("HOMOGENEOUS_DIRICHLET"));
  parser_m.AddOptionNoDefault(
      "BCManager/bc_list",
      "A nested JSON object describing the boundary conditions. The keys "
      "inside bc_list should be the tag of the boundary condition being set. "
      "Inside this tag their should be {key,value} pairs. Valid {key,value} "
      "pairs are: {type, HOMOGENEOUS_DIRICHLET || HOMOGENEOUS_NEUMANN || "
      "DIRICHLET || NEUMANN}, {value, double}, {precice, true || false}",
      false);
}

void BoundaryConditionManager::SetBoundaryConditionsFromInput(void) {
  const auto default_bc = this->BoundaryConditionNameToEnum(
      parser_m["BCManager/bc_default"].get<std::string>());
  assert(default_bc == BoundaryConditionType::HOMOGENEOUS_DIRICHLET ||
         default_bc == BoundaryConditionType::HOMOGENEOUS_NEUMANN);

  // Set default everywhere, then update (if needed) when calls to
  // SetBoundaryConditionType happen.
  std::fill(boundary_conditions_m.begin(), boundary_conditions_m.end(),
            BoundaryCondition(default_bc, false, false));
  if (default_bc == BoundaryConditionType::HOMOGENEOUS_DIRICHLET) {
    boundary_condition_counts_m[details::BCType::H_DIRICHLET] +=
        boundary_conditions_m.size();
  } else if (default_bc == BoundaryConditionType::HOMOGENEOUS_NEUMANN) {
    boundary_condition_counts_m[details::BCType::H_NEUMANN] +=
        boundary_conditions_m.size();
  }

  if (!parser_m.OptionSet("BCManager/bc_list")) {
    // No list provided, use only default
    return;
  }

  const auto& bc_list = parser_m["BCManager/bc_list"];

  for (const auto& el : bc_list.items()) {
    const auto& key = el.key();
    const auto& value = el.value();
    const int tag = std::stoi(key);
    assert(tag > 0);
    assert(tag <= this->GetNumberOfBoundaryConditions());
    assert(value.contains("type"));
    const auto bc_type =
        this->BoundaryConditionNameToEnum(value["type"].get<std::string>());
    this->SetBoundaryConditionType(tag, bc_type, false, false);

    double bc_value = 0.0;
    if (bc_type == BoundaryConditionType::DIRICHLET ||
        bc_type == BoundaryConditionType::NEUMANN) {
      assert(value.contains("value"));
      bc_value = value["value"].get<double>();
      this->SetBoundaryConditionValues(tag, bc_value);
    }
    bool precice_condition = false;
    if (value.contains("precice")) {
      precice_condition = value["precice"].get<bool>();
      if (precice_condition) {
        this->SetBoundaryConditionAsPrecice(tag);
      }
    }
  }
}

bool BoundaryConditionManager::HasBeenInitialized(void) const {
  return mesh_m != nullptr;
}

BoundaryConditionType BoundaryConditionManager::BoundaryConditionNameToEnum(
    const std::string a_bc_name) {
  if (a_bc_name == "HOMOGENEOUS_DIRICHLET") {
    return BoundaryConditionType::HOMOGENEOUS_DIRICHLET;
  } else if (a_bc_name == "HOMOGENEOUS_NEUMANN") {
    return BoundaryConditionType::HOMOGENEOUS_NEUMANN;
  } else if (a_bc_name == "DIRICHLET") {
    return BoundaryConditionType::DIRICHLET;
  } else if (a_bc_name == "HOMOGENEOUS_DIRICHLET") {
    return BoundaryConditionType::NEUMANN;
  } else {
    // FIXME: Make proper error
    std::cout << "Unknown boundary condition name in input file of: "
              << a_bc_name << std::endl;
  }
  // NOTE: Will never get this far.
  return BoundaryConditionType::HOMOGENEOUS_DIRICHLET;
}

int BoundaryConditionManager::GetConditionCountType(
    const BoundaryCondition& a_bc) {
  switch (a_bc.GetBCType()) {
    case BoundaryConditionType::HOMOGENEOUS_DIRICHLET: {
      return details::BCType::H_DIRICHLET;
      break;
    }
    case BoundaryConditionType::HOMOGENEOUS_NEUMANN: {
      return details::BCType::H_NEUMANN;
      break;
    }
    case BoundaryConditionType::DIRICHLET: {
      return a_bc.IsTimeVarying() ? details::BCType::TV_DIRICHLET
                                  : details::BCType::DIRICHLET;
      break;
    }
    case BoundaryConditionType::NEUMANN: {
      return a_bc.IsTimeVarying() ? details::BCType::TV_NEUMANN
                                  : details::BCType::NEUMANN;
      break;
    }
    default:
      // FIXME: Make actual error handler.
      std::cout << "Unknown BC type in BoundaryConditionManager of: "
                << static_cast<int>(a_bc.GetBCType()) << std::endl;
      std::exit(-1);
  }
}

}  // namespace chyps
