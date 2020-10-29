// This file is part of the Coupled Hypersonic Protected System (CHyPS)
// Simulator
//
//
// Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef CHYPS_BOUNDARY_CONDITION_MANAGER_HPP_
#define CHYPS_BOUNDARY_CONDITION_MANAGER_HPP_

#include <utility>

#include "chyps/boundary_condition.hpp"
#include "chyps/input_parser.hpp"

namespace chyps {

// Forward declare Mesh to prevent cyclical reference
class Mesh;

class BoundaryConditionManager {
 public:
  /// \brief Constructor for boundary condition manager with a specified default
  /// BC type.
  ///
  BoundaryConditionManager(InputParser& a_parser);

  /// \brief Initialize the boundary condition manager. Must be called after
  /// initializing the associated a_mesh.
  ///
  /// Only BoundaryConditionType::HOMOGENEOUS_DIRICHLET and
  /// BoundaryConditionType::HOMOGENEOUS_NEUMANN are valid default types.
  ///
  /// NOTE: A reference is held to the provided mesh. Do not invalidate it.
  void Initialize(const Mesh& a_mesh);

  /// \brief Set the type of boundary condition for a_tag. The value of a_tag
  /// must have corresponding boundary elements in a_mesh.
  void SetBoundaryConditionType(const int a_tag,
                                const BoundaryConditionType& a_bc_type,
                                const bool a_spatially_varying,
                                const bool a_time_varying);

  /// \brief For a_tag, return the list of corresponding boundary vertex
  /// positions and the index from the supplied a_mesh. The positions are
  /// ordered as X1/Y1/Z1, X2/Y2/Z2, etc. Only existing dimensions in a_mesh are
  /// used, meaning the std::vector<double>::size will be dim*NVert long.
  ///
  /// NOTE: Before using this method, the boundary with a_tag must first be set
  /// as a spatially varying DIRICHLET or NEUMANN condition type with the
  /// SetBoundaryConditionType method.
  std::pair<const std::vector<double>*, const std::vector<int>*>
  GetBoundaryVertices(const int a_tag) const;

  /// \brief Set the values for a_tag. Length of a_values must be equal to the
  /// number of vertices that exist for that tag if a spatially varying BC.
  /// Otherwise, a_values should be of size 1.
  void SetBoundaryConditionValues(const int a_tag,
                                  const std::vector<double>& a_values);

  /// \brief Set the values for a_tag that is a non-spatially varying
  /// inhomogeneous boundary condition.
  void SetBoundaryConditionValues(const int a_tag, const double a_value);

  /// \brief Mark as a boundary condition that will be set through preCICE
  /// coupling to another solver.
  void SetBoundaryConditionAsPrecice(const int a_tag);

  /// \brief Return pointer to the data buffer for the values of the boundary
  /// condition marked with a_tag. This is only allowed for preCICE boundary
  /// conditions (which implies it is also either of DIRICHLET or NEUMANN type).
  double* GetDataBuffer(const int a_tag);

  /// \brief Return number of boundary conditions.
  int GetNumberOfBoundaryConditions(void) const;

  /// \brief Return boundary condition for boundary element tag a_tag.
  const BoundaryCondition& GetBoundaryCondition(const int a_tag) const;

  /// \brief Return number of homogeneous Dirichlet conditions set.
  int GetNumberOfHomogeneousDirichletConditions(void) const;

  /// \brief Return number of inhomogeneous Dirichlet conditions set.
  int GetNumberOfDirichletConditions(void) const;

  /// \brief Return number of time-varying Dirichlet conditions.
  ///
  /// Note: This must be <= the number of inhomogeneous Dirichlet conditions.
  int GetNumberOfTimeVaryingDirichletConditions(void) const;

  /// \brief Return number of homogeneous Neumann conditions set.
  int GetNumberOfHomogeneousNeumannConditions(void) const;

  /// \brief Return number of inhomogeneous Neumann conditions set.
  int GetNumberOfNeumannConditions(void) const;

  /// \brief Return number of time-varying Neumann conditions.
  ///
  /// Note: This must be <= the number of inhomogeneous Neumann conditions.
  int GetNumberOfTimeVaryingNeumannConditions(void) const;

 private:
  void GatherOptions(void);
  bool AllOptionsSupplied(void) const;
  void SetBoundaryConditionsFromInput(void);
  bool HasBeenInitialized(void) const;
  static BoundaryConditionType BoundaryConditionNameToEnum(
      const std::string a_bc_name);
  static int GetConditionCountType(const BoundaryCondition& a_bc);

  InputParser& parser_m;
  const Mesh* mesh_m;
  std::vector<int> boundary_condition_counts_m;
  std::vector<BoundaryCondition> boundary_conditions_m;
  std::vector<bool> precice_condition_m;
  std::vector<std::vector<double>> vertex_positions_m;
  std::vector<std::vector<int>> indices_m;
  std::vector<std::vector<double>> values_m;
};

}  // namespace chyps

#endif  // CHYPS_BOUNDARY_CONDITION_MANAGER_HPP_
