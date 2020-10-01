# This file is part of the Coupled Hypersonic Protected System (CHyPS) Simulator.
#
# Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

#[=======================================================================[.rst:
FindPetsc
-------

Finds the Petsc library.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

``Petsc::Petsc``
The Petsc library

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``Petsc_FOUND``
True if the system has the Petsc library.
``Petsc_VERSION``
The version of the Petsc library which was found.
``Petsc_INCLUDE_DIRS``
Include directories needed to use Petsc.
``Petsc_LIBRARIES``
Libraries needed to link to Petsc.

Cache Variables
^^^^^^^^^^^^^^^

The following cache variables may also be set:

``Petsc_INCLUDE_DIR``
The directory containing ``Petsc.h``.
``Petsc_LIBRARY``
The path to the Petsc library.

#]=======================================================================]

# Try to grab information from PkgConfig
find_package(PkgConfig)
pkg_check_modules(PC_Petsc QUIET Petsc)

find_path(Petsc_INCLUDE_DIR
NAMES petsc.h
PATHS ${PC_Petsc_INCLUDE_DIRS} ${Petsc_DIR}/include
PATH_SUFFIXES Petsc
)
find_library(Petsc_LIBRARY
NAMES petsc
PATHS ${PC_Petsc_LIBRARY_DIRS} ${Petsc_DIR}/lib
)

set(Petsc_VERSION ${PC_Petsc_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Petsc
FOUND_VAR Petsc_FOUND
REQUIRED_VARS
Petsc_LIBRARY
Petsc_INCLUDE_DIR
VERSION_VAR Petsc_VERSION
)

if(NOT Petsc_FOUND)
message("Petsc library not found. Define Petsc_DIR with path to Petsc.")
return()
endif()

if(Petsc_FOUND)
set(Petsc_LIBRARIES ${Petsc_LIBRARY})
set(Petsc_INCLUDE_DIRS ${Petsc_INCLUDE_DIR})
set(Petsc_DEFINITIONS ${PC_Petsc_CFLAGS_OTHER})
endif()

if(Petsc_FOUND AND NOT TARGET Petsc::Petsc)
add_library(Petsc::Petsc UNKNOWN IMPORTED)
set_target_properties(Petsc::Petsc PROPERTIES
IMPORTED_LOCATION "${Petsc_LIBRARY}"
INTERFACE_COMPILE_OPTIONS "${PC_Petsc_CFLAGS_OTHER}"
INTERFACE_INCLUDE_DIRECTORIES "${Petsc_INCLUDE_DIR}"
)
endif()

mark_as_advanced(
Petsc_INCLUDE_DIR
Petsc_LIBRARY
)
