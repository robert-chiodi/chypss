# This file is part of the Coupled Hypersonic Protected System (CHyPS) Simulator.
#
# Copyright (C) 2020 Robert Chiodi <chiodi@illinois.edu>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

#[=======================================================================[.rst:
FindSpdlog
-------

Finds the Spdlog library.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

``Spdlog::Spdlog``
The Spdlog library

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``Spdlog_FOUND``
True if the system has the Spdlog library.
``Spdlog_VERSION``
The version of the Spdlog library which was found.
``Spdlog_INCLUDE_DIRS``
Include directories needed to use Spdlog.
``Spdlog_LIBRARIES``
Libraries needed to link to Spdlog.

Cache Variables
^^^^^^^^^^^^^^^

The following cache variables may also be set:

``Spdlog_INCLUDE_DIR``
The directory containing ``Spdlog.h``.
``Spdlog_LIBRARY``
The path to the Spdlog library.

#]=======================================================================]

# Try to grab information from PkgConfig
find_package(PkgConfig)
pkg_check_modules(PC_Spdlog QUIET Spdlog)

if(NOT Spdlog_DIR)
  find_path(Spdlog_INCLUDE_DIR
  NAMES spdlog.hpp
  PATHS ${PC_Spdlog_INCLUDE_DIRS}
  PATH_SUFFIXES spdlog
  )
  find_library(Spdlog_LIBRARY
  NAMES spdlog
  PATHS ${PC_Spdlog_LIBRARY_DIRS}
  )
else()
  find_path(Spdlog_INCLUDE_DIR
  NAMES spdlog/spdlog.h
  PATHS ${Spdlog_DIR}/include
  PATH_SUFFIXES spdlog
  NO_DEFAULT_PATH
  )
  find_library(Spdlog_LIBRARY
  NAMES spdlog
  PATHS ${Spdlog_DIR}/lib
  NO_DEFAULT_PATH
  )
endif()

set(Spdlog_VERSION ${PC_Spdlog_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Spdlog
FOUND_VAR Spdlog_FOUND
REQUIRED_VARS
Spdlog_LIBRARY
Spdlog_INCLUDE_DIR
VERSION_VAR Spdlog_VERSION
)

if(NOT Spdlog_FOUND)
message("Spdlog library not found. Define Spdlog_DIR with path to Spdlog.")
return()
endif()

if(Spdlog_FOUND)
set(Spdlog_LIBRARIES ${Spdlog_LIBRARY})
set(Spdlog_INCLUDE_DIRS ${Spdlog_INCLUDE_DIR})
set(Spdlog_DEFINITIONS ${PC_Spdlog_CFLAGS_OTHER})
endif()

if(Spdlog_FOUND AND NOT TARGET Spdlog::Spdlog)
add_library(Spdlog::Spdlog UNKNOWN IMPORTED)
set_target_properties(Spdlog::Spdlog PROPERTIES
IMPORTED_LOCATION "${Spdlog_LIBRARY}"
INTERFACE_COMPILE_OPTIONS "${PC_Spdlog_CFLAGS_OTHER}"
INTERFACE_INCLUDE_DIRECTORIES "${Spdlog_INCLUDE_DIR}"
)
endif()

mark_as_advanced(
Spdlog_INCLUDE_DIR
Spdlog_LIBRARY
)
