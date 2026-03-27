#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#                   Unified External Dependencies (all platforms)
#
# This file replaces the separate OpenSeesDependenciesUnix.cmake and
# OpenSeesDependenciesWin.cmake with a single cross-platform dependency
# discovery module.  It is included from the platform blocks in the main
# CMakeLists.txt when Conan is NOT being used.
#==============================================================================

# Guard: nothing to do when Conan already resolved everything.
if(OPS_USING_CONAN)
  return()
endif()

# --------------------------------------------------------------------------
# TCL
# --------------------------------------------------------------------------
# Normalize variable names coming from different providers:
#   CMake FindTCL  -> TCL_LIBRARY, TCL_INCLUDE_PATH
#   vcpkg / Conan  -> TCL_LIBRARIES (plural), TCL_INCLUDE_DIR
#   User flags     -> either form
#
# Strategy:
#   1. Map alternate names so pre-set values are recognised.
#   2. Run find_package only if still missing.
#   3. Map again (find_package may set the "other" names).
#   4. Fatal error if still unresolved.

if(NOT TCL_LIBRARY AND DEFINED TCL_LIBRARIES)
  set(TCL_LIBRARY "${TCL_LIBRARIES}")
endif()
if(NOT TCL_INCLUDE_PATH AND DEFINED TCL_INCLUDE_DIR)
  set(TCL_INCLUDE_PATH "${TCL_INCLUDE_DIR}")
endif()

if(NOT TCL_LIBRARY OR NOT TCL_INCLUDE_PATH)
  find_package(TCL REQUIRED)
endif()

# Post-find normalization (find_package may set only one naming convention).
if(NOT TCL_INCLUDE_PATH AND DEFINED TCL_INCLUDE_DIR)
  set(TCL_INCLUDE_PATH "${TCL_INCLUDE_DIR}")
endif()
if(NOT TCL_LIBRARY AND DEFINED TCL_LIBRARIES)
  set(TCL_LIBRARY "${TCL_LIBRARIES}")
endif()

if(NOT TCL_LIBRARY OR NOT TCL_INCLUDE_PATH)
  message(FATAL_ERROR
    "Tcl discovery failed.  Set -DTCL_LIBRARY=<path> and "
    "-DTCL_INCLUDE_PATH=<path> explicitly.")
endif()

# Publish the canonical names used by the rest of the build.
include_directories(${TCL_INCLUDE_PATH})
set(TCL_LIBRARIES ${TCL_LIBRARY})

message(STATUS "TCL_LIBRARY      = ${TCL_LIBRARY}")
message(STATUS "TCL_INCLUDE_PATH = ${TCL_INCLUDE_PATH}")

# --------------------------------------------------------------------------
# MySQL (optional on all platforms)
# --------------------------------------------------------------------------
find_package(MySQL)

# --------------------------------------------------------------------------
# HDF5 (prefer static; not REQUIRED — the build handles HDF5_FOUND=FALSE)
# --------------------------------------------------------------------------
set(HDF5_USE_STATIC_LIBRARIES ON)
find_package(HDF5)
if(HDF5_FOUND)
  message(STATUS "HDF5 found version: ${HDF5_VERSION}")
  message(STATUS "HDF5_LIBRARIES = ${HDF5_LIBRARIES}")
  if(HDF5_VERSION VERSION_LESS 1.12)
    message(STATUS "HDF5 VERSION OLD: ${HDF5_VERSION}")
  endif()
endif()
