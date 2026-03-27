#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#                            External Libraries
#
# DEPRECATED: This file is superseded by OpenSeesDependencies.cmake
# (unified cross-platform).  Kept for backward compatibility with
# external build scripts that may include it directly.
#==============================================================================
message(AUTHOR_WARNING
  "OpenSeesDependenciesWin.cmake is deprecated.  "
  "Use OpenSeesDependencies.cmake instead.")


if (NOT TCL_LIBRARY AND DEFINED TCL_LIBRARIES)
  set(TCL_LIBRARY "${TCL_LIBRARIES}")
endif()
if (NOT TCL_INCLUDE_PATH AND DEFINED TCL_INCLUDE_DIR)
  set(TCL_INCLUDE_PATH "${TCL_INCLUDE_DIR}")
endif()
if (NOT TCL_LIBRARY OR NOT TCL_INCLUDE_PATH)
  find_package(TCL REQUIRED)
endif()
find_package(MySQL ) # Not required

#
#sudo apt-get install libhdf5-serial-dev
#  - installed version 1.10
#  - we need version 1.12
#  - hdf5-1.12.1-linux-centos7-x86_64-gcc485-static.tar.gz

set (HDF5_USE_STATIC_LIBRARIES ON)
find_package(HDF5)
if (HDF5_FOUND)
     message(STATUS "HDF5 found version: ${HDF5_VERSION}")
     message(STATUS "HDF5_CXX_DEFINITIONS = ${HDF5_CXX_DEFINITIONS}")
     message(STATUS "HDF5_LIBRARIES = ${HDF5_LIBRARIES}" )
  if (HDF_VERSION VERSION_LESS 1.12)
     message(STATUS "HDF5 VERSION OLD: ${HDF5_VERSION}" )
  endif()
endif()


