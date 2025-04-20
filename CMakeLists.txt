# Main CMakeLists.txt structure
cmake_minimum_required(VERSION 3.16)
project(OpenSeesFramework LANGUAGES C CXX Fortran)

# Setup basic properties
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(OPS_EXTERN_SOURCE_DIR "${PROJECT_SOURCE_DIR}/OTHER/")
set(OPS_EXTERNALS_DIR "${PROJECT_SOURCE_DIR}/OTHER/")
set(OPS_BUNDLED_DIR "${PROJECT_SOURCE_DIR}/OTHER/")
set(OPS_SRC_DIR "${PROJECT_SOURCE_DIR}/SRC/")

# Set up the module path to find our custom CMake modules
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake/modules")

# Include user config
include(${PROJECT_SOURCE_DIR}/Conf.cmake)

# Include functionality modules
include(cmake/modules/OpenSeesFunctions.cmake)
include(cmake/modules/ConanDetection.cmake)
include(cmake/modules/CompilerSettings.cmake)
include(cmake/modules/DependencyHandling.cmake)
include(cmake/modules/GlobalIncludes.cmake)

# Include library definitions
include(cmake/libraries/NumericalLibraries.cmake)
include(cmake/libraries/OpenSeesLibraries.cmake)

# Include executable definitions
include(cmake/executables/OpenSeesExe.cmake)
include(cmake/executables/OpenSeesSP.cmake)  # Only if MPI found
include(cmake/executables/OpenSeesMP.cmake)  # Only if MPI found
include(cmake/executables/OpenSeesPy.cmake)  # Only if Python found

# Apply build options
include(cmake/options/ApplyOptions.cmake)

# Install rules
include(cmake/InstallRules.cmake)

# Add sources to targets
add_subdirectory(${OPS_SRC_DIR})