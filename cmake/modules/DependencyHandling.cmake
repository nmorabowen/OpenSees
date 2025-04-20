# DependencyHandling.cmake
# Handles external dependencies based on Conan version or user-supplied paths

# Create interface libraries for dependencies
add_library(OPS_External_packages INTERFACE)
add_library(OPS_OS_Specific_libs INTERFACE)

# Handle different Conan versions or manual dependency setup
if(USING_CONAN2)
  # With Conan 2, we use find_package for dependencies
  if(DEFINED CONAN_GENERATORS_DIR)
    # Add generators directory to CMAKE_PREFIX_PATH and CMAKE_MODULE_PATH
    list(APPEND CMAKE_PREFIX_PATH "${CONAN_GENERATORS_DIR}")
    list(APPEND CMAKE_MODULE_PATH "${CONAN_GENERATORS_DIR}")
    
    MESSAGE("Using Conan generators directory: ${CONAN_GENERATORS_DIR}")
    MESSAGE("CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")
    MESSAGE("CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
    
    # Try to include the toolchain file
    if(EXISTS "${CONAN_GENERATORS_DIR}/conan_toolchain.cmake")
      include("${CONAN_GENERATORS_DIR}/conan_toolchain.cmake")
    endif()
  endif()
  
  # Find each package explicitly
  find_package(Eigen3 CONFIG REQUIRED)
  find_package(TCL NAMES TCL CONFIG REQUIRED)
  find_package(ZLIB CONFIG REQUIRED)
  find_package(HDF5 CONFIG REQUIRED)
  
  # Set variables for rest of the build
  set(Eigen3_FOUND TRUE)
  set(HDF5_FOUND TRUE)
  set(HDF5_LIBRARIES hdf5::hdf5 zlib::zlib)
  set(HDF5_VERSION "1.14.0")
  
  # TCL settings
  include_directories(${TCL_INCLUDE_DIRS})
  set(TCL_INCLUDE_PATH ${TCL_INCLUDE_DIRS})
  set(TCL_LIBRARY tcl::tcl)
  set(TCL_LIBRARIES tcl::tcl)
  
  # Get paths for copying TCL libraries if possible
  if(TARGET tcl::tcl)
    get_target_property(TCL_LIB_LOCATION tcl::tcl LOCATION)
    if(TCL_LIB_LOCATION)
      get_filename_component(TCL_LIB_DIR "${TCL_LIB_LOCATION}" DIRECTORY)
    endif()
  endif()

elseif(USING_CONAN1)
  # Using Conan 1
  include("${CMAKE_BINARY_DIR}/conanbuildinfo.cmake")
  conan_basic_setup()
  set(Eigen3_FOUND TRUE)
  set(HDF5_FOUND TRUE)
  set(HDF5_LIBRARIES ${CONAN_LIBS_HDF5} ${CONAN_LIBS_ZLIB})
  set(HDF5_VERSION "1.12.0")
  set(CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR} ${CMAKE_MODULE_PATH})
  set(CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR} ${CMAKE_PREFIX_PATH})
  set(TCL_LIBRARIES ${TCL_LIBRARY})
  include_directories(${TCL_INCLUDE_PATH})
else()
  # Not using Conan - find dependencies manually
  include(OpenSeesFunctions)
  set(NOT_USING_CONAN TRUE)
  set(CONAN_LIBS)
  find_package(HDF5 REQUIRED)
  find_package(TCL REQUIRED)
  find_package(Eigen3 REQUIRED)
  include_directories(${TCL_INCLUDE_DIR})
  set(TCL_INCLUDE_PATH ${TCL_INCLUDE_DIR})
  set(TCL_LIBRARY ${TCL_LIBRARIES})
endif()

# Rest of the file remains the same
# OS-specific configuration
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  add_compile_definitions(_LINUX _UNIX _TCL85)
  if(NOT_USING_CONAN)
    include(OpenSeesDependenciesUnix)
  endif()
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  set(OPS_Use_Graphics_Option None)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  add_compile_definitions(_LINUX _UNIX _TCL85)
  if(NOT_USING_CONAN)
    include(OpenSeesDependenciesUnix)
  endif()  
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
  target_link_libraries(OPS_OS_Specific_libs INTERFACE wsock32 ws2_32)
  add_compile_definitions(_WIN32 _TCL85)
  if(NOT_USING_CONAN)
    include(OpenSeesDependenciesWin)
  endif()      
endif()

# Handle other packages
if(NOT DEFINED LAPACK_LIBRARIES) 
  find_package(LAPACK)
else()
  set(LAPACK_FOUND TRUE CACHE BOOL "providing own lapack lib")
endif()

if(NOT DEFINED Python_LIBRARIES) 
  find_package(Python COMPONENTS Interpreter Development)
else()
  set(PYTHON_FOUND TRUE CACHE BOOL "providing own python lib")
endif()

find_package(MPI)
find_package(MKL)

# Set up LAPACK if not found
if(LAPACK_FOUND)
  message(STATUS "LAPACK was found.")
  message(STATUS "LAPACK_LINKER_FLAGS = ${LAPACK_LINKER_FLAGS}")
  message(STATUS "LAPACK_LIBRARIES = ${LAPACK_LIBRARIES}")
else()
  add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/BLAS") 
  add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/LAPACK")
  set(LAPACK_LIBRARIES LAPACK)
endif()

# Log Python status
if(PYTHON_FOUND)
  message("Python_FOUND:${Python_FOUND}")
  message("Python_LIBRARIES:${Python_LIBRARIES}")
  message("Python_INCLUDES:${Python_INCLUDE_DIRS}")  
else()
  message(STATUS "PYTHON NOT FOUND")
endif()

# Handle MPI-dependent packages
if(MPI_FOUND)
  message(STATUS "MPI was found.")
  message(STATUS "MPI was found .. path added ${MPI_C_INCLUDE_DIRS} ${MPI_FOUND}")
  include_directories(SYSTEM ${MPI_C_INCLUDE_DIRS})
  
  # MUMPS handling
  if(NOT DEFINED MUMPS_DIR)
    set(MUMPS_FLAG -D_NOMUMPS)
    set(MUMPS_LIBRARIES "")
  else()
    set(MUMPS_FLAG -D_MUMPS)  
    if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
      set(MUMPS_LIBRARIES "${MUMPS_DIR}/dmumps;${MUMPS_DIR}/mumps_common;${MUMPS_DIR}/pord")
    else()
      set(MUMPS_LIBRARIES "${MUMPS_DIR}/libdmumps.a;${MUMPS_DIR}/libmumps_common.a;${MUMPS_DIR}/libpord.a")
    endif()
  endif()

  add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/SuperLU_DIST_4.3/SRC")
  add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/METIS")    
else()
  message(STATUS "MPI was NOT found.")
endif()

# MKL handling
if(MKL_FOUND)
  message(STATUS "MKL was found.")
  if(NOT DEFINED SCALAPACK_LIBRARIES)
    set(SCALAPACK_LIBRARIES ${MKL_LIBRARIES})
  endif()
  if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(MKL_LPATH ${MKL_ROOT}/lib/intel64)
    set(SCALAPACK_LIBRARIES "${MKL_LPATH}/mkl_scalapack_ilp64.lib;${MKL_LPATH}/mkl_intel_ilp64.lib;${MKL_LPATH}/mkl_sequential.lib;${MKL_LPATH}/mkl_core.lib;${MKL_LPATH}/mkl_blacs_intelmpi_ilp64.lib")
  endif() 
else()
  message(STATUS "MKL NOT found .. user to provide -DSCALAPACK_LIBRARIES=")
endif()

message(STATUS "SCALAPACK_LIBRARIES=${SCALAPACK_LIBRARIES}")