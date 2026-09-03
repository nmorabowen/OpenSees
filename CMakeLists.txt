#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#     (c) Copyright 1999-2021 The Regents of the University of California
#                             All Rights Reserved
# (Copyright and Disclaimer @ http://www.berkeley.edu/OpenSees/copyright.html)
#
#------------------------------------------------------------------------------
#
# Major sections:
# - General Setup
# - External Libreries (Conan or User supplied paths)
# - Libraries
# - Executables
# - Option Application
#
#==============================================================================

cmake_minimum_required(VERSION 3.16)
project(OpenSeesFramework LANGUAGES C CXX Fortran)


set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(OPS_EXTERN_SOURCE_DIR "${PROJECT_SOURCE_DIR}/OTHER/")
set(OPS_EXTERNALS_DIR "${PROJECT_SOURCE_DIR}/OTHER/")
set(OPS_BUNDLED_DIR "${PROJECT_SOURCE_DIR}/OTHER/")
set(OPS_SRC_DIR "${PROJECT_SOURCE_DIR}/SRC/")

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake/cmake")
include(OpenSeesFunctions)

# include user config
include(${PROJECT_SOURCE_DIR}/Conf.cmake)

#
# for ubuntu this has to be before the Compiler SPECIFIC stuff below .. go figure!
#

find_package(MPI)



#cmake -DCMAKE_RULE_MESSAGES:BOOL=OFF -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON .

#
# Compiler SPECIFIC
#



if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
 MESSAGE("STATUS" "COMPILER: Clang (macOS Configuration)")

 set(CMAKE_OSX_DEPLOYMENT_TARGET "13.0")
 set(BUILD_SHARED_LIBS OFF)

 # 1. Ask gfortran for the absolute path to the static (.a) files
 execute_process(COMMAND gfortran --print-file-name=libgfortran.a OUTPUT_VARIABLE GF_STATIC OUTPUT_STRIP_TRAILING_WHITESPACE)
 execute_process(COMMAND gfortran --print-file-name=libquadmath.a OUTPUT_VARIABLE QM_STATIC OUTPUT_STRIP_TRAILING_WHITESPACE)
 execute_process(COMMAND gfortran --print-file-name=libgcc.a OUTPUT_VARIABLE GCC_STATIC OUTPUT_STRIP_TRAILING_WHITESPACE)

 # 2. CLEAR THE IMPLICIT DYNAMIC FLAGS
 # This prevents CMake from adding '-lgfortran -lquadmath' to the end of the command
 set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "")

 # 3. Force the linker to use these specific archive files
 # link_libraries applies to all targets created after this line
 link_libraries("${GF_STATIC}" "${QM_STATIC}" "${GCC_STATIC}")
    
 # Optional: ensure we search for static libs first for other dependencies
 set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.dylib")
 # set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
 set(BUILD_SHARED_LIBS OFF)
 # set(CMAKE_EXE_LINKER_FLAGS "-static")

elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")

 # set(CMAKE_FIND_LIBRARY_SUFFIXES ".so")
 set(BUILD_SHARED_LIBS OFF)
 # set(CMAKE_EXE_LINKER_FLAGS "-static")
 add_compile_options(-fPIC)
 set(CMAKE_POSITION_INDEPENDENT_CODE ON)
 set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
 MESSAGE("COMPILER: GNU")

elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
 MESSAGE("COMPILER: Intel")
 add_compile_options(-fPIC)

elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
 MESSAGE("COMPILER: MSVC")
 set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
 set(BUILD_SHARED_LIBS OFF)
 set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded")
 # ASDPlasticMaterial3D's OPS_AllASDPlasticMaterial3Ds.cpp instantiates a very
 # large product of templated material specializations (every EL x YF x PF x
 # hardening combo, including StiffSoil + HoekBrown). MSVC's COFF format caps
 # sections at ~64K; this file overflows. /bigobj lifts it to ~4G. The flag is
 # harmless for files that don't need it, so set it globally for MSVC.
 add_compile_options(/bigobj)
else()
  MESSAGE("COMPILER: UKNOWN COMPILER ${CMAKE_CXX_COMPILER_ID}")
endif()


#
# external packages
#

add_library(OPS_External_packages INTERFACE)
add_library(OPS_OS_Specific_libs INTERFACE)

#
# Conan (manages external dependencies for a typical build)
#

MESSAGE("${CMAKE_BINARY_DIR}")


if(CONAN_EXPORTED)
    message(STATUS "USING CONAN (via Transparent Integration)")
    set(USING_CONAN TRUE)
else()
    message(STATUS "NOT USING CONAN (System Search)")
    include(OpenSeesFunctions) # Your custom logic
endif()

find_package(ZLIB REQUIRED)
find_package(HDF5 REQUIRED)
find_package(TCL REQUIRED)

if(TARGET tcl::tcl)
	  
    message(STATUS "Found tcl::tcl target")
    include_directories(${TCL_INCLUDE_DIR})
    if (TCL_INCLUDE_DIR)
        find_file(TCL_INIT_FILE init.tcl
    	    PATHS 
            "${TCL_INCLUDE_DIR}/../lib"
            "${TCL_INCLUDE_DIR}/../lib/tcl8.6"
            NO_DEFAULT_PATH
    	)
     endif()

    if(TCL_INIT_FILE)
        message(STATUS "Confirmed TCL_INIT_FILE at: ${TCL_INIT_FILE}")
    else()
        message(WARNING "init.tcl not found automatically. Check your Tcl installation structure.")
    endif()

  else()
  
    # This is where Stampede3 usually lands. 
    # We manually create the lowercase target so the linker is happy.
    if(TCL_FOUND AND NOT TARGET tcl::tcl)
        add_library(tcl::tcl UNKNOWN IMPORTED)
        set_target_properties(tcl::tcl PROPERTIES
            IMPORTED_LOCATION "${TCL_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${TCL_INCLUDE_PATH}"
        )
    endif()

    set(TCL_LIBRARIES tcl::tcl)
    include_directories(${TCL_INCLUDE_PATH})
    
    # Use the actual library path found to find init.tcl
    get_filename_component(TCL_LIB_DIR "${TCL_LIBRARY}" DIRECTORY)
    get_filename_component(TCL_INIT_FILE "${TCL_LIB_DIR}/tcl8.6/init.tcl" ABSOLUTE)

    message(STATUS "2 TCL_INIT_FILE at ${TCL_INIT_FILE}")
endif()

message(STATUS "3 TCL_INIT_FILE at ${TCL_INIT_FILE}")


find_package(Eigen3 REQUIRED)

if(TARGET Eigen3::Eigen)
    # This pulls the include directory directly from the target object
    get_target_property(EIGEN_INC Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
    include_directories(${EIGEN_INC})
endif()
 
#if(Eigen3_FOUND)
#$   message(STATUS "OPS >>> FOUND EIGEN at ${Eigen3_INCLUDE_DIR}")  
#   include_directories(${Eigen3_INCLUDE_DIR})
#   add_compile_definitions(_EIGEN3)
#else()
#   message(STATUS "OPS >>> Could not find Eigen3")
#endif()


#
# add git log tag
#

find_package(Git QUIET)

if(GIT_FOUND)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} log -1 --format=%H
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
else()
    set(GIT_VERSION "unknown")
endif()

# Pass this to your code as a preprocessor definition.
# Ladruno: scope to C/C++ ONLY. Passing the quoted 40-char git-SHA define to Fortran
# (ifx /fpp) via a response file overruns an internal ifx buffer and crashes xfortcom
# (error #10273, 0xC0000005 / 0xC0000409). Fortran never uses OPENSEES_VERSION. The
# crash only surfaces when response files are forced (CMAKE_NINJA_FORCE_RESPONSE_FILE)
# and/or long build paths inflate the rsp. See Ladruno_implementation/06_profiler.md.
add_compile_definitions($<$<COMPILE_LANGUAGE:C,CXX>:OPENSEES_VERSION="${GIT_VERSION}">)


#
# define user-selectable options for end target
#

set_property(CACHE OPS_FINAL_TARGET PROPERTY STRINGS
    G3 OpenSees OpenSeesMP OpenSeesSP OpenSeesPy OpenSeesPyMP)

# Ladruno ADR 1000: the hierarchical CMS eigensolver is an opt-in library.
# OFF must leave every stock/serial target byte-identical at the source-list level.
option(LADRUNO_CMS "Build the independent two-level LadrunoCMS eigensolver" OFF)
option(LADRUNO_CMS_BUILD_TESTS "Build standalone LadrunoCMS C++ numerical checks" OFF)

# Ladruno ADR-87 D2: mutation gates. A mutation build SABOTAGES one physics
# family so CI can prove that family's tests actually test the physics (they
# must turn RED). NEVER ship or benchmark a binary built this way -- ask it
# what it is with the `ladrunoMutation` verb, which reports "none" for a
# normal build. Driver: Ladruno_scripts/mutation_gate.py.
#
#   cmake -DLADRUNO_MUTATE_FAMILY=CONTINUUM -DLADRUNO_MUTATE_MODE=ZERO ...
#
# Empty family (the default) adds no definition at all, so a normal build is
# untouched -- the header then defaults every family to LADRUNO_MUT_NONE.
set(LADRUNO_MUTATE_FAMILY "" CACHE STRING
    "ADR-87 D2 mutation gate: physics family to sabotage (CONTINUUM, UP, CONTACT, EXPLICIT, SANISAND). Empty = normal build.")
set(LADRUNO_MUTATE_MODE "ZERO" CACHE STRING
    "ADR-87 D2 mutation mode: ZERO (force->0), SCALE (force*1.5), IDENT (tangent->I).")
set_property(CACHE LADRUNO_MUTATE_FAMILY PROPERTY STRINGS "" CONTINUUM UP CONTACT EXPLICIT SANISAND)
set_property(CACHE LADRUNO_MUTATE_MODE PROPERTY STRINGS ZERO SCALE IDENT)

if(NOT LADRUNO_MUTATE_FAMILY STREQUAL "")
  set(_ladruno_valid_families CONTINUUM UP CONTACT EXPLICIT SANISAND)
  if(NOT LADRUNO_MUTATE_FAMILY IN_LIST _ladruno_valid_families)
    message(FATAL_ERROR
      "LADRUNO_MUTATE_FAMILY='${LADRUNO_MUTATE_FAMILY}' is not a known family. "
      "Known: ${_ladruno_valid_families}. A typo must not silently produce an "
      "UNMUTATED binary -- that would make the gate report 'tests cannot detect "
      "deleted physics' when nothing was ever deleted.")
  endif()
  if(LADRUNO_MUTATE_MODE STREQUAL "ZERO")
    set(_ladruno_mut_code 1)
  elseif(LADRUNO_MUTATE_MODE STREQUAL "SCALE")
    set(_ladruno_mut_code 2)
  elseif(LADRUNO_MUTATE_MODE STREQUAL "IDENT")
    set(_ladruno_mut_code 3)
  else()
    message(FATAL_ERROR "LADRUNO_MUTATE_MODE='${LADRUNO_MUTATE_MODE}' unknown; use ZERO, SCALE or IDENT.")
  endif()
  add_compile_definitions(LADRUNO_MUTATE_${LADRUNO_MUTATE_FAMILY}=${_ladruno_mut_code})
  message(WARNING
    "ADR-87 D2 MUTATION BUILD: ${LADRUNO_MUTATE_FAMILY}=${LADRUNO_MUTATE_MODE}. "
    "This binary computes DELIBERATELY WRONG physics. Do not ship, install or "
    "benchmark it. `ladrunoMutation` will report it.")
endif()

#
# OS Configuration:
#  defines: compile flags & user supplied paths and includes if not using conan
#

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  add_compile_definitions(_LINUX _UNIX _TCL85)
  if (NOT_USING_CONAN)
    include(OpenSeesDependenciesUnix)
  endif()
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  set(OPS_Use_Graphics_Option None)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  add_compile_definitions(_LINUX _UNIX _TCL85)
  if (NOT_USING_CONAN)
    include(OpenSeesDependenciesUnix)
  endif()
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    target_link_libraries(OPS_OS_Specific_libs INTERFACE wsock32 ws2_32)
    add_compile_definitions(_WIN32 _TCL85)
  if (NOT_USING_CONAN)
    include(OpenSeesDependenciesWin)
  endif()
endif()

#==============================================================================
#                            General Setup
#
#==============================================================================
#----------------------------------------------------------------
# Compilers
#----------------------------------------------------------------

# Fortran
#--------------------------------------
enable_language(Fortran)

# C++
#--------------------------------------

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD 11)
set(CMAKE_FORTRAN_STANDARD 08)

set(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> cqls <TARGET> <LINK_FLAGS> <OBJECTS>")
#add_compile_options(-g)

# Warnings
#opensees_add_cxx_flag(
#    GNU  -Wall
#    MSVC /W0
#)

# Floating-point
#opensees_add_cxx_flag(
#    GNU  -ffloat-store
#    MSVC /fp:precise
#)

#----------------------------------------------------------------
# Global Includes
#----------------------------------------------------------------
# include paths to main abstract classes

# NOTE BeamIntegration and MatrixUtil need to be removed from element/forceBEamColumn

include_directories(
  ${OPS_SRC_DIR}
  ${OPS_SRC_DIR}/matrix
  ${OPS_SRC_DIR}/handler
  ${OPS_SRC_DIR}/database
  ${OPS_SRC_DIR}/element
  ${OPS_SRC_DIR}/element/forceBeamColumn
  ${OPS_SRC_DIR}/element/nonlinearBeamColumn/matrixutil
  ${OPS_SRC_DIR}/coordTransformation
  ${OPS_SRC_DIR}/damping
  ${OPS_SRC_DIR}/tagged
  ${OPS_SRC_DIR}/tagged/storage
  ${OPS_SRC_DIR}/recorder
  ${OPS_SRC_DIR}/renderer
  ${OPS_SRC_DIR}/damage
  ${OPS_SRC_DIR}/recorder/response
  ${OPS_SRC_DIR}/material
  ${OPS_SRC_DIR}/material/section
  ${OPS_SRC_DIR}/material/uniaxial
  ${OPS_SRC_DIR}/material/nD
  ${OPS_SRC_DIR}/graph/graph
  ${OPS_SRC_DIR}/graph/numberer
  ${OPS_SRC_DIR}/graph/partitioner
  ${OPS_SRC_DIR}/domain/component
  ${OPS_SRC_DIR}/domain/domain
  ${OPS_SRC_DIR}/domain/subdomain
  ${OPS_SRC_DIR}/domain/load
  ${OPS_SRC_DIR}/domain/loadBalancer
  ${OPS_SRC_DIR}/domain/pattern
  ${OPS_SRC_DIR}/domain/groundMotion
  ${OPS_SRC_DIR}/domain/node
  ${OPS_SRC_DIR}/domain/constraints
  ${OPS_SRC_DIR}/domain/region
  ${OPS_SRC_DIR}/analysis/algorithm
  ${OPS_SRC_DIR}/analysis/dof_grp
  ${OPS_SRC_DIR}/analysis/fe_ele
  ${OPS_SRC_DIR}/analysis/algorithm/equiSolnAlgo
  ${OPS_SRC_DIR}/analysis/algorithm/eigenAlgo
  ${OPS_SRC_DIR}/analysis/algorithm/domainDecompAlgo
  ${OPS_SRC_DIR}/analysis/analysis
  ${OPS_SRC_DIR}/analysis/integrator
  ${OPS_SRC_DIR}/analysis/handler
  ${OPS_SRC_DIR}/analysis/numberer
  ${OPS_SRC_DIR}/analysis/model
  ${OPS_SRC_DIR}/convergenceTest
  ${OPS_SRC_DIR}/modelbuilder
  ${OPS_SRC_DIR}/system_of_eqn
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE
  ${OPS_SRC_DIR}/system_of_eqn/eigenSOE
  ${OPS_SRC_DIR}/actor/actor
  ${OPS_SRC_DIR}/actor/channel
  ${OPS_SRC_DIR}/actor/objectBroker
  ${OPS_SRC_DIR}/actor/message
  ${OPS_BUNDLED_DIR}/eigenAPI
)

# TODO: Temporary fix; include all directories with .h files through glob.
# The recursive glob produces ~200 directories, which makes the per-target
# include flag list >18 KB. The cmd.exe /C wrapper that Ninja uses for the
# Fortran preprocess rule has a hard ~8 KB command-line limit, so Fortran
# builds explode. Fortran sources don't need any C/C++ headers, so we wrap
# the include set in a generator expression that restricts it to C and CXX.
file(GLOB_RECURSE ops_include_files CONFIGURE_DEPENDS
    ${OPS_SRC_DIR}/*.h
    ${OPS_BUNDLED_DIR}/*.h
)
set(OPS_INCLUDE_DIRS "")
foreach(filepath ${ops_include_files})
    get_filename_component(dir_path ${filepath} PATH)
    set(OPS_INCLUDE_DIRS ${OPS_INCLUDE_DIRS} ${dir_path})
endforeach()
list(REMOVE_DUPLICATES OPS_INCLUDE_DIRS)

foreach(_d IN LISTS OPS_INCLUDE_DIRS)
    include_directories($<$<COMPILE_LANGUAGE:C,CXX>:${_d}>)
endforeach()
# end TODO

# ##############################################################
#
# NUMERICAL and OpenSees LIBRARIES
#
# ##############################################################

#
# NUMERICAL LIBRARIES
#

add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/SuperLU_5.1.1")
set (SUPERLU_LIBRARIES SUPERLU)
add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/UMFPACK")
set (UMFPACK_LIBRARIES UMFPACK)
add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/AMD")
set (AMD_LIBRARIES AMD)
add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/ARPACK")
set (ARPACK_LIBRARIES ARPACK)
add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/CSPARSE")
set (CSPARSE_LIBRARIES CSPARSE)
add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/ITPACK")
set (ITPACK_LIBRARIES ITPACK)
add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/tetgen1.4.3")
set (TETGEN_LIBRARIES tet)
add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/Triangle")
set (TRIANGLE_LIBRARIES triangle)
# add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/eigenAPI")
# set (EIGENAPI_LIBRARIES EIGENAPI)


if (NOT DEFINED LAPACK_LIBRARIES)
   find_package(LAPACK)
else()
   set (LAPACK_FOUND TRUE CACHE BOOL "providing own lapack lib")
endif()

if (NOT DEFINED Python_LIBRARIES)
   find_package(Python COMPONENTS Interpreter Development)
else()
   set (PYTHON_FOUND TRUE CACHE BOOL "providing own python lib")
endif()


find_package(MPI)
find_package(MKL)

# --- Workaround: on Windows, MPI/MKL/compiler install under "C:\Program Files (x86)\..."
# CMake's Ninja generator emits broken response-file quoting for paths with spaces+parens
# when used as Fortran include flags. Convert MPI include dirs to 8.3 short names.
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    foreach(_var MPI_C_INCLUDE_DIRS MPI_CXX_INCLUDE_DIRS MPI_Fortran_INCLUDE_DIRS
                 MPI_C_INCLUDE_PATH MPI_CXX_INCLUDE_PATH MPI_Fortran_INCLUDE_PATH)
        if(DEFINED ${_var})
            set(_short_list "")
            foreach(_p IN LISTS ${_var})
                if(EXISTS "${_p}")
                    file(TO_NATIVE_PATH "${_p}" _native)
                    execute_process(
                        COMMAND cmd /c for %I in ("${_native}") do @echo %~sI
                        OUTPUT_VARIABLE _short
                        OUTPUT_STRIP_TRAILING_WHITESPACE)
                    if(_short AND EXISTS "${_short}")
                        file(TO_CMAKE_PATH "${_short}" _short)
                        list(APPEND _short_list "${_short}")
                    else()
                        list(APPEND _short_list "${_p}")
                    endif()
                else()
                    list(APPEND _short_list "${_p}")
                endif()
            endforeach()
            set(${_var} "${_short_list}")
        endif()
    endforeach()
endif()

#set (LAPACK_FOUND FALSE)
#set (BLAS_FOUND FALSE)

if(LAPACK_FOUND)
  message(STATUS "LAPACK was found.")
  message(STATUS "LAPACK_LINKER_FLAGS = ${LAPACK_LINKER_FLAGS}")
  message(STATUS "LAPACK_LIBRARIES = ${LAPACK_LIBRARIES}" )
else()
  add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/BLAS")
  add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/LAPACK")
  set(LAPACK_LIBRARIES LAPACK)
endif()

if(PYTHON_FOUND)
  message("Python_FOUND:${Python_FOUND}")
  #message("Python_VERSION:${Python_VERSION}")
  #message("Python_Development_FOUND:${Python_Development_FOUND}")
  message("Python_LIBRARIES:${Python_LIBRARIES}")
  message("Python_INCLUDES:${Python_INCLUDE_DIRS}")
else()
  message(STATUS "PYTHON NOT FOUND")
endif()


if(MPI_FOUND)
  message(STATUS "MPI was found.")
  message(STATUS "MPI was found .. path added ${MPI_C_INCLUDE_DIRS} ${MPI_FOUND}")
  include_directories(SYSTEM ${MPI_C_INCLUDE_DIRS})
  #add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/MUMPS_5.4.1")
  #if(NOT CMAKE_SYSTEM_NAME STREQUAL "Windows")
  if (NOT DEFINED MUMPS_DIR)
  #    include(ExternalProject)
  #    ExternalProject_Add(mumps
  #       GIT_REPOSITORY    https://github.com/scivision/mumps.git
  #       GIT_TAG           v5.5.1.5
  #       SOURCE_DIR        "${CMAKE_BINARY_DIR}/mumps-src"
  #       BINARY_DIR        "${CMAKE_BINARY_DIR}"
  #       INSTALL_COMMAND   ""
  #     )
  #
  #    set (MUMPS_LIBRARIES "${CMAKE_BINARY_DIR}/libdmumps.lib;${CMAKE_BINARY_DIR}/libmumps_common.lib;${CMAKE_BINARY_DIR}/libpord.lib")
     set (MUMPS_FLAG -D_NOMUMPS)
     set (MUMPS_LIBRARIES "")
  else()
     set (MUMPS_FLAG -D_MUMPS)
     if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
       set (MUMPS_LIBRARIES "${MUMPS_DIR}/dmumps.lib;${MUMPS_DIR}/mumps_common.lib;${MUMPS_DIR}/pord.lib")
     else()
       set (MUMPS_LIBRARIES "${MUMPS_DIR}/libdmumps.a;${MUMPS_DIR}/libmumps_common.a;${MUMPS_DIR}/libpord.a")
     endif()

  endif()

  add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/SuperLU_DIST_4.3/SRC")
  add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/METIS")
else()
  message(STATUS "MPI was NOT found.")
endif()


if(MKL_FOUND)
  message(STATUS "MKL was found.")
  if (NOT DEFINED SCALAPACK_LIBRARIES)
     set (SCALAPACK_LIBRARIES ${MKL_LIBRARIES})
  endif()
  if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
     if(EXISTS "${MKL_ROOT}/lib/mkl_core.lib")
      set(MKL_LPATH "${MKL_ROOT}/lib")
    elseif(EXISTS "${MKL_ROOT}/lib/intel64/mkl_core.lib")
      set(MKL_LPATH "${MKL_ROOT}/lib/intel64")
    else()
      message(FATAL_ERROR "MKL libraries not found in expected locations")
    endif()

    # LP64 (32-bit MPI integers) to match a default-built MUMPS, which uses
    # the standard 32-bit Fortran INTEGER and Intel MPI's default LP64 layer.
    # If you rebuild MUMPS with -Dintsize64=ON, switch every "lp64" back to
    # "ilp64" and ensure Intel MPI is in ILP64 mode.
    set (SCALAPACK_LIBRARIES "${MKL_LPATH}/mkl_scalapack_lp64.lib;${MKL_LPATH}/mkl_intel_lp64.lib;${MKL_LPATH}/mkl_sequential.lib;${MKL_LPATH}/mkl_core.lib;${MKL_LPATH}/mkl_blacs_intelmpi_lp64.lib")
  endif()

else()
  message(STATUS "MKL NOT found .. user to provide -DSCALAPACK_LIBRARIES=")

endif()

message(STATUS "SCALAPACK_LIBRARIES=${SCALAPACK_LIBRARIES}")

# ---------------------------------------------------------------------------
# Ladruno ADR-75 P1b: MKL PARDISO for the DESKTOP (serial) targets.
#
# ADR-75 measured the gap: the serial build has no threaded sparse-direct
# solver at all (Lane-B baseline -- UmfPack 22.711 s, single-threaded; `system
# Pardiso` and `system Mumps` both report "unknown system type"). MUMPS stays
# the CLUSTER solver and is untouched here.
#
# THREADING LAYER (the reason this is a separate variable from
# SCALAPACK_LIBRARIES): SCALAPACK_LIBRARIES deliberately uses mkl_sequential,
# because the parallel targets get their parallelism from MPI. Desktop PARDISO
# is the opposite case -- its entire value is shared-memory threading -- so the
# serial targets link mkl_intel_thread + libiomp5md instead. Mixing the two
# layers inside ONE target/process is unsupported by MKL, so they are kept on
# disjoint targets: threaded here (serial only), sequential there (parallel).
# ---------------------------------------------------------------------------
set(PARDISO_FLAG "")
set(PARDISO_MKL_LIBRARIES "")
if(MKL_FOUND AND DEFINED MKL_LPATH)
  set(PARDISO_FLAG -D_PARDISO)
  set(PARDISO_MKL_LIBRARIES
      "${MKL_LPATH}/mkl_intel_lp64.lib;${MKL_LPATH}/mkl_intel_thread.lib;${MKL_LPATH}/mkl_core.lib")
  # OpenMP runtime for mkl_intel_thread (oneAPI ships libiomp5md next to the
  # compiler); harmless if the linker already resolves it.
  find_library(IOMP5_LIB NAMES libiomp5md iomp5md iomp5
               HINTS "${MKL_ROOT}/../compiler/latest/lib"
                     "${MKL_ROOT}/../compiler/latest/windows/compiler/lib/intel64_win"
                     ENV LIB)
  if(IOMP5_LIB)
    list(APPEND PARDISO_MKL_LIBRARIES ${IOMP5_LIB})
  endif()
  message(STATUS "PARDISO (desktop, threaded MKL) ENABLED: ${PARDISO_MKL_LIBRARIES}")
else()
  message(STATUS "PARDISO DISABLED (no MKL) - serial builds fall back to UmfPack")
endif()

# ---------------------------------------------------------------------------
# Ladruno ADR43 L2-profile: opt-in MKL FEAST on a non-Windows host (esmeralda).
# The FEAST/PARDISO paths in FeastEigenSolver.cpp + LadrunoBlockZKernel.cpp are
# historically gated on _WIN32 (MKL was Windows/oneAPI-only on this fork). With
# LADRUNO_MKL_FEAST_LINUX=ON on a non-Windows host that HAS MKL, define
# _LADRUNO_MKL_FEAST (flips those guards on) and link the explicit LP64/
# sequential/core 3-layer (carrying dfeast_srci/dfeast_scsrgv/feastinit/pardiso)
# via LADRUNO_MKL_FEAST_LIBS for the parallel targets — NOT the single libmkl_rt
# dispatcher, whose runtime kernel selection is fragile (see the find block).
# Default OFF => the Windows/oneAPI and Ubuntu Zone-A reference-LAPACK builds
# are byte-identical.
# ---------------------------------------------------------------------------
option(LADRUNO_MKL_FEAST_LINUX "Enable MKL FEAST on a non-Windows host" OFF)
set(LADRUNO_MKL_FEAST_LIBS "")
if(LADRUNO_MKL_FEAST_LINUX AND NOT WIN32)
  # Prefer the explicit LP64/sequential/core 3-layer (the canonical robust MKL
  # link) over the single libmkl_rt dispatcher, which does fragile runtime
  # kernel selection (observed on esmeralda: apt libmkl_rt fails to load its
  # libmkl_def kernel). MKL_RT_HINT points at the MKL lib dir (e.g. a spack
  # oneAPI mkl/<ver>/lib). LP64 => MKL_INT==int, matching the FEAST code.
  find_library(LADRUNO_MKL_LP64 NAMES mkl_intel_lp64
    HINTS ${MKL_RT_HINT} ENV MKLROOT PATHS /usr/lib/x86_64-linux-gnu
    PATH_SUFFIXES lib lib/intel64)
  find_library(LADRUNO_MKL_SEQ NAMES mkl_sequential
    HINTS ${MKL_RT_HINT} ENV MKLROOT PATHS /usr/lib/x86_64-linux-gnu
    PATH_SUFFIXES lib lib/intel64)
  find_library(LADRUNO_MKL_CORE NAMES mkl_core
    HINTS ${MKL_RT_HINT} ENV MKLROOT PATHS /usr/lib/x86_64-linux-gnu
    PATH_SUFFIXES lib lib/intel64)
  if(LADRUNO_MKL_LP64 AND LADRUNO_MKL_SEQ AND LADRUNO_MKL_CORE)
    # Macro + libs are attached to the OPS_SysOfEqn object library once it is
    # created (below), NOT globally: the FEAST/PARDISO sources
    # (FeastEigenSolver.cpp, LadrunoBlockZKernel.cpp) live ONLY in OPS_SysOfEqn,
    # so a dir-global add_compile_definitions would make every target compile the
    # MKL path while only the parallel exes carried the libs -> the serial/Py/PyMP
    # targets fail to link (adversarial-review CRITICAL). PUBLIC libs on the object
    # library propagate MKL to every consumer (OpenSeesLIB -> all exes) in lockstep.
    set(LADRUNO_MKL_FEAST_ON TRUE)
    # --start-group resolves the circular deps among the three MKL layers.
    set(LADRUNO_MKL_FEAST_LIBS
        -Wl,--start-group ${LADRUNO_MKL_LP64} ${LADRUNO_MKL_SEQ} ${LADRUNO_MKL_CORE} -Wl,--end-group
        pthread m dl)
    message(STATUS "Ladruno: MKL FEAST on Linux ENABLED (3-layer: ${LADRUNO_MKL_LP64})")
  else()
    message(FATAL_ERROR "LADRUNO_MKL_FEAST_LINUX=ON but the MKL 3-layer "
      "(mkl_intel_lp64/mkl_sequential/mkl_core) was not all found; "
      "pass -DMKL_RT_HINT=<MKL lib dir>")
  endif()
endif()

#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#  message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()



add_library(OPS_Numerics INTERFACE)

target_link_libraries(OPS_Numerics INTERFACE
       ${ARPACK_LIBRARIES}
       ${CSPARSE_LIBRARIES}
       ${SUPERLU_LIBRARIES}
       ${UMFPACK_LIBRARIES}
       ${TETGEN_LIBRARIES}
       ${TRIANGLE_LIBRARIES}
       ${AMD_LIBRARIES}
       ${ITPACK_LIBRARIES}
       ${LAPACK_LIBRARIES}
       ${EIGENAPI_LIBRARIES}
)

set(TCL_LIBRARIES ${TCL_LIBRARY})


#
# OpenSees OBJECT LIBRARIES
#

add_library(OPS_Matrix             OBJECT)
add_library(OPS_Actor              OBJECT)
add_library(OPS_ObjectBroker       OBJECT)
add_library(OPS_Handler            OBJECT)
add_library(OPS_Recorder           OBJECT)
add_library(OPS_Reliability        OBJECT)
add_library(OPS_Tagged             OBJECT)
add_library(OPS_Utilities          OBJECT)
add_library(OPS_ModelBuilder       OBJECT)
add_library(OPS_Domain             OBJECT)
add_library(OPS_SysOfEqn           OBJECT)
add_library(OPS_SysOfEqnPython     OBJECT)

# Ladruno ADR 1000: declare the target before final MP targets are linked;
# SRC/system_of_eqn/ladrunoCMS supplies its sources later in the configure pass.
if(LADRUNO_CMS)
  if(NOT MPI_FOUND)
    message(FATAL_ERROR "LADRUNO_CMS=ON requires MPI")
  endif()
  if(NOT DEFINED MUMPS_DIR)
    message(FATAL_ERROR "LADRUNO_CMS=ON requires MUMPS_DIR")
  endif()
  add_library(OPS_LadrunoCMS STATIC)
  target_compile_definitions(OPS_LadrunoCMS PUBLIC _LADRUNO_CMS)
endif()

# Ladruno ADR43 L2: attach the MKL FEAST macro + libs to the ONE object library
# that compiles the FEAST/PARDISO sources (FeastEigenSolver.cpp,
# LadrunoBlockZKernel.cpp live only in OPS_SysOfEqn). PRIVATE macro => only these
# sources see the MKL path; PUBLIC libs => every consumer (OpenSeesLIB -> all
# exes + the openseesmp module) links MKL in lockstep, so there is no
# macro-global/libs-partial skew (adversarial-review CRITICAL). Gated OFF by
# default (LADRUNO_MKL_FEAST_ON set only when the opt-in flag + MKL are present).
if(LADRUNO_MKL_FEAST_ON)
  target_compile_definitions(OPS_SysOfEqn PRIVATE _LADRUNO_MKL_FEAST)
  target_link_libraries(OPS_SysOfEqn PUBLIC ${LADRUNO_MKL_FEAST_LIBS})
endif()

# ---------------------------------------------------------------------------
# Ladruno (TIMs, 2026-09-03): opt-in MKL PARDISO on a non-Windows host.
#
# ADR-75 P1b enables `system Pardiso` only where `MKL_LPATH` is set, and
# `MKL_LPATH` is assigned only inside the Windows branch from `.lib` file
# names -- so no MKL on Linux could ever turn it on, and the serial Linux
# build fell back to UmfPack. Measured on esmeralda (build e95a1c74f, TIMs
# act work/ape/response-curve-matrix): UmfPack's 32-bit index path returns
# `numeric analysis returns -1` on a 93 246-DOF LadrunoUP model with 28 GB
# free, its 64-bit path (`useLongIndices`) never runs the symbolic analysis
# in this build, and SuperLU took >18x UmfPack's whole gravity stage on a
# 49 626-DOF mesh without finishing one solve. PARDISO handles the same
# sizes on the Windows desktop. This block mirrors the LADRUNO_MKL_FEAST_LINUX
# pattern above: explicit LP64/sequential/core 3-layer found through
# MKL_RT_HINT, never the libmkl_rt dispatcher.
#
# SEQUENTIAL layer on purpose: the Linux targets already link the sequential
# MKL layer for BLAS/LAPACK/ScaLAPACK, and one process must not carry two
# threading layers (MKL's own rule, restated at ADR-75 P1b above). A threaded
# Linux PARDISO is a later, separate decision that must move the whole link
# to mkl_gnu_thread/mkl_intel_thread at once. Default OFF => every existing
# build is byte-identical.
# ---------------------------------------------------------------------------
option(LADRUNO_MKL_PARDISO_LINUX "Enable MKL PARDISO (sequential layer) on a non-Windows host" OFF)
if(LADRUNO_MKL_PARDISO_LINUX AND NOT WIN32)
  find_library(LADRUNO_MKL_LP64 NAMES mkl_intel_lp64
    HINTS ${MKL_RT_HINT} ENV MKLROOT PATHS /usr/lib/x86_64-linux-gnu
    PATH_SUFFIXES lib lib/intel64)
  find_library(LADRUNO_MKL_SEQ NAMES mkl_sequential
    HINTS ${MKL_RT_HINT} ENV MKLROOT PATHS /usr/lib/x86_64-linux-gnu
    PATH_SUFFIXES lib lib/intel64)
  find_library(LADRUNO_MKL_CORE NAMES mkl_core
    HINTS ${MKL_RT_HINT} ENV MKLROOT PATHS /usr/lib/x86_64-linux-gnu
    PATH_SUFFIXES lib lib/intel64)
  # mkl_pardiso.h / mkl_types.h: the spack oneAPI layout is <mkl>/lib beside
  # <mkl>/include; the apt layout is /usr/include/mkl.
  find_path(LADRUNO_MKL_INCLUDE NAMES mkl_pardiso.h
    HINTS ${MKL_RT_HINT}/../include ${MKL_RT_HINT}/../../include ENV MKLROOT
    PATHS /usr/include/mkl /usr/include
    PATH_SUFFIXES include)
  if(LADRUNO_MKL_LP64 AND LADRUNO_MKL_SEQ AND LADRUNO_MKL_CORE AND LADRUNO_MKL_INCLUDE)
    set(LADRUNO_MKL_PARDISO_ON TRUE)
    set(PARDISO_FLAG -D_PARDISO)
    set(PARDISO_MKL_LIBRARIES
        -Wl,--start-group ${LADRUNO_MKL_LP64} ${LADRUNO_MKL_SEQ} ${LADRUNO_MKL_CORE} -Wl,--end-group
        pthread m dl)
    # The PARDISO SOE/solver sources compile in OPS_SysOfEqn (see
    # SRC/system_of_eqn/linearSOE/CMakeLists.txt); they need the MKL headers,
    # and PUBLIC libs on the object library reach every consumer in lockstep,
    # exactly as the FEAST block does.
    target_include_directories(OPS_SysOfEqn PRIVATE ${LADRUNO_MKL_INCLUDE})
    target_link_libraries(OPS_SysOfEqn PUBLIC ${PARDISO_MKL_LIBRARIES})
    message(STATUS "Ladruno: MKL PARDISO on Linux ENABLED (sequential 3-layer: ${LADRUNO_MKL_LP64}; headers ${LADRUNO_MKL_INCLUDE})")
  else()
    message(FATAL_ERROR "LADRUNO_MKL_PARDISO_LINUX=ON but the MKL 3-layer "
      "(mkl_intel_lp64/mkl_sequential/mkl_core) or mkl_pardiso.h was not found; "
      "pass -DMKL_RT_HINT=<MKL lib dir> (headers are looked for beside it)")
  endif()
endif()
#add_library(OPS_SysOfEqnSP        OBJECT)
add_library(OPS_Analysis           OBJECT)
add_library(OPS_ConvergenceTest    OBJECT)
add_library(OPS_Thermal            OBJECT)
add_library(OPS_Element            OBJECT)
add_library(OPS_ElementFortran     OBJECT)
#add_library(OPS_ElementSP         OBJECT)
add_library(OPS_Material           OBJECT)
add_library(OPS_MaterialFortran    OBJECT)
add_library(OPS_Damage             OBJECT)
add_library(OPS_Database           OBJECT)
add_library(OPS_INTERPRETER        OBJECT)



#
# Optional Extensions
#

add_library(OPS_Paraview           OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Renderer           OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Graphics           OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Graphics_Default   OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Graphics_GL        OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_PFEM               OBJECT EXCLUDE_FROM_ALL)

#
# Tcl Interpreter Library & special include and link directives
#

# Two source files inside OPS_InterpTcl declare `theDomain` as either
# `Domain` or `PartitionedDomain` depending on _PARALLEL_PROCESSING /
# _PARALLEL_INTERPRETERS, and another file (commands.cpp) defines the
# matching global. A single STATIC OPS_InterpTcl compiles them once with
# the wrong type for SP/MP and produces unresolved-symbol link errors.
# Solution: keep OPS_InterpTcl STATIC for everything else, but pull the
# two parallel-sensitive sources out and let each exe compile them with
# its own compile definitions (see OPS_INTERP_TCL_PER_TARGET_SOURCES).
add_library(OPS_InterpTcl          STATIC)
set(OPS_INTERP_TCL_PER_TARGET_SOURCES
    ${OPS_SRC_DIR}/modelbuilder/tcl/myCommands.cpp
    ${OPS_SRC_DIR}/element/UWelements/Tcl_generateInterfacePoints.cpp
)

# ---------------------------------------------------------------------------
# Ladruno ADR-78: MPI-sensitive sources that must be compiled PER TARGET
#
# Same shape as Patch 8 above, and the same reason. Anything in the
# OPS_Analysis OBJECT library is compiled ONCE with no parallel define and
# folded into every target -- so an `#ifdef _PARALLEL_*` block written there
# compiles to nothing in EVERY build, OpenSeesMP included. Measured, not
# assumed: a guarded MPI_Abort placed in LadrunoContactHandler.cpp vanished
# silently and the mutation deck hung exactly as it had before the change.
#
# So each parallel-sensitive function is isolated in its own TU and compiled
# per-target with that target's defines. Keep this list trivial -- it is
# duplicated into every binary -- and keep the OPS_Analysis-side file free of
# `#ifdef _PARALLEL_*` entirely, so the trap cannot silently reopen.
#
#   LadrunoContactAbort.cpp     ADR-78 P1 -- the contact fatal exit, which must
#                               MPI_Abort rather than let a failing rank strand
#                               its peers in the next collective.
#   LadrunoAutoPenaltyReduce.cpp  ADR-78 follow-up -- the `constraints Auto`
#                               global-stiffness allreduce. It had sat dead in
#                               AutoConstraintHandler.cpp since that handler was
#                               contributed, so under MPI every rank sized its
#                               penalty from rank-local stiffness and the ranks
#                               sharing an interface constraint disagreed on it.
set(OPS_MPI_PER_TARGET_SOURCES
    ${OPS_SRC_DIR}/analysis/handler/LadrunoContactAbort.cpp
    ${OPS_SRC_DIR}/analysis/handler/LadrunoAutoPenaltyReduce.cpp
)

target_compile_definitions(OPS_InterpTcl PUBLIC _TCL85)
if(NOT EXISTS "${CMAKE_BINARY_DIR}/conanbuildinfo.cmake")
  target_include_directories(OPS_InterpTcl PUBLIC ${TCL_INCLUDE_PATH})
  target_link_libraries(OPS_InterpTcl PRIVATE ${TCL_LIBRARIES})
endif()

# ---------------------------------------------------------------------------
# Ladruno Patch 9: per-define Python command-engine static libs
#
# OpenSeesCommands.cpp / OpenSeesMiscCommands.cpp are the ONLY sources in
# SRC/interpreter carrying _PARALLEL_INTERPRETERS / _PARALLEL_PROCESSING
# guards (MPI_MachineBroker construction, parallel SOE selection, the
# getPID/getNP/barrier/send/recv commands). They normally live in the single
# OPS_INTERPRETER OBJECT lib -- compiled once with NO parallel define and
# folded into OpenSeesLIB -- so those parallel paths are compiled out of
# every target. A define on a consuming target cannot retroactively
# recompile a built OBJECT lib, and compiling these TUs both seq (shared)
# and parallel (per-target) into one link would collide (LNK2005). So,
# exactly as Patch 8 (OPS_INTERP_TCL_PER_TARGET_SOURCES) does for the Tcl
# side, the two sources are removed from OPS_INTERPRETER (see
# SRC/interpreter/CMakeLists.txt) and rebuilt here into two STATIC libs,
# one per define. Each final target links exactly one:
#   OPS_InterpPyCmds      no parallel define  -> Tcl OpenSees/SP/MP, G3,
#                          sequential OpenSeesPy (byte-identical to before)
#   OPS_InterpPyCmds_MP   _PARALLEL_INTERPRETERS -> OpenSeesPyMP
#
# A _PARALLEL_PROCESSING variant for a Python "openseessp" was prototyped
# then REMOVED: the Python/DL command engine has no SP domain-decomposition
# subsystem (SP is Tcl-only). See the definitive write-up in
# Ladruno_implementation/02_openseespymp.md / compilation journal Patch 9.
#
# The Tcl SP/MP exes deliberately link the SEQ variant: they drive
# parallelism through commands.cpp, never the Python command engine, and
# must stay byte-identical (linking _MP/_SP would double the
# MPI_MachineBroker setup). OpenSeesCommands.h declares theMachineBroker /
# theChannels / numChannels UNCONDITIONALLY, so the class ABI is
# define-independent -- OpenSeesLIB's seq-compiled TUs link cleanly against
# the parallel-compiled OpenSeesCommands.o. These libs are intentionally
# kept OUT of OpenSeesLIB so OpenSeesPyMP (which links OpenSeesLIB) does not
# also pull the seq variant.
set(OPS_INTERP_PY_PARALLEL_SOURCES
    ${OPS_SRC_DIR}/interpreter/OpenSeesCommands.cpp
    ${OPS_SRC_DIR}/interpreter/OpenSeesMiscCommands.cpp
)

add_library(OPS_InterpPyCmds STATIC ${OPS_INTERP_PY_PARALLEL_SOURCES})

# Ladruno ADR-75 P1b: `system Pardiso` registration lives in OpenSeesCommands.cpp
# under #ifdef _PARDISO. This lib is shared by the Tcl OpenSees/SP/MP exes AND
# sequential OpenSeesPy, so the flag is applied to every consumer -- the
# PARDISOGenLin* objects are in OpenSeesLIB (all of them link it) and the MKL
# PARDISO() implementation resolves from each target's own MKL layer (threaded
# for serial, the existing sequential ScaLAPACK layer for parallel).
if(PARDISO_FLAG)
    target_compile_definitions(OPS_InterpPyCmds PRIVATE ${PARDISO_FLAG})
    target_include_directories(OPS_InterpPyCmds PRIVATE
        ${OPS_SRC_DIR}/system_of_eqn/linearSOE/pardiso)
endif()

if(MPI_FOUND)
    add_library(OPS_InterpPyCmds_MP STATIC ${OPS_INTERP_PY_PARALLEL_SOURCES})

    foreach(_pylib OPS_InterpPyCmds_MP)
        if(DEFINED MUMPS_INCLUDE_DIR)
            target_include_directories(${_pylib} PRIVATE ${MUMPS_INCLUDE_DIR} ${MPI_CXX_INCLUDE_DIRS})
        else()
            target_include_directories(${_pylib} PRIVATE ${MUMPS_DIR}/_deps/mumps-src/include ${MPI_CXX_INCLUDE_DIRS})
        endif()
        target_compile_options(${_pylib} PRIVATE ${MPI_CXX_COMPILE_FLAGS})
        # Under the parallel define OpenSeesMiscCommands.cpp #includes <metis.h>
        # and OpenSeesCommands.cpp pulls the distributed SuperLU wrappers.
        # PRIVATE-link the same targets the Tcl SP/MP exes do so their
        # INTERFACE include dirs reach these TUs at compile time (consumers
        # such as OpenSeesPyMP link METIS/SUPERLU_DIST in their own right).
        target_link_libraries(${_pylib} PRIVATE METIS SUPERLU_DIST)
    endforeach()

    if (DEFINED OPENMPI)
        target_compile_definitions(OPS_InterpPyCmds_MP PRIVATE _PARALLEL_INTERPRETERS ${MUMPS_FLAG} _OPENMPI)
    else()
        target_compile_definitions(OPS_InterpPyCmds_MP PRIVATE _PARALLEL_INTERPRETERS ${MUMPS_FLAG})
    endif()
    if(LADRUNO_CMS)
        target_compile_definitions(OPS_InterpPyCmds_MP PRIVATE _LADRUNO_CMS)
    endif()
endif()

#
# OpenSees G3 Library - a slimmed down OpenSees
#

add_library(G3)

target_link_libraries(G3
    coordTransformation
    damping
    OPS_Matrix
    OPS_Analysis
    OPS_ModelBuilder
    OPS_Domain
    OPS_ConvergenceTest
    OPS_Element
    OPS_Material
    OPS_Recorder
    OPS_Handler
    OPS_SysOfEqn
    OPS_Tagged
    OPS_Utilities
    graph
    OPS_Actor
    OPS_ObjectBroker
    OPS_Numerics
    OPS_INTERPRETER
    OPS_InterpPyCmds   # Patch 9: parallel-sensitive Python cmds (seq variant)
)

#
# OpenSees Libray
#

add_library(OpenSeesLIB  EXCLUDE_FROM_ALL)

#
# For UMFPACK long int DLONG version -- Gustavo, Michael
#
target_compile_definitions(OpenSeesLIB PUBLIC _UMFPACK_DLONG)

target_link_libraries(OpenSeesLIB
    ${OPS_Extension_List}
    OPS_Actor
    OPS_Analysis
    OPS_ConvergenceTest
    OPS_Damage
    OPS_Database
    OPS_Domain
    graph
    OPS_Element
    OPS_ElementFortran
    OPS_Handler
    OPS_INTERPRETER
    OPS_Material
    OPS_MaterialFortran
    OPS_Material_YieldSurface
    OPS_Matrix
    OPS_ModelBuilder
    OPS_ObjectBroker
    OPS_PFEM
    OPS_Recorder
    coordTransformation
    damping
    OPS_Renderer
    OPS_Section_Repres
    OPS_Section_YieldSurface
    OPS_SysOfEqn
    OPS_Thermal
    OPS_OS_Specific_libs
    OPS_Tagged
    OPS_Utilities
)


# ########################################################################
#
# EXECUTABLES
#
# ########################################################################

# Adding EXCLUDE_FROM_ALL prevents these from being made by default.

#
# OpenSees Tcl Interpreter
#



add_executable(OpenSees EXCLUDE_FROM_ALL
  ${OPS_SRC_DIR}/tcl/tclAppInit.cpp
  ${OPS_SRC_DIR}/tcl/tclMain.cpp
  ${OPS_SRC_DIR}/tcl/commands.cpp
  ${OPS_SRC_DIR}/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp
  ${OPS_INTERP_TCL_PER_TARGET_SOURCES}
  ${OPS_MPI_PER_TARGET_SOURCES}   # Ladruno ADR-78 (contact abort + auto-penalty allreduce)
)

if(TARGET TCL::TCL)
    # --- CONAN 2 LOGIC ---
    message(STATUS "Conan TCL found at: ${tcl_PACKAGE_FOLDER}")
    
    # In Conan 2, the library files are typically in the 'lib' subfolder of the package
    set(CONAN_TCL_LIB_DIR "${tcl_PACKAGE_FOLDER}/lib")

    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib/tcl8.6)
    
    add_custom_command(
        TARGET OpenSees POST_BUILD
        COMMENT "Copying init.tcl from Conan TCL to binary directory"
        COMMAND ${CMAKE_COMMAND} -E copy
        "${CONAN_TCL_LIB_DIR}/tcl8.6/init.tcl"
        "${CMAKE_CURRENT_BINARY_DIR}/lib/tcl8.6"
    )
    
    # Note: You don't usually need target_include_directories here 
    # because target_link_libraries(OpenSees PRIVATE TCL::TCL) handles it.

else()
    # --- SYSTEM / NON-CONAN LOGIC ---
    message(STATUS "Using System TCL")
    
    target_include_directories(OpenSees PUBLIC ${TCL_INCLUDE_PATH})

    file(COPY ${TCL_INCLUDE_PATH}/../lib/tcl8.6  # Adjust path based on your system layout
         DESTINATION ${PROJECT_BINARY_DIR}/lib/tcl8.6
         FILES_MATCHING PATTERN "*.tcl")
endif()


target_link_libraries(OpenSees
  OPS_InterpTcl
  OPS_InterpPyCmds
  coordTransformation
  OpenSeesLIB
  OPS_Reliability
  OPS_ReliabilityTcl
  OPS_Numerics
  OPS_Recorder
  ${CMAKE_DL_LIBS}
  HDF5::HDF5
  tcl::tcl
  ZLIB::ZLIB
  Eigen3::Eigen
  ${PARDISO_MKL_LIBRARIES}   # Ladruno ADR-75 P1b (threaded MKL; empty w/o MKL)
)

# Ladruno ADR-75 P1d: `system Pardiso` in the TCL chain (SRC/tcl/commands.cpp,
# which is compiled straight into this exe -- it is NOT part of a shared lib).
#
# SERIAL exe only. CORRECTED RATIONALE (adversarial review): an earlier version
# of this comment said SP/MP are excluded because "mixing the threaded and
# sequential MKL layers in one process is unsupported". That is a non-sequitur
# -- OPS_InterpPyCmds is already compiled WITH ${PARDISO_FLAG} and linked into
# OpenSeesSP/MP, so those binaries already contain the PARDISO SOE/solver
# objects, each resolving against exactly one MKL layer. Nothing is mixed.
#
# The real reasons are: (1) SP/MP link the SEQUENTIAL MKL layer, so PARDISO's
# entire value -- shared-memory threading -- is absent there; and (2) a serial
# CSR SOE has no defined meaning against a PartitionedDomain, where the solve
# must be distributed. Cluster runs use `system Mumps`. commands.cpp carries an
# explicit #else that REFUSES `system Pardiso` in those builds rather than
# letting it fall through to a stale global SOE.
if(PARDISO_FLAG)
    target_compile_definitions(OpenSees PRIVATE ${PARDISO_FLAG})
    target_include_directories(OpenSees PRIVATE
        ${OPS_SRC_DIR}/system_of_eqn/linearSOE/pardiso)
endif()


#
# OpenSeesSP Tcl Parallel Interpreter
#


if(MPI_FOUND)

    add_executable(OpenSeesSP EXCLUDE_FROM_ALL
      ${OPS_SRC_DIR}/tcl/mpiMain.cpp
      ${OPS_SRC_DIR}/tcl/tclMain.cpp
      ${OPS_SRC_DIR}/tcl/commands.cpp
      ${OPS_SRC_DIR}/domain/domain/partitioned/PartitionedDomain.cpp
      ${OPS_SRC_DIR}/domain/domain/partitioned/PartitionedDomainEleIter.cpp
      ${OPS_SRC_DIR}/domain/domain/partitioned/PartitionedDomainSubIter.cpp
      ${OPS_SRC_DIR}/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp
      ${OPS_SRC_DIR}/actor/machineBroker/MPI_MachineBroker.cpp
      ${OPS_SRC_DIR}/system_of_eqn/linearSOE/diagonal/MPIDiagonalSOE.cpp
      ${OPS_SRC_DIR}/system_of_eqn/linearSOE/diagonal/MPIDiagonalSolver.cpp
      ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsSOE.cpp
      ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsSolver.cpp
      ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsParallelSOE.cpp
      ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.cpp
      ${OPS_SRC_DIR}/system_of_eqn/linearSOE/sparseGEN/DistributedSuperLU.cpp
      ${OPS_SRC_DIR}/system_of_eqn/linearSOE/sparseGEN/DistributedSparseGenColLinSOE.cpp
      ${OPS_SRC_DIR}/system_of_eqn/linearSOE/sparseGEN/DistributedSparseGenRowLinSOE.cpp
      ${OPS_SRC_DIR}/domain/subdomain/ActorSubdomain.cpp
      ${OPS_SRC_DIR}/domain/subdomain/ShadowSubdomain.cpp
      ${OPS_INTERP_TCL_PER_TARGET_SOURCES}
  ${OPS_MPI_PER_TARGET_SOURCES}   # Ladruno ADR-78 (contact abort + auto-penalty allreduce)
    )


    if(DEFINED MUMPS_INCLUDE_DIR)
        target_include_directories(OpenSeesSP PRIVATE ${MUMPS_INCLUDE_DIR} ${MPI_CXX_INCLUDE_DIRS})
    else()
        target_include_directories(OpenSeesSP PRIVATE ${MUMPS_DIR}/_deps/mumps-src/include ${MPI_CXX_INCLUDE_DIRS})
    endif()
    target_compile_options(OpenSeesSP PRIVATE ${MPI_CXX_COMPILE_FLAGS})

    if (DEFINED OPENMPI)
        target_compile_definitions(OpenSeesSP
        PUBLIC _PARALLEL_PROCESSING ${MUMPS_FLAG} _OPENMPI)
    else()
        target_compile_definitions(OpenSeesSP
        PUBLIC _PARALLEL_PROCESSING ${MUMPS_FLAG})
    endif()

    target_link_libraries(OpenSeesSP
       OPS_InterpTcl
       OPS_InterpPyCmds
       OpenSeesLIB
       OPS_Reliability
       OPS_ReliabilityTcl
       OPS_Recorder
       METIS
       SUPERLU_DIST
       OPS_Numerics
       ${MUMPS_LIBRARIES}
       ${CMAKE_DL_LIBS}
       HDF5::HDF5
       tcl::tcl
       ZLIB::ZLIB
       Eigen3::Eigen
       ${MPI_CXX_LIBRARIES}
       ${SCALAPACK_LIBRARIES}
       ${MPI_Fortran_LIBRARIES}
       ${MPI_CXX_LINK_FLAGS}
    )

endif()


#
# OpenSeesMP Tcl MPI Interpreter
#

if(MPI_FOUND)

    add_executable(OpenSeesMP EXCLUDE_FROM_ALL
      ${OPS_SRC_DIR}/tcl/mpiParameterMain.cpp
      ${OPS_SRC_DIR}/tcl/tclMain.cpp
      ${OPS_SRC_DIR}/tcl/commands.cpp
      ${OPS_SRC_DIR}/domain/domain/partitioned/PartitionedDomain.cpp
      ${OPS_SRC_DIR}/domain/domain/partitioned/PartitionedDomainEleIter.cpp
      ${OPS_SRC_DIR}/domain/domain/partitioned/PartitionedDomainSubIter.cpp
      ${OPS_SRC_DIR}/actor/machineBroker/MPI_MachineBroker.cpp
      ${OPS_SRC_DIR}/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp
      ${OPS_SRC_DIR}/system_of_eqn/linearSOE/diagonal/MPIDiagonalSOE.cpp
      ${OPS_SRC_DIR}/system_of_eqn/linearSOE/diagonal/MPIDiagonalSolver.cpp
      ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsSOE.cpp
      ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsSolver.cpp
      ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsParallelSOE.cpp
      ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.cpp
      # Ladruno ADR43 P3c-MPI: distributed FEAST inner solve (self-registers
      # the LadrunoFeastInnerSolve factory; _PARALLEL_INTERPRETERS targets only)
      ${OPS_SRC_DIR}/system_of_eqn/eigenSOE/LadrunoDistBlockZKernel.cpp
      ${OPS_SRC_DIR}/domain/subdomain/ActorSubdomain.cpp
      ${OPS_SRC_DIR}/domain/subdomain/ShadowSubdomain.cpp
      ${OPS_INTERP_TCL_PER_TARGET_SOURCES}
  ${OPS_MPI_PER_TARGET_SOURCES}   # Ladruno ADR-78 (contact abort + auto-penalty allreduce)
    )

    if(DEFINED MUMPS_INCLUDE_DIR)
        target_include_directories(OpenSeesMP PRIVATE ${MUMPS_INCLUDE_DIR} ${MPI_CXX_INCLUDE_DIRS})
    else()
        target_include_directories(OpenSeesMP PRIVATE ${MUMPS_DIR}/_deps/mumps-src/include ${MPI_CXX_INCLUDE_DIRS})
    endif()

    target_compile_options(OpenSeesMP PRIVATE ${MPI_CXX_COMPILE_FLAGS})

    if (NOT DEFINED MUMPS_DIR AND TARGET mumps)
        add_dependencies(OpenSeesMP mumps)
    endif()

    if (DEFINED OPENMPI)
        target_compile_definitions(OpenSeesMP
        PUBLIC _PARALLEL_INTERPRETERS ${MUMPS_FLAG} _OPENMPI)
    else()
        target_compile_definitions(OpenSeesMP
        PUBLIC _PARALLEL_INTERPRETERS ${MUMPS_FLAG})
    endif()

    target_link_libraries(OpenSeesMP
       OPS_InterpTcl
       OPS_InterpPyCmds
       $<$<BOOL:${LADRUNO_CMS}>:OPS_LadrunoCMS>
       OpenSeesLIB
       OPS_Reliability
       OPS_ReliabilityTcl
       OPS_Recorder
       METIS
       SUPERLU_DIST
       OPS_Numerics
       ${MUMPS_LIBRARIES}
       ${CMAKE_DL_LIBS}
       HDF5::HDF5
       tcl::tcl	
       ZLIB::ZLIB		
       Eigen3::Eigen
       ${MPI_CXX_LIBRARIES}
       ${SCALAPACK_LIBRARIES}
       ${MPI_Fortran_LIBRARIES}
       ${MPI_CXX_LINK_FLAGS}
    )

    if(LADRUNO_CMS)
      target_compile_definitions(OpenSeesMP PRIVATE _LADRUNO_CMS)
    endif()

endif()


#
# OpenSeesPy Module
#

add_library(OpenSeesPy SHARED EXCLUDE_FROM_ALL
    ${OPS_SRC_DIR}/interpreter/PythonModule.cpp
    ${OPS_SRC_DIR}/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp
    ${OPS_MPI_PER_TARGET_SOURCES}   # Ladruno ADR-78 (contact abort + auto-penalty allreduce)
)

set_target_properties(OpenSeesPy PROPERTIES PREFIX "")

target_include_directories(OpenSeesPy PUBLIC ${Python_INCLUDE_DIRS})

target_link_libraries(OpenSeesPy
   OPS_InterpPyCmds
   OpenSeesLIB
   OPS_SysOfEqnPython
   OPS_Reliability
   OPS_Recorder
   OPS_Numerics
   HDF5::HDF5
   ZLIB::ZLIB
   Eigen3::Eigen
   ${Python_LIBRARIES}
   ${PARDISO_MKL_LIBRARIES}  # Ladruno ADR-75 P1b (threaded MKL; empty w/o MKL)
)


#
# OpenSeesPyMP Module (MPI-enabled Python module, imports as `openseesmp`)
#
# Ladruno Patch 9: the Python analogue of the Tcl OpenSeesMP exe. Same
# parallel source set as OpenSeesMP (minus the Tcl mains/commands.cpp), the
# Python entry trio added via SRC/interpreter/CMakeLists.txt, and the
# parallel-compiled command engine supplied by OPS_InterpPyCmds_MP. NOT in
# OpenSeesLIB's link set, so no seq/parallel OpenSeesCommands.o collision.
#

if(MPI_FOUND)

add_library(OpenSeesPyMP SHARED EXCLUDE_FROM_ALL
  ${OPS_MPI_PER_TARGET_SOURCES}   # Ladruno ADR-78 (contact abort + auto-penalty allreduce)
  ${OPS_SRC_DIR}/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp
  ${OPS_SRC_DIR}/domain/domain/partitioned/PartitionedDomain.cpp
  ${OPS_SRC_DIR}/domain/domain/partitioned/PartitionedDomainEleIter.cpp
  ${OPS_SRC_DIR}/domain/domain/partitioned/PartitionedDomainSubIter.cpp
  ${OPS_SRC_DIR}/actor/machineBroker/MPI_MachineBroker.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/diagonal/MPIDiagonalSOE.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/diagonal/MPIDiagonalSolver.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsSOE.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsSolver.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsParallelSOE.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.cpp
  # Ladruno ADR43 P3c-MPI: distributed FEAST inner solve (self-registers the
  # LadrunoFeastInnerSolve factory; _PARALLEL_INTERPRETERS targets only)
  ${OPS_SRC_DIR}/system_of_eqn/eigenSOE/LadrunoDistBlockZKernel.cpp
  ${OPS_SRC_DIR}/domain/subdomain/ActorSubdomain.cpp
  ${OPS_SRC_DIR}/domain/subdomain/ShadowSubdomain.cpp
)

set_target_properties(OpenSeesPyMP PROPERTIES
  PREFIX ""
  OUTPUT_NAME "openseesmp")

target_include_directories(OpenSeesPyMP PUBLIC ${Python_INCLUDE_DIRS})
if(DEFINED MUMPS_INCLUDE_DIR)
    target_include_directories(OpenSeesPyMP PRIVATE ${MUMPS_INCLUDE_DIR} ${MPI_CXX_INCLUDE_DIRS})
else()
    target_include_directories(OpenSeesPyMP PRIVATE ${MUMPS_DIR}/_deps/mumps-src/include ${MPI_CXX_INCLUDE_DIRS})
endif()

target_compile_options(OpenSeesPyMP PRIVATE ${MPI_CXX_COMPILE_FLAGS})

if (NOT DEFINED MUMPS_DIR AND TARGET mumps)
    add_dependencies(OpenSeesPyMP mumps)
endif()

if (DEFINED OPENMPI)
    target_compile_definitions(OpenSeesPyMP
    PUBLIC _PARALLEL_INTERPRETERS ${MUMPS_FLAG} _OPENMPI)
else()
    target_compile_definitions(OpenSeesPyMP
    PUBLIC _PARALLEL_INTERPRETERS ${MUMPS_FLAG})
endif()

target_link_libraries(OpenSeesPyMP
   OPS_InterpPyCmds_MP
   $<$<BOOL:${LADRUNO_CMS}>:OPS_LadrunoCMS>
   OpenSeesLIB
   OPS_SysOfEqnPython
   OPS_Reliability
   OPS_Recorder
   OPS_Numerics
   METIS
   SUPERLU_DIST
   ${MUMPS_LIBRARIES}
   ${CMAKE_DL_LIBS}
   HDF5::HDF5
   ZLIB::ZLIB
   Eigen3::Eigen
   ${Python_LIBRARIES}
   ${MPI_CXX_LIBRARIES}
   ${SCALAPACK_LIBRARIES}
   ${MPI_Fortran_LIBRARIES}
   ${MPI_CXX_LINK_FLAGS}
)

if(LADRUNO_CMS)
  target_compile_definitions(OpenSeesPyMP PRIVATE _LADRUNO_CMS)
endif()

endif()


#
# (A Python "openseessp" / OpenSeesPySP SHARED target was prototyped here
#  and REMOVED. The Python/DL command engine has no _PARALLEL_PROCESSING
#  domain-decomposition subsystem — SP is Tcl-only via OpenSeesSP.exe.
#  Definitive analysis + decision: Ladruno_implementation/02_openseespymp.md
#  and compilation journal Patch 9. For parallel Python use OpenSeesPyMP.)
#


#
# INSTALL
#

install (TARGETS OpenSees DESTINATION bin)
if (TCL_INIT_FILE)
   install(FILES "${TCL_INIT_FILE}"
       DESTINATION lib/tcl8.6 
   )
endif()
install(PROGRAMS $<TARGET_FILE:OpenSeesPy> 
        DESTINATION bin 
        RENAME opensees.so)


# Set selected executable to be built by default
set_target_properties(${OPS_FINAL_TARGET} PROPERTIES EXCLUDE_FROM_ALL OFF)

# Add sources to targets
add_subdirectory(${OPS_SRC_DIR})




#==============================================================================
#                            Configure targets
#
#==============================================================================

#target_link_libraries(OPS_Numerics     INTERFACE ${OPS_Numlib_List})
#message("Numlibs: ${OPS_Numlib_List}")
#target_link_libraries(OPS_Element      PRIVATE ${OPS_Element_List})
#target_link_libraries(OPS_ObjectBroker PRIVATE ${OPS_Element_List})
#target_link_libraries(OPS_InterpTcl    PRIVATE ${OPS_Element_List})



#==============================================================================
#                            Apply Options
#
#==============================================================================
#----------------------------
# FMK
#----------------------------
if(FMK)
    add_compile_definitions(
        _HAVE_Damage2p
        _HAVE_PSUMAT
        _HAVE_PML
        _FILIP_LHNMYS
    )
endif()

#----------------------------
# Extensions
#----------------------------
message("OPS >>> Configuring OpenSees extensions")
foreach(extension IN LISTS OPS_SysOfEqn_List OPS_Element_List OPS_Extension_List)
    string(TOUPPER "${extension}" ext_flag)
    string(REGEX REPLACE "^OPS_" "OPSDEF_" ext_flag "${ext_flag}")
    add_compile_definitions(${ext_flag})
endforeach()
foreach(extension IN LISTS OPS_Exclude_List)
    string(TOUPPER "${extension}" ext_flag)
    string(REGEX REPLACE "^OPS_" "OPS_EXCLUDE_" ext_flag "${ext_flag}")
    message("    Adding macro definition '${ext_flag}'")
    add_compile_definitions(${ext_flag})
endforeach()


# Renderer
#----------------------------
if (${OPS_Use_Graphics_Option} STREQUAL "Base")
  target_link_libraries(${OPS_FINAL_TARGET} OPS_Graphics_Default OPS_Graphics)

elseif (${OPS_Use_Graphics_Option} STREQUAL "OpenGL")
  message("OPS >>> Including OpenGL graphics option")
  set_source_files_properties(
      "${OPS_SRC_DIR}/recorder/AlgorithmIncrements.cpp"
      "${OPS_SRC_DIR}/recorder/FilePlotter.cpp"
      "${OPS_SRC_DIR}/renderer/main.cpp"
      "${OPS_SRC_DIR}/renderer/OpenGlDevice.cpp"
      "${OPS_SRC_DIR}/renderer/OpenGlRenderer.cpp"
      "${OPS_SRC_DIR}/tcl/TclFeViewer.cpp"
      "${OPS_SRC_DIR}/tcl/TclVideoPlayer.cpp"

    PROPERTIES COMPILE_DEFINITIONS _AGL
  )
  target_link_libraries(${OPS_FINAL_TARGET} OPS_Graphics_GL OPS_Graphics)
else()
  add_compile_definitions(_NOGRAPHICS)
endif()


# Reliability
#----------------------------
#if ("OPS_Reliability" IN_LIST OPS_Extension_List)
    add_compile_definitions(_RELIABILITY)
    #target_link_libraries(${OPS_FINAL_TARGET} OPS_Reliability)
#endif()


#----------------------------
# HDF5
#----------------------------

if(HDF5_FOUND)
   include_directories(${HDF5_INCLUDE_DIR})
   set(_hdf5_libs hdf5 hdf5_cpp)
   if (HDF5_VERSION VERSION_GREATER_EQUAL 1.12.0)
         add_compile_definitions(_H5DRM)
         add_compile_definitions(_HDF5)
         # Ladruno (DRM study): add_compile_definitions only reaches targets
         # created AFTER this point, but OPS_INTERPRETER is defined ~line 594
         # (well before here), so OpenSeesPatternCommands.cpp was compiled
         # WITHOUT _H5DRM -> ops.pattern('H5DRM',...) returned "unknown pattern
         # type" from openseespy even though the class is linked. Apply the
         # defines to that already-created target explicitly.
         if (TARGET OPS_INTERPRETER)
            target_compile_definitions(OPS_INTERPRETER PRIVATE _H5DRM _HDF5)
         endif()
         message(STATUS "OPS >>> Have HDF5 and VERSION >= 1.12.0")
   endif()
else()
   message(STATUS "OPS >>> Could not find HDF5")
endif()

#----------------------------
# Eigen3
#----------------------------
if(Eigen3_FOUND)
   include_directories(${Eigen3_INCLUDE_DIR})
   add_compile_definitions(_EIGEN3)
else()
   message(STATUS "OPS >>> Could not find Eigen3")
endif()

if (OPS_Use_Dev_Directories)
  add_subdirectory("${PROJECT_SOURCE_DIR}/DEVELOPER/")
endif()
