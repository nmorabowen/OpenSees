# CompilerSettings.cmake
# Configure compiler settings for OpenSees

# Find MPI before compiler-specific settings
find_package(MPI)

# Compiler-specific settings
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  MESSAGE("COMPILER: Clang")
  # set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(BUILD_SHARED_LIBS OFF)
  # set(CMAKE_EXE_LINKER_FLAGS "-static")
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".so")
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
else()
  MESSAGE("COMPILER: UNKNOWN COMPILER ${CMAKE_CXX_COMPILER_ID}")
endif()

if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    message(STATUS "Detected IntelLLVM Fortran Compiler")

    # Safe default flags for Intel oneAPI ifx
    set(CMAKE_Fortran_FLAGS_RELEASE "/O2 /nologo /QxHost /fp:precise /fpp /Qlocation,link,\"${CMAKE_LINKER}\"" CACHE STRING "" FORCE)

    # Prevent path splitting issues
    set(CMAKE_Fortran_PREPROCESS_SOURCE TRUE)
    set(CMAKE_Fortran_PREPROCESS_FLAG "/fpp")

    # Explicitly quote any flags involving paths with spaces
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /fpp /Qm64")
endif()

# Enable language support
enable_language(Fortran)

# Set language standards
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD 11)
set(CMAKE_FORTRAN_STANDARD 08)

set(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> cqls <TARGET> <LINK_FLAGS> <OBJECTS>")