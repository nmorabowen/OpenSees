# OpenSeesFunctions.cmake
# Utility functions for OpenSees build

# Function to add compiler flags based on compiler
function(opensees_add_cxx_flag)
  foreach(flag ${ARGN})
    if("${flag}" MATCHES "^([A-Z]+):(.+)$")
      set(compiler "${CMAKE_MATCH_1}")
      set(compiler_flag "${CMAKE_MATCH_2}")
      if("${CMAKE_CXX_COMPILER_ID}" MATCHES "${compiler}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${compiler_flag}" PARENT_SCOPE)
      endif()
    endif()
  endforeach()
endfunction()

# Function to determine architecture prefix
function(get_arch_prefix prefix)
  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(${prefix} "x64" PARENT_SCOPE)
  else()
    set(${prefix} "x86" PARENT_SCOPE)
  endif()
endfunction()

# Function to extract library information
function(get_lib_info lib_path name path)
  get_filename_component(${path} ${lib_path} DIRECTORY PARENT_SCOPE)
  get_filename_component(${name} ${lib_path} NAME_WE PARENT_SCOPE)
endfunction()