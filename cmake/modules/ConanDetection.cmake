# ConanDetection.cmake
# Detects Conan version and sets up appropriate variables

# Define a function to check if we're using Conan 2
function(detect_conan_version)
  # Check if user provided a CONAN_GENERATORS_DIR
  if(DEFINED CONAN_GENERATORS_DIR)
    if(EXISTS "${CONAN_GENERATORS_DIR}/conan_toolchain.cmake")
      set(USING_CONAN2 TRUE PARENT_SCOPE)
      set(USING_CONAN TRUE PARENT_SCOPE)
      MESSAGE("USING CONAN 2 with user-defined generators directory: ${CONAN_GENERATORS_DIR}")
      return()
    endif()
  endif()

  # Check for generators subdirectory with Conan 2 files
  if(EXISTS "${CMAKE_BINARY_DIR}/generators/conan_toolchain.cmake")
    set(USING_CONAN2 TRUE PARENT_SCOPE)
    set(USING_CONAN TRUE PARENT_SCOPE)
    set(CONAN_GENERATORS_DIR "${CMAKE_BINARY_DIR}/generators" PARENT_SCOPE)
    MESSAGE("USING CONAN 2 (generators subdirectory)")
    return()
  endif()
  
  # Check main build directory for Conan 2 files
  if(EXISTS "${CMAKE_BINARY_DIR}/conan_toolchain.cmake")
    set(USING_CONAN2 TRUE PARENT_SCOPE)
    set(USING_CONAN TRUE PARENT_SCOPE)
    set(CONAN_GENERATORS_DIR "${CMAKE_BINARY_DIR}" PARENT_SCOPE)
    MESSAGE("USING CONAN 2 (main build directory)")
    return()
  endif()
  
  # Check for Conan 1
  if(EXISTS "${CMAKE_BINARY_DIR}/conanbuildinfo.cmake")
    set(USING_CONAN1 TRUE PARENT_SCOPE)
    set(USING_CONAN TRUE PARENT_SCOPE)
    MESSAGE("USING CONAN 1")
    return()
  endif()
  
  # If we get here, we're not using Conan
  set(USING_CONAN FALSE PARENT_SCOPE)
  MESSAGE("NOT USING CONAN - Could not find Conan files in expected locations")
  MESSAGE("  Checked: ${CMAKE_BINARY_DIR}/generators/conan_toolchain.cmake")
  MESSAGE("  Checked: ${CMAKE_BINARY_DIR}/conan_toolchain.cmake")
  MESSAGE("  Checked: ${CMAKE_BINARY_DIR}/conanbuildinfo.cmake")
endfunction()

# Call the function to detect Conan version
detect_conan_version()