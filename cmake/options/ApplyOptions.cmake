# ApplyOptions.cmake
# Apply build options to the OpenSees build

# Set selected executable to be built by default
set_target_properties(${OPS_FINAL_TARGET} PROPERTIES EXCLUDE_FROM_ALL OFF)

# FMK options
if(FMK)
    add_compile_definitions(
        _HAVE_Damage2p    
        _HAVE_PSUMAT
        _HAVE_PML
        _FILIP_LHNMYS
    )    
endif()

# Configure OpenSees extensions
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

# Graphics options
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
add_compile_definitions(_RELIABILITY)

# HDF5
if(HDF5_FOUND)
   include_directories(${HDF5_INCLUDE_DIR})
   set(_hdf5_libs hdf5 hdf5_cpp)
   if (HDF5_VERSION VERSION_GREATER_EQUAL 1.12.0)
         add_compile_definitions(_H5DRM)
         add_compile_definitions(_HDF5)
         message(STATUS "OPS >>> Have HDF5 and VERSION >= 1.12.0")    
   endif()
else()
   message(STATUS "OPS >>> Could not find HDF5")
endif()

# Eigen3
if(Eigen3_FOUND)
   include_directories(${Eigen3_INCLUDE_DIR})
   add_compile_definitions(_EIGEN3)
else()
   message(STATUS "OPS >>> Could not find Eigen3")
endif()

# Developer directories (if enabled)
if (OPS_Use_Dev_Directories)
  add_subdirectory("${PROJECT_SOURCE_DIR}/DEVELOPER/")
endif()