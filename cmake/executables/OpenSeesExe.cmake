# OpenSeesExe.cmake
# Defines the main OpenSees executable

# OpenSees Tcl Interpreter
add_executable(OpenSees EXCLUDE_FROM_ALL 
  ${OPS_SRC_DIR}/tcl/tclAppInit.cpp
  ${OPS_SRC_DIR}/tcl/tclMain.cpp
  ${OPS_SRC_DIR}/tcl/commands.cpp
  ${OPS_SRC_DIR}/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp
)

# Copy TCL files based on Conan version
if(USING_CONAN2)
  # Conan 2 handling
  FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib/tcl8.6)
  add_custom_command(
    TARGET OpenSees POST_BUILD
    COMMENT "Copying TCL files to ${CMAKE_CURRENT_BINARY_DIR}"
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${TCL_LIB_DIR}/tcl8.6
    ${CMAKE_CURRENT_BINARY_DIR}/lib/tcl8.6
  )
elseif(USING_CONAN1)
  # Original Conan 1 handling
  FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib/tcl8.6)
  add_custom_command(
    TARGET OpenSees POST_BUILD
    COMMENT "Copying init.tcl to ${CMAKE_CURRENT_BINARY_DIR}"
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CONAN_LIB_DIRS_TCL}/tcl8.6/init.tcl
    ${CMAKE_CURRENT_BINARY_DIR}/lib/tcl8.6
  )
else()
  # Original non-Conan handling
  FILE(COPY ${TCL_LIBRARY_PATH}/tcl8.6
      DESTINATION ${PROJECT_BINARY_DIR}/lib/tcl8.6
      FILES_MATCHING PATTERN *.tcl)
endif()

# Link libraries based on Conan version
if(USING_CONAN2)
  # Conan 2 explicit dependencies
  target_link_libraries(OpenSees
    OPS_InterpTcl 
    coordTransformation
    OpenSeesLIB
    OPS_Reliability
    OPS_ReliabilityTcl  
    OPS_Numerics
    OPS_Recorder
    ${CMAKE_DL_LIBS} 
    hdf5::hdf5
    tcl::tcl
    ZLIB::ZLIB
    Eigen3::Eigen
  )
elseif(USING_CONAN1)
  # Original Conan 1 linking
  target_link_libraries(OpenSees
    OPS_InterpTcl 
    coordTransformation
    OpenSeesLIB
    OPS_Reliability
    OPS_ReliabilityTcl  
    OPS_Numerics
    OPS_Recorder
    ${CMAKE_DL_LIBS} 
    ${HDF5_LIBRARIES} 
    ${CONAN_LIBS}
  )
else()
  # Non-Conan linking
  target_link_libraries(OpenSees
    OPS_InterpTcl 
    coordTransformation
    OpenSeesLIB
    OPS_Reliability
    OPS_ReliabilityTcl  
    OPS_Numerics
    OPS_Recorder
    ${CMAKE_DL_LIBS} 
    ${HDF5_LIBRARIES}
  )
endif()