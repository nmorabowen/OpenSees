# InstallRules.cmake
# Installation rules for OpenSees

# Install the OpenSees executable
install(TARGETS OpenSees DESTINATION bin)

# Install TCL initialization files
if(USING_CONAN2)
  # For Conan 2, we need to get the TCL initialization file location
  set(TCL_INIT_FILE "${TCL_LIB_DIR}/tcl8.6/init.tcl")
elseif(USING_CONAN1)
  # For Conan 1, we use the CONAN_LIB_DIRS_TCL variable
  set(TCL_INIT_FILE "${CONAN_LIB_DIRS_TCL}/tcl8.6/init.tcl")
else()
  # Without Conan, we use the TCL_INCLUDE_PATH to find the lib directory
  set(TCL_INIT_FILE "${TCL_INCLUDE_PATH}/../lib/tcl8.6/init.tcl")
endif()

# Install the init.tcl file
install(FILES ${TCL_INIT_FILE} DESTINATION lib/tcl8.6)