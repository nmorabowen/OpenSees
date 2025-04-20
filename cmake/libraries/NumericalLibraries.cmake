# NumericalLibraries.cmake
# Set up numerical libraries for OpenSees

# Add numerical subdirectories
add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/SuperLU_5.1.1")
set(SUPERLU_LIBRARIES SUPERLU)
add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/UMFPACK")
set(UMFPACK_LIBRARIES UMFPACK)
add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/AMD")
set(AMD_LIBRARIES AMD)
add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/ARPACK")
set(ARPACK_LIBRARIES ARPACK)
add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/CSPARSE")
set(CSPARSE_LIBRARIES CSPARSE)
add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/tetgen1.4.3")
set(TETGEN_LIBRARIES tet)
add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/Triangle")
set(TRIANGLE_LIBRARIES triangle)
# add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/eigenAPI")
# set(EIGENAPI_LIBRARIES EIGENAPI)

# Create interface library for numerical packages
add_library(OPS_Numerics INTERFACE)

target_link_libraries(OPS_Numerics INTERFACE
  ${ARPACK_LIBRARIES}
  ${CSPARSE_LIBRARIES} 
  ${SUPERLU_LIBRARIES}
  ${UMFPACK_LIBRARIES}
  ${TETGEN_LIBRARIES}
  ${TRIANGLE_LIBRARIES}                     
  ${AMD_LIBRARIES}
  ${AMD_LIBRARIES}       
  ${LAPACK_LIBRARIES}
  ${EIGENAPI_LIBRARIES}
)