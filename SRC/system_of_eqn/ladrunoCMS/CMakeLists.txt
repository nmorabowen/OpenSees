# Ladruno ADR 1000: independent two-level component-mode synthesis library.

target_sources(OPS_LadrunoCMS
  PRIVATE
    LadrunoCMSOptions.cpp
    LadrunoCMSMumps.cpp
    LadrunoCMSLocalLanczos.cpp
    LadrunoCMSHierarchy.cpp
    LadrunoCMSSubspaceRefiner.cpp
    LadrunoCMSEigenSOE.cpp
    LadrunoCMSEigenSolver.cpp
  PUBLIC
    LadrunoCMSOptions.h
    LadrunoCMSMumps.h
    LadrunoCMSLocalLanczos.h
    LadrunoCMSHierarchy.h
    LadrunoCMSSubspaceRefiner.h
    LadrunoCMSEigenSOE.h
    LadrunoCMSEigenSolver.h
)

target_include_directories(OPS_LadrunoCMS
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
    ${CMAKE_CURRENT_LIST_DIR}/../eigenSOE
  PRIVATE
    ${MPI_CXX_INCLUDE_DIRS}
)

if(DEFINED MUMPS_INCLUDE_DIR)
  target_include_directories(OPS_LadrunoCMS PRIVATE ${MUMPS_INCLUDE_DIR})
else()
  target_include_directories(OPS_LadrunoCMS PRIVATE ${MUMPS_DIR}/_deps/mumps-src/include)
endif()

target_compile_options(OPS_LadrunoCMS PRIVATE ${MPI_CXX_COMPILE_FLAGS})
target_compile_options(OPS_LadrunoCMS PRIVATE
  $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>:-Wall>
  $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>:-Wextra>
  $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>:-Wpedantic>
  $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>:-Wconversion>
  $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>:-Wshadow>
)
target_link_libraries(OPS_LadrunoCMS
  PRIVATE
    METIS
  PUBLIC
    MPI::MPI_CXX
    ${MUMPS_LIBRARIES}
    ${SCALAPACK_LIBRARIES}
    ${MPI_Fortran_LIBRARIES}
    ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}
)

# BLAS/LAPACK layer. When MKL is present, SCALAPACK_LIBRARIES above ALREADY
# carries the full static sequential MKL stack (scalapack_lp64 + intel_lp64 +
# sequential + core + blacs) -- exactly what OpenSeesSP/OpenSeesMP link. Adding
# ${LAPACK_LIBRARIES} on top pulls in the *dynamic* MKL interface that
# find_package(LAPACK) resolves to (mkl_intel_lp64_dll + mkl_intel_thread_dll +
# mkl_core_dll + libiomp5md); mixing the static and dynamic layers in one link
# is unsupported by MKL and fails LNK2005 (mkl_serv_verbose_mode /
# mkl_serv_xerbla multiply defined), taking down every ladruno_cms_*_check and
# any MP target built with LADRUNO_CMS=ON. Only fall back to LAPACK_LIBRARIES on
# a non-MKL box, where SCALAPACK_LIBRARIES is user-supplied and carries no BLAS.
if(NOT MKL_FOUND)
  target_link_libraries(OPS_LadrunoCMS PUBLIC ${LAPACK_LIBRARIES})
endif()

if(LADRUNO_CMS_BUILD_TESTS)
  add_executable(ladruno_cms_mumps_check EXCLUDE_FROM_ALL
    ${PROJECT_SOURCE_DIR}/tests/ladruno_cms_mumps_check.cpp
  )
  target_link_libraries(ladruno_cms_mumps_check PRIVATE OPS_LadrunoCMS)

  add_executable(ladruno_cms_lanczos_check EXCLUDE_FROM_ALL
    ${PROJECT_SOURCE_DIR}/tests/ladruno_cms_lanczos_check.cpp
  )
  target_link_libraries(ladruno_cms_lanczos_check PRIVATE OPS_LadrunoCMS)

  add_executable(ladruno_cms_hierarchy_check EXCLUDE_FROM_ALL
    ${PROJECT_SOURCE_DIR}/tests/ladruno_cms_hierarchy_check.cpp
  )
  target_link_libraries(ladruno_cms_hierarchy_check PRIVATE OPS_LadrunoCMS)

  add_executable(ladruno_cms_subspace_check EXCLUDE_FROM_ALL
    ${PROJECT_SOURCE_DIR}/tests/ladruno_cms_subspace_check.cpp
  )
  target_link_libraries(ladruno_cms_subspace_check PRIVATE OPS_LadrunoCMS)

  add_executable(ladruno_cms_assembly_check EXCLUDE_FROM_ALL
    ${PROJECT_SOURCE_DIR}/tests/ladruno_cms_assembly_check.cpp
  )
  target_link_libraries(ladruno_cms_assembly_check PRIVATE OPS_LadrunoCMS)

  add_executable(ladruno_cms_topology_check EXCLUDE_FROM_ALL
    ${PROJECT_SOURCE_DIR}/tests/ladruno_cms_topology_check.cpp
  )
  target_link_libraries(ladruno_cms_topology_check PRIVATE
    OPS_LadrunoCMS OpenSeesLIB)
endif()
