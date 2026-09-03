#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_SysOfEqn
  PRIVATE
    DomainSolver.cpp
    LinearSOE.cpp
    LinearSOESolver.cpp
  PUBLIC
    DomainSolver.h
    LinearSOE.h
    LinearSOESolver.h
)

target_include_directories(OPS_SysOfEqn PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(bandGEN)
add_subdirectory(bandSPD)
add_subdirectory(diagonal)
add_subdirectory(fullGEN)
add_subdirectory(sparseGEN)
add_subdirectory(sparseSYM)
add_subdirectory(umfGEN)
add_subdirectory(sparsePython)

add_subdirectory(profileSPD)
#add_subdirectory(cg)
#add_subdirectory(petsc)
#add_subdirectory(mumps)
add_subdirectory(itpack)

# Ladruno ADR-75 P1b: MKL PARDISO — the desktop/shared-memory sparse-direct
# solver. Gated on MKL because it IS MKL; a build without MKL simply has no
# `system Pardiso` (and falls back to UmfPack). Unlike mumps/ above, the sources
# can live in the shared OPS_SysOfEqn because the MKL *link layer* (threaded for
# serial, sequential for the ScaLAPACK-carrying parallel targets) is selected
# per target, not here.
# Ladruno (TIMs, 2026-09-03): also when the Linux opt-in found the MKL
# 3-layer + headers explicitly (LADRUNO_MKL_PARDISO_LINUX in the top-level
# CMakeLists), which is the case on a host where find_package(MKL) either
# fails or lands on an apt MKL whose kernels do not load.
if(MKL_FOUND OR LADRUNO_MKL_PARDISO_ON)
    add_subdirectory(pardiso)
endif()

