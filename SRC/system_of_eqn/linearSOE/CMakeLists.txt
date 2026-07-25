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
if(MKL_FOUND)
    add_subdirectory(pardiso)
endif()

