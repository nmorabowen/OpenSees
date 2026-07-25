#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#
# Ladruno ADR-75 P1b: build the MKL PARDISO SOE/solver pair.
#
# The parent only add_subdirectory()s this when MKL_FOUND, so a non-MKL build
# (e.g. the Zone-A Ubuntu reference-LAPACK runner) never compiles these and
# `system Pardiso` stays unregistered -- users there fall back to UmfPack.
#
# The objects land in OPS_SysOfEqn (=> OpenSeesLIB => every target), but the MKL
# LINK layer is chosen per target on purpose: the serial targets take the
# THREADED layer (that is the whole point of desktop PARDISO) while the parallel
# targets keep the sequential layer they already use for ScaLAPACK. These
# translation units are layer-agnostic -- they only call PARDISO().

target_sources(OPS_SysOfEqn
    PRIVATE
        PARDISOGenLinSOE.cpp
        PARDISOGenLinSolver.cpp

    PUBLIC
        PARDISOGenLinSOE.h
        PARDISOGenLinSolver.h
)

target_include_directories(OPS_SysOfEqn PUBLIC ${CMAKE_CURRENT_LIST_DIR})

# mkl_pardiso.h / mkl_types.h
if(DEFINED MKL_ROOT)
    target_include_directories(OPS_SysOfEqn PUBLIC ${MKL_ROOT}/include)
endif()
