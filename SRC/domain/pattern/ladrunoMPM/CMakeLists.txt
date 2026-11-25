#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# LadrunoMPMPattern — the explicit MPM engine host (Ladruno fork, WP-232 mode F).
#   PATTERN_TAG_LadrunoMPM = 33027. A LoadPattern subclass owning the
#   background grid, the particles and one NDMaterial copy per particle; the
#   `mpm` command (one LadrunoCommandTable.h row) drives it. The numerics are the
#   header-only, OpenSees-free LadrunoMPMCore.h, checked standalone by
#   tests/mpm_engine_check.cpp. See Ladruno_implementation/232_ladruno_mpm_guide.md.
target_sources(OPS_Domain
    PRIVATE
    LadrunoMPMPattern.cpp           # Ladruno (WP-232): the host LoadPattern
    OPS_LadrunoMPM.cpp              # Ladruno (WP-232): the `mpm` command
    LadrunoMPMIO.cpp                # Ladruno (WP-232 c2): contract 1 reader, contract 2 writer
    LadrunoMPMImplicit.cpp          # Ladruno (WP-232f): mode A, particles as LadrunoMPMPoint elements
    LadrunoMPMSoil.cpp              # Ladruno (WP-232k): soil particles -- stage walk, guard, queries, mode-F gridvel
    LadrunoMPMTwoPhase.cpp          # Ladruno (WP-232l): two-phase particles -- route E, pfix, plate, traction, debug hooks
    PUBLIC
    LadrunoMPMPattern.h
    LadrunoMPMCore.h
    LadrunoMPMIO.h
)

target_include_directories(OPS_Domain PUBLIC ${CMAKE_CURRENT_LIST_DIR})

# Ladruno WP-253: the threaded mode-F step is bit-identical to 1 thread only if its gathers and
# the serial scatters round alike, so floating-point contraction is pinned OFF for these
# translation units on GCC / Clang (a -march with FMA, e.g. a spack wrapper's, would otherwise
# let the compiler fuse a*b+c differently in the two loops). Without -march, x86-64 GCC emits no
# FMA at all (checked: g++ 15 -O2/-O3 on the core, zero vfmadd), so no existing Linux record
# changes. MSVC needs nothing: /fp:precise (the default) does not contract without /fp:contract
# (VS 2022), and the build sets no /arch:AVX2 that would provide FMA.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  set(_ladruno_mpm_sources
      LadrunoMPMPattern.cpp OPS_LadrunoMPM.cpp LadrunoMPMIO.cpp LadrunoMPMImplicit.cpp
      LadrunoMPMSoil.cpp LadrunoMPMTwoPhase.cpp)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.18)
    set_property(SOURCE ${_ladruno_mpm_sources} TARGET_DIRECTORY OPS_Domain
                 APPEND PROPERTY COMPILE_OPTIONS "-ffp-contract=off")
  else()
    message(WARNING "Ladruno WP-253: CMake < 3.18 cannot set -ffp-contract=off on the MPM sources;"
                    " the threaded MPM step's bit-identity then relies on no -march with FMA")
  endif()
endif()
