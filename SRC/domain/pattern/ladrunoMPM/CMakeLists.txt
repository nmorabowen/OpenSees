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
    PUBLIC
    LadrunoMPMPattern.h
    LadrunoMPMCore.h
    LadrunoMPMIO.h
)

target_include_directories(OPS_Domain PUBLIC ${CMAKE_CURRENT_LIST_DIR})
