#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# LadrunoDispBeamColumn — regularized displacement-based fiber frame element
# (Ladruno, ADR 32). Stage 0+1: 2D element with the Tier-1 per-IP characteristic
# length (lch) channel for mesh-objective softening.
target_sources(OPS_Element
    PRIVATE
    LadrunoDispBeamColumn2d.cpp
    LadrunoDispBeamColumn3d.cpp
    OPS_LadrunoDispBeamColumn.cpp
    PUBLIC
    LadrunoDispBeamColumn2d.h
    LadrunoDispBeamColumn3d.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
