#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# LadrunoForceBeamColumn — regularized force-based (flexibility) fiber frame
# element with Tier-1 per-IP/element/user characteristic length (lch) and
# Tier-2 embedded/series cohesive fracture hinges (Ladruno WP-200).
target_sources(OPS_Element
    PRIVATE
    LadrunoForceBeamColumn2d.cpp
    LadrunoForceBeamColumn3d.cpp
    OPS_LadrunoForceBeamColumn.cpp
    PUBLIC
    LadrunoForceBeamColumn2d.h
    LadrunoForceBeamColumn3d.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
