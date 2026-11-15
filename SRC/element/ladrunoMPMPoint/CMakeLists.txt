#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# LadrunoMPMPoint — the mode-A MPM particle element (Ladruno fork, WP-232f).
#   ELE_TAG_LadrunoMPMPoint = 33026. One integration point and one NDMaterial
#   copy per material point; its nodes are the particle's active GIMP support,
#   re-pointed every step by the MPM host (SRC/domain/pattern/ladrunoMPM). Not
#   user-creatable (no parser). See Ladruno_implementation/232_ladruno_mpm_guide.md.
target_sources(OPS_Element
    PRIVATE
    LadrunoMPMPoint.cpp
    PUBLIC
    LadrunoMPMPoint.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
