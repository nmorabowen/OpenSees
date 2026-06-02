#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# LadrunoIMKBeam — 3D concentrated-plasticity IMK beam macro element (Ladruno).
target_sources(OPS_Element
    PRIVATE
    LadrunoIMKBeam.cpp
    OPS_LadrunoIMKBeam.cpp
    PUBLIC
    LadrunoIMKBeam.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
