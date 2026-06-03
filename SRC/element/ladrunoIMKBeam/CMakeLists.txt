#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# LadrunoIMKBeam — 2D/3D concentrated-plasticity IMK beam macro elements
# (Ladruno). Both dimensions share the per-axis hinge kernel LadrunoIMKHinge.h.
target_sources(OPS_Element
    PRIVATE
    LadrunoIMKBeam.cpp
    OPS_LadrunoIMKBeam.cpp
    LadrunoIMKBeam2d.cpp
    OPS_LadrunoIMKBeam2d.cpp
    PUBLIC
    LadrunoIMKBeam.h
    LadrunoIMKBeam2d.h
    LadrunoIMKHinge.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
