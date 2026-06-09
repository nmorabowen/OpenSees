#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# LadrunoKinematicCoupling — RBE2 / kinematic-coupling element (Ladruno fork, ADR 29).
# Reuses the scalar bipenalty helpers from LadrunoEmbeddedKernel.h (on the include path);
# assembles its own B/K (no kernel dependency for the rotation/transport blocks).
target_sources(OPS_Element
    PRIVATE
    LadrunoKinematicCoupling.cpp
    OPS_LadrunoKinematicCoupling.cpp
    PUBLIC
    LadrunoKinematicCoupling.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
