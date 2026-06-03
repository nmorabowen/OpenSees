#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# LadrunoEmbeddedRebar — embedded-reinforcement coupling element (Ladruno fork).
target_sources(OPS_Element
    PRIVATE
    LadrunoEmbeddedRebar.cpp
    OPS_LadrunoEmbeddedRebar.cpp
    PUBLIC
    LadrunoEmbeddedRebar.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
