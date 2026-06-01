#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# LadrunoBrick — unified 8-node hexahedral solid element (Ladruno fork).
# Reuses shp3d (trilinear shape functions) from ../brick.
target_sources(OPS_Element
    PRIVATE
    LadrunoBrick.cpp
    OPS_LadrunoBrick.cpp
    PUBLIC
    LadrunoBrick.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
