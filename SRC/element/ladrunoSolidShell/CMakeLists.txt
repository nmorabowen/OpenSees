#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# LadrunoSolidShell — 8-node ANS/EAS solid-shell element (Ladruno fork, ADR 64).
target_sources(OPS_Element
    PRIVATE
    LadrunoSolidShell.cpp
    OPS_LadrunoSolidShell.cpp
    PUBLIC
    LadrunoSolidShell.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
