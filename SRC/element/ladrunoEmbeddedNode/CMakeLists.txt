#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# LadrunoEmbeddedNode — general node-to-host coupling element (Ladruno fork, ADR 23).
# Reuses LadrunoEmbeddedKernel.{cpp,h} from ../ladrunoEmbeddedRebar (on the include path).
target_sources(OPS_Element
    PRIVATE
    LadrunoEmbeddedNode.cpp
    OPS_LadrunoEmbeddedNode.cpp
    PUBLIC
    LadrunoEmbeddedNode.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
