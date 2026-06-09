#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# LadrunoDistributingCoupling — RBE3 / distributing-coupling element (Ladruno fork, ADR 28).
# Reuses LadrunoEmbeddedKernel.{cpp,h} from ../ladrunoEmbeddedRebar (on the include path).
target_sources(OPS_Element
    PRIVATE
    LadrunoDistributingCoupling.cpp
    OPS_LadrunoDistributingCoupling.cpp
    PUBLIC
    LadrunoDistributingCoupling.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
