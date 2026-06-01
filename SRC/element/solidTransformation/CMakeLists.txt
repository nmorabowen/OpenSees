#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# SolidTransformation — geometry-method layer for solid elements (Ladruno fork).
# The solid-element analogue of the beam CrdTransf/CorotCrdTransf abstraction.
# v1 ships the identity method (SolidTransformationLinear); corot/finite add
# source files here as they land. Element-agnostic — no classTags.h entry.
target_sources(OPS_Element
    PRIVATE
    SolidTransformationLinear.cpp
    SolidTransformationFinite.cpp
    PUBLIC
    SolidTransformation.h
    SolidTransformationLinear.h
    SolidTransformationFinite.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
