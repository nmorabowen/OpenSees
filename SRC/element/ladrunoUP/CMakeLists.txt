#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# LadrunoUP — Biot u-p saturated-porous continuum family (Ladruno fork, ADR 71).
#   ELE_TAG_LadrunoUP = 33017. One unified element over internal stateless shape
#   providers (P1 lanes: T3/Q4/H8 equal-order; P3: Bézier Taylor–Hood).
# See Ladruno_implementation/71_ladruno_up_family_adr.md and
# Ladruno_implementation/_adr71_implementation_plan.md.
target_sources(OPS_Element
    PRIVATE
    LadrunoUP.cpp        # Ladruno (ADR 71 P1): the unified u-p element
    OPS_LadrunoUP.cpp    # Ladruno (ADR 71 P1): parser/factory (both interpreters)
    PUBLIC
    LadrunoUP.h
    LadrunoUPKernel.h    # Ladruno (ADR 71): shared Biot u-p block kernel (header-only)
    LadrunoUPShapes.h    # Ladruno (ADR 71): stateless shape providers (header-only)
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
