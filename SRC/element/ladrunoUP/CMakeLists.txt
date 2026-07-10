#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# LadrunoUP — Biot u-p saturated-porous continuum family (Ladruno fork, ADR 71).
#   ELE_TAG_LadrunoUP = 33017 (reserved; #define + element sources land in P1).
#
# P0 STUB: header-only shared kernel + shape providers. No compiled sources yet
# (the LadrunoUP element / OPS_ parser arrive in P1), so this adds ZERO build
# impact — it only publishes the two pure headers and their include directory so
# the P0 standalone harness (tests/ladruno_up_kernel_check.cpp) and the future
# P1 element can find them. See Ladruno_implementation/71_ladruno_up_family_adr.md
# and Ladruno_implementation/_adr71_implementation_plan.md.
target_sources(OPS_Element
    PUBLIC
    LadrunoUPKernel.h    # Ladruno (ADR 71): shared Biot u-p block kernel (header-only)
    LadrunoUPShapes.h    # Ladruno (ADR 71): stateless shape providers (header-only)
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
