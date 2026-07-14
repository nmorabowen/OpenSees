#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# LadrunoPorousOverlay — persistent-fluid staggered u-p overlay (Ladruno fork,
#   ADR 73). PATTERN_TAG_LadrunoPorousOverlay = 33022. A LoadPattern subclass
#   that owns a pore-pressure field OUTSIDE the Domain DOF graph over a geometry
#   snapshot of a named region of plain solid cells (T3/Q4/H8), coupling to the
#   skeleton through +Q·p (applyLoad) and advancing its own SPD fluid system at
#   commit (fixed-stress RATE form). Reuses SRC/element/ladrunoUP/LadrunoUPKernel.h
#   + LadrunoUPShapes.h verbatim (header-only, relative include).
# See Ladruno_implementation/73_ladruno_porous_overlay_adr.md and
# Ladruno_implementation/_adr73_implementation_plan.md (WP1 pins).
target_sources(OPS_Domain
    PRIVATE
    LadrunoPorousOverlay.cpp        # Ladruno (ADR 73 P1): the overlay engine
    OPS_LadrunoPorousOverlay.cpp    # Ladruno (ADR 73 P1): pattern parser/factory
    PUBLIC
    LadrunoPorousOverlay.h
)

target_include_directories(OPS_Domain PUBLIC ${CMAKE_CURRENT_LIST_DIR})
