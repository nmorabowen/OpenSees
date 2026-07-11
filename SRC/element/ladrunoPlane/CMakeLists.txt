#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# LadrunoPlane — unified 2D continuum elements (Ladruno fork):
#   LadrunoQuad (4-node, 33007) + LadrunoCST (3-node, 33008) +
#   LadrunoLST (6-node T6, 33016; ADR 70 P3).
# See Ladruno_implementation/25_ladruno_plane_elements_adr.md and
# Ladruno_implementation/70_ladruno_plane_finite_triangles_adr.md
target_sources(OPS_Element
    PRIVATE
    LadrunoQuad.cpp
    OPS_LadrunoQuad.cpp
    LadrunoCST.cpp
    OPS_LadrunoCST.cpp
    LadrunoLST.cpp
    OPS_LadrunoLST.cpp
    PUBLIC
    LadrunoQuad.h
    LadrunoCST.h
    LadrunoLST.h
    LadrunoFiniteStrain2DKernel.h   # Ladruno (ADR 70): shared 2D finite-strain kernel (header-only)
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
