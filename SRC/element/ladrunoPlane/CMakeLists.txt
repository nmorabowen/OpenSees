#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# LadrunoPlane — unified 2D continuum elements (Ladruno fork):
#   LadrunoQuad (4-node, classTag 33006) + LadrunoCST (3-node, classTag 33007).
# See Ladruno_implementation/25_ladruno_plane_elements_adr.md
target_sources(OPS_Element
    PRIVATE
    LadrunoQuad.cpp
    OPS_LadrunoQuad.cpp
    LadrunoCST.cpp
    OPS_LadrunoCST.cpp
    PUBLIC
    LadrunoQuad.h
    LadrunoCST.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
