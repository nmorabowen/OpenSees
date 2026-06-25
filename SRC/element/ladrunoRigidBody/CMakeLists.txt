#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# LadrunoRigidBody — 6-DOF rigid-body element (Ladruno fork, ADR 58).
# Zero-stiffness Element owning a private internal CoM Node + rigid-link MP slaving
# (Joint3D pattern). P1 = translation/ballistic; SO(3) rotation self-integrated at P2.
target_sources(OPS_Element
    PRIVATE
    LadrunoRigidBody.cpp
    OPS_LadrunoRigidBody.cpp
    PUBLIC
    LadrunoRigidBody.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
