#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Element
    PRIVATE
        LadrunoMembrane.cpp
        OPS_LadrunoMembrane.cpp
    PUBLIC
        LadrunoMembrane.h
        LadrunoMembraneKernel.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

option(LADRUNO_MUTATE_MEMBRANE "Enable ADR-87 mutation gate for LadrunoMembrane" OFF)
if(LADRUNO_MUTATE_MEMBRANE)
    target_compile_definitions(OPS_Element PRIVATE LADRUNO_MUTATE_MEMBRANE=1)
endif()
