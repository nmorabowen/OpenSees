#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#
#   Ladruno fork: ADR-39 ContactDomain subsystem (Domain-side state).
#==============================================================================
target_sources(OPS_Domain
    PRIVATE
      LadrunoContactDomain.cpp
      LadrunoContactSurface.cpp
    PUBLIC
      LadrunoContactDomain.h
      LadrunoContactSurface.h
)

target_include_directories(OPS_Domain PUBLIC ${CMAKE_CURRENT_LIST_DIR})
