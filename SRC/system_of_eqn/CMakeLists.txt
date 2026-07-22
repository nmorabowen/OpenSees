#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_SysOfEqn
    PRIVATE
    PUBLIC
    Solver.h
    SystemOfEqn.h
)

target_include_directories(OPS_SysOfEqn PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_sources(OPS_SysOfEqnPython
    PRIVATE
    PUBLIC
)

target_include_directories(OPS_SysOfEqn PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(linearSOE)
add_subdirectory(eigenSOE)

# Ladruno ADR 1000: independent, opt-in CMS library for MP targets only.
if(LADRUNO_CMS)
  add_subdirectory(ladrunoCMS)
endif()
