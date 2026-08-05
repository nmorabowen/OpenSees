#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Analysis
    PRIVATE
      DOF_Numberer.cpp
      PlainNumberer.cpp
      ParallelNumberer.cpp
      LadrunoParallelNumberer.cpp
    PUBLIC
      DOF_Numberer.h
      PlainNumberer.h
      ParallelNumberer.h
      LadrunoParallelNumberer.h
)

#target_include_directories(OPS_Analysis PUBLIC ${CMAKE_CURRENT_LIST_DIR})
