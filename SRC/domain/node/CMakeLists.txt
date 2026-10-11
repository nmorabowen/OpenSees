#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Domain
  PRIVATE
    Node.cpp
    NodalLoad.cpp
    RelativeMassDamping.cpp
  PUBLIC
    Node.h
    NodalLoad.h
    RelativeMassDamping.h
)

target_include_directories(OPS_Domain PUBLIC ${CMAKE_CURRENT_LIST_DIR})
