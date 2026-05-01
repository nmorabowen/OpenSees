#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

# myCommands.cpp is added per-target (see OPS_INTERP_TCL_PER_TARGET_SOURCES
# in the top-level CMakeLists.txt) because it declares `theDomain` as
# either Domain or PartitionedDomain depending on _PARALLEL_PROCESSING /
# _PARALLEL_INTERPRETERS, which only the consuming exe knows.

target_sources(OPS_INTERPRETER
    PRIVATE
        Block2D.cpp
        Block3D.cpp
    PUBLIC
        Block2D.h
        Block3D.h
)

target_include_directories(OPS_InterpTcl PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(OPS_INTERPRETER PUBLIC ${CMAKE_CURRENT_LIST_DIR})

