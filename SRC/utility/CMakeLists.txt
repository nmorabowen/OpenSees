#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Utilities
    PRIVATE
    Timer.cpp 
    FileIter.cpp 
    File.cpp 
    SimulationInformation.cpp 
    StringContainer.cpp
    PeerNGA.cpp
    LadrunoThreads.cpp
    # NB LadrunoParallelBuild.cpp is deliberately NOT here -- it must be
    # compiled once per TARGET to see _PARALLEL_*; it is listed in the root
    # CMakeLists' OPS_MPI_PER_TARGET_SOURCES (Ladruno WP-107 red-team B2).
    PUBLIC
    Timer.h
    FileIter.h
    File.h
    SimulationInformation.h
    StringContainer.h
    LadrunoThreads.h
)

target_include_directories(OPS_Utilities PUBLIC ${CMAKE_CURRENT_LIST_DIR})

# Ladruno Stack Profiler (timing core always; HDF5 writer guarded inside).
add_subdirectory(profiler)
