#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Actor
    PRIVATE
      ChannelAddress.cpp
    PUBLIC
      ChannelAddress.h
)

if(MPI_FOUND)

target_sources(OpenSeesMP
  PRIVATE
    MPI_ChannelAddress.cpp
  PUBLIC
    MPI_ChannelAddress.h
)

target_sources(OpenSeesSP
  PRIVATE
    MPI_ChannelAddress.cpp
  PUBLIC
    MPI_ChannelAddress.h
)

# Ladruno Patch 9: OpenSeesPyMP needs the same MPI channel address as the
# Tcl parallel exes (MPI_Channel pulls MPI_ChannelAddress).
target_sources(OpenSeesPyMP
  PRIVATE
    MPI_ChannelAddress.cpp
  PUBLIC
    MPI_ChannelAddress.h
)

endif()

target_include_directories(OPS_Actor PUBLIC ${CMAKE_CURRENT_LIST_DIR})

