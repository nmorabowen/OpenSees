#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Actor
    PRIVATE
      Channel.cpp
      HTTP.cpp
      Socket.cpp
      TCP_Socket.cpp
      UDP_Socket.cpp      
    PUBLIC
      Channel.h
      Socket.h
      TCP_Socket.h
      UDP_Socket.h      
)

if(MPI_FOUND)

target_sources(OpenSeesMP
    PRIVATE
      MPI_Channel.cpp
    PUBLIC
      MPI_Channel.h
)

target_sources(OpenSeesSP
    PRIVATE
      MPI_Channel.cpp
    PUBLIC
      MPI_Channel.h
)

# Ladruno Patch 9: OpenSeesPyMP needs the same MPI channel as the Tcl
# parallel exes (MPI_MachineBroker references MPI_Channel::MPI_Channel(int)).
target_sources(OpenSeesPyMP
    PRIVATE
      MPI_Channel.cpp
    PUBLIC
      MPI_Channel.h
)

endif()



target_include_directories(OPS_Actor PUBLIC ${CMAKE_CURRENT_LIST_DIR})
