#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Element
    PRIVATE
    FourNodeTetrahedron.cpp
    TenNodeTetrahedron.cpp
    TenNodeTetrahedronThermal.cpp
    ThermalBoundaryTri6.cpp
    PUBLIC
    FourNodeTetrahedron.h
    TenNodeTetrahedron.h
    TenNodeTetrahedronThermal.h
    ThermalBoundaryTri6.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
