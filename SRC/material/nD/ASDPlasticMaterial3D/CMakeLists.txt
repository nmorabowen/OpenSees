#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Material
    PRIVATE
      # ASDPlasticMaterial3DGlobals.cpp
      OPS_AllASDPlasticMaterial3Ds.cpp
    PUBLIC
      AllASDPlasticMaterial3Ds.h
      ASDPlasticMaterial3D.h
      ElasticityBase.h
      YieldFunctionBase.h
      ASDPlasticMaterial3DGlobals.h
      ASDPlasticMaterial3DTraits.h
      PlasticFlowBase.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

# /bigobj is set globally for MSVC at the top-level CMakeLists — see the note
# there. Needed by OPS_AllASDPlasticMaterial3Ds.cpp which instantiates >64K
# COFF sections worth of templated material specializations.

add_subdirectory(ElasticityModels)
add_subdirectory(PlasticFlowDirections)
add_subdirectory(YieldFunctions)
