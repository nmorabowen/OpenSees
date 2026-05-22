
# Ladruno Patch 9: OpenSeesCommands.cpp and OpenSeesMiscCommands.cpp are
# deliberately NOT listed here. They carry the only _PARALLEL_* guards in
# this directory and are built per-define into the OPS_InterpPyCmds{,_MP,_SP}
# STATIC libs (see the "Patch 9" block in the top-level CMakeLists.txt).
# Every consumer of OPS_INTERPRETER links exactly one of those libs.
target_sources(OPS_INTERPRETER
  PRIVATE
    DL_Interpreter.cpp
    OpenSeesUniaxialMaterialCommands.cpp
    OpenSeesElementCommands.cpp
    OpenSeesTimeSeriesCommands.cpp
    OpenSeesPatternCommands.cpp
    OpenSeesSectionCommands.cpp
    OpenSeesOutputCommands.cpp
    OpenSeesCrdTransfCommands.cpp
    OpenSeesDampingCommands.cpp
    OpenSeesBeamIntegrationCommands.cpp
    OpenSeesNDMaterialCommands.cpp
    OpenSeesParameterCommands.cpp
    OpenSeesFrictionModelCommands.cpp
    OpenSeesReliabilityCommands.cpp
    OpenSeesNDTestCommands.cpp
    OpenSeesIGACommands.cpp
)


target_sources(OpenSeesPy
   PRIVATE
     PythonModule.cpp
     PythonWrapper.cpp
     PythonStream.cpp
)

if(MPI_FOUND)
# Ladruno Patch 9: OpenSeesPyMP Python entry sources. PythonMPIModule.cpp is
# a thin shim that #includes PythonModule.cpp with OPS_PY_MODULE_NAME
# overridden to "openseesmp"; PythonModule.cpp itself is NOT compiled as its
# own TU for this target (only via the #include) to avoid a duplicate
# PyInit_/moduledef. The OpenSeesPyMP target itself defines
# _PARALLEL_INTERPRETERS.
target_sources(OpenSeesPyMP
   PRIVATE
     PythonMPIModule.cpp
     PythonWrapper.cpp
     PythonStream.cpp
)
# (OpenSeesPySP entry sources were removed — Python has no SP subsystem;
#  see Ladruno Patch 9 docs.)
endif()
