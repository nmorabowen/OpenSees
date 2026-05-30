#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Recorder
    PRIVATE
      DamageRecorder.cpp
      DatastoreRecorder.cpp
      DriftRecorder.cpp
      ElementRecorder.cpp
      ElementRecorderRMS.cpp
      EnergyBalanceRecorder.cpp
      EnvelopeDriftRecorder.cpp
      EnvelopeElementRecorder.cpp
      EnvelopeNodeRecorder.cpp
      GSA_Recorder.cpp
      GmshRecorder.cpp
      MaxNodeDispRecorder.cpp
      NodeRecorder.cpp
      NodeRecorderRMS.cpp
      NormElementRecorder.cpp
      NormEnvelopeElementRecorder.cpp
      PatternRecorder.cpp
      PVDRecorder.cpp      
      Recorder.cpp
      RemoveRecorder.cpp
      VTK_Recorder.cpp
    PUBLIC
      DamageRecorder.h
      DatastoreRecorder.h
      DriftRecorder.h
      ElementRecorder.h
      ElementRecorderRMS.h
      EnergyBalanceRecorder.h
      EnvelopeDriftRecorder.h
      EnvelopeElementRecorder.h
      EnvelopeNodeRecorder.h
      GSA_Recorder.h
      GmshRecorder.h
      MaxNodeDispRecorder.h
      NodeRecorder.h
      NodeRecorderRMS.h
      NormElementRecorder.h
      NormEnvelopeElementRecorder.h
      PatternRecorder.h
      PVDRecorder.h      
      Recorder.h
      RemoveRecorder.h
      VTK_Recorder.h
)


if(HDF5_FOUND)
   if (HDF5_VERSION VERSION_GREATER_EQUAL 1.12.0)
      target_sources(OPS_Recorder
        PRIVATE
          MPCORecorder.cpp
          MPCORecorderLadruno.cpp
          MPCOL_NodeResults.cpp
          MPCOL_ElementResults.cpp
          MPCOL_DomainResults.cpp
          MPCOL_Sinks.cpp
          VTKHDF_Recorder.cpp
        PUBLIC
         MPCORecorder.h
         MPCORecorderLadruno.h
         MPCOL_Types.h
         MPCOL_Hdf5.h
         MPCOL_ResultIO.h
         MPCOL_NodeResults.h
         MPCOL_ElementResults.h
         MPCOL_DomainResults.h
         MPCOL_Sinks.h
         EnergyBalanceKernel.h
         VTKHDF_Recorder.h
      )
   endif()
endif()

target_sources(OPS_Paraview
    PRIVATE
      PVDRecorder.cpp
    PUBLIC
      PVDRecorder.h
)

target_sources(OPS_Graphics
    PRIVATE
      AlgorithmIncrements.cpp
      FilePlotter.cpp
    PUBLIC
      AlgorithmIncrements.h
      FilePlotter.h
)

add_subdirectory(response)

