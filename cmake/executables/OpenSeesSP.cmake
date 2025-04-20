# OpenSeesSP.cmake
# Configuration for OpenSeesSP (Parallel) executable

# Only configure if MPI was found
if(MPI_FOUND)
    add_executable(OpenSeesSP EXCLUDE_FROM_ALL 
      ${OPS_SRC_DIR}/tcl/mpiMain.cpp
      ${OPS_SRC_DIR}/tcl/tclMain.cpp
      ${OPS_SRC_DIR}/tcl/commands.cpp
      ${OPS_SRC_DIR}/domain/domain/partitioned/PartitionedDomain.cpp
      ${OPS_SRC_DIR}/domain/domain/partitioned/PartitionedDomainEleIter.cpp
      ${OPS_SRC_DIR}/domain/domain/partitioned/PartitionedDomainSubIter.cpp
      ${OPS_SRC_DIR}/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp
      ${OPS_SRC_DIR}/actor/machineBroker/MPI_MachineBroker.cpp
      ${OPS_SRC_DIR}/system_of_eqn/linearSOE/diagonal/MPIDiagonalSOE.cpp
      ${OPS_SRC_DIR}/system_of_eqn/linearSOE/diagonal/MPIDiagonalSolver.cpp
      ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsSOE.cpp
      ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsSolver.cpp
      ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsParallelSOE.cpp
      ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.cpp
      ${OPS_SRC_DIR}/system_of_eqn/linearSOE/sparseGEN/DistributedSuperLU.cpp
      ${OPS_SRC_DIR}/system_of_eqn/linearSOE/sparseGEN/DistributedSparseGenColLinSOE.cpp
      ${OPS_SRC_DIR}/system_of_eqn/linearSOE/sparseGEN/DistributedSparseGenRowLinSOE.cpp
      ${OPS_SRC_DIR}/domain/subdomain/ActorSubdomain.cpp
      ${OPS_SRC_DIR}/domain/subdomain/ShadowSubdomain.cpp
    )

    target_include_directories(OpenSeesSP PRIVATE ${MUMPS_DIR}/_deps/mumps-src/include ${MPI_CXX_INCLUDE_DIRS})
    target_compile_options(OpenSeesSP PRIVATE ${MPI_CXX_COMPILE_FLAGS})

    if (DEFINED OPENMPI)
        target_compile_definitions(OpenSeesSP 
        PUBLIC _PARALLEL_PROCESSING ${MUMPS_FLAG} _OPENMPI)
    else()
        target_compile_definitions(OpenSeesSP 
        PUBLIC _PARALLEL_PROCESSING ${MUMPS_FLAG})
    endif()

    # Link libraries based on Conan version
    if(USING_CONAN2)
        target_link_libraries(OpenSeesSP 
           OPS_InterpTcl 
           OpenSeesLIB
           OPS_Reliability
           OPS_ReliabilityTcl
           OPS_Recorder
           METIS
           SUPERLU_DIST
           OPS_Numerics
           ${MUMPS_LIBRARIES}
           ${CMAKE_DL_LIBS} 
           hdf5::hdf5
           tcl::tcl
           ZLIB::ZLIB
           Eigen3::Eigen
           ${MPI_CXX_LIBRARIES}
           ${SCALAPACK_LIBRARIES}
           ${MPI_Fortran_LIBRARIES}       
           ${MPI_CXX_LINK_FLAGS}
        )
    else()
        target_link_libraries(OpenSeesSP 
           OPS_InterpTcl 
           OpenSeesLIB
           OPS_Reliability
           OPS_ReliabilityTcl
           OPS_Recorder      
           METIS
           SUPERLU_DIST
           OPS_Numerics
           ${MUMPS_LIBRARIES}
           ${CMAKE_DL_LIBS} 
           ${HDF5_LIBRARIES} 
           ${CONAN_LIBS}
           ${MPI_CXX_LIBRARIES}
           ${SCALAPACK_LIBRARIES}
           ${MPI_Fortran_LIBRARIES}       
           ${MPI_CXX_LINK_FLAGS}
        )
    endif()
endif()