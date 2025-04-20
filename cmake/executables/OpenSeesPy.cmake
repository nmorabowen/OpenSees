# OpenSeesPy.cmake
# Configuration for OpenSeesPy Python module

# Only configure if Python was found
if(Python_FOUND)
    # OpenSeesPy Module
    add_library(OpenSeesPy SHARED EXCLUDE_FROM_ALL
        ${OPS_SRC_DIR}/interpreter/PythonModule.cpp
        ${OPS_SRC_DIR}/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp
    )

    set_target_properties(OpenSeesPy PROPERTIES PREFIX "")

    target_include_directories(OpenSeesPy PUBLIC ${Python_INCLUDE_DIRS})

    # Link with appropriate libraries
    if(USING_CONAN2)
        target_link_libraries(OpenSeesPy
            OpenSeesLIB 
            OPS_Reliability
            OPS_Recorder
            OPS_Numerics
            hdf5::hdf5
            ZLIB::ZLIB
            Eigen3::Eigen 
            ${Python_LIBRARIES}
        )
    else()
        target_link_libraries(OpenSeesPy
            OpenSeesLIB 
            OPS_Reliability
            OPS_Recorder
            OPS_Numerics 
            ${HDF5_LIBRARIES} 
            ${CONAN_LIBS} 
            ${Python_LIBRARIES}
        )
    endif()
endif()