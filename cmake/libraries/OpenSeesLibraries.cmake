
# OpenSeesLibraries.cmake
# Define OpenSees libraries

# Set TCL_LIBRARIES
set(TCL_LIBRARIES ${TCL_LIBRARY})

# OpenSees OBJECT LIBRARIES
add_library(OPS_Matrix             OBJECT)
add_library(OPS_Actor              OBJECT)
add_library(OPS_ObjectBroker       OBJECT)
add_library(OPS_Handler            OBJECT)
add_library(OPS_Recorder           OBJECT)
add_library(OPS_Reliability        OBJECT)
add_library(OPS_Tagged             OBJECT)
add_library(OPS_Utilities          OBJECT)
add_library(OPS_ModelBuilder       OBJECT)
add_library(OPS_Domain             OBJECT)
add_library(OPS_SysOfEqn           OBJECT)
add_library(OPS_Analysis           OBJECT)
add_library(OPS_ConvergenceTest    OBJECT)
add_library(OPS_Thermal            OBJECT)
add_library(OPS_Element            OBJECT)
add_library(OPS_ElementFortran     OBJECT)
add_library(OPS_Material           OBJECT)
add_library(OPS_MaterialFortran    OBJECT)
add_library(OPS_Damage             OBJECT)
add_library(OPS_Database           OBJECT)
add_library(OPS_INTERPRETER        OBJECT)

# Optional Extensions
add_library(OPS_Paraview           OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Renderer           OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Graphics           OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Graphics_Default   OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Graphics_GL        OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_PFEM               OBJECT EXCLUDE_FROM_ALL)

# Tcl Interpreter Library & special include and link directives
add_library(OPS_InterpTcl          STATIC)

target_compile_definitions(OPS_InterpTcl PUBLIC _TCL85)
if(NOT USING_CONAN)
  target_include_directories(OPS_InterpTcl PUBLIC ${TCL_INCLUDE_PATH})
  target_link_libraries(OPS_InterpTcl PRIVATE ${TCL_LIBRARIES})
endif()

# OpenSees G3 Library - a slimmed down OpenSees
add_library(G3) 

target_link_libraries(G3
    coordTransformation
    damping
    OPS_Matrix
    OPS_Analysis 
    OPS_ModelBuilder
    OPS_Domain
    OPS_ConvergenceTest
    OPS_Element
    OPS_Material
    OPS_Recorder
    OPS_Handler
    OPS_SysOfEqn
    OPS_Tagged 
    OPS_Utilities
    graph
    OPS_Actor 
    OPS_ObjectBroker
    OPS_Numerics
    OPS_INTERPRETER
)

# OpenSees Libray
add_library(OpenSeesLIB  EXCLUDE_FROM_ALL)

target_link_libraries(OpenSeesLIB
    ${OPS_Extension_List}
    OPS_Actor 
    OPS_Analysis 
    OPS_ConvergenceTest
    OPS_Damage
    OPS_Database
    OPS_Domain
    graph
    OPS_Element
    OPS_ElementFortran    
    OPS_Handler
    OPS_INTERPRETER
    OPS_Material
    OPS_MaterialFortran    
    OPS_Material_YieldSurface
    OPS_Matrix
    OPS_ModelBuilder
    OPS_ObjectBroker
    OPS_PFEM
    OPS_Recorder
    coordTransformation
    damping
    OPS_Renderer
    OPS_Section_Repres
    OPS_Section_YieldSurface
    OPS_SysOfEqn
    OPS_Thermal 
    OPS_OS_Specific_libs
    OPS_Tagged 
    OPS_Utilities
)