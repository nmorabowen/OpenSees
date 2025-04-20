# GlobalIncludes.cmake
# Set up global include directories

# Include paths to main abstract classes
include_directories(
  ${OPS_SRC_DIR}
  ${OPS_SRC_DIR}/matrix
  ${OPS_SRC_DIR}/handler
  ${OPS_SRC_DIR}/database
  ${OPS_SRC_DIR}/element
  ${OPS_SRC_DIR}/element/forceBeamColumn
  ${OPS_SRC_DIR}/element/nonlinearBeamColumn/matrixutil
  ${OPS_SRC_DIR}/coordTransformation
  ${OPS_SRC_DIR}/damping
  ${OPS_SRC_DIR}/tagged
  ${OPS_SRC_DIR}/tagged/storage
  ${OPS_SRC_DIR}/recorder
  ${OPS_SRC_DIR}/renderer
  ${OPS_SRC_DIR}/damage
  ${OPS_SRC_DIR}/recorder/response
  ${OPS_SRC_DIR}/material
  ${OPS_SRC_DIR}/material/section
  ${OPS_SRC_DIR}/material/uniaxial
  ${OPS_SRC_DIR}/material/nD
  ${OPS_SRC_DIR}/graph/graph
  ${OPS_SRC_DIR}/graph/numberer
  ${OPS_SRC_DIR}/graph/partitioner
  ${OPS_SRC_DIR}/domain/component
  ${OPS_SRC_DIR}/domain/domain
  ${OPS_SRC_DIR}/domain/subdomain
  ${OPS_SRC_DIR}/domain/load
  ${OPS_SRC_DIR}/domain/loadBalancer
  ${OPS_SRC_DIR}/domain/pattern
  ${OPS_SRC_DIR}/domain/groundMotion
  ${OPS_SRC_DIR}/domain/node
  ${OPS_SRC_DIR}/domain/constraints
  ${OPS_SRC_DIR}/domain/region
  ${OPS_SRC_DIR}/analysis/algorithm
  ${OPS_SRC_DIR}/analysis/dof_grp
  ${OPS_SRC_DIR}/analysis/fe_ele
  ${OPS_SRC_DIR}/analysis/algorithm/equiSolnAlgo
  ${OPS_SRC_DIR}/analysis/algorithm/eigenAlgo
  ${OPS_SRC_DIR}/analysis/algorithm/domainDecompAlgo
  ${OPS_SRC_DIR}/analysis/analysis
  ${OPS_SRC_DIR}/analysis/integrator
  ${OPS_SRC_DIR}/analysis/handler
  ${OPS_SRC_DIR}/analysis/numberer
  ${OPS_SRC_DIR}/analysis/model
  ${OPS_SRC_DIR}/convergenceTest
  ${OPS_SRC_DIR}/modelbuilder
  ${OPS_SRC_DIR}/system_of_eqn
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE
  ${OPS_SRC_DIR}/system_of_eqn/eigenSOE
  ${OPS_SRC_DIR}/actor/actor
  ${OPS_SRC_DIR}/actor/channel
  ${OPS_SRC_DIR}/actor/objectBroker
  ${OPS_SRC_DIR}/actor/message
  ${OPS_BUNDLED_DIR}/eigenAPI
)

# TODO: Temporary fix; include all directories with .h files through glob
file(GLOB_RECURSE ops_include_files CONFIGURE_DEPENDS 
    ${OPS_SRC_DIR}/*.h
    ${OPS_BUNDLED_DIR}/*.h
)
set(OPS_INCLUDE_DIRS "")
foreach(filepath ${ops_include_files})
    get_filename_component(dir_path ${filepath} PATH)
    set(OPS_INCLUDE_DIRS ${OPS_INCLUDE_DIRS} ${dir_path})
endforeach()
list(REMOVE_DUPLICATES OPS_INCLUDE_DIRS)

include_directories(${OPS_INCLUDE_DIRS})