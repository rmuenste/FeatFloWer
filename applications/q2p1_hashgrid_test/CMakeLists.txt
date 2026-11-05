include(${CMAKE_SOURCE_DIR}/cmake/modules/CustomTargets.cmake)

# q2p1_hashgrid_test - the Layer-2 FBM classification gate.
# Real partitioned mesh, real Q2 DOF sweep, pe's accelerated point query checked against the
# exhaustive one on every DOF of every rank. No fluid solve. See README.md.
set(src_q2p1_hashgrid_test
hashgrid_gate_driver.f90
app_init.f90
q2p1_hashgrid_test.f90
)

add_executable(q2p1_hashgrid_test ${src_q2p1_hashgrid_test})

target_link_libraries(q2p1_hashgrid_test
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(q2p1_hashgrid_test PUBLIC ${FF_APPLICATION_INCLUDE_PATH} ${CMAKE_CURRENT_BINARY_DIR})

target_compile_options(q2p1_hashgrid_test PRIVATE ${Fortran_FLAGS})

set_target_properties(q2p1_hashgrid_test PROPERTIES LINKER_LANGUAGE Fortran)

add_dependencies(q2p1_hashgrid_test metis)

IF(${OUT_OF_SOURCE_BUILD})
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

# Stage the run directory. The mesh is the smallest committed partitioned one, QBOX9Z3
# (unit cube, 9^3 hexes, 3 z-slab subdomains -> mpirun -np 4), taken from the EL tier-2
# terminal-velocity case rather than duplicated here. Its coarse project lands beside the
# binary and its partitions become _mesh/, which is what the deck's MeshFolder expects.
# No libmetis.so is staged: the mesh is pre-partitioned, so the solver never calls the partitioner.
set(GATE_MESH_SRC ${CMAKE_SOURCE_DIR}/applications/q2p1_el_pipeflow/tier2_cases/terminal_velocity/mesh)
add_custom_target(q2p1_hashgrid_test_stage
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/_data
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/_mesh
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          ${CMAKE_SOURCE_DIR}/_data/MG.dat
          ${CMAKE_CURRENT_BINARY_DIR}/_data/MG.dat
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          ${CMAKE_CURRENT_SOURCE_DIR}/_data/q2p1_param.dat
          ${CMAKE_CURRENT_BINARY_DIR}/_data/q2p1_param.dat
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          ${CMAKE_CURRENT_SOURCE_DIR}/pe_user_config.json
          ${CMAKE_CURRENT_BINARY_DIR}/pe_user_config.json
  COMMAND ${CMAKE_COMMAND} -E copy_directory
          ${GATE_MESH_SRC}/quiescent_box_9
          ${CMAKE_CURRENT_BINARY_DIR}/quiescent_box_9
  COMMAND ${CMAKE_COMMAND} -E copy_directory
          ${GATE_MESH_SRC}/partitions
          ${CMAKE_CURRENT_BINARY_DIR}/_mesh
  # FullC0ntact's init_fc_rigid_body reads start/sampleRigidBody.xml and start/data.TXT
  # unconditionally at bootstrap (their values are superseded, the FILES must exist). The
  # copies here come from q2p1_dns_drag, the certified particle lineage; the heat app's set
  # points FullC0ntact at a CGAL config file and a Reiloy mesh and aborts the master.
  COMMAND ${CMAKE_COMMAND} -E copy_directory
          ${CMAKE_CURRENT_SOURCE_DIR}/start
          ${CMAKE_CURRENT_BINARY_DIR}/start
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS q2p1_hashgrid_test metis
  COMMENT "Staging the q2p1_hashgrid_test gate run directory (QBOX9Z3 mesh, deck, pe config)"
)
