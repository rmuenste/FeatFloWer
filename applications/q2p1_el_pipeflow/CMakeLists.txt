include(${CMAKE_SOURCE_DIR}/cmake/modules/CustomTargets.cmake)

set(src_q2p1_el_pipeflow
  ${CMAKE_CURRENT_SOURCE_DIR}/app_init.f90
  q2p1_el_pipeflow.f90
)

add_executable(q2p1_el_pipeflow ${src_q2p1_el_pipeflow})

target_link_libraries(q2p1_el_pipeflow
  ${FF_APPLICATION_LIBS}
  inshape3dcore
)
target_include_directories(q2p1_el_pipeflow PUBLIC ${FF_APPLICATION_INCLUDE_PATH})
target_compile_options(q2p1_el_pipeflow PRIVATE ${Fortran_FLAGS})
target_compile_definitions(q2p1_el_pipeflow PRIVATE BUILD_Q2P1_EL_PIPEFLOW)

set_target_properties(q2p1_el_pipeflow PROPERTIES LINKER_LANGUAGE Fortran)

if(${OUT_OF_SOURCE_BUILD})
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

  if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/example.json")
    file(COPY
      "${CMAKE_SOURCE_DIR}/applications/q2p1_el_frozen_trace/example.json"
      DESTINATION "${CMAKE_CURRENT_BINARY_DIR}"
    )
  endif()

  if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/unit_cube_27_case")
    file(COPY
      "${CMAKE_SOURCE_DIR}/applications/q2p1_el_frozen_trace/unit_cube_27_case"
      DESTINATION "${CMAKE_CURRENT_BINARY_DIR}"
    )
  endif()

  if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/_data/q2p1_param.dat")
    file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/_data/q2p1_param.dat"
      DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/_data"
    )
  endif()
endif()

add_custom_target(q2p1_el_pipeflow_stage
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          ${CMAKE_BINARY_DIR}/extern/libraries/metis-4.0.3/Lib/libmetis.so
          ${CMAKE_CURRENT_BINARY_DIR}/libmetis.so
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          ${CMAKE_SOURCE_DIR}/_data/MG.dat
          ${CMAKE_CURRENT_BINARY_DIR}/_data/MG.dat
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          ${CMAKE_CURRENT_SOURCE_DIR}/_data/q2p1_param.dat
          ${CMAKE_CURRENT_BINARY_DIR}/_data/q2p1_param.dat
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          ${CMAKE_SOURCE_DIR}/applications/q2p1_el_frozen_trace/example.json
          ${CMAKE_CURRENT_BINARY_DIR}/example.json
  COMMAND ${CMAKE_COMMAND} -E copy_directory
          ${CMAKE_SOURCE_DIR}/applications/q2p1_el_frozen_trace/unit_cube_27_case
          ${CMAKE_CURRENT_BINARY_DIR}/unit_cube_27_case
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS q2p1_el_pipeflow metis
  COMMENT "Staging q2p1_el_pipeflow runtime files in build directory"
)

foreach(_el_tier2_case
    momentum_conservation
    momentum_conservation_long
    momentum_conservation_semi
    terminal_velocity
    saffman_lift
    straddling_conservation)
  set(_el_tier2_case_src ${CMAKE_CURRENT_SOURCE_DIR}/tier2_cases/${_el_tier2_case})
  if("${_el_tier2_case}" STREQUAL "straddling_conservation" OR
     "${_el_tier2_case}" STREQUAL "momentum_conservation_long" OR
     "${_el_tier2_case}" STREQUAL "momentum_conservation_semi" OR
     "${_el_tier2_case}" STREQUAL "saffman_lift")
    set(_el_tier2_mesh_src ${CMAKE_CURRENT_SOURCE_DIR}/tier2_cases/momentum_conservation/mesh)
  else()
    set(_el_tier2_mesh_src ${_el_tier2_case_src}/mesh)
  endif()
  # Self-contained staging: depends only on the binary (not q2p1_el_pipeflow_stage,
  # which stages the frozen-trace mesh and copies a build-specific libmetis that is
  # vestigial here — the EL cases use pre-partitioned meshes and are not linked to
  # metis). Runtime directories (_data, _mesh, ...) are created at configure time by
  # createDefaultDirectories above.
  add_custom_target(q2p1_el_pipeflow_tier2_${_el_tier2_case}_stage
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_SOURCE_DIR}/_data/MG.dat
            ${CMAKE_CURRENT_BINARY_DIR}/_data/MG.dat
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${_el_tier2_case_src}/q2p1_param.dat
            ${CMAKE_CURRENT_BINARY_DIR}/_data/q2p1_param.dat
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${_el_tier2_case_src}/example.json
            ${CMAKE_CURRENT_BINARY_DIR}/example.json
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${_el_tier2_case_src}/particles.xyz
            ${CMAKE_CURRENT_BINARY_DIR}/particles.xyz
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS q2p1_el_pipeflow
    COMMENT "Staging q2p1_el_pipeflow Tier 2 case ${_el_tier2_case}"
  )

  # Self-contained mesh: each case carries its own coarse project (mesh/quiescent_box_9)
  # and pre-partitioned submeshes (mesh/partitions/<MeshFolder>). The runtime
  # PartitionReader reads _mesh/<MeshFolder>/sub<rank>/GRID0001.tri directly (no
  # ./partitioner step), so pre-placing _mesh is sufficient. The committed source dir is
  # named "partitions" (not "_mesh") so it is not caught by the _mesh/* gitignore rule.
  # Only cases that have a mesh/ dir get this staging; cases still being set up stay
  # config-only.
  if(EXISTS "${_el_tier2_mesh_src}")
    add_custom_command(TARGET q2p1_el_pipeflow_tier2_${_el_tier2_case}_stage POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_directory
              ${_el_tier2_mesh_src}/quiescent_box_9
              ${CMAKE_CURRENT_BINARY_DIR}/quiescent_box_9
      COMMAND ${CMAKE_COMMAND} -E copy_directory
              ${_el_tier2_mesh_src}/partitions
              ${CMAKE_CURRENT_BINARY_DIR}/_mesh
      COMMENT "Staging self-contained mesh for ${_el_tier2_case}"
    )
  endif()
endforeach()

# Validation-campaign cases (validation_cases/<case>): same self-contained staging
# pattern as the Tier-2 cases. Each case carries its own coarse project folder
# (mesh/<project>) and pre-partitioned submeshes (mesh/partitions/<MeshFolder>);
# a case may share another case's mesh via a "mesh" symlink (v3_ss_frozen ->
# pipe_hp_check).
foreach(_el_val_case
    pipe_hp_check
    v3_ss_frozen
    v1b_tencate_settling
    v2_rz_settling
    v3_ss_coupled
    unit_periodic_contact
    fourway_shear)
  set(_el_val_case_src ${CMAKE_CURRENT_SOURCE_DIR}/validation_cases/${_el_val_case})
  add_custom_target(q2p1_el_pipeflow_val_${_el_val_case}_stage
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_SOURCE_DIR}/_data/MG.dat
            ${CMAKE_CURRENT_BINARY_DIR}/_data/MG.dat
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${_el_val_case_src}/q2p1_param.dat
            ${CMAKE_CURRENT_BINARY_DIR}/_data/q2p1_param.dat
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${_el_val_case_src}/example.json
            ${CMAKE_CURRENT_BINARY_DIR}/example.json
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${_el_val_case_src}/particles.xyz
            ${CMAKE_CURRENT_BINARY_DIR}/particles.xyz
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS q2p1_el_pipeflow
    COMMENT "Staging q2p1_el_pipeflow validation case ${_el_val_case}"
  )
  if(EXISTS "${_el_val_case_src}/mesh")
    file(GLOB _el_val_projects RELATIVE ${_el_val_case_src}/mesh
         ${_el_val_case_src}/mesh/*)
    foreach(_el_val_proj ${_el_val_projects})
      if(NOT "${_el_val_proj}" STREQUAL "partitions")
        add_custom_command(TARGET q2p1_el_pipeflow_val_${_el_val_case}_stage POST_BUILD
          COMMAND ${CMAKE_COMMAND} -E copy_directory
                  ${_el_val_case_src}/mesh/${_el_val_proj}
                  ${CMAKE_CURRENT_BINARY_DIR}/${_el_val_proj}
          COMMENT "Staging coarse project ${_el_val_proj} for ${_el_val_case}"
        )
      endif()
    endforeach()
    add_custom_command(TARGET q2p1_el_pipeflow_val_${_el_val_case}_stage POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_directory
              ${_el_val_case_src}/mesh/partitions
              ${CMAKE_CURRENT_BINARY_DIR}/_mesh
      COMMENT "Staging partitioned mesh for ${_el_val_case}"
    )
  endif()
endforeach()
