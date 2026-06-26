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
    terminal_velocity
    saffman_lift
    straddling_conservation)
  add_custom_target(q2p1_el_pipeflow_tier2_${_el_tier2_case}_stage
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_CURRENT_SOURCE_DIR}/tier2_cases/${_el_tier2_case}/q2p1_param.dat
            ${CMAKE_CURRENT_BINARY_DIR}/_data/q2p1_param.dat
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_CURRENT_SOURCE_DIR}/tier2_cases/${_el_tier2_case}/example.json
            ${CMAKE_CURRENT_BINARY_DIR}/example.json
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_CURRENT_SOURCE_DIR}/tier2_cases/${_el_tier2_case}/particles.xyz
            ${CMAKE_CURRENT_BINARY_DIR}/particles.xyz
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS q2p1_el_pipeflow_stage
    COMMENT "Staging q2p1_el_pipeflow Tier 2 case ${_el_tier2_case}"
  )
endforeach()
