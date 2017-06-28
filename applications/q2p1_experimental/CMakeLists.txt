# CMakeLists of q2p1_fc_ext application
set(src_quadLS_ext
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_solver.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_laplace.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_diff.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_massrho.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_proj.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_conv.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_shear.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_stress.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_surftens.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_force.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_gravity.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_BMatrix.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_mpi.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_main.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_var.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_user.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_mg.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_def.f90
)

# CMakeLists of q2p1_fc_ext application
set(src_cinterface
${CMAKE_SOURCE_DIR}/source/cinterface/cinterface.f90
)

# source files for standad q2p1 
set(src_q2p1_experimental
${Elements} 
${src_q2p1} 
${src_pp3d} 
${src_mpi} 
${src_PLin} 
${src_LinSc} 
${src_quadLS_ext} 
${src_visco} 
${src_mesh}
${src_cinterface}
app_init.f90
postprocessing.f90
q2p1_experimental.f90
)


# Add the application executable
add_executable(q2p1_experimental ${src_q2p1_experimental}) 

source_group(src_quadLS FILES ${src_q2p1})
source_group(src_pp3d FILES ${src_pp3d})
source_group(src_mpi FILES ${src_mpi})
source_group(src_PLin FILES ${src_PLin})
source_group(src_LinSc FILES ${src_LinSc})
source_group(src_quadLS FILES ${src_quadLS_ext})
source_group(src_visco FILES ${src_visco})
source_group(src_mesh FILES ${src_mesh})
source_group(src_cinterface FILES ${src_cinterface})
source_group(Elements FILES ${Elements})


target_link_libraries(q2p1_experimental 
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(q2p1_experimental PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(q2p1_experimental PRIVATE ${Fortran_FLAGS})

set_target_properties(q2p1_experimental PROPERTIES LINKER_LANGUAGE Fortran)

#message(STATUS "q2p1_fc_ext source dir: ${CMAKE_CURRENT_SOURCE_DIR}")
#message(STATUS "q2p1_fc_ext bin dir: ${CMAKE_CURRENT_BINARY_DIR}")

IF(${OUT_OF_SOURCE_BUILD})
  # check whether the $ENV{Q2P1_MESH_DIR} variable is set
  IF(NOT $ENV{Q2P1_MESH_DIR} STREQUAL "")
    set(DIRECTORYLINKS meshes)
    IF(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/_adc)
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        $ENV{Q2P1_MESH_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}/_adc)
    ENDIF()
  ELSE(NOT $ENV{Q2P1_MESH_DIR} STREQUAL "")
    set(DIRECTORYLINKS _adc meshes)
  ENDIF(NOT $ENV{Q2P1_MESH_DIR} STREQUAL "")

  set(FILELINKS partitioner)
  set(DIRECTORYCOPIES _data start)
  set(DEF_DIRECTORIES _dump _gmv _mesh _ns solution testresults _vtk)

  # establish the file links
  FOREACH(flink ${FILELINKS})
    IF(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${flink})
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink 
        ${CMAKE_CURRENT_SOURCE_DIR}/${flink}
        ${CMAKE_CURRENT_BINARY_DIR}/${flink})
    ENDIF ()    
  ENDFOREACH()

  # create the default directories if they are missing
  FOREACH(dir ${DEF_DIRECTORIES})
    IF(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${dir})
      file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir})
    ENDIF()
  ENDFOREACH()

  # copy the data directories
  FOREACH(dir ${DIRECTORYCOPIES})
    IF(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${dir})
      execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/${dir} ${CMAKE_CURRENT_BINARY_DIR}/${dir})
    ENDIF()
  ENDFOREACH()
ENDIF()
