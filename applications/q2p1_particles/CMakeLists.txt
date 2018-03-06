# CMakeLists of q2p1_fc_ext application
set(src_cinterface
${CMAKE_SOURCE_DIR}/source/cinterface/cinterface.f90
)

# source files for standad q2p1 
set(src_q2p1_fc_ext
${Elements} 
${src_q2p1} 
${src_pp3d} 
${src_mpi} 
${src_PLin} 
${src_LinSc} 
${src_quadLS_app} 
${src_visco} 
${src_mesh}
${src_cinterface}
${postprocessing}
app_init.f90
q2p1_particles.f90
)

# Add the application executable
add_executable(q2p1_particles ${src_q2p1_fc_ext}) 

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

target_link_libraries(q2p1_particles 
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(q2p1_particles PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(q2p1_particles PRIVATE ${Fortran_FLAGS})

set_target_properties(q2p1_particles PROPERTIES LINKER_LANGUAGE Fortran)

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

  set(DIRECTORYCOPIES _data start)
  set(DEF_DIRECTORIES _dump _gmv _mesh _ns solution testresults _vtk)

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
