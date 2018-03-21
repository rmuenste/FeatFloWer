IF(${OUT_OF_SOURCE_BUILD})
  # check whether the $ENV{Q2P1_MESH_DIR} variable is set
  if(NOT WIN32)
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
  endif(NOT WIN32)

  set(FILELINKS partitioner)
  set(DIRECTORYCOPIES _data start)
  set(DEF_DIRECTORIES _dump _gmv _mesh _ns solution testresults _vtk)

  if(NOT WIN32)
    # establish the file links
    FOREACH(flink ${FILELINKS})
      IF(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${flink})
        execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink 
          ${CMAKE_SOURCE_DIR}/${flink}
          ${CMAKE_CURRENT_BINARY_DIR}/${flink})
      ENDIF ()    
    ENDFOREACH()
  endif(NOT WIN32)

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
