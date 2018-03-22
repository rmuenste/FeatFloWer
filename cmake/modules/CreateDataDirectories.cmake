function(createDefaultDirectories bdir sdir)
  # check whether the $ENV{Q2P1_MESH_DIR} variable is set
  if(NOT WIN32)
    IF(NOT $ENV{Q2P1_MESH_DIR} STREQUAL "")
      set(DIRECTORYLINKS meshes)
      IF(NOT EXISTS ${bdir}/_adc)
        execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
          $ENV{Q2P1_MESH_DIR}
          ${bdir}/_adc)
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
      IF(NOT EXISTS ${bdir}/${flink})
        execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink 
          ${CMAKE_SOURCE_DIR}/${flink}
          ${bdir}/${flink})
      ENDIF ()    
    ENDFOREACH()
  endif(NOT WIN32)

  # create the default directories if they are missing
  FOREACH(dir ${DEF_DIRECTORIES})
    IF(NOT EXISTS ${bdir}/${dir})
      file(MAKE_DIRECTORY ${bdir}/${dir})
    ENDIF()
  ENDFOREACH()

  # copy the data directories
  FOREACH(dir ${DIRECTORYCOPIES})
    IF(NOT EXISTS ${bdir}/${dir})
      execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${sdir}/${dir} ${bdir}/${dir})
    ENDIF()
  ENDFOREACH()
endfunction(createDefaultDirectories bdir sdir)
