IF(NOT ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  set(OUT_OF_SOURCE_BUILD True)
  message(STATUS "Configuring out-of source build")
else()
  set(OUT_OF_SOURCE_BUILD False)
ENDIF()

IF(${OUT_OF_SOURCE_BUILD})
  # check whether the $ENV{Q2P1_MESH_DIR} variable is set
  IF(NOT $ENV{Q2P1_MESH_DIR} STREQUAL "")
    message(STATUS "Found mesh directory : $ENV{Q2P1_MESH_DIR}")
  ELSE(NOT $ENV{Q2P1_MESH_DIR} STREQUAL "")
    message(AUTHOR_WARNING "Mesh directory variable not found.")
  ENDIF(NOT $ENV{Q2P1_MESH_DIR} STREQUAL "")
endif()


IF(${OUT_OF_SOURCE_BUILD})

  # check whether the $ENV{Q2P1_MESH_DIR} variable is set
  IF(NOT $ENV{Q2P1_MESH_DIR} STREQUAL "")
    set(DIRECTORYLINKS meshes)
    IF(NOT EXISTS ${CMAKE_BINARY_DIR}/_adc)
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        $ENV{Q2P1_MESH_DIR}
        ${CMAKE_BINARY_DIR}/_adc)
    ENDIF()
  ELSE(NOT $ENV{Q2P1_MESH_DIR} STREQUAL "")
    set(DIRECTORYLINKS _adc meshes)
  ENDIF(NOT $ENV{Q2P1_MESH_DIR} STREQUAL "")

  set(FILELINKS partitioner)
  set(DIRECTORYCOPIES _data start)
  set(DEF_DIRECTORIES _dump _gmv _mesh _ns solution testresults _vtk)

  file(COPY ${CMAKE_SOURCE_DIR}/_data DESTINATION ${CMAKE_BINARY_DIR})
  file(COPY ${CMAKE_SOURCE_DIR}/start DESTINATION ${CMAKE_BINARY_DIR})

  file(COPY ${CMAKE_SOURCE_DIR}/testcases DESTINATION ${CMAKE_BINARY_DIR})

  file(COPY ${CMAKE_SOURCE_DIR}/testcases/fac/test1 DESTINATION ${CMAKE_BINARY_DIR})
  file(COPY ${CMAKE_SOURCE_DIR}/testcases/fac_nnewt/test2 DESTINATION ${CMAKE_BINARY_DIR})
  file(COPY ${CMAKE_SOURCE_DIR}/testcases/fallingparticle/test3 DESTINATION ${CMAKE_BINARY_DIR})
  file(COPY ${CMAKE_SOURCE_DIR}/testcases/fac_visco/test4 DESTINATION ${CMAKE_BINARY_DIR})
  file(COPY ${CMAKE_SOURCE_DIR}/testcases/fac_visco_elastic/test5 DESTINATION ${CMAKE_BINARY_DIR})

  # establish the file links
  FOREACH(flink ${FILELINKS})
    IF(NOT EXISTS ${CMAKE_BINARY_DIR}/${flink})
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink 
        ${CMAKE_SOURCE_DIR}/${flink}
        ${CMAKE_BINARY_DIR}/${flink})
    ENDIF ()    
  ENDFOREACH()

  # create the default directories if they are missing
  FOREACH(dir ${DEF_DIRECTORIES})
    IF(NOT EXISTS ${CMAKE_BINARY_DIR}/${dir})
      file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${dir})
    ENDIF()
  ENDFOREACH()

  # copy the data directories
  FOREACH(dir ${DIRECTORYCOPIES})
    IF(NOT EXISTS ${CMAKE_BINARY_DIR}/${dir})
      message(STATUS "${CMAKE_BINARY_DIR}/${dir} does not exist, copying directory")
      execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/${dir} ${CMAKE_BINARY_DIR}/${dir})
    ENDIF()
  ENDFOREACH()

ENDIF()

# continue to create the neccessary directory structure:
# we know the _dump directory exists now, we need to create the subdirectories 
set(SUBDIRS 01 02 03 04 05 06 07 08 09 10)
FOREACH(sd ${SUBDIRS})
  IF(NOT EXISTS ${CMAKE_BINARY_DIR}/_dump/${sd})
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/_dump/${sd})
  ENDIF()  
ENDFOREACH()
