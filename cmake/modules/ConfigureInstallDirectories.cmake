install(TARGETS q2p1
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

set(DEF_DIRS _mesh _vtk _data _dump _adc start testresults testcases)

FOREACH(dir ${DEF_DIRS})
  list(APPEND mydirs ${CMAKE_BINARY_DIR}/${dir})
ENDFOREACH()

#message(installdirs "${installdirs}")
#message(FATAL_ERROR "${mydirs}")

install(DIRECTORY ${mydirs} DESTINATION bin)
install(PROGRAMS partitioner DESTINATION bin)

# create the default directories if they are missing
#FOREACH(dir ${DEF_DIRECTORIES})
#  IF(NOT EXISTS ${CMAKE_INSTALL_PREFIX}/bin/${dir})
#    file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/bin/${dir})
#  ENDIF()
#ENDFOREACH()
#
## copy the data directories
#FOREACH(dir ${DIRECTORYCOPIES})
#  IF(NOT EXISTS ${CMAKE_INSTALL_PREFIX}/bin/${dir})
#    file(COPY ${CMAKE_SOURCE_DIR}/${dir} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
#  ENDIF()
#ENDFOREACH()
#
#IF(NOT $ENV{Q2P1_MESH_DIR} STREQUAL "")
#  IF(NOT EXISTS ${CMAKE_INSTALL_PREFIX}/bin/_adc)
#    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink 
#      $ENV{Q2P1_MESH_DIR}
#      ${CMAKE_INSTALL_PREFIX}/bin/_adc)
#  ENDIF()    
#ENDIF()
