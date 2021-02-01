include(${CMAKE_MODULE_PATH}/CustomTargets.cmake)  


# source files for standad q2p1 
set(src_q2p1_sse
# f90getopt.f90
app_init.f90
q2p1_sse.f90
)

set(CMAKE_SKIP_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "${ORIGIN}/../lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(SSE_TARGET_NAME "q2p1_sse")

# Add the application executable
add_executable(${SSE_TARGET_NAME} ${src_q2p1_sse}) 

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

if(GIT_FOUND)
  writeCommitHash(${SSE_TARGET_NAME})
endif(GIT_FOUND)

if(MODULE_COMMAND)
  addCompilerInfo("${SSE_TARGET_NAME}")
endif(MODULE_COMMAND)

if(Q1P1_SSE_LINK_STATIC_FORTRAN)
  target_link_libraries(${SSE_TARGET_NAME} 
    ${FF_APPLICATION_LIBS} -static-libgfortran
    )
else(Q1P1_SSE_LINK_STATIC_FORTRAN)
  target_link_libraries(${SSE_TARGET_NAME} 
    ${FF_APPLICATION_LIBS}
    )
endif(Q1P1_SSE_LINK_STATIC_FORTRAN)

target_include_directories(${SSE_TARGET_NAME} PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(${SSE_TARGET_NAME} PRIVATE ${Fortran_FLAGS})

set_target_properties(${SSE_TARGET_NAME} PROPERTIES LINKER_LANGUAGE Fortran)

add_dependencies(${SSE_TARGET_NAME} s3d_mesher)

IF(${OUT_OF_SOURCE_BUILD})
  include(${CMAKE_MODULE_PATH}/CreateDataDirectories.cmake)  
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

  # create the sse specific directories 
  IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/_1D")
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/_1D")
  ENDIF()

  IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/_hist")
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/_hist")
  ENDIF()

  IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/_prot0")
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/_prot0")
  ENDIF()

  IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/_prot1")
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/_prot1")
  ENDIF()

  IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/_RTD")
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/_RTD")
  ENDIF()

  # copy the data directories
  IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/_data_BU")
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/_data_BU" "${CMAKE_CURRENT_BINARY_DIR}/_data_BU" )
  ENDIF()

  IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/_ianus")
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/_ianus" "${CMAKE_CURRENT_BINARY_DIR}/_ianus" )
  ENDIF()

  IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/q2p1_sse_start.py")
    file(COPY q2p1_sse_start.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  ENDIF()

  file(COPY "${CMAKE_SOURCE_DIR}/tools/e3d_scripts/e3d_start.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
 
  file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/conv_check.sh" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

ENDIF()

#=========================================================================
#                       INSTALL CONFIGURATION
#=========================================================================
install(TARGETS ${SSE_TARGET_NAME}
    RUNTIME DESTINATION bin/${SSE_TARGET_NAME}
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

set(DEF_DIRS _mesh _vtk _data _dump start _1D _hist _prot0 _prot1 _RTD _data_BU _ianus)

FOREACH(dir ${DEF_DIRS})
  list(APPEND mydirs ${CMAKE_CURRENT_BINARY_DIR}/${dir})
ENDFOREACH()

install(DIRECTORY ${mydirs} DESTINATION bin/${SSE_TARGET_NAME})

#=========================================================================
#                       INSTALL CONFIGURATION PYTHON
#=========================================================================
# install the partitioner directory so the python sources can be found by
# the driver script
install(DIRECTORY ${CMAKE_SOURCE_DIR}/tools/partitioner DESTINATION bin/${SSE_TARGET_NAME})

# install the driver script
set(E3D_SCRIPTS "${CMAKE_SOURCE_DIR}/tools/e3d_scripts/e3d_start.py")
list(APPEND E3D_SCRIPTS "${CMAKE_CURRENT_SOURCE_DIR}/conv_check.sh")

install(FILES ${E3D_SCRIPTS} DESTINATION bin/${SSE_TARGET_NAME}
    PERMISSIONS WORLD_EXECUTE OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE)

add_test(q2p1-sse python3 ./e3d_start.py -n 16 -f _ianus/TSE/Conv -x -a 0)



#=========================================================================
#                       INSTALL CONFIGURATION Q2P1_GENDIE
#=========================================================================
install(TARGETS ${SSE_TARGET_NAME}
    RUNTIME DESTINATION bin/q2p1_gendie
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

set(DEF_DIRS _mesh _vtk _data _dump start _1D _hist _prot0 _prot1 _RTD _ianus)

FOREACH(dir ${DEF_DIRS})
  list(APPEND mydirs ${CMAKE_CURRENT_BINARY_DIR}/${dir})
ENDFOREACH()
list(APPEND mydirs ${CMAKE_CURRENT_SOURCE_DIR}/_data_BU)

install(DIRECTORY ${mydirs} DESTINATION bin/q2p1_gendie)

#=========================================================================
#                       INSTALL CONFIGURATION Q2P1_GENDIE PYTHON
#=========================================================================
# install the partitioner directory so the python sources can be found by
# the driver script
install(DIRECTORY ${CMAKE_SOURCE_DIR}/tools/partitioner DESTINATION bin/q2p1_gendie)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/testcases/gendie/VEKA_S DESTINATION bin/q2p1_gendie)

# install the driver script
set(E3D_SCRIPTS "${CMAKE_SOURCE_DIR}/tools/e3d_scripts/e3d_start.py")
list(APPEND E3D_SCRIPTS "${CMAKE_CURRENT_SOURCE_DIR}/conv_check.sh")
list(APPEND E3D_SCRIPTS "${CMAKE_SOURCE_DIR}/tools/computeAreas.py")
list(APPEND E3D_SCRIPTS "${CMAKE_SOURCE_DIR}/applications/mesh_ref/RunnerScript.sh")
list(APPEND E3D_SCRIPTS "${CMAKE_BINARY_DIR}/applications/mesh_ref/meshref")
list(APPEND E3D_SCRIPTS "${CMAKE_BINARY_DIR}/applications/mesh_ref/param_BU.cfg")

install(FILES ${E3D_SCRIPTS} DESTINATION bin/q2p1_gendie
    PERMISSIONS WORLD_EXECUTE OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE)

install(DIRECTORY "${CMAKE_SOURCE_DIR}/applications/mesh_ref/PATCHES" DESTINATION bin/q2p1_gendie)
