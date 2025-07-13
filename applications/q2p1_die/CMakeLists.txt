include(CustomTargets)  


# source files for standad q2p1 
set(src_q2p1_die
# f90getopt.f90
app_init.f90
q2p1_die.f90
)

set(CMAKE_SKIP_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "${ORIGIN}/../lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(SSE_DIE_TARGET_NAME "q2p1_die")

# Add the application executable
add_executable(${SSE_DIE_TARGET_NAME} ${src_q2p1_die}) 

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
  writeCommitHash(${SSE_DIE_TARGET_NAME})
endif(GIT_FOUND)

if(MODULE_COMMAND)
  addCompilerInfo("${SSE_DIE_TARGET_NAME}")
endif(MODULE_COMMAND)

if(Q1P1_SSE_DIE_LINK_STATIC_FORTRAN)
  target_link_libraries(${SSE_DIE_TARGET_NAME} 
    ${FF_APPLICATION_LIBS} -static-libgfortran
    )
else(Q1P1_SSE_DIE_LINK_STATIC_FORTRAN)
  target_link_libraries(${SSE_DIE_TARGET_NAME} 
    ${FF_APPLICATION_LIBS}
    )
endif(Q1P1_SSE_DIE_LINK_STATIC_FORTRAN)

target_include_directories(${SSE_DIE_TARGET_NAME} PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(${SSE_DIE_TARGET_NAME} PRIVATE ${Fortran_FLAGS})

set_target_properties(${SSE_DIE_TARGET_NAME} PROPERTIES LINKER_LANGUAGE Fortran)

add_dependencies(${SSE_DIE_TARGET_NAME} s3d_mesher)

IF(${OUT_OF_SOURCE_BUILD})
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

  # create the die specific directories 
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

  IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/e3d_input")
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/e3d_input" "${CMAKE_CURRENT_BINARY_DIR}/e3d_input" )
  ENDIF()
  
  IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/q2p1_die_start.py")
    file(COPY q2p1_die_start.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  ENDIF()

  file(COPY "${CMAKE_SOURCE_DIR}/tools/e3d_scripts/e3d_start_die.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
 
  file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/conv_check_DIE.sh" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/plot_conv_DIE.gp" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  
ENDIF()

#=========================================================================
#                       INSTALL CONFIGURATION
#=========================================================================
install(TARGETS ${SSE_DIE_TARGET_NAME}
    RUNTIME DESTINATION bin/${SSE_DIE_TARGET_NAME}
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

set(DEF_DIRS _mesh _vtk _data _dump start _1D _hist _prot0 _prot1 _RTD _data_BU e3d_input)

FOREACH(dir ${DEF_DIRS})
  list(APPEND mydirs ${CMAKE_CURRENT_BINARY_DIR}/${dir})
ENDFOREACH()

install(DIRECTORY ${mydirs} DESTINATION bin/${SSE_DIE_TARGET_NAME})

#=========================================================================
#                       INSTALL CONFIGURATION
#=========================================================================
# install the partitioner directory so the python sources can be found by
# the driver script
install(DIRECTORY ${CMAKE_SOURCE_DIR}/tools/partitioner DESTINATION bin/${SSE_DIE_TARGET_NAME})

# install the driver script
set(E3D_SCRIPTS "${CMAKE_SOURCE_DIR}/tools/e3d_scripts/e3d_start_die.py")
list(APPEND E3D_SCRIPTS "${CMAKE_CURRENT_SOURCE_DIR}/conv_check_DIE.sh")
list(APPEND E3D_SCRIPTS "${CMAKE_CURRENT_SOURCE_DIR}/plot_conv_DIE.gp")

install(FILES ${E3D_SCRIPTS} DESTINATION bin/${SSE_DIE_TARGET_NAME}
    PERMISSIONS WORLD_EXECUTE OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE)

#add_test(q2p1-die python3 ./e3d_start_die.py -n 16 -f _ianus/TSE/Conv -x -a 0)
