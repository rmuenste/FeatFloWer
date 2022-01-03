# CMakeLists of q1_scalar application
set(SSE_TARGET_NAME "q1_scalar")

# source files for standad q2p1 
set(src_q1_scalar
app_init.f90
q1_scalar.f90
)

# Add the application executable
add_executable(q1_scalar ${src_q1_scalar}) 

source_group(src_quadLS FILES ${src_q2p1})
source_group(src_pp3d FILES ${src_pp3d})
source_group(src_mpi FILES ${src_mpi})
source_group(src_PLin FILES ${src_PLin})
source_group(src_LinSc FILES ${src_LinSc})
source_group(src_quadLS FILES ${src_quadLS_q1_scalar})
source_group(src_visco FILES ${src_visco})
source_group(src_mesh FILES ${src_mesh})
source_group(src_cinterface FILES ${src_cinterface})
source_group(Elements FILES ${Elements})

target_link_libraries(q1_scalar 
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(q1_scalar PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(q1_scalar PRIVATE ${Fortran_FLAGS})

set_target_properties(q1_scalar PROPERTIES LINKER_LANGUAGE Fortran)

#message(STATUS "q2p1_fc_ext source dir: ${CMAKE_CURRENT_SOURCE_DIR}")
#message(STATUS "q2p1_fc_ext bin dir: ${CMAKE_CURRENT_BINARY_DIR}")

IF(${OUT_OF_SOURCE_BUILD})
  include(${CMAKE_MODULE_PATH}/CreateDataDirectories.cmake)  
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/q1_scalar_start.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

#=========================================================================
#                       INSTALL CONFIGURATION
#=========================================================================
install(TARGETS ${SSE_TARGET_NAME}
    RUNTIME DESTINATION bin/${SSE_TARGET_NAME}
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

set(DEF_DIRS _data _dump _gmv _mesh _ns solution start testresults _vtk)

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
set(E3D_SCRIPTS "${CMAKE_CURRENT_SOURCE_DIR}/q1_scalar_start.py")
install(FILES ${E3D_SCRIPTS} DESTINATION bin/${SSE_TARGET_NAME}
    PERMISSIONS WORLD_EXECUTE OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE)
