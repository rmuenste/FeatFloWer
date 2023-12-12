# CMakeLists of autoparam application

# source files for standard q2p1 
set(src_autoparam
meshproc_def.f90
meshproc_main.f90
)

# Add the application executable
add_executable(autoparam ${src_autoparam}) 

target_link_libraries(autoparam 
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(autoparam PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(autoparam PRIVATE ${Fortran_FLAGS})

set_target_properties(autoparam PROPERTIES LINKER_LANGUAGE Fortran)

IF(${OUT_OF_SOURCE_BUILD})
  include(${FF_MODULE_PATH}/CreateDataDirectories.cmake)  
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Mesh)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Mesh/PAR)
file(COPY Mesh/mesh.tri DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Mesh)
file(COPY param_vtu_mesher.cfg DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

#=========================================================================
#                       INSTALL CONFIGURATION HEAT PYTHON
#=========================================================================
install(TARGETS autoparam
    RUNTIME DESTINATION bin/heat
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
    
set(AUTOPARAM "${CMAKE_CURRENT_SOURCE_DIR}/param_vtu_mesher.cfg")

install(FILES ${AUTOPARAM} DESTINATION bin/heat
    PERMISSIONS WORLD_EXECUTE OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE)

#=========================================================================
#                       INSTALL CONFIGURATION Q2P1_GENDIE PYTHON
#=========================================================================
install(TARGETS autoparam
    RUNTIME DESTINATION bin/q2p1_gendie
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
    
set(AUTOPARAM "${CMAKE_CURRENT_SOURCE_DIR}/param_vtu_mesher.cfg")

install(FILES ${AUTOPARAM} DESTINATION bin/q2p1_gendie
    PERMISSIONS WORLD_EXECUTE OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE)
