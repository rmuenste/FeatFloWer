# CMakeLists of stream application

# source files for standard q2p1 
set(src_stream
stream.f90
)

# Add the application executable
add_executable(stream ${src_stream}) 

target_link_libraries(stream 
#  ${FF_APPLICATION_LIBS}
  )

target_include_directories(stream PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(stream PRIVATE ${Fortran_FLAGS})

set_target_properties(stream PROPERTIES LINKER_LANGUAGE Fortran)

IF(${OUT_OF_SOURCE_BUILD})
  include(${FF_MODULE_PATH}/CreateDataDirectories.cmake)  
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()


#=========================================================================
#                       INSTALL CONFIGURATION Q2P1_GENDIE PYTHON
#=========================================================================
install(TARGETS stream
    RUNTIME DESTINATION bin/q2p1_gendie
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
    
#set(stream "${CMAKE_CURRENT_SOURCE_DIR}/param.cfg")
# list(APPEND E3D_SCRIPTS "${CMAKE_CURRENT_SOURCE_DIR}/conv_check.sh")

install(FILES ${stream} DESTINATION bin/q2p1_gendie
    PERMISSIONS WORLD_EXECUTE OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE)
