set(src_quadLS_mesher
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_Sigma_User.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_var.f90
)

# source files for standard q2p1 
set(src_s3d_mesher
${src_quadLS_mesher} 
s3d_mesher.f90
)

set(CMAKE_SKIP_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "${ORIGIN}/../lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Add the application executable
add_executable(s3d_mesher ${src_s3d_mesher}) 

  target_link_libraries(s3d_mesher 
    ff_util
    )

target_include_directories(s3d_mesher PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(s3d_mesher PRIVATE ${Fortran_FLAGS})

set_target_properties(s3d_mesher PROPERTIES LINKER_LANGUAGE Fortran)

add_custom_command(
  TARGET s3d_mesher 
  POST_BUILD 
  COMMAND ${CMAKE_COMMAND} -E copy 
 ${CMAKE_CURRENT_BINARY_DIR}/s3d_mesher
 ${CMAKE_BINARY_DIR}/applications/q2p1_sse
 ${CMAKE_BINARY_DIR}/applications/q2p1_die
 ${CMAKE_BINARY_DIR}/applications/heat
)

#=========================================================================
#                       INSTALL CONFIGURATION
#=========================================================================
install(TARGETS s3d_mesher
    RUNTIME DESTINATION bin/q2p1_sse
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

install(TARGETS s3d_mesher
    RUNTIME DESTINATION bin/q2p1_die
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

install(TARGETS s3d_mesher
    RUNTIME DESTINATION bin/heat
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

