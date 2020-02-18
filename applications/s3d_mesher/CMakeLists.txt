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

#=========================================================================
#                       INSTALL CONFIGURATION
#=========================================================================
install(TARGETS s3d_mesher
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

