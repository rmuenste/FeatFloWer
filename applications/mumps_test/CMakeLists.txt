
# source files for standad q2p1 
set(src_mumps_test
mumps_test.f90
)

# Add the application executable
add_executable(mumps_test ${src_mumps_test}) 

target_link_libraries(mumps_test 
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(mumps_test PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(mumps_test PRIVATE ${Fortran_FLAGS})

set_target_properties(mumps_test PROPERTIES LINKER_LANGUAGE Fortran)

IF(${OUT_OF_SOURCE_BUILD})
  include(${CMAKE_MODULE_PATH}/CreateDataDirectories.cmake)  
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()
