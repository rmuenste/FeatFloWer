# CMakeLists of mumps application
set(src_cinterface
${CMAKE_SOURCE_DIR}/source/cinterface/cinterface.f90
)

# source files for standad q2p1 
set(src_mumps_test
${Elements} 
${src_q2p1} 
${src_pp3d} 
${src_mpi} 
${src_PLin} 
${src_LinSc} 
${src_quadLS_app} 
${src_visco} 
${src_mesh}
${src_cinterface}
${postprocessing}
app_init.f90
postprocessing.f90
mumps_test.f90
)
   

# Add the application executable
add_executable(mumps_test ${src_mumps_test}) 

source_group(src_quadLS FILES ${src_q2p1})
source_group(src_pp3d FILES ${src_pp3d})
source_group(src_mpi FILES ${src_mpi})
source_group(src_PLin FILES ${src_PLin})
source_group(src_LinSc FILES ${src_LinSc})
source_group(src_quadLS FILES ${src_quadLS_mumps_test})
source_group(src_visco FILES ${src_visco})
source_group(src_mesh FILES ${src_mesh})
source_group(src_cinterface FILES ${src_cinterface})
source_group(Elements FILES ${Elements})

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
