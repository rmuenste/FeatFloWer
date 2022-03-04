# CMakeLists of q1_scalar_partfil application

# source files for standad q2p1 
set(src_q1_scalar_partfil
app_init.f90
q1_scalar_partfil.f90
)

set(Q1SSE_TARGET_NAME "q1_scalar_partfil")

# Add the application executable
add_executable(${Q1SSE_TARGET_NAME} ${src_q1_scalar_partfil}) 

source_group(src_quadLS FILES ${src_q2p1})
source_group(src_pp3d FILES ${src_pp3d})
source_group(src_mpi FILES ${src_mpi})
source_group(src_PLin FILES ${src_PLin})
source_group(src_LinSc FILES ${src_LinSc})
source_group(src_quadLS FILES ${src_quadLS_q1_scalar_partfil})
source_group(src_visco FILES ${src_visco})
source_group(src_mesh FILES ${src_mesh})
source_group(src_cinterface FILES ${src_cinterface})
source_group(Elements FILES ${Elements})

target_link_libraries(${Q1SSE_TARGET_NAME} 
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(${Q1SSE_TARGET_NAME} PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(${Q1SSE_TARGET_NAME} PRIVATE ${Fortran_FLAGS})

set_target_properties(${Q1SSE_TARGET_NAME} PROPERTIES LINKER_LANGUAGE Fortran)

IF(${OUT_OF_SOURCE_BUILD})
  include(${CMAKE_MODULE_PATH}/CreateDataDirectories.cmake)  
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

#=========================================================================
#                       INSTALL CONFIGURATION
#=========================================================================
install(TARGETS ${Q1SSE_TARGET_NAME}
    RUNTIME DESTINATION bin/q2p1_sse_partfil
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
