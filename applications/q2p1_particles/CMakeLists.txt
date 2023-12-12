# CMakeLists of q2p1_fc_ext application

# source files for standad q2p1 
set(src_q2p1_fc_ext
app_init.f90
q2p1_particles.f90
)

# Add the application executable
add_executable(q2p1_particles ${src_q2p1_fc_ext}) 

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

target_link_libraries(q2p1_particles 
  ${FF_APPLICATION_LIBS} pthread
  )

target_include_directories(q2p1_particles PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(q2p1_particles PRIVATE ${Fortran_FLAGS})

set_target_properties(q2p1_particles PROPERTIES LINKER_LANGUAGE Fortran)

#message(STATUS "q2p1_fc_ext source dir: ${CMAKE_CURRENT_SOURCE_DIR}")
#message(STATUS "q2p1_fc_ext bin dir: ${CMAKE_CURRENT_BINARY_DIR}")

IF(${OUT_OF_SOURCE_BUILD})
  include(${FF_MODULE_PATH}/CreateDataDirectories.cmake)  
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()
