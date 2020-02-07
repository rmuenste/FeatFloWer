# CMakeLists of q2p1_devel application

# source files for standad q2p1 
set(src_q2p1_devel
app_init.f90
q2p1_devel.f90
)


# Add the application executable
add_executable(q2p1_devel ${src_q2p1_devel}) 

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

target_link_libraries(q2p1_devel 
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(q2p1_devel PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(q2p1_devel PRIVATE ${Fortran_FLAGS})

set_target_properties(q2p1_devel PROPERTIES LINKER_LANGUAGE Fortran)

IF(${OUT_OF_SOURCE_BUILD})
  include(${CMAKE_MODULE_PATH}/CreateDataDirectories.cmake)  
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

file(COPY q2p1_devel.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
