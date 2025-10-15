# source files for standad q2p1 
set(src_q2p1_fc_ext
app_init.f90
q2p1_geom_test.f90
)

# Add the application executable
add_executable(q2p1_geom_test ${src_q2p1_fc_ext}) 

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

target_link_libraries(q2p1_geom_test 
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(q2p1_geom_test PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(q2p1_geom_test PRIVATE ${Fortran_FLAGS})

target_compile_options(q2p1_geom_test PRIVATE -DSKIP_DYNAMICS=1)

set_target_properties(q2p1_geom_test PROPERTIES LINKER_LANGUAGE Fortran)

#add_dependencies(q2p1_geom_test CGAL_PRJ)

IF(${OUT_OF_SOURCE_BUILD})
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()
