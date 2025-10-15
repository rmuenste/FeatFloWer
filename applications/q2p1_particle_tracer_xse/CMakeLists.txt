

# source files for standad q2p1 
set(src_q2p1_particle_tracer_xse
app_init.f90
q2p1_particle_tracer.f90
)


# Add the application executable
add_executable(q2p1_particle_tracer_xse ${src_q2p1_particle_tracer_xse}) 

source_group(src_quadLS FILES ${src_q2p1})
source_group(src_pp3d FILES ${src_pp3d})
source_group(src_mpi FILES ${src_mpi})
source_group(src_PLin FILES ${src_PLin})
source_group(src_LinSc FILES ${src_LinSc})
source_group(src_quadLS FILES ${src_quadLS_ext})
source_group(src_visco FILES ${src_visco})
source_group(src_mesh FILES ${src_mesh})
source_group(src_particles FILES ${src_particles})
source_group(src_cinterface FILES ${src_cinterface})
source_group(Elements FILES ${Elements})

target_link_libraries(q2p1_particle_tracer_xse 
  ${FF_APPLICATION_LIBS}
  ff_particles
  )

target_include_directories(q2p1_particle_tracer_xse PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(q2p1_particle_tracer_xse PRIVATE ${Fortran_FLAGS})

set_target_properties(q2p1_particle_tracer_xse PROPERTIES LINKER_LANGUAGE Fortran)

IF(${OUT_OF_SOURCE_BUILD})
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

#file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
#file(COPY ${CMAKE_SOURCE_DIR}/tools/dashboard/q2p1_ctest_start.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

#add_test(q2p1-fac-newt python ./q2p1_ctest_start.py)

#=========================================================================
#                       INSTALL CONFIGURATION Q2P1_GENDIE PYTHON
#=========================================================================
install(TARGETS q2p1_particle_tracer_xse
    RUNTIME DESTINATION bin/q2p1_gendie
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
    
set(E3D_SCRIPTS "${CMAKE_SOURCE_DIR}/tools/e3d_scripts/config_particles.dat")

install(FILES ${E3D_SCRIPTS} DESTINATION bin/q2p1_gendie/_data
    PERMISSIONS WORLD_EXECUTE OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE)

