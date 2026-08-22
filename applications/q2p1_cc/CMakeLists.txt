
# CC-specific sources. Core FEATFLOW sources are provided by
# FF_APPLICATION_LIBS; compiling them again here creates duplicate COMMON
# blocks and module state in the executable.
set(src_q2p1_cc
app_init.f90
assemblies_cc.f
iso_assemblies.f
postprocessing.f90
lin_transport_cc.f90
q2p1_def_cc.f90
q2p1_mg_cc.f90
q2p1_cc_umfpacksolver.f90
q2p1_cc.f90
q2p1_transport_cc.f90
q2p1_var_newton.f90
)

# Add the application executable
add_executable(q2p1_cc ${src_q2p1_cc}) 
add_dependencies(q2p1_cc metis)
add_custom_command(TARGET q2p1_cc POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          $<TARGET_FILE:metis> ${CMAKE_CURRENT_BINARY_DIR}/libmetis.so)

source_group(src_quadLS FILES ${src_q2p1})
source_group(src_pp3d FILES ${src_pp3d})
source_group(src_mpi FILES ${src_mpi})
source_group(src_PLin FILES ${src_PLin})
source_group(src_LinSc FILES ${src_LinSc})
source_group(src_quadLS FILES ${src_quadLS_ext})
source_group(src_visco FILES ${src_visco})
source_group(src_mesh FILES ${src_mesh})
source_group(Elements FILES ${Elements})

target_link_libraries(q2p1_cc 
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(q2p1_cc PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(q2p1_cc PRIVATE ${Fortran_FLAGS})

set_target_properties(q2p1_cc PROPERTIES LINKER_LANGUAGE Fortran)

IF(${OUT_OF_SOURCE_BUILD})
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/_data/q2p1_param.dat
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/_data)
file(COPY ${CMAKE_SOURCE_DIR}/_data/MG.dat
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/_data)
file(COPY ${CMAKE_SOURCE_DIR}/tools/partitioner
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/tools/PyPartitioner.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY q2p1_cc.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
