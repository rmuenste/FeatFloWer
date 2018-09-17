
# CMakeLists of q2p1_fc_ext application
set(src_cinterface
${CMAKE_SOURCE_DIR}/source/cinterface/cinterface.f90
)

# source files for standad q2p1 
set(src_q2p1_cc
${Elements} 
${src_util}
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

target_link_libraries(q2p1_cc 
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(q2p1_cc PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(q2p1_cc PRIVATE ${Fortran_FLAGS})

set_target_properties(q2p1_cc PROPERTIES LINKER_LANGUAGE Fortran)

IF(${OUT_OF_SOURCE_BUILD})
  include(${CMAKE_MODULE_PATH}/CreateDataDirectories.cmake)  
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

file(COPY q2p1_cc.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set(NPROCS 0)
if(NOT $ENV{NUM_PROCS} STREQUAL "")
  set(NPROCS $ENV{NUM_PROCS})
endif(NOT $ENV{NUM_PROCS} STREQUAL "")

add_test(q2p1-fac3D-CC python ./q2p1_cc.py -n ${PROCS})
