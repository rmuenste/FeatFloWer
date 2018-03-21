# CMakeLists of q2p1_fc_ext application
set(src_cinterface
${CMAKE_SOURCE_DIR}/source/cinterface/cinterface.f90
)

# source files for standad q2p1 
set(src_q2p1_optic_tweezers
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
q2p1_optic_tweezers.f90
)

# Add the application executable
add_executable(q2p1_optic_tweezers ${src_q2p1_optic_tweezers}) 

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

target_link_libraries(q2p1_optic_tweezers 
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(q2p1_optic_tweezers PUBLIC 
${FF_APPLICATION_INCLUDE_PATH} 
)

target_compile_options(q2p1_optic_tweezers PRIVATE ${Fortran_FLAGS})

target_compile_definitions(q2p1_optic_tweezers PUBLIC OPTIC_FORCES=1)

include(${CMAKE_MODULE_PATH}/CreateDataDirectories.cmake)  
