# CMakeLists of s3d_mesher application
set(src_quadLS_ext
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_solver.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_proj.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_force.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_torque.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_Sigma_User.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_mpi.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_main.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_var.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_user.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_mg.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_def.f90
)

# CMakeLists of q2p1_fc_ext application
set(src_cinterface
${CMAKE_SOURCE_DIR}/source/cinterface/cinterface.f90
)

# source files for standard q2p1 
set(src_s3d_mesher
${Elements} 
${src_assemblies} 
${src_q2p1} 
${src_pp3d} 
${src_mpi} 
${src_ini} 
${src_util} 
${src_PLin} 
${src_LinSc} 
${src_quadLS_ext} 
${src_visco} 
${src_mesh}
${src_cinterface}
${postprocessing}
s3d_mesher.f90
)

# Add the application executable
add_executable(s3d_mesher ${src_s3d_mesher}) 

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

target_link_libraries(s3d_mesher 
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(s3d_mesher PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(s3d_mesher PRIVATE ${Fortran_FLAGS})

set_target_properties(s3d_mesher PROPERTIES LINKER_LANGUAGE Fortran)

