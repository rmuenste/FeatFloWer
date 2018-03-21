# CMakeLists of laplace application
set(src_quadLS_mumps_test
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_solver.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_laplace.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_diff.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_massrho.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_proj.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_conv.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_shear.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_stress.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_surftens.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_force.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_gravity.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_BMatrix.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_mpi.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_main.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_var.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_user.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_mg.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_def.f90
)

# CMakeLists of laplace application
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
${src_quadLS_mumps_test} 
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

include(${CMAKE_MODULE_PATH}/CreateDataDirectories.cmake)  
