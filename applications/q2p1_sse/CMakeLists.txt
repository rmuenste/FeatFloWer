# CMakeLists of q2p1_fc_ext application
set(src_quadLS_ext
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
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_torque.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_Sigma_User.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_gravity.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_BMatrix.f
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

# source files for standad q2p1 
set(src_q2p1_sse
${Elements} 
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
app_init.f90
q2p1_sse.f90
)

# Add the application executable
add_executable(q2p1_sse ${src_q2p1_sse}) 

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

target_link_libraries(q2p1_sse 
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(q2p1_sse PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(q2p1_sse PRIVATE ${Fortran_FLAGS})

set_target_properties(q2p1_sse PROPERTIES LINKER_LANGUAGE Fortran)

IF(${OUT_OF_SOURCE_BUILD})
  include(${CMAKE_MODULE_PATH}/CreateDataDirectories.cmake)  
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

  # create the sse specific directories 
    IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/_1D")
      file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/_1D")
    ENDIF()

    IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/_hist")
      file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/_hist")
    ENDIF()

    IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/q2p1_sse_start.py")
      file(COPY q2p1_sse_start.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    ENDIF()

ENDIF()

