# CMakeLists of q2p1_fc_ext application

# source files for standad q2p1 
set(src_q2p1_fc_ext
../../source/Elements/e013.f
../../source/Elements/e012.f
../../source/src_pp3d/projma.f
../../source/src_pp3d/coeff.f
../../source/src_pp3d/indat3d.f
../../source/src_pp3d/gupwd.f
../../source/src_pp3d/nsdef.f
../../source/src_pp3d/bndry.f
../../source/src_pp3d/conv.f
../../source/src_pp3d/util.f
../../source/src_pp3d/avs3d.f
../../source/src_pp3d/period.f
../../source/src_pp3d/parq3d.f
../../source/src_pp3d/xmrout.f
../../source/src_pp3d/rdparm.f
../../source/src_pp3d/trsort.f
../../source/src_pp3d/orsc.f
../../source/src_pp3d/mgrout.f
../../source/src_pp3d/inout.f
../../source/src_pp3d/user.f
../../source/src_pp3d/optcnl.f
../../source/src_pp3d/dfkt.f
../../source/src_pp3d/error.f
../../source/src_pp3d/bmul.f
../../source/src_pp3d/diff.f
../../source/src_pp3d/nu_turb.f
../../source/src_pp3d/prostp.f
../../source/src_quadLS/QuadSc_solver.f
../../source/src_quadLS/QuadSc_laplace.f
../../source/src_quadLS/QuadSc_diff.f
../../source/src_quadLS/QuadSc_massrho.f
../../source/src_quadLS/QuadSc_proj.f
../../source/src_quadLS/QuadSc_conv.f
../../source/src_quadLS/QuadSc_shear.f
../../source/src_quadLS/QuadSc_stress.f
../../source/src_quadLS/QuadSc_force.f90
../../source/src_quadLS/QuadSc_gravity.f
../../source/src_quadLS/QuadSc_BMatrix.f
../../source/src_mpi/smooth_mpi.f
../../source/src_mpi/master_mpi.f
../../source/src_LinSc/LinSc_conv.f
../../source/src_LinSc/LinSc_diff.f
../../source/src_LinSc/LinSc_solver.f
../../source/src_pp3d/def_feat.f90
../../source/OutputProfiles.f90
../../source/src_quadLS/QuadSc_mpi.f90
../../source/src_quadLS/QuadSc_main.f90
../../source/src_quadLS/QuadSc_var.f90
../../source/src_quadLS/QuadSc_user.f90
../../source/src_quadLS/QuadSc_mg.f90
../../source/src_quadLS/QuadSc_def.f90
../../source/tetra.f90
../../source/Init.f90
../../source/Statistics.f90
../../source/Umbrella.f90
../../source/src_mpi/pp3d_mpi.f90
../../source/src_LinSc/LinSc_def.f90
../../source/src_LinSc/LinSc_dfkt.f90
../../source/src_LinSc/LinSc_main.f90
../../source/src_PLin/PLinSc_def.f90
../../source/src_PLin/PLinSc_main.f90
../../source/src_PLin/PLinSc_matstr.f
../../source/src_PLin/PLinSc_mass.f
../../source/src_PLin/PLinSc_user.f90
../../source/src_PLin/PLinSc_conv.f
../../source/src_PLin/PLinSc_norm.f
../../source/src_PLin/PLinSc_flux.f
../../source/src_PLin/PLinSc_intphase.f
../../source/src_PLin/PLinSc_intpol.f
../../source/src_PLin/PLinSc_limiter.f
../../source/src_PLin/PLinSc_mpi.f
../../source/src_PLin/PLinSc_density.f90
../../source/ProcCtrl.f90
../../source/UmfpackSolver.f90
../../source/inverse.f90
../../source/Parametrization.f90
../../source/distance.f90
../../source/3x3EigenV.f
../../source/src_visco/Visco_integration.f
../../source/src_visco/Visco_solver.f
../../source/src_visco/Visco_main.f90
../../source/src_visco/Visco_def.f90
../../source/src_mesh/mesh_refine.f90
../../source/cinterface/cinterface.f90
app_init.f90
postprocessing.f90
q2p1_fc_ext.f90
)


# Add the application executable
add_executable(q2p1_fc_ext ${src_q2p1_fc_ext}) 

target_link_libraries(q2p1_fc_ext 
  amd
  umfpack4
  feat2d
  feat3d
  sysutils 
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  inshape3dcore
  Utility
  Math
  OpenMeshCore 
  OpenMeshTools
  ${LIBRT_LIBRARY}
  ${MPI_Fortran_LIBRARIES}
  )

target_include_directories(q2p1_fc_ext PUBLIC ${MPI_Fortran_INCLUDE_PATH})

target_compile_options(q2p1_fc_ext PRIVATE ${Fortran_FLAGS})

#message(STATUS "q2p1_fc_ext source dir: ${CMAKE_CURRENT_SOURCE_DIR}")
#message(STATUS "q2p1_fc_ext bin dir: ${CMAKE_CURRENT_BINARY_DIR}")

IF(${OUT_OF_SOURCE_BUILD})
  # check whether the $ENV{Q2P1_MESH_DIR} variable is set
  IF(NOT $ENV{Q2P1_MESH_DIR} STREQUAL "")
    set(DIRECTORYLINKS meshes)
    IF(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/_adc)
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        $ENV{Q2P1_MESH_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}/_adc)
    ENDIF()
  ELSE(NOT $ENV{Q2P1_MESH_DIR} STREQUAL "")
    set(DIRECTORYLINKS _adc meshes)
  ENDIF(NOT $ENV{Q2P1_MESH_DIR} STREQUAL "")

  set(FILELINKS partitioner)
  set(DIRECTORYCOPIES _data start)
  set(DEF_DIRECTORIES _dump _gmv _mesh _ns solution testresults _vtk)

  # establish the file links
  FOREACH(flink ${FILELINKS})
    IF(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${flink})
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink 
        ${CMAKE_SOURCE_DIR}/${flink}
        ${CMAKE_CURRENT_BINARY_DIR}/${flink})
    ENDIF ()    
  ENDFOREACH()

  # create the default directories if they are missing
  FOREACH(dir ${DEF_DIRECTORIES})
    IF(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${dir})
      file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir})
    ENDIF()
  ENDFOREACH()

  # copy the data directories
  FOREACH(dir ${DIRECTORYCOPIES})
    IF(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${dir})
      execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/${dir} ${CMAKE_CURRENT_BINARY_DIR}/${dir})
    ENDIF()
  ENDFOREACH()
ENDIF()
