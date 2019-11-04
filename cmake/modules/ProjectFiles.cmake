
#=========================================================================
#                         project directories
#=========================================================================
set(src_pp3d
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/projma.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/coeff.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/indat3d.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/gupwd.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/nsdef.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/bndry.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/conv.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/util.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/avs3d.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/period.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/parq3d.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/xmrout.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/rdparm.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/trsort.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/orsc.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/mgrout.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/inout.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/user.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/optcnl.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/dfkt.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/error.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/bmul.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/diff.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/nu_turb.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/prostp.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/def_feat.f90
  #source/src_pproc/pproc_main.f90
  )

set(postprocessing
  ${CMAKE_SOURCE_DIR}/source/postprocessing/solution_io.f90
  ${CMAKE_SOURCE_DIR}/source/postprocessing/post_utils.f90
  ${CMAKE_SOURCE_DIR}/source/postprocessing/visualization_output.f90
  )

set(Elements
  ${CMAKE_SOURCE_DIR}/source/Elements/e013.f
  ${CMAKE_SOURCE_DIR}/source/Elements/e012.f
  )

set(src_util
  ${CMAKE_SOURCE_DIR}/source/src_util/ReadExtrud3DParameters.f90
  ${CMAKE_SOURCE_DIR}/source/src_util/types.f90
  ${CMAKE_SOURCE_DIR}/source/src_util/OctTreeSearch.f90
  )

set(src_fbm
  ${CMAKE_SOURCE_DIR}/source/src_fbm/fbm_aux.f90
  ${CMAKE_SOURCE_DIR}/source/src_fbm/fbm_main.f90
)

set(src_assemblies
${CMAKE_SOURCE_DIR}/source/assemblies/QuadSc_laplace.f
${CMAKE_SOURCE_DIR}/source/assemblies/QuadSc_diff.f  
${CMAKE_SOURCE_DIR}/source/assemblies/QuadSc_massrho.f
${CMAKE_SOURCE_DIR}/source/assemblies/QuadSc_conv.f
${CMAKE_SOURCE_DIR}/source/assemblies/QuadSc_shear.f
${CMAKE_SOURCE_DIR}/source/assemblies/QuadSc_stress.f
${CMAKE_SOURCE_DIR}/source/assemblies/QuadSc_surftens.f
${CMAKE_SOURCE_DIR}/source/assemblies/QuadSc_gravity.f
${CMAKE_SOURCE_DIR}/source/assemblies/QuadSc_BMatrix.f
${CMAKE_SOURCE_DIR}/source/assemblies/QuadSc_barMmatrix.f
)

set(src_quadLS_app
${src_util}
${src_fbm}
${src_assemblies}
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_solver.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_proj.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_force.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_cylforce.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_torque.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_Sigma_User.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_mpi.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_main.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_var.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_user.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_mg.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_def.f90
${CMAKE_SOURCE_DIR}/source/initialization/app_initialization.f90
)


set(src_mpi
  ${CMAKE_SOURCE_DIR}/source/src_mpi/smooth_mpi.f
  ${CMAKE_SOURCE_DIR}/source/src_mpi/master_mpi.f
  )

set(src_LinSc
  ${CMAKE_SOURCE_DIR}/source/src_LinSc/LinSc_conv.f
  ${CMAKE_SOURCE_DIR}/source/src_LinSc/LinSc_diff.f
  ${CMAKE_SOURCE_DIR}/source/src_LinSc/LinSc_solver.f
  ${CMAKE_SOURCE_DIR}/source/src_LinSc/LinSc_def.f90
  ${CMAKE_SOURCE_DIR}/source/src_LinSc/LinSc_dfkt.f90
  ${CMAKE_SOURCE_DIR}/source/src_LinSc/LinSc_main.f90
  ${CMAKE_SOURCE_DIR}/source/src_LinSc/LinSc_mg.f90
  ${CMAKE_SOURCE_DIR}/source/src_LinSc/LinSc_mpi.f90
  )

set(src_PLin
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_def.f90
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_main.f90
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_matstr.f
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_mass.f
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_user.f90
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_conv.f
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_norm.f
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_flux.f
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_intphase.f
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_intpol.f
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_limiter.f
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_mpi.f
  ${CMAKE_SOURCE_DIR}/source/src_PLin/PLinSc_density.f90
  )

set(src_visco
  ${CMAKE_SOURCE_DIR}/source/src_visco/Visco_integration.f
  ${CMAKE_SOURCE_DIR}/source/src_visco/Visco_solver.f
  ${CMAKE_SOURCE_DIR}/source/src_visco/Visco_main.f90
  ${CMAKE_SOURCE_DIR}/source/src_visco/Visco_def.f90
  )

set(src_mesh
  ${CMAKE_SOURCE_DIR}/source/src_mesh/mesh_refine.f90
  ${CMAKE_SOURCE_DIR}/source/src_mesh/umbrella_smoother.f90
  ${CMAKE_SOURCE_DIR}/source/src_mesh/geometry_processing.f90
  )

# source files for standard q2p1 
set(src_q2p1
  ${CMAKE_SOURCE_DIR}/source/OutputProfiles.f90
  ${CMAKE_SOURCE_DIR}/source/tetra.f90
  ${CMAKE_SOURCE_DIR}/source/Init.f90
  ${CMAKE_SOURCE_DIR}/source/Statistics.f90
  ${CMAKE_SOURCE_DIR}/source/Umbrella.f90
  ${CMAKE_SOURCE_DIR}/source/src_mpi/pp3d_mpi.f90
  ${CMAKE_SOURCE_DIR}/source/ProcCtrl.f90
  ${CMAKE_SOURCE_DIR}/source/UmfpackSolver.f90
  ${CMAKE_SOURCE_DIR}/source/inverse.f90
  ${CMAKE_SOURCE_DIR}/source/Parametrization.f90
  ${CMAKE_SOURCE_DIR}/source/src_ini/iniparser.f90
  ${CMAKE_SOURCE_DIR}/source/distance.f90
  ${CMAKE_SOURCE_DIR}/source/3x3EigenV.f
  )


set(src_ini_c
  ${CMAKE_SOURCE_DIR}/source/src_ini/isdirectory.c
  ${CMAKE_SOURCE_DIR}/source/src_ini/mkdir_recursive.c
  )

set(src_particles
  ${CMAKE_SOURCE_DIR}/source/src_particles/part_input.f90
  ${CMAKE_SOURCE_DIR}/source/src_particles/part_step.f90
  ${CMAKE_SOURCE_DIR}/source/src_particles/part_tracer.f90
  )
#=========================================================================
#         This will create source folders in Visual Studio
#=========================================================================
source_group(src_quadLS FILES ${src_q2p1})
source_group(src_pp3d FILES ${src_pp3d})
source_group(src_mpi FILES ${src_mpi})
source_group(src_PLin FILES ${src_PLin})
source_group(src_LinSc FILES ${src_LinSc})
source_group(src_quadLS FILES ${src_quadLS})
source_group(src_visco FILES ${src_visco})
source_group(src_mesh FILES ${src_mesh})
source_group(src_particles FILES ${src_particles})
source_group(Elements FILES ${Elements})
source_group(postprocessing FILES ${postprocessing})