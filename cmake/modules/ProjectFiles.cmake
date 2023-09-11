#=========================================================================
#                         project libraries
#=========================================================================

#*************************************************************************
set(src_cinterface
${CMAKE_SOURCE_DIR}/source/cinterface/cinterface.f90
)

#=========================================================================
#                     CInterface Library Source
#=========================================================================
add_library(ff_cinterface ${src_cinterface})
target_link_libraries(ff_cinterface ff_util ff_mesh)
target_include_directories(ff_cinterface PUBLIC ${FF_APPLICATION_INCLUDE_PATH})
target_compile_options(ff_cinterface PUBLIC ${Fortran_FLAGS})

set(postprocessing
  ${CMAKE_SOURCE_DIR}/source/postprocessing/solution_io.f90
  ${CMAKE_SOURCE_DIR}/source/postprocessing/post_utils.f90
  ${CMAKE_SOURCE_DIR}/source/postprocessing/visualization_output.f90
  )

#=========================================================================
#                    Postprocessing Library Source
#=========================================================================
add_library(ff_postprocessing ${postprocessing})
target_link_libraries(ff_postprocessing ff_cinterface ff_util ff_mesh)
target_include_directories(ff_postprocessing PUBLIC ${FF_APPLICATION_INCLUDE_PATH})
target_compile_options(ff_postprocessing PUBLIC ${Fortran_FLAGS})

set(src_ini_aux
  ${CMAKE_SOURCE_DIR}/source/src_util/getpid_wrapper.c
  ${CMAKE_SOURCE_DIR}/source/src_ini/isdirectory.c
  ${CMAKE_SOURCE_DIR}/source/src_ini/mkdir_recursive.c
)

set(src_ini_c
  ${CMAKE_SOURCE_DIR}/source/src_ini/iniparser.f90
)

#=========================================================================
#                         Ini Library Source
#=========================================================================
add_library(ff_ini_aux ${src_ini_aux})
add_library(ff_ini_c ${src_ini_c})
target_link_libraries(ff_ini_c ff_ini_aux)
target_include_directories(ff_ini_c PUBLIC ${FF_APPLICATION_INCLUDE_PATH})
target_compile_options(ff_ini_c PUBLIC ${Fortran_FLAGS})

set(src_pp3d
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/projma.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/coeff.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/indat3d.f
  ${CMAKE_SOURCE_DIR}/source/src_pp3d/inout.f
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

#=========================================================================
#                       PP3D Library Source
#=========================================================================


set(src_mpi
  ${CMAKE_SOURCE_DIR}/source/src_mpi/smooth_mpi.f
  ${CMAKE_SOURCE_DIR}/source/src_mpi/master_mpi.f
  )

set(src_util
  ${src_pp3d}
  ${CMAKE_SOURCE_DIR}/source/Parametrization.f90
  ${src_mpi}
  ${CMAKE_SOURCE_DIR}/source/src_mpi/pp3d_mpi.f90
  ${CMAKE_SOURCE_DIR}/source/src_util/min_sphere.f90
  ${CMAKE_SOURCE_DIR}/source/src_util/OutputMumpsMatrices.f90
  ${CMAKE_SOURCE_DIR}/source/src_mpi/parentcomm.f90
  ${CMAKE_SOURCE_DIR}/source/src_util/FolderManagement.f90
  ${CMAKE_SOURCE_DIR}/source/src_util/MPI_DumpOutputProfiles.f90
  ${CMAKE_SOURCE_DIR}/source/src_util/ReadExtrud3DParameters.f90
  ${CMAKE_SOURCE_DIR}/source/src_util/types.f90
  ${CMAKE_SOURCE_DIR}/source/src_util/OctTreeSearch.f90
  ${CMAKE_SOURCE_DIR}/source/src_util/get_pid.f90
  ${CMAKE_SOURCE_DIR}/source/src_util/f90getopt.f90
  ${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_Sigma_User.f90
  ${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_var.f90
  )

#=========================================================================
#                      Util Library Source
#=========================================================================
add_library(ff_util ${src_util})
target_link_libraries(ff_util ff_ini_c ${FF_DEFAULT_LIBS})
target_include_directories(ff_util PUBLIC ${FF_APPLICATION_INCLUDE_PATH})
target_compile_options(ff_util PUBLIC ${Fortran_FLAGS})

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
${CMAKE_SOURCE_DIR}/source/assemblies/QuadSc_Interface.f90
)

#=========================================================================
#                      Assemblies Library Source
#=========================================================================
add_library(ff_assemblies ${src_assemblies})
target_link_libraries(ff_assemblies ff_util ff_elements ff_mesh ${FF_DEFAULT_LIBS})
target_include_directories(ff_assemblies PUBLIC ${FF_APPLICATION_INCLUDE_PATH})
target_compile_options(ff_assemblies PUBLIC ${Fortran_FLAGS})

set(src_mesh
  ${CMAKE_SOURCE_DIR}/source/src_mesh/mesh_refine.f90
  ${CMAKE_SOURCE_DIR}/source/src_mesh/umbrella_smoother.f90
  ${CMAKE_SOURCE_DIR}/source/src_mesh/geometry_processing.f90
  )

#=========================================================================
#                       Mesh Library Source
#=========================================================================
add_library(ff_mesh ${src_mesh})
target_link_libraries(ff_mesh ff_util ${FF_DEFAULT_LIBS})
target_include_directories(ff_mesh PUBLIC ${FF_APPLICATION_INCLUDE_PATH})
target_compile_options(ff_mesh PUBLIC ${Fortran_FLAGS})

set(Elements
  ${CMAKE_SOURCE_DIR}/source/Elements/e013.f
  ${CMAKE_SOURCE_DIR}/source/Elements/e012.f
  )

#=========================================================================
#                      Elements Library Source
#=========================================================================
add_library(ff_elements ${Elements})
target_include_directories(ff_elements PUBLIC ${FF_APPLICATION_INCLUDE_PATH})
target_compile_options(ff_elements PUBLIC ${Fortran_FLAGS})

set(le_solvers 
  ${CMAKE_SOURCE_DIR}/source/UmfpackSolver.f90
)

if(USE_MUMPS)
  list(APPEND le_solvers ${CMAKE_SOURCE_DIR}/source/MumpsSolver.f90)
endif(USE_MUMPS)

if(USE_HYPRE)
  list(APPEND le_solvers ${CMAKE_SOURCE_DIR}/source/HypreSolver.f90)
endif(USE_HYPRE)

add_library(ff_le_solvers ${le_solvers})
target_link_libraries(ff_le_solvers ff_util ff_mesh ${FF_DEFAULT_LIBS})
target_include_directories(ff_le_solvers PUBLIC ${FF_APPLICATION_INCLUDE_PATH})
target_compile_options(ff_le_solvers PUBLIC ${Fortran_FLAGS})

set(src_fbm
  ${CMAKE_SOURCE_DIR}/source/src_fbm/fbm_aux.f90
  ${CMAKE_SOURCE_DIR}/source/src_fbm/fbm_main.f90
)

#=========================================================================
#                         FBM Library Source
#=========================================================================
add_library(ff_fbm ${src_fbm})
target_link_libraries(ff_fbm ff_util ff_mesh ff_cinterface ${FF_DEFAULT_LIBS})
target_include_directories(ff_fbm PUBLIC ${FF_APPLICATION_INCLUDE_PATH})
target_compile_options(ff_fbm PUBLIC ${Fortran_FLAGS})

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

#=========================================================================
#                       LinSc Library Source
#=========================================================================
add_library(ff_LinSc ${src_LinSc})
target_link_libraries(ff_LinSc ${FF_DEFAULT_LIBS} ff_util ff_mesh ff_assemblies ff_elements ff_le_solvers)
target_include_directories(ff_LinSc PUBLIC ${FF_APPLICATION_INCLUDE_PATH})
target_compile_options(ff_LinSc PUBLIC ${Fortran_FLAGS})

set(src_particles
  ${CMAKE_SOURCE_DIR}/source/src_particles/part_input.f90
  ${CMAKE_SOURCE_DIR}/source/src_particles/part_step.f90
  ${CMAKE_SOURCE_DIR}/source/src_particles/part_tracer.f90
  )

#=========================================================================
#                      Particles Library Source
#=========================================================================
add_library(ff_particles ${src_particles})
target_link_libraries(ff_particles ff_util ff_mesh ff_fbm ff_le_solvers ff_postprocessing ff_quadLS_app ${FF_DEFAULT_LIBS})
target_include_directories(ff_particles PUBLIC ${FF_APPLICATION_INCLUDE_PATH})
target_compile_options(ff_particles PUBLIC ${Fortran_FLAGS})

# source files for standard q2p1 
set(src_q2p1
  ${CMAKE_SOURCE_DIR}/source/OutputProfiles.f90
  ${CMAKE_SOURCE_DIR}/source/PID.f90
  ${CMAKE_SOURCE_DIR}/source/tetra.f90
  ${CMAKE_SOURCE_DIR}/source/Statistics.f90
  ${CMAKE_SOURCE_DIR}/source/inverse.f90
  ${CMAKE_SOURCE_DIR}/source/distance.f90
  ${CMAKE_SOURCE_DIR}/source/3x3EigenV.f
  )

#=========================================================================
#                       Q2P1 Library Source
#=========================================================================
add_library(ff_q2p1 ${src_q2p1})
target_link_libraries(ff_q2p1 ff_util ff_mesh ff_fbm ff_postprocessing ${FF_DEFAULT_LIBS})
target_include_directories(ff_q2p1 PUBLIC ${FF_APPLICATION_INCLUDE_PATH})
target_compile_options(ff_q2p1 PUBLIC ${Fortran_FLAGS})


#*************************************************************************

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

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

set(src_quadLS_app_only
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_solver.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_proj.f
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_force.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_cylforce.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_torque.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_Sigma_User.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_mpi.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_main.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_user.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_mg.f90
${CMAKE_SOURCE_DIR}/source/src_quadLS/QuadSc_def.f90
${CMAKE_SOURCE_DIR}/source/initialization/app_initialization.f90
${CMAKE_SOURCE_DIR}/source/Init.f90
${CMAKE_SOURCE_DIR}/source/Umbrella.f90
${CMAKE_SOURCE_DIR}/source/ProcCtrl.f90
)

set(src_visco
  ${CMAKE_SOURCE_DIR}/source/src_visco/Visco_integration.f
  ${CMAKE_SOURCE_DIR}/source/src_visco/Visco_solver.f
  ${CMAKE_SOURCE_DIR}/source/src_visco/Visco_main.f90
  ${CMAKE_SOURCE_DIR}/source/src_visco/Visco_def.f90
  )

#=========================================================================
#                        QuadLS Library Source
#=========================================================================
add_library(ff_quadLS_app ${src_quadLS_app_only} ${src_visco})
target_link_libraries(ff_quadLS_app ff_util ff_mesh ff_assemblies ff_elements ff_le_solvers ff_fbm ff_LinSc ff_q2p1 ${FF_DEFAULT_LIBS})
target_include_directories(ff_quadLS_app PUBLIC ${FF_APPLICATION_INCLUDE_PATH})
target_compile_options(ff_quadLS_app PUBLIC ${Fortran_FLAGS})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


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

##=========================================================================
##                       PLinSc Library Source
##=========================================================================
#add_library(lib_PLin ${src_PLin})


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
