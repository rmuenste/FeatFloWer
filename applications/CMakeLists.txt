
# Add the subdirectories of the applications
add_subdirectory(q2p1_fc_ext)

add_subdirectory(q2p1_fac_nnewt)

add_subdirectory(q2p1_fac_visco)

add_subdirectory(q2p1_fbm_spring)

add_subdirectory(q2p1_bench_sedimentation)

add_subdirectory(laplace)

add_subdirectory(q2p1_particles)

add_subdirectory(q2p1_experimental)

add_subdirectory(q2p1_euler_fibers)

add_subdirectory(q2p1_particle_tracer)

if(USE_CGAL)
  add_subdirectory(q2p1_geom_test)
  add_subdirectory(q2p1_bloodflow)
  add_subdirectory(q2p1_sse)
endif(USE_CGAL)

if(USE_OPTICALTWEEZERS)
  add_subdirectory(q2p1_optic_tweezers)
endif(USE_OPTICALTWEEZERS)

if(USE_MUMPS)
  if(${Q2P1_BUILD_ID} MATCHES "intel")
    add_subdirectory(mumps_test)
    add_subdirectory(q2p1_cc)
  else(${Q2P1_BUILD_ID} MATCHES "intel")
    message(FATAL_ERROR "The MUMPS library requires the intel compiler. Choose a build-id with the intel compiler to resolve this issue.")
  endif(${Q2P1_BUILD_ID} MATCHES "intel")
endif(USE_MUMPS)

