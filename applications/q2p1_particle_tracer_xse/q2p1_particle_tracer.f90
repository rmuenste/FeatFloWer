PROGRAM Q2P1_PARTICLE_TRACER

  include 'defs_include.h'

  use app_initialization, only: init_q2p1_particle_tracer
  use particle

  integer            :: ufile,ilog
  real               :: tt0 = 0.0

  call init_q2p1_particle_tracer(ufile)

  call Init_Particle_xse(ufile)

  call Transport_Particle_xse(ufile)
  
  call Finalize_Particles(tt0,ufile)

END PROGRAM Q2P1_PARTICLE_TRACER
