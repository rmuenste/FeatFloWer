PROGRAM Q2P1_PARTICLE_TRACER

  include 'defs_include.h'

  use solution_io, only: postprocessing_app

  use app_initialization, only: init_q2p1_particle_tracer

  use post_utils,  only: handle_statistics,&
                         print_time,&
                         sim_finalize

  use fbmaux                         
  use particle

  integer            :: iOGMV,iTout
  character(len=200) :: command
  character(len=60)  :: CPP3D
  real               :: dout = 0.0
  integer            :: ufile,ilog
  real               :: tt0 = 0.0
  real               :: dtt0 = 0.0
  real*8, dimension(8) :: x,y,z
  !real*8 :: xx,yy,zz,DPARX,DPARY,DPARZ
  integer :: ielem = 0

  call init_q2p1_particle_tracer(ufile)

  call Init_Particle(ufile)

  call Transport_Particle(ufile)

  call sim_finalize(tt0,ufile)


END PROGRAM Q2P1_PARTICLE_TRACER
