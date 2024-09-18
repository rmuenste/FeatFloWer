PROGRAM Q2P1_FC_EXT

  include 'defs_include.h'

  use solution_io, only: postprocessing_app

  use app_initialization, only: init_q2p1_app

  use post_utils,  only: handle_statistics,&
                         print_time,&
                         sim_finalize
  USE var_QuadScalar, ONLY :  myTimer

  integer            :: iOGMV,iTout
  character(len=200) :: command
  character(len=60)  :: CPP3D
  real               :: dout = 0.0
  integer            :: ufile,ilog
  real               :: tt0 = 0.0
  real               :: dtt0 = 0.0
  real               :: dtt10 = 0.0

  !-------INIT PHASE-------

  call init_q2p1_app(ufile)

  CALL ZTIME(tt0)
  call ztime(dtt0)

  dout = Real(INT(timens/dtgmv)+1)*dtgmv

  !-------MAIN LOOP-------

  DO itns=1,nitns

  itnsr=0
  timnsh=timens
  dt=tstep
  timens=timens+dt

  ! Solve Navier-Stokes (add discretization in name + equation or quantity)
  CALL Transport_q2p1_UxyzP_fc_ext(ufile,inonln_u,itns)

  IF (bTracer) THEN
    ! Solve transport equation for linear scalar
    CALL Transport_LinScalar(ufile,inonln_t)
  ELSE
    inonln_t = 2
  END IF

  call postprocessing_app(dout, inonln_u, inonln_t,ufile)

  call print_time(timens, timemx, tstep, itns, nitns, ufile, uterm)

  call handle_statistics(tt0,itns)

  istep_ns = istep_ns + 1
  ! Exit if done
  IF (timemx.LE.(timens+1D-10)) EXIT

  END DO

!   write(*,'(I0,A,6I5,6ES12.4)') myid," : ",myTimer%n,myTimer%t
  CALL Output_MPI_Timings()
  
  call sim_finalize(tt0,ufile)


END PROGRAM Q2P1_FC_EXT
!
!-----------------------------------------------------------------------------------
!
SUBROUTINE Output_MPI_Timings()
use mpi
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier,MPI_COMM_SUBS,NUMNODES,subNODES
USE var_QuadScalar, ONLY : myStat,myCommTimer,myTimer
implicit none
integer :: n1,n2,i,j,ierr
REAL*8,allocatable :: d1(:),d2(:)
character :: cFMT*(256)


if (myid.eq.master) return

n1 = 14
allocate(d1(n1))

if (myid.eq.1) then
 n2 = n1*subnodes
 allocate(d2(n2))
 d2=0d0
end if

!d1 = 1d0
d1(1:6) = DBLE(myTimer(1:6)%t)
d1(7:14) = DBLE(myCommTimer(1:8)%t)

CALL MPI_GATHER(d1, n1, MPI_DOUBLE_PRECISION, d2, n1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SUBS, ierr)

if (myid.eq.1) then
 open(file='_data/CommStats.txt',UNIT=5465)
 cFMT = ' '
 write(cFMT(1:),'(A,I0,A,A)') "(",subnodes-1,"((ES12.4,(','))",",ES12.4))"
 write(*,*) adjustl(trim(cFMT))
 DO i=1,n1
  write(5465,cFMT) (d2(i + (j-1)*n1),j=1,subnodes)
 end do
 close(5465)
!  write(*,*) d2
END IF

deallocate(d1)
if (myid.eq.1) then
 deallocate(d2)
end if

END SUBROUTINE Output_MPI_Timings