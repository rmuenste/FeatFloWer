module post_utils
  USE var_QuadScalar,ONLY:knvt,knet,knat,knel
  !-------------------------------------------------------------------------------------------------
  ! A module for saving the solution values to 
  ! a file. This output(dump) is mainly done
  ! to a binary file.
  !-------------------------------------------------------------------------------------------------

  ! a variable for counting the outputs
  integer :: ifile = 0

contains
!
!-------------------------------------------------------------------------------------------------
! A simple time stepping routine
!-------------------------------------------------------------------------------------------------
! @param dt The current time step 
! @param inlU   
! @param inlT 
! @param filehandle Unit of the output file
!
subroutine TimeStepCtrl(dt,inlU,inlT, filehandle)

  USE PP3D_MPI,only :myid,ShowID

  INTEGER IADTIM

  REAL*8  TIMEMX,DTMIN,DTMAX

  COMMON /NSADAT/ TIMEMX,DTMIN,DTMAX,IADTIM

  integer, intent(in) :: filehandle

  integer :: inlU,inlT
  integer :: iMem,nMEm=2
  real*8  :: dt, dt_old
  character(len=9) :: char_dt
  data iMem/0/

  IF (IADTIM.EQ.0) RETURN

  iMem = iMem + 1
  dt_old = dt
  IF (((inlU.GT.3).OR. (inlT.GT.5)).AND.iMem.GE.nMem) THEN
    dt=MAX(dt/1.1d0,DTMIN)
    WRITE(char_dt,'(D9.2)') dt
    READ (char_dt,'(D9.2)') dt
  END IF
  IF (((inlU.LT.3).AND.(inlT.LT.4)).AND.iMem.GE.nMem) THEN
    dt=MIN(1.1d0*dt,DTMAX)
    WRITE(char_dt,'(D9.2)') dt
    READ (char_dt,'(D9.2)') dt
  END IF

  IF (dt.NE.dt_old.AND.myid.eq.ShowID) THEN
    WRITE(MTERM,1) dt_old,dt
    WRITE(filehandle,1) dt_old,dt
  END IF

  IF (dt.NE.dt_old) iMem = 0

  1  FORMAT('Time step change from ',D9.2,' to ',D9.2)

END SUBROUTINE TimeStepCtrl
!
!-------------------------------------------------------------------------------------------------
! A routine for an application to start the statistics output
!-------------------------------------------------------------------------------------------------
! @param dttt0 
! @param istepns   
!
subroutine handle_statistics(dttt0, istepns)
  include 'defs_include.h'

  implicit none

  real, intent(in) :: dttt0
  integer, intent(in) :: istepns
  real :: dttx = 0.0

  ! Statistics reset
  IF (istepns.eq.1) THEN
    CALL ResetTimer()
    CALL ZTIME(dttt0)
  END IF

  IF (MOD(istepns,10).EQ.5) THEN
    CALL ZTIME(dttx)
    CALL StatOut(dttx-dttt0,0)
  END IF

end subroutine handle_statistics
!
!-------------------------------------------------------------------------------------------------
! An output routine for printing time to the screen 
!-------------------------------------------------------------------------------------------------
! @param dtimens Simulation time 
! @param dtimemx Maximum simulation time 
! @param dt Time step of the simulation
! @param istepns Discrete time step of the simulation   
! @param istepmaxns Maximum discrete time step   
! @param ufile Handle to the log file  
! @param uterm Handle to the terminal  
!
subroutine print_time(dtimens, dtimemx, dt, istepns, istepmaxns, ufile,uterm)

  include 'defs_include.h'
  USE PP3D_MPI,only :myid,ShowID

  implicit none

  real*8, intent(in)  :: dtimens, dtimemx, dt
  integer, intent(in) :: istepns , istepmaxns
  integer, intent(in) :: ufile, uterm

  IF (myid.eq.showid) THEN
    write(uTERM,5)
    write(uFILE,5)
    write(uTERM,'(A5,D12.4,A1,D12.4,A8,I8,A1,I8,A7,D12.4)')&
      "time:", dtimens,"/",dtimemx," | itns:",istepns,"/",istepmaxns,&
      " | dt:",dt
    write(uFILE,'(A5,D12.4,A1,D12.4,A8,I8,A1,I8,A7,D12.4)')&
      "time:", dtimens,"/",dtimemx," | itns:",istepns,"/",istepmaxns,&
      " | dt:",dt
    write(uTERM,5)
    write(uFILE,5)
  END IF

  5 FORMAT(104('='))

end subroutine print_time
!
!-------------------------------------------------------------------------------------------------
! A wrapper for mpi finalize 
!-------------------------------------------------------------------------------------------------
! @param dttt0 Simulation start time 
! @param filehandle protocol file handle 
! ----------------------------------------------
subroutine sim_finalize(dttt0, filehandle)
  use Mesh_Structures, only: release_mesh
  USE PP3D_MPI, ONLY : myid,master,showid,Barrier_myMPI
  use var_QuadScalar, only: istep_ns,mg_mesh

  real, intent(inout) :: dttt0
  integer, intent(in) :: filehandle

  integer :: ierr
  integer :: terminal = 6
  real :: time,time_passed

  CALL ZTIME(time)

  time_passed = time - dttt0
  CALL StatOut(time_passed,0)

  CALL StatOut(time_passed,terminal)

  ! Save the final solution vector in unformatted form
  istep_ns = istep_ns - 1
  CALL SolToFile(-1)

  IF (myid.eq.showid) THEN
    WRITE(*,*) "PP3D_LES has successfully finished. "
    WRITE(filehandle,*) "PP3D_LES has successfully finished. "
  END IF

  call release_mesh(mg_mesh)
  CALL Barrier_myMPI()
  CALL MPI_Finalize(ierr)

end subroutine sim_finalize

end module post_utils
