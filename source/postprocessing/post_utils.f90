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
! A routine for an application to start the statistics output
!-------------------------------------------------------------------------------------------------
! @param dttt0 
! @param istepns   
!
subroutine handle_statistics(dttt0, istepns)

  USE var_QuadScalar, ONLY: myStat
  use def_Feat

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
! @param term_out Handle to the terminal  
!
subroutine print_time_q1(dtimens, dtimemx, dt,dttMIN,dttMAX, istepns, istepmaxns, ufile, term_out)

  use def_Feat
  USE PP3D_MPI,only :myid,ShowID

  implicit none

  real*8, intent(in)  :: dtimens, dtimemx, dt,dttMIN,dttMAX
  integer, intent(in) :: istepns , istepmaxns
  integer, intent(in) :: ufile, term_out

  IF (myid.eq.showid) THEN
    write(term_out,5)
    write(uFILE,5)
    write(term_out,'(A5,D12.4,A1,D12.4,A8,I8,A1,I8,A7,3D12.4)')&
      "time:", dtimens,"/",dtimemx," | itns:",istepns,"/",istepmaxns,&
      " | dt:",dt,dttMIN,dttMAX
    write(uFILE,'(A5,D12.4,A1,D12.4,A8,I8,A1,I8,A7,3D12.4)')&
      "time:", dtimens,"/",dtimemx," | itns:",istepns,"/",istepmaxns,&
      " | dt:",dt,dttMIN,dttMAX
    write(term_out,5)
    write(uFILE,5)
  END IF

  5 FORMAT(104('='))

end subroutine print_time_q1

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
! @param term_out Handle to the terminal  
!
subroutine print_time(dtimens, dtimemx, dt, istepns, istepmaxns, ufile, term_out)

  use def_Feat
  USE PP3D_MPI,only :myid,ShowID

  implicit none

  real*8, intent(in)  :: dtimens, dtimemx, dt
  integer, intent(in) :: istepns , istepmaxns
  integer, intent(in) :: ufile, term_out

  IF (myid.eq.showid) THEN
    write(term_out,5)
    write(uFILE,5)
    write(term_out,'(A5,D12.4,A1,D12.4,A8,I8,A1,I8,A7,D12.4)')&
      "time:", dtimens,"/",dtimemx," | itns:",istepns,"/",istepmaxns,&
      " | dt:",dt
    write(uFILE,'(A5,D12.4,A1,D12.4,A8,I8,A1,I8,A7,D12.4)')&
      "time:", dtimens,"/",dtimemx," | itns:",istepns,"/",istepmaxns,&
      " | dt:",dt
    write(term_out,5)
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
  use def_Feat
  use Mesh_Structures, only: release_mesh
  USE PP3D_MPI, ONLY : myid,master,showid,Barrier_myMPI
  use var_QuadScalar, only: istep_ns,mg_mesh
  use solution_io, only: write_sol_to_file

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
  !CALL SolToFile(-1)
  call write_sol_to_file(insavn, timens)

  IF (myid.eq.showid) THEN
    WRITE(*,*) "PP3D_LES has successfully finished. "
    WRITE(filehandle,*) "PP3D_LES has successfully finished. "
  END IF

  call release_mesh(mg_mesh)
  CALL Barrier_myMPI()
  CALL MPI_Finalize(ierr)

end subroutine sim_finalize
!
!-------------------------------------------------------------------------------------------------
! A wrapper for mpi finalize 
!-------------------------------------------------------------------------------------------------
! @param dttt0 Simulation start time 
! @param filehandle protocol file handle 
! ----------------------------------------------
subroutine sim_finalize_sse(dttt0, filehandle)
  use def_Feat
  USE PP3D_MPI, ONLY : myid,master,showid,Barrier_myMPI,MPI_COMM_WORLD
  use var_QuadScalar, only: istep_ns,mg_mesh,DivergedSolution,myErrorCode
  use Mesh_Structures, only: release_mesh
  use solution_io, only: write_sol_to_file,write_sse_1d_sol
  use Sigma_User, only: myProcess

  real, intent(inout) :: dttt0
  integer, intent(in) :: filehandle

  integer :: ierr
  integer :: terminal = 6
  real :: time,time_passed

  CALL ZTIME(time)

  time_passed = time - dttt0
  CALL StatOut(time_passed,0)

  CALL StatOut(time_passed,terminal)
  
  if (.not.DivergedSolution) then

  ! Save the final solution vector in unformatted form
  istep_ns = istep_ns - 1
  !CALL SolToFile(-1)
!  CALL Release_ListFiles_SSE(int(myProcess%Angle))
!  call write_sol_to_file(insavn, timens)
!  call write_sse_1d_sol()
  
   IF (myid.eq.showid) THEN
     WRITE(*,*) "Q2P1_SSE has successfully finished. "
     WRITE(filehandle,*) "Q2P1_SSE has successfully finished. "
   END IF
  ELSE
   IF (myid.eq.showid) THEN
     WRITE(*,*) "Q2P1_SSE has been stopped due to divergence ..."
     WRITE(filehandle,*) "Q2P1_SSE has been stopped due to divergence ..."
   END IF
  end if


  call release_mesh(mg_mesh)
  CALL Barrier_myMPI()
  CALL MPI_Finalize(ierr)

end subroutine sim_finalize_sse
!
!-------------------------------------------------------------------------------------------------
! A wrapper for mpi finalize 
!-------------------------------------------------------------------------------------------------
! @param dttt0 Simulation start time 
! @param filehandle protocol file handle 
! ----------------------------------------------
subroutine sim_finalize_sse_mesh(dttt0, filehandle)
  use def_Feat
  USE PP3D_MPI, ONLY : myid,master,showid,MPI_COMM_WORLD
  use var_QuadScalar, only: istep_ns,mg_mesh,DivergedSolution,myErrorCode
  use Mesh_Structures, only: release_mesh
  use solution_io, only: write_sol_to_file,write_sse_1d_sol
  use Sigma_User, only: myProcess

  real, intent(inout) :: dttt0
  integer, intent(in) :: filehandle

  integer :: ierr
  
  IF (myid.eq.showid) THEN
    WRITE(*,*) "Q2P1_SSE_MESH has successfully finished. "
    WRITE(filehandle,*) "Q2P1_SSE_MESH has successfully finished. "
  END IF

  call release_mesh(mg_mesh)
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
  CALL MPI_Finalize(ierr)

end subroutine sim_finalize_sse_mesh

end module post_utils
