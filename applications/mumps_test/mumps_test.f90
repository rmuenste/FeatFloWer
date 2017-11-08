PROGRAM MUMPS_TEST

  include 'defs_include.h'
  use var_QuadScalar, only: istep_ns
  use solution_io, only: postprocessing_app
  use post_utils,  only: handle_statistics,&
    print_time,&
    sim_finalize

  integer            :: iOGMV,iTout
  character(len=200) :: command
  character(len=60)  :: CPP3D
  real               :: dout = 0.0
  integer            :: ufile=111, uterm,ilog
  real               :: tt0 = 0.0
  real               :: dtt0 = 0.0

  !-------INIT PHASE-------

  include 'mpif.h'
  include 'dmumps_struc.h'
  type (dmumps_struc) mumps_par
  integer ierr, i
  character cFile*(20)
  integer :: myFile = 55
  integer :: istat

  integer :: iloop

  call MPI_INIT(IERR)
  ! Define a communicator for the package.
  mumps_par%COMM = MPI_COMM_WORLD

  do iloop = 1,2

    !  Initialize an instance of the package
    !  for L U factorization (sym = 0, with working host)
    mumps_par%JOB = -1
    mumps_par%SYM = 0
    mumps_par%PAR = 0
    call DMUMPS(mumps_par)

    if (mumps_par%INFOG(1).LT.0) then
      write(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",&
        "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1),&
        "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
      goto 500
    end if

    write(cFile,'(A,I2.2,A)')  'input_',mumps_par%MYID,'.txt'

    open(file=adjustl(trim(cFile)),unit=myFile, &
      action="read",iostat=istat)

    !  Define problem on the host (processor 0)
    if ( mumps_par%myid .eq. 0 ) then
      read(myFile,*) mumps_par%N
      read(myFile,*) mumps_par%NZ
      allocate( mumps_par%IRN ( mumps_par%NZ ) )
      allocate( mumps_par%JCN ( mumps_par%NZ ) )

      !         ALLOCATE( mumps_par%A   ( mumps_par%NZ ) )
      !mumps_par%A  = 0d0

      do I = 1, mumps_par%NZ
      read(myFile,*) mumps_par%IRN(I),mumps_par%JCN(I)
      end do
      allocate( mumps_par%RHS ( mumps_par%N  ) )
      do I = 1, mumps_par%N
      read(myFile,*) mumps_par%RHS(I)
      !          WRITE(*,*) mumps_par%rhs(I)
      end do
    else

      read(myFile,*) mumps_par%NZ_loc
      allocate( mumps_par%IRN_loc ( mumps_par%NZ_loc) )
      allocate( mumps_par%JCN_loc ( mumps_par%NZ_loc) )
      allocate( mumps_par%A_loc( mumps_par%NZ_loc) )
      do I = 1, mumps_par%NZ_loc
      read(myFile,*) mumps_par%IRN_loc(I),mumps_par%JCN_loc(I),&
        mumps_par%A_loc(I)
      end do

    end if

    close(myFile)

    !  Call package for solution
    MUMPS_PAR%icntl(6)  = 7
    !     pivot order (automatic)
    MUMPS_PAR%icntl(7)  = 7
    !     scaling (automatic)
    MUMPS_PAR%icntl(8)  = 7
    !     no transpose
    MUMPS_PAR%icntl(9)  = 1
    !     max steps for iterative refinement
    MUMPS_PAR%icntl(10) = 0
    !     statistics info
    MUMPS_PAR%icntl(11) = 0 
    !     controls parallelism
    MUMPS_PAR%icntl(12) = 0
    !     use ScaLAPACK for root node
    MUMPS_PAR%icntl(13) = 0
    !     percentage increase in estimated workspace
    MUMPS_PAR%icntl(14) = 100

    mumps_par%ICNTL(5)=0
    mumps_par%ICNTL(18)=3
    !     mumps_par%WRITE_PROBLEM='matrix'

    mumps_par%JOB = 6
    CALL DMUMPS(mumps_par)

    !       mumps_par%JOB = 2
    !       CALL DMUMPS(mumps_par)
    ! 
    !       mumps_par%JOB = 3
    !       CALL DMUMPS(mumps_par)

    if (mumps_par%INFOG(1).lt.0) then
      write(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",&
        "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1),& 
        "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
      !    goto 500
    end if

    !  !  Solution has been assembled on the host
    !  if ( mumps_par%MYID .eq. 0 ) then
    !    write( 6, * ) ' Solution is '
    !    do i=1,mumps_par%N
    !    write( 6, * ) mumps_par%RHS(I)
    !    end do
    !  end if

    !  Deallocate user data
!    if ( mumps_par%MYID .eq. 0 )then
      !         DEALLOCATE( mumps_par%A )
!    end if

    if(allocated( mumps_par%RHS)   ) deallocate( mumps_par%RHS)
    if(allocated( mumps_par%IRN)   ) deallocate( mumps_par%IRN)
    if(allocated( mumps_par%JCN)   ) deallocate( mumps_par%JCN)
                                   
    if(allocated( mumps_par%IRN_loc)) deallocate( mumps_par%IRN_loc)
    if(allocated( mumps_par%JCN_loc)) deallocate( mumps_par%JCN_loc)
    if(allocated( mumps_par%A_loc) ) deallocate( mumps_par%A_loc)

    !  Destroy the instance (deallocate internal data structures)
    mumps_par%JOB = -2
    call DMUMPS(mumps_par)
    if (mumps_par%INFOG(1).LT.0) then
      write(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",&
        "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1),& 
        "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
      goto 500
    end if

  end do

  500  call mpi_finalize(ierr)

  stop

  !       write(*,*) mumps_par%MYID,"... is done!"
  !       pause

END PROGRAM MUMPS_TEST
