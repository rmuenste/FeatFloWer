PROGRAM MUMPS_TEST

!   include 'defs_include.h'

  integer            :: iOGMV,iTout
  character(len=200) :: command
  character(len=60)  :: CPP3D
  real               :: dout = 0.0
  integer            :: ufile=111, ilog
  real               :: tt0 = 0.0
  real               :: dtt0 = 0.0

  !-------INIT PHASE-------

  include 'mpif.h'
  include 'dmumps_struc.h'
  type (dmumps_struc) mumps_par
  integer ierr, i
  character cFile*(100)
  integer :: myFile = 55
  integer :: istat,numnodes,subnodes
  integer*4, allocatable, dimension(:), target :: myJCN
  integer*4, allocatable, dimension(:), target :: myIRN

  integer :: iloop

  call MPI_INIT(IERR)
  ! Define a communicator for the package.
  mumps_par%COMM = MPI_COMM_WORLD
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,numnodes,IERR)
  subnodes=numnodes-1

  do iloop = 1,1

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

    write(cFile,'(A,I2.2,A,I2.2,A)')  'input_',subnodes,'/input_',mumps_par%MYID,'.txt'
!     write(*,*) adjustl(trim(cFile)),subnodes
!     pause
    open(file=adjustl(trim(cFile)),unit=myFile, &
      action="read",iostat=istat)

    !  Define problem on the host (processor 0)
    if ( mumps_par%myid .eq. 0 ) then
      read(myFile,*) mumps_par%N
      read(myFile,*) mumps_par%NZ

      IF (associated(mumps_par%IRN)) THEN
        deallocate(mumps_par%IRN)
        write(*,*)'Deallocated-------------------------- IRN:',iloop
        nullify(mumps_par%IRN)
      endif

      if(.not.allocated(myJCN))then
        allocate(myJCN(mumps_par%NZ))
      end if
        mumps_par%JCN => myJCN

      if(.not.allocated(myIRN))then
        allocate(myIRN(mumps_par%NZ))
      end if
        mumps_par%IRN => myIRN

      IF (associated(mumps_par%A)) THEN
        deallocate(mumps_par%A)
        write(*,*)'Deallocated-------------------------- A:',iloop
        nullify(mumps_par%A)
      endif

      !         ALLOCATE( mumps_par%A   ( mumps_par%NZ ) )
      !mumps_par%A  = 0d0

      do I = 1, mumps_par%NZ
      read(myFile,*) mumps_par%IRN(I),mumps_par%JCN(I)
      end do

      IF (associated(mumps_par%RHS)) THEN
        deallocate(mumps_par%RHS)
        nullify(mumps_par%RHS)
      endif

      allocate( mumps_par%RHS ( mumps_par%N  ) )

      do I = 1, mumps_par%N
      read(myFile,*) mumps_par%RHS(I)
      end do
    else

      read(myFile,*) mumps_par%NZ_loc

      IF (associated(mumps_par%IRN_loc)) THEN
        deallocate(mumps_par%IRN_loc)
        nullify(mumps_par%IRN_loc)
      endif

      allocate( mumps_par%IRN_loc ( mumps_par%NZ_loc) )

      IF (associated(mumps_par%JCN_loc)) THEN
        deallocate(mumps_par%JCN_loc)
        nullify(mumps_par%JCN_loc)
      endif

      allocate( mumps_par%JCN_loc ( mumps_par%NZ_loc) )

      if (associated(mumps_par%A_loc)) THEN
        deallocate(mumps_par%A_loc)
        nullify(mumps_par%A_loc)
      endif

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
    MUMPS_PAR%icntl(13) = 2
    !     percentage increase in estimated workspace
    MUMPS_PAR%icntl(14) = 2

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

    !     !  Solution has been assembled on the host
    !     if ( mumps_par%MYID .eq. 0 ) then
    !       write( 6, * ) ' Solution is '
    !       do i=1,mumps_par%N
    !       write( 6, * ) mumps_par%RHS(I)
    !       end do
    !     end if

    !  Deallocate user data
    !    if ( mumps_par%MYID .eq. 0 )then
    !         DEALLOCATE( mumps_par%A )
    !    end if

    IF (associated(mumps_par%IRN)) THEN
      deallocate(mumps_par%IRN)
      nullify(mumps_par%IRN)
    endif

    IF (associated(mumps_par%JCN)) THEN
      deallocate(mumps_par%JCN)
      write(*,*)'Deallocated-------------------------- End End JCN:',iloop
      nullify(mumps_par%JCN)
    endif

    IF (associated(mumps_par%RHS)) THEN
      write(cFile,'(A,I2.2,A,I2.2,A)')  'input_',subnodes,'/solution_',mumps_par%MYID,'.txt'
      open(file=adjustl(trim(cFile)),unit=myFile, action="write",iostat=istat)
      write(myFile,*) mumps_par%RHS
      close(myFile)
      deallocate(mumps_par%RHS)
      nullify(mumps_par%RHS)
    endif

    IF (associated(mumps_par%IRN_loc)) THEN
      deallocate(mumps_par%IRN_loc)
      nullify(mumps_par%IRN_loc)
    endif

    IF (associated(mumps_par%JCN_loc)) THEN
      deallocate(mumps_par%JCN_loc)
      nullify(mumps_par%JCN_loc)
    endif

    IF (associated(mumps_par%A_loc)) THEN
      deallocate(mumps_par%A_loc)
      nullify(mumps_par%A_loc)
    endif

    !  Destroy the instance (deallocate internal data structures)
    mumps_par%JOB = -2
    call DMUMPS(mumps_par)
    if (mumps_par%INFOG(1).LT.0) then
      write(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",&
        "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1),&
        "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2)
      goto 500
    end if

    IF (associated(mumps_par%future_niv2)) THEN
      deallocate(mumps_par%future_niv2)
      nullify(mumps_par%future_niv2)
    endif

    !------------------------------------------------

      IF (associated(mumps_par%IRN)) THEN
        deallocate(mumps_par%IRN)
        write(*,*)'Deallocated-------------------------- IRN:',iloop
        nullify(mumps_par%IRN)
      endif

      IF (associated(mumps_par%JCN)) THEN
        deallocate(mumps_par%JCN)
        write(*,*)'Deallocated-------------------------- JCN:',iloop
        nullify(mumps_par%JCN)
      endif

      IF (associated(mumps_par%A)) THEN
        deallocate(mumps_par%A)
        write(*,*)'Deallocated-------------------------- A:',iloop
        nullify(mumps_par%A)
      endif

      IF (associated(mumps_par%RHS)) THEN
        deallocate(mumps_par%RHS)
        nullify(mumps_par%RHS)
      endif
    !----------------------------------------------------------------------k

  end do


    500  call mpi_finalize(ierr)

    stop

    !       write(*,*) mumps_par%MYID,"... is done!"
    !       pause

  END PROGRAM MUMPS_TEST
