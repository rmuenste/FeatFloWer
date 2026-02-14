MODULE ProcCtrl_mod
!===============================================================================
! MODULE: ProcCtrl_mod
!
! DESCRIPTION:
!   Runtime process control and finalization routines for the FeatFloWer solver.
!   Provides mechanisms for:
!   - Runtime parameter reloading via ProcCtrl.txt file
!   - Solution dumps and visualization output
!   - Proper MPI finalization
!
! PUBLIC SUBROUTINES:
!   - ProcessControl: Main runtime control loop for parameter reloading
!   - Finalize: Standard finalization routine
!   - Finalize_Particles: Finalization routine for particle simulations
!
! PRIVATE SUBROUTINES:
!   - ReloadParameterFile: Helper routine for parameter file reloading
!===============================================================================
 IMPLICIT NONE
 PRIVATE

 ! Public interface
 PUBLIC :: ProcessControl, Finalize, Finalize_Particles

 ! Named constants for file handling
 INTEGER, PARAMETER :: PROCCTRL_UNIT = 989      ! File unit for ProcCtrl.txt
 INTEGER, PARAMETER :: DUMP_FILE_MIN = 0        ! Min file number for Dump_Out
 INTEGER, PARAMETER :: DUMP_FILE_MAX = 99       ! Max file number for Dump_Out
 INTEGER, PARAMETER :: GMV_FILE_MIN  = 0        ! Min file number for GMV_Out
 INTEGER, PARAMETER :: GMV_FILE_MAX  = 9999     ! Max file number for GMV_Out

CONTAINS

!===============================================================================
SUBROUTINE ProcessControl(MFILE,MTERM)
!-------------------------------------------------------------------------------
! DESCRIPTION:
!   Monitors ProcCtrl.txt for runtime commands and executes them.
!   Supports: Dump_Out, Dump_In, GMV_Out, Reload_Velo, Reload_Pres,
!             Reload_SimPar, Reload_Prop
!-------------------------------------------------------------------------------
 USE PP3D_MPI, ONLY : myid,master,showid,ShareValueK_myMPI,&
                      Barrier_myMPI,ShareValueC_myMPI
 USE Transport_Q2P1, ONLY : myDataFile,QuadSc,LinSc,Properties
 USE param_parser, ONLY : GDATNEW,GetVeloParameters,GetPresParameters,GetPhysiclaParameters
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: MFILE, MTERM
 CHARACTER(len=200) :: string
 CHARACTER(len=100) :: command, cvalue
 CHARACTER(len=7)   :: CSimPar
 INTEGER :: iEnd(1), iCmnd(1), iExist(1), iPos, iFile, istat
 LOGICAL :: bExist

 ! Check if control file exists
 iExist = 0
 IF (myid == showid) THEN
  INQUIRE(FILE='_data/ProcCtrl.txt',EXIST=bExist)
  IF (bExist) iExist = 1
 END IF

 CALL ShareValueK_myMPI(iExist,1,showid)

 IF (iExist(1) == 0) RETURN

 ! Open the control file with error handling to prevent MPI deadlock
 IF (myid == showid) THEN
   OPEN(UNIT=PROCCTRL_UNIT,FILE='_data/ProcCtrl.txt',STATUS='OLD',ACTION='READWRITE',IOSTAT=istat)
   IF (istat /= 0) THEN
     WRITE(MTERM,'(A)') 'ERROR: Cannot open _data/ProcCtrl.txt'
     WRITE(MTERM,'(A,I0)') '  IOSTAT = ', istat
     iExist(1) = 0  ! Signal open failure
   END IF
 END IF

 ! Re-broadcast to handle open failure across all processes
 CALL ShareValueK_myMPI(iExist,1,showid)
 IF (iExist(1) == 0) RETURN

 ! Main control loop: read and process commands
 DO

   iCmnd = 0

   IF (myid == showid) THEN
    READ (UNIT=PROCCTRL_UNIT,FMT='(A)',IOSTAT=iEnd(1)) string
   END IF

   CALL ShareValueK_myMPI(iEnd,1,showid)
   IF (iEnd(1) == -1) EXIT

   CALL ShareValueC_myMPI(string,200,showid)

   ! Parse command and value from the control string
   iPos = INDEX(string,"=")
   IF (iPos /= 0) THEN
    READ(string(1:iPos-1),'(A)') command
    READ(string(iPos+1:),'(A)') cvalue
    cvalue = TRIM(ADJUSTL(cvalue))
   ELSE
    command = "Nothing to be done ... "
    cvalue = ""  ! Initialize to prevent undefined behavior
    IF (myid == showid) THEN
      WRITE(MTERM,'(A,A)') 'WARNING: Malformed line in ProcCtrl.txt: ', TRIM(string)
    END IF
   END IF

!================================================================================
   SELECT CASE (TRIM(ADJUSTL(command)))
! -------------------------------------------------------------------------------
   CASE ("Dump_Out")

    ! Read and validate file number with error handling
    READ(cvalue,'(I10)',IOSTAT=istat) iFile
    IF (istat /= 0) THEN
     IF (myid == showid) THEN
       WRITE(MTERM,'(A,A)') 'ERROR: Invalid file number for Dump_Out: ', TRIM(cvalue)
     END IF
    ELSE IF (iFile < DUMP_FILE_MIN .or. iFile > DUMP_FILE_MAX) THEN
     IF (myid == showid) THEN
       WRITE(MTERM,'(A,I0,A,I0,A,I0,A)') 'ERROR: File number out of range: ', iFile, &
                                          ' (valid: ', DUMP_FILE_MIN, '-', DUMP_FILE_MAX, ')'
     END IF
    ELSE
     CALL SolToFile(iFile)
     CALL FBM_ToFile()
    END IF

! -------------------------------------------------------------------------------
   CASE ("Dump_In")

    IF (myid == showid) WRITE(*,*) cvalue
    CALL SolFromFile(cvalue,1)
    CALL FBM_FromFile()

! -------------------------------------------------------------------------------
   CASE ("GMV_Out")

    ! Read and validate file number with error handling
    READ(cvalue,'(I10)',IOSTAT=istat) iFile
    IF (istat /= 0) THEN
     IF (myid == showid) THEN
       WRITE(MTERM,'(A,A)') 'ERROR: Invalid file number for GMV_Out: ', TRIM(cvalue)
     END IF
    ELSE IF (iFile < GMV_FILE_MIN .or. iFile > GMV_FILE_MAX) THEN
     IF (myid == showid) THEN
       WRITE(MTERM,'(A,I0,A,I0,A,I0,A)') 'ERROR: File number out of range: ', iFile, &
                                          ' (valid: ', GMV_FILE_MIN, '-', GMV_FILE_MAX, ')'
     END IF
    ELSE
     CALL Output_Profiles(iFile)
    END IF

! -------------------------------------------------------------------------------
   CASE ("Reload_Velo")
    CALL ReloadParameterFile(MFILE, MTERM, cvalue, "Velocity", bExist)
    IF (bExist) CALL GetVeloParameters(QuadSc%prm, QuadSc%cName, MFILE)

! -------------------------------------------------------------------------------
   CASE ("Reload_Pres")
    CALL ReloadParameterFile(MFILE, MTERM, cvalue, "Pressure", bExist)
    IF (bExist) CALL GetPresParameters(LinSc%prm, LinSc%cName, MFILE)

! -------------------------------------------------------------------------------
   CASE ("Reload_SimPar")
    CALL ReloadParameterFile(MFILE, MTERM, cvalue, "Simulation", bExist)
    IF (bExist) THEN
      CSimPar = "SimPar"
      CALL GDATNEW(CSimPar, 1)
    END IF

! -------------------------------------------------------------------------------
   CASE ("Reload_Prop")
    CALL ReloadParameterFile(MFILE, MTERM, cvalue, "Physical", bExist)
    IF (bExist) CALL GetPhysiclaParameters(Properties, Properties%cName, MFILE)

   END SELECT
!================================================================================

   CALL Barrier_myMPI()

  END DO

  ! Reset control file to default state
  IF (myid == showid) REWIND(PROCCTRL_UNIT)
  IF (myid == showid) WRITE(PROCCTRL_UNIT,*) "Nothing to be done ... "
  IF (myid == showid) CLOSE(PROCCTRL_UNIT)

END SUBROUTINE ProcessControl


!===============================================================================
SUBROUTINE ReloadParameterFile(MFILE, MTERM, filename, param_type, bExist)
!-------------------------------------------------------------------------------
! DESCRIPTION:
!   Helper routine for parameter file reloading operations.
!   Handles common tasks: logging, file existence checking, and setting myDataFile.
!   Reduces code duplication for Reload_Velo, Reload_Pres, Reload_SimPar,
!   and Reload_Prop commands.
!
! ARGUMENTS:
!   MFILE      - File unit for output logging
!   MTERM      - Terminal unit for output logging
!   filename   - Name of parameter file to reload
!   param_type - Descriptive name of parameter type (e.g., "Velocity")
!   bExist     - Returns .TRUE. if file exists and reload should proceed
!-------------------------------------------------------------------------------
 USE PP3D_MPI, ONLY : myid, showid
 USE Transport_Q2P1, ONLY : myDataFile
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: MFILE, MTERM
 CHARACTER(len=*), INTENT(IN) :: filename, param_type
 LOGICAL, INTENT(OUT) :: bExist

 ! Log reload action to terminal and file
 IF (myid == showid) THEN
   WRITE(MTERM,'(A,A,A,A)') "Reloading ", TRIM(param_type), " parameters from file: ", TRIM(ADJUSTL(filename))
   WRITE(MFILE,'(A,A,A,A)') "Reloading ", TRIM(param_type), " parameters from file: ", TRIM(ADJUSTL(filename))
 END IF

 ! Check if parameter file exists
 INQUIRE(FILE=TRIM(ADJUSTL(filename)), EXIST=bExist)
 IF (bExist) THEN
   myDataFile = TRIM(ADJUSTL(filename))
 ELSE
   IF (myid == showid) THEN
     WRITE(MTERM,'(A,A)') "WARNING: File does not exist: ", TRIM(ADJUSTL(filename))
   END IF
 END IF

END SUBROUTINE ReloadParameterFile


!===============================================================================
SUBROUTINE Finalize(MFILE,MTERM)
!-------------------------------------------------------------------------------
! DESCRIPTION:
!   Standard finalization routine for the solver.
!   Saves the final solution and terminates MPI.
!-------------------------------------------------------------------------------
 USE PP3D_MPI, ONLY : myid,master,showid,Barrier_myMPI
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: MFILE, MTERM
 INTEGER :: iERR

 ! Save the final solution vector in unformatted form
 CALL SolToFile(-1)

 IF (myid == showid) THEN
   WRITE(MTERM,*) "PP3D_LES has successfully finished. "
   WRITE(MFILE,*) "PP3D_LES has successfully finished. "
 END IF

 CALL Barrier_myMPI()
 CALL MPI_Finalize(ierr)

END SUBROUTINE Finalize


!===============================================================================
SUBROUTINE Finalize_Particles(MFILE,MTERM)
!-------------------------------------------------------------------------------
! DESCRIPTION:
!   Finalization routine for particle simulations.
!   Similar to Finalize but without saving the solution vector.
!-------------------------------------------------------------------------------
 USE PP3D_MPI, ONLY : myid,master,showid,Barrier_myMPI
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: MFILE, MTERM
 INTEGER :: iERR

!  ! Save the final solution vector in unformatted form
!  CALL SolToFile(-1)

 IF (myid == showid) THEN
   WRITE(MTERM,*) "PP3D_LES has successfully finished. "
   WRITE(MFILE,*) "PP3D_LES has successfully finished. "
 END IF

 CALL Barrier_myMPI()
 CALL MPI_Finalize(ierr)

END SUBROUTINE Finalize_Particles

END MODULE ProcCtrl_mod
