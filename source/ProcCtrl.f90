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
 USE pco_coarse_operator, ONLY : PCO_PRESSURE_PARAMETERS_RELOADED, PCO_PRESSURE_RELOAD_ADMIT
 USE pressure_param_agreement, ONLY : PRESSURE_PARAMETERS_AGREE
 USE types, ONLY : tParamP
 USE PP3D_MPI, ONLY : myid,master,showid,ShareValueK_myMPI,&
                      Barrier_myMPI,ShareValueC_myMPI
 USE Transport_Q2P1, ONLY : myDataFile,QuadSc,LinSc,Properties
 USE param_parser, ONLY : GDATNEW,GetVeloParameters,GetPresParameters,GetPhysiclaParameters,&
                          ValidateSolverTypes
 USE timestep_control, ONLY : SetSimulationTimeStep
 USE solution_io, ONLY : write_sol_to_file
 USE EL_CONFIG, ONLY : EL_INLET_RELOAD_NOTICE
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: MFILE, MTERM
 CHARACTER(len=200) :: string, line
 CHARACTER(len=100) :: command, cvalue
 CHARACTER(len=7)   :: CSimPar
 INTEGER :: iEnd(1), iCmnd(1), iExist(1), iPos, iFile, istat
 LOGICAL :: bExist, bAgreed
 TYPE(tParamP) :: PresPrmOld
 DOUBLE PRECISION :: TSTEP, THETA, THSTEP, TIMENS, EPSNS
 INTEGER :: NITNS, ITNS
 COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
 INTEGER :: INSAV, INSAVN
 COMMON /NSSAV/  INSAV, INSAVN

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

   ! Parse command and value from the control string.
   ! Empty lines and comment lines are ignored silently.
   line = TRIM(ADJUSTL(string))
   command = ""
   cvalue  = ""

   IF (LEN_TRIM(line) > 0) THEN
    IF (line(1:1) /= "#" .AND. line(1:1) /= "!" .AND. line(1:1) /= ";") THEN
      iPos = INDEX(line,"=")
      IF (iPos /= 0) THEN
        READ(line(1:iPos-1),'(A)') command
        READ(line(iPos+1:),'(A)') cvalue
        cvalue = TRIM(ADJUSTL(cvalue))
      ELSE
        IF (myid == showid) THEN
          WRITE(MTERM,'(A,A)') 'WARNING: Malformed line in ProcCtrl.txt: ', TRIM(line)
        END IF
      END IF
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
     ! Complete restart dump (write_sol_to_file also writes the
     ! MaterialDistribution field; the legacy SolToFile dump is not
     ! restart-complete - v24f seg-2 lesson, 2026-09-04).
     CALL write_sol_to_file(INSAVN, TIMENS, iFile)
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
    IF (bExist) THEN
     CALL GetVeloParameters(QuadSc%prm, QuadSc%cName, MFILE)
     ! A runtime reload must not be able to smuggle in an unsupported type.
     CALL ValidateSolverTypes(QuadSc%prm%MGprmIn%CrsSolverType,&
                              LinSc%prm%MGprmIn%CrsSolverType,&
                              LinSc%prm%MGprmIn%MinLev,&
                              LinSc%prm%MGprmIn%MedLev,MFILE)
    END IF

! -------------------------------------------------------------------------------
   CASE ("Reload_Pres")
    CALL ReloadParameterFile(MFILE, MTERM, cvalue, "Pressure", bExist)
    IF (bExist) THEN
     PresPrmOld = LinSc%prm
     CALL GetPresParameters(LinSc%prm, LinSc%cName, MFILE)
     ! every rank parsed the file on its own: the ranks must run with the same settings (the
     ! number of smoothing steps or cycles decides the number of exchanges of each rank)
     CALL PRESSURE_PARAMETERS_AGREE(LinSc%prm, PresPrmOld, bAgreed)
     IF (.NOT.bAgreed .AND. myid == showid) THEN
       WRITE(MTERM,'(A)') 'WARNING: the ranks read different pressure parameters from the reload file; '// &
                          'previous values kept on every rank'
       WRITE(MFILE,'(A)') 'WARNING: the ranks read different pressure parameters from the reload file; '// &
                          'previous values kept on every rank'
     END IF
     ! a periodic coarse operator undoes what it cannot take (before the generic validation,
     ! which would abort on it)
     CALL PCO_PRESSURE_RELOAD_ADMIT()
     ! A runtime reload must not be able to smuggle in an unsupported type.
     CALL ValidateSolverTypes(QuadSc%prm%MGprmIn%CrsSolverType,&
                              LinSc%prm%MGprmIn%CrsSolverType,&
                              LinSc%prm%MGprmIn%MinLev,&
                              LinSc%prm%MGprmIn%MedLev,MFILE)
     ! a periodic coarse operator rebuilds its solver copy for a new type 1-4
     CALL PCO_PRESSURE_PARAMETERS_RELOADED()
    END IF

! -------------------------------------------------------------------------------
   CASE ("Reload_SimPar")
    CALL ReloadParameterFile(MFILE, MTERM, cvalue, "Simulation", bExist)
    IF (bExist) THEN
      CSimPar = "SimPar"
      CALL GDATNEW(CSimPar, 1)
      CALL SetSimulationTimeStep(TSTEP)
      ! the SimPar@EL* keys are start-up settings (the generic parser reads none of them): an
      ! inlet definition the reload file would change is reported, the run keeps its own
      CALL EL_INLET_RELOAD_NOTICE(myDataFile, MFILE, MTERM, myid == showid)
    END IF

! -------------------------------------------------------------------------------
   CASE ("Reload_Prop")
    CALL ReloadParameterFile(MFILE, MTERM, cvalue, "Physical", bExist)
    IF (bExist) CALL GetPhysiclaParameters(Properties, Properties%cName, MFILE)

   END SELECT
!================================================================================

   CALL Barrier_myMPI()

  END DO

  ! Reset control file by truncating it to zero length.
  IF (myid == showid) THEN
    REWIND(PROCCTRL_UNIT)
    ENDFILE(PROCCTRL_UNIT)
    CLOSE(PROCCTRL_UNIT)
  END IF

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
!   bExist     - Returns .TRUE. if file exists and reload should proceed; the same on
!                every rank (the file must be visible to all of them, agreed by an
!                allreduce), so the reload branches and their collectives are entered by
!                all ranks or by none. Collective: every rank calls this for the command.
!-------------------------------------------------------------------------------
 USE PP3D_MPI, ONLY : myid, showid, MPI_COMM_WORLD, MPI_LOGICAL, MPI_LAND
 USE Transport_Q2P1, ONLY : myDataFile
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: MFILE, MTERM
 CHARACTER(len=*), INTENT(IN) :: filename, param_type
 LOGICAL, INTENT(OUT) :: bExist
 LOGICAL :: bHere
 INTEGER :: ierr

 ! Log reload action to terminal and file
 IF (myid == showid) THEN
   WRITE(MTERM,'(A,A,A,A)') "Reloading ", TRIM(param_type), " parameters from file: ", TRIM(ADJUSTL(filename))
   WRITE(MFILE,'(A,A,A,A)') "Reloading ", TRIM(param_type), " parameters from file: ", TRIM(ADJUSTL(filename))
 END IF

 ! Check if parameter file exists (on every rank)
 INQUIRE(FILE=TRIM(ADJUSTL(filename)), EXIST=bHere)
 CALL MPI_Allreduce(bHere, bExist, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
 IF (bExist) THEN
   myDataFile = TRIM(ADJUSTL(filename))
 ELSE
   IF (myid == showid) THEN
     WRITE(MTERM,'(A,A)') "WARNING: File does not exist (on every rank): ", TRIM(ADJUSTL(filename))
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

SUBROUTINE ProcessControl_External(MFILE,MTERM)
 USE ProcCtrl_mod, ONLY : ProcessControl
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: MFILE, MTERM

 CALL ProcessControl(MFILE,MTERM)
END SUBROUTINE ProcessControl_External
