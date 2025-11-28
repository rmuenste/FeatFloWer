MODULE param_parser
!-------------------------------------------------------------------------------------------------
! Parameter parsing module for q2p1_param.dat file
!
! This module contains routines to parse velocity (Velo@), pressure (Pres@),
! physical property (Prop@), and simulation (SimPar@) parameters from the parameter file.
!
! Centralized parameter parsing: GDATNEW, GetVeloParameters, GetPresParameters, GetPhysiclaParameters
!-------------------------------------------------------------------------------------------------
USE PP3D_MPI, ONLY: myid, showID, master
USE iniparser
USE var_QuadScalar, ONLY: myDataFile, GAMMA, iCommSwitch, BaSynch, &
  myMatrixRenewal, bNonNewtonian, cGridFileName, nSubCoarseMesh, cFBM_File, &
  bTracer, cProjectFile, bMeshAdaptation, myExport, cAdaptedMeshFile, &
  nUmbrellaSteps, nInitUmbrellaSteps, bNoOutflow, bViscoElastic, bRefFrame, &
  bSteadyState, Properties, dCGALtoRealFactor, nUmbrellaStepsLvl, &
  nMainUmbrellaSteps, bBoundaryCheck, Transform, postParams, &
  ProlongationDirection, bNS_Stabilization, b2DViscoBench, b3DViscoBench, &
  SSE_HAS_ANGLE, extruder_angle, ApplicationString, VersionString, &
  MaxLevelKnownToMaster, GammaDot, AlphaRelax, RadParticle
USE types, ONLY: tParamV, tParamP, tProperties

IMPLICIT NONE

! Parameters
INTEGER, PARAMETER :: NNLEV = 9

! COMMON block variables - /OUTPUT/
INTEGER :: M, MT, MKEYB, MTERM, MERR, MPROT, MSYS, MTRC, IRECL8

! COMMON block variables - /MGPAR/
INTEGER :: ILEV, NLEV, NLMIN, NLMAX, ICYCLE
INTEGER :: KPRSM(NNLEV), KPOSM(NNLEV)

! COMMON block variables - /NSPAR/
DOUBLE PRECISION :: TSTEP, THETA, THSTEP, TIMENS, EPSNS
INTEGER :: NITNS, ITNS

! COMMON block variables - /NSADAT/
DOUBLE PRECISION :: TIMEMX, DTMIN, DTMAX, DTFACT, TIMEIN
DOUBLE PRECISION :: EPSADI, EPSADL, EPSADU
DOUBLE PRECISION :: PRDIF1, PRDIF2
INTEGER :: IEPSAD, IADIN, IREPIT, IADTIM

! COMMON block variables - /NSSAV/
INTEGER :: INSAV, INSAVN

! COMMON block variables - /NSSAVF/
DOUBLE PRECISION :: DTFILM, DTFILO, DTAVS, DTAVSO, DTGMV, DTGMVO
INTEGER :: IFUSAV, IFPSAV, IFXSAV, IGID, IGMV, IFINIT

! COMMON block variables - /FILES/
INTEGER :: IMESH1, MMESH1, MFILE1
INTEGER :: ISTART, MSTART
INTEGER :: ISOL, MSOL
CHARACTER(len=60) :: CPARM1, CMESH1, CFILE1, CSTART, CSOL

! COMMON block declarations
COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,ICYCLE,KPRSM,KPOSM
COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
COMMON /NSADAT/ TIMEMX,DTMIN,DTMAX,DTFACT,TIMEIN,EPSADI,EPSADL,&
                EPSADU,IEPSAD,IADIN,IREPIT,IADTIM,PRDIF1,PRDIF2
COMMON /NSSAV/  INSAV,INSAVN
COMMON /NSSAVF/ DTFILM,DTFILO,DTAVS,DTAVSO,DTGMV,DTGMVO,&
                IFUSAV,IFPSAV,IFXSAV,IGID,IGMV,IFINIT
COMMON /FILES/ IMESH1,MMESH1,CPARM1,CMESH1,MFILE1,CFILE1,&
               ISTART,MSTART,CSTART,ISOL,MSOL,CSOL

PRIVATE
PUBLIC :: GDATNEW
PUBLIC :: GetVeloParameters, GetPresParameters, GetPhysiclaParameters
PUBLIC :: write_param_real, write_param_int

CONTAINS

!-------------------------------------------------------------------------------------------------
! Shared helper subroutine to parse parameter file lines
! Finds positions of '@' and '=' characters in a string
!-------------------------------------------------------------------------------------------------
SUBROUTINE StrStuct(string, iAt, iEq, bOK)
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: string
  INTEGER, INTENT(OUT) :: iAt, iEq
  LOGICAL, INTENT(OUT) :: bOK
  INTEGER :: i, n

  n = LEN(string)
  iAt = 0
  iEq = 0
  DO i = 1, n
    IF (string(i:i) == '@') iAt = i
    IF (string(i:i) == '=') iEq = i
  END DO

  bOk = .false.
  IF (iAt /= 0 .and. iEq /= 0) bOk = .true.

END SUBROUTINE StrStuct

!-------------------------------------------------------------------------------------------------
! Parse velocity solver parameters (Velo@ section)
!-------------------------------------------------------------------------------------------------
SUBROUTINE GetVeloParameters(myParam, cName, mfile)
  IMPLICIT NONE
  TYPE(tParamV), INTENT(INOUT) :: myParam
  CHARACTER(len=7), INTENT(IN) :: cName
  INTEGER, INTENT(IN) :: mfile

  ! Local variables
  INTEGER :: iEnd, iAt, iEq, istat
  INTEGER :: myFile = 90909090
  CHARACTER(len=100) :: string
  CHARACTER(len=50) :: param, cPar
  CHARACTER(len=7) :: cVar
  CHARACTER(len=20) :: out_string
  LOGICAL :: bOK

  OPEN (UNIT=myFile, FILE=TRIM(ADJUSTL(myDataFile)), action='read', iostat=istat)
  IF (istat /= 0) THEN
    WRITE(*,*) "Could not open data file: ", myDataFile
    STOP
  END IF

  ! Initialize default values
  myParam%MGprmIn%RLX = 0.66d0
  myParam%MGprmIn%CrsSolverType = 1

  DO
    READ (UNIT=myFile, FMT='(A100)', IOSTAT=iEnd) string
    IF (iEnd == -1) EXIT
    CALL StrStuct(string, iAt, iEq, bOK)

    IF (bOK) THEN
      READ(string(1:iAt-1), *) cVar

      IF (TRIM(ADJUSTL(cVar)) == TRIM(ADJUSTL(cName))) THEN
        READ(string(iAt+1:iEq-1), *) cPar

        SELECT CASE (TRIM(ADJUSTL(cPar)))

          ! Note: gfortran requires explicit format length for integers
          CASE ("iMass")
            READ(string(iEq+1:), *) myParam%iMass
            WRITE(out_string, '(I20)') myParam%iMass
            CALL write_param_int(mfile, cVar, cPar, out_string, myParam%iMass)

          CASE ("SolvType")
            READ(string(iEq+1:), *) param
            IF (TRIM(ADJUSTL(param)) == "Jacobi") myParam%SolvType = 1
            IF (TRIM(ADJUSTL(param)) == "BiCGSTAB") myParam%SolvType = 2
            CALL write_param_int(mfile, cVar, cPar, out_string, myParam%SolvType)

          CASE ("defCrit")
            READ(string(iEq+1:), *) myParam%defCrit
            CALL write_param_real(mfile, cVar, cPar, out_string, myParam%defCrit)

          CASE ("Alpha")
            READ(string(iEq+1:), *) myParam%Alpha
            CALL write_param_real(mfile, cVar, cPar, out_string, myParam%Alpha)

          CASE ("MinDef")
            READ(string(iEq+1:), *) myParam%MinDef
            CALL write_param_real(mfile, cVar, cPar, out_string, myParam%MinDef)

          CASE ("NLmin")
            READ(string(iEq+1:), *) myParam%NLmin
            CALL write_param_int(mfile, cVar, cPar, out_string, myParam%NLmin)

          CASE ("NLmax")
            READ(string(iEq+1:), *) myParam%NLmax
            CALL write_param_int(mfile, cVar, cPar, out_string, myParam%NLmax)

          CASE ("MGMinLev")
            READ(string(iEq+1:), *) myParam%MGprmIn%MinLev
            CALL write_param_int(mfile, cVar, cPar, out_string, myParam%MGprmIn%MinLev)

          CASE ("MGMedLev")
            READ(string(iEq+1:), *) myParam%MGprmIn%MedLev
            CALL write_param_int(mfile, cVar, cPar, out_string, myParam%MGprmIn%MedLev)

          CASE ("MGMinIterCyc")
            READ(string(iEq+1:), *) myParam%MGprmIn%MinIterCycle
            CALL write_param_int(mfile, cVar, cPar, out_string, myParam%MGprmIn%MinIterCycle)

          CASE ("MGMaxIterCyc")
            READ(string(iEq+1:), *) myParam%MGprmIn%MaxIterCycle
            CALL write_param_int(mfile, cVar, cPar, out_string, myParam%MGprmIn%MaxIterCycle)

          CASE ("MGSmoothSteps")
            READ(string(iEq+1:), *) myParam%MGprmIn%nSmootherSteps
            CALL write_param_int(mfile, cVar, cPar, out_string, myParam%MGprmIn%nSmootherSteps)

          CASE ("MGIterCoarse")
            READ(string(iEq+1:), *) myParam%MGprmIn%nIterCoarse
            CALL write_param_int(mfile, cVar, cPar, out_string, myParam%MGprmIn%nIterCoarse)

          CASE ("MG_VANKA")
            READ(string(iEq+1:), *) myParam%MGprmIn%VANKA
            CALL write_param_int(mfile, cVar, cPar, out_string, myParam%MGprmIn%VANKA)

          CASE ("MGCrsRelaxPrm")
            READ(string(iEq+1:), *) myParam%MGprmIn%CrsRelaxPrm
            CALL write_param_real(mfile, cVar, cPar, out_string, myParam%MGprmIn%CrsRelaxPrm)

          CASE ("MGCrsRelaxParPrm")
            READ(string(iEq+1:), *) myParam%MGprmIn%CrsRelaxParPrm
            CALL write_param_real(mfile, cVar, cPar, out_string, myParam%MGprmIn%CrsRelaxParPrm)

          CASE ("MGCrsSolverType")
            READ(string(iEq+1:), *) myParam%MGprmIn%CrsSolverType
            CALL write_param_int(mfile, cVar, cPar, out_string, myParam%MGprmIn%CrsSolverType)

          CASE ("MGDefImprCoarse")
            READ(string(iEq+1:), *) myParam%MGprmIn%DefImprCoarse
            CALL write_param_real(mfile, cVar, cPar, out_string, myParam%MGprmIn%DefImprCoarse)

          CASE ("MGCriterion1")
            READ(string(iEq+1:), *) myParam%MGprmIn%Criterion1
            CALL write_param_real(mfile, cVar, cPar, out_string, myParam%MGprmIn%Criterion1)

          CASE ("MGCriterion2")
            READ(string(iEq+1:), *) myParam%MGprmIn%Criterion2
            CALL write_param_real(mfile, cVar, cPar, out_string, myParam%MGprmIn%Criterion2)

          CASE ("MGCycType")
            READ(string(iEq+1:), *) param
            myParam%MGprmIn%CycleType = TRIM(ADJUSTL(param))
            IF (myid == showid) WRITE(mterm, '(A,A)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", myParam%MGprmIn%CycleType
            IF (myid == showid) WRITE(mfile, '(A,A)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", myParam%MGprmIn%CycleType

          CASE ("MGRelaxPrm")
            READ(string(iEq+1:), *) myParam%MGprmIn%RLX
            CALL write_param_real(mfile, cVar, cPar, out_string, myParam%MGprmIn%RLX)
        END SELECT

      END IF
    END IF
  END DO

  CLOSE (myFile)

  ! Set maximum level and write final parameter
  myParam%MGprmIn%MaxLev = NLMAX
  out_string = ""
  WRITE(out_string, '(I20)') myParam%iMass
  CALL write_param_int(mfile, "Velo", "MGMaxLev", out_string, myParam%MGprmIn%MaxLev)

END SUBROUTINE GetVeloParameters

!-------------------------------------------------------------------------------------------------
! Parse pressure solver parameters (Pres@ section)
!-------------------------------------------------------------------------------------------------
SUBROUTINE GetPresParameters(myParam, cName, mfile)
  IMPLICIT NONE
  TYPE(tParamP), INTENT(INOUT) :: myParam
  CHARACTER(len=7), INTENT(IN) :: cName
  INTEGER, INTENT(IN) :: mfile

  ! Local variables
  INTEGER :: iEnd, iAt, iEq, istat
  INTEGER :: myFile = 90909090
  CHARACTER(len=100) :: string
  CHARACTER(len=50) :: param, cPar
  CHARACTER(len=7) :: cVar
  CHARACTER(len=20) :: out_string
  LOGICAL :: bOK

  OPEN (UNIT=myFile, FILE=TRIM(ADJUSTL(myDataFile)), action='read', iostat=istat)
  IF (istat /= 0) THEN
    WRITE(*,*) "Could not open data file: ", myDataFile
    STOP
  END IF

  ! Initialize default values
  myParam%MGprmIn%RLX = 0.66d0

  DO
    READ (UNIT=myFile, FMT='(A100)', IOSTAT=iEnd) string
    IF (iEnd == -1) EXIT
    CALL StrStuct(string, iAt, iEq, bOK)

    IF (bOK) THEN
      READ(string(1:iAt-1), *) cVar

      IF (TRIM(ADJUSTL(cVar)) == TRIM(ADJUSTL(cName))) THEN
        READ(string(iAt+1:iEq-1), *) cPar

        SELECT CASE (TRIM(ADJUSTL(cPar)))

          CASE ("MGMinLev")
            READ(string(iEq+1:), *) myParam%MGprmIn%MinLev
            CALL write_param_int(mfile, cVar, cPar, out_string, myParam%MGprmIn%MinLev)

          CASE ("MGMedLev")
            READ(string(iEq+1:), *) myParam%MGprmIn%MedLev
            CALL write_param_int(mfile, cVar, cPar, out_string, myParam%MGprmIn%MedLev)

          CASE ("MGMinIterCyc")
            READ(string(iEq+1:), *) myParam%MGprmIn%MinIterCycle
            CALL write_param_int(mfile, cVar, cPar, out_string, myParam%MGprmIn%MinIterCycle)

          CASE ("MGMaxIterCyc")
            READ(string(iEq+1:), *) myParam%MGprmIn%MaxIterCycle
            CALL write_param_int(mfile, cVar, cPar, out_string, myParam%MGprmIn%MaxIterCycle)

          CASE ("MGSmoothSteps")
            READ(string(iEq+1:), *) myParam%MGprmIn%nSmootherSteps
            CALL write_param_int(mfile, cVar, cPar, out_string, myParam%MGprmIn%nSmootherSteps)

          CASE ("MGIterCoarse")
            READ(string(iEq+1:), *) myParam%MGprmIn%nIterCoarse
            CALL write_param_int(mfile, cVar, cPar, out_string, myParam%MGprmIn%nIterCoarse)

          CASE ("MGDefImprCoarse")
            READ(string(iEq+1:), *) myParam%MGprmIn%DefImprCoarse
            CALL write_param_real(mfile, cVar, cPar, out_string, myParam%MGprmIn%DefImprCoarse)

          CASE ("MGCrsRelaxPrm")
            READ(string(iEq+1:), *) myParam%MGprmIn%CrsRelaxPrm
            CALL write_param_real(mfile, cVar, cPar, out_string, myParam%MGprmIn%CrsRelaxPrm)

          CASE ("MGCriterion1")
            READ(string(iEq+1:), *) myParam%MGprmIn%Criterion1
            CALL write_param_real(mfile, cVar, cPar, out_string, myParam%MGprmIn%Criterion1)

          CASE ("MGCriterion2")
            READ(string(iEq+1:), *) myParam%MGprmIn%Criterion2
            CALL write_param_real(mfile, cVar, cPar, out_string, myParam%MGprmIn%Criterion2)

          CASE ("MGCycType")
            READ(string(iEq+1:), *) param
            myParam%MGprmIn%CycleType = TRIM(ADJUSTL(param))
            IF (myid == showid) WRITE(mterm, '(A,A)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", myParam%MGprmIn%CycleType
            IF (myid == showid) WRITE(mfile, '(A,A)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", myParam%MGprmIn%CycleType

          CASE ("MGCrsSolverType")
            READ(string(iEq+1:), *) myParam%MGprmIn%CrsSolverType
            CALL write_param_int(mfile, cVar, cPar, out_string, myParam%MGprmIn%CrsSolverType)

          CASE ("MGRelaxPrm")
            READ(string(iEq+1:), *) myParam%MGprmIn%RLX
            CALL write_param_real(mfile, cVar, cPar, out_string, myParam%MGprmIn%RLX)

        END SELECT

      END IF
    END IF
  END DO

  CLOSE (myFile)

END SUBROUTINE GetPresParameters

!-------------------------------------------------------------------------------------------------
! Parse physical properties (Prop@ section)
! Note: Function name has typo "Physicla" instead of "Physical" - preserved for compatibility
!-------------------------------------------------------------------------------------------------
SUBROUTINE GetPhysiclaParameters(Props, cName, mfile)
  IMPLICIT NONE
  TYPE(tProperties), INTENT(INOUT) :: Props
  CHARACTER(len=7), INTENT(IN) :: cName
  INTEGER, INTENT(IN) :: mfile

  ! Local variables
  INTEGER :: iEnd, iAt, iEq, istat
  INTEGER :: myFile = 90909090
  CHARACTER(len=100) :: string
  CHARACTER(len=50) :: param, cPar
  CHARACTER(len=7) :: cVar
  LOGICAL :: bOK

  OPEN (UNIT=myFile, FILE=TRIM(ADJUSTL(myDataFile)), action='read', iostat=istat)
  IF (istat /= 0) THEN
    WRITE(*,*) "Could not open data file: ", myDataFile
    STOP
  END IF

  IF (myid == showid) WRITE(mfile, '(47("-"),A10,47("-"))') TRIM(ADJUSTL(cName))
  IF (myid == showid) WRITE(mterm, '(47("-"),A10,47("-"))') TRIM(ADJUSTL(cName))

  DO
    READ (UNIT=myFile, FMT='(A100)', IOSTAT=iEnd) string
    IF (iEnd == -1) EXIT
    CALL StrStuct(string, iAt, iEq, bOK)

    IF (bOK) THEN
      READ(string(1:iAt-1), *) cVar

      IF (TRIM(ADJUSTL(cVar)) == TRIM(ADJUSTL(cName))) THEN
        READ(string(iAt+1:iEq-1), *) cPar

        SELECT CASE (TRIM(ADJUSTL(cPar)))

          CASE ("nTPSubSteps")
            READ(string(iEq+1:), *) Props%nTPSubSteps
            IF (myid == showid) WRITE(mterm, '(A,I5)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%nTPSubSteps
            IF (myid == showid) WRITE(mfile, '(A,I5)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%nTPSubSteps

          CASE ("nTPFSubSteps")
            READ(string(iEq+1:), *) Props%nTPFSubSteps
            IF (myid == showid) WRITE(mterm, '(A,I5)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%nTPFSubSteps
            IF (myid == showid) WRITE(mfile, '(A,I5)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%nTPFSubSteps

          CASE ("nTPIterations")
            READ(string(iEq+1:), *) Props%nTPIterations
            IF (myid == showid) WRITE(mterm, '(A,I5)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%nTPIterations
            IF (myid == showid) WRITE(mfile, '(A,I5)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%nTPIterations

          CASE ("Material")
            READ(string(iEq+1:), *) Props%Material
            IF (myid == showid) WRITE(mterm, '(A,A)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%Material
            IF (myid == showid) WRITE(mfile, '(A,A)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%Material

          CASE ("Gravity")
            READ(string(iEq+1:), *) Props%Gravity
            IF (myid == showid) WRITE(mterm, '(A,3E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%Gravity
            IF (myid == showid) WRITE(mfile, '(A,3E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%Gravity

          CASE ("Density")
            READ(string(iEq+1:), *) Props%Density
            IF (myid == showid) WRITE(mterm, '(A,2E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%Density
            IF (myid == showid) WRITE(mfile, '(A,2E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%Density

          CASE ("Viscosity")
            READ(string(iEq+1:), *) Props%Viscosity
            IF (myid == showid) WRITE(mterm, '(A,2E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%Viscosity
            IF (myid == showid) WRITE(mfile, '(A,2E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%Viscosity

          CASE ("DiffCoeff")
            READ(string(iEq+1:), *) Props%DiffCoeff
            IF (myid == showid) WRITE(mterm, '(A,2E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%DiffCoeff
            IF (myid == showid) WRITE(mfile, '(A,2E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%DiffCoeff

          CASE ("Sigma")
            READ(string(iEq+1:), *) Props%Sigma
            IF (myid == showid) WRITE(mterm, '(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%Sigma
            IF (myid == showid) WRITE(mfile, '(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%Sigma

          CASE ("DiracEps")
            READ(string(iEq+1:), *) Props%DiracEps
            IF (myid == showid) WRITE(mterm, '(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%DiracEps
            IF (myid == showid) WRITE(mfile, '(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%DiracEps

          CASE ("PowerLawExp")
            READ(string(iEq+1:), *) Props%PowerLawExp
            IF (myid == showid) WRITE(mterm, '(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%PowerLawExp
            IF (myid == showid) WRITE(mfile, '(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%PowerLawExp

          CASE ("GAMMA")
            READ(string(iEq+1:), *) GAMMA
            IF (myid == showid) WRITE(mterm, '(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", GAMMA
            IF (myid == showid) WRITE(mfile, '(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", GAMMA

          CASE ("ViscoLambda")
            READ(string(iEq+1:), *) Props%ViscoLambda
            IF (myid == showid) WRITE(mterm, '(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%ViscoLambda
            IF (myid == showid) WRITE(mfile, '(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%ViscoLambda

          CASE ("ViscoAlphaImp")
            READ(string(iEq+1:), *) Props%ViscoAlphaImp
            IF (myid == showid) WRITE(mterm, '(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%ViscoAlphaImp
            IF (myid == showid) WRITE(mfile, '(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%ViscoAlphaImp

          CASE ("ViscoAlphaExp")
            READ(string(iEq+1:), *) Props%ViscoAlphaExp
            IF (myid == showid) WRITE(mterm, '(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%ViscoAlphaExp
            IF (myid == showid) WRITE(mfile, '(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%ViscoAlphaExp

          CASE ("NS_StabAlpha_Imp")
            READ(string(iEq+1:), *) Props%NS_StabAlpha_Imp
            IF (myid == showid) WRITE(mterm, '(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%NS_StabAlpha_Imp
            IF (myid == showid) WRITE(mfile, '(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%NS_StabAlpha_Imp

          CASE ("NS_StabAlpha_Exp")
            READ(string(iEq+1:), *) Props%NS_StabAlpha_Exp
            IF (myid == showid) WRITE(mterm, '(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%NS_StabAlpha_Exp
            IF (myid == showid) WRITE(mfile, '(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%NS_StabAlpha_Exp

          CASE ("ViscoModel")
            READ(string(iEq+1:), *) Props%ViscoModel
            IF (myid == showid) WRITE(mterm, '(A,I2)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%ViscoModel
            IF (myid == showid) WRITE(mfile, '(A,I2)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%ViscoModel

          CASE ("ForceScale")
            READ(string(iEq+1:), *) Props%ForceScale
            IF (myid == showid) WRITE(mterm, '(A,6E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%ForceScale
            IF (myid == showid) WRITE(mfile, '(A,6E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", Props%ForceScale
        END SELECT

      END IF
    END IF
  END DO

  IF (myid == showid) WRITE(mfile, '(47("-"),A10,47("-"))') TRIM(ADJUSTL(cName))
  IF (myid == showid) WRITE(mterm, '(47("-"),A10,47("-"))') TRIM(ADJUSTL(cName))

  CLOSE (myFile)

END SUBROUTINE GetPhysiclaParameters

!-------------------------------------------------------------------------------------------------
! Helper subroutine to write real parameter to protocol files
!-------------------------------------------------------------------------------------------------
SUBROUTINE write_param_real(unitid, cVar, cPar, out_string, val)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: unitid
  CHARACTER(len=*), INTENT(IN) :: cVar
  CHARACTER(len=*), INTENT(IN) :: cPar
  CHARACTER(len=*), INTENT(INOUT) :: out_string
  DOUBLE PRECISION, INTENT(IN) :: val

  out_string = " "
  WRITE(out_string, '(E16.8)') val
  IF (myid == showid) WRITE(mterm, '(A,A)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", ADJUSTL(out_string)
  IF (myid == showid) WRITE(unitid, '(A,A)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", ADJUSTL(out_string)

END SUBROUTINE write_param_real

!-------------------------------------------------------------------------------------------------
! Helper subroutine to write integer parameter to protocol files
!-------------------------------------------------------------------------------------------------
SUBROUTINE write_param_int(unitid, cVar, cPar, out_string, val)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: unitid
  CHARACTER(len=*), INTENT(IN) :: cVar
  CHARACTER(len=*), INTENT(IN) :: cPar
  CHARACTER(len=*), INTENT(INOUT) :: out_string
  INTEGER, INTENT(IN) :: val

  out_string = " "
  WRITE(out_string, '(I20)') val
  IF (myid == showid) WRITE(mterm, '(A,A)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", ADJUSTL(out_string)
  IF (myid == showid) WRITE(unitid, '(A,A)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ", ADJUSTL(out_string)

END SUBROUTINE write_param_int

!-------------------------------------------------------------------------------------------------
! Parse simulation parameters (SimPar@ section)
! Modernized version extracted from Init.f90
!-------------------------------------------------------------------------------------------------
SUBROUTINE GDATNEW (cName,iCurrentStatus)
  ! File unit constants
  INTEGER, PARAMETER :: PARAM_FILE_UNIT = 888
  INTEGER, PARAMETER :: MESH_FILE_UNIT = 61
  INTEGER, PARAMETER :: PROTOCOL_FILE_UNIT = 62
  INTEGER, PARAMETER :: START_FILE_UNIT = 63
  INTEGER, PARAMETER :: SOL_FILE_UNIT = 64

  ! Subroutine arguments
  CHARACTER(len=7) :: cName
  INTEGER :: iCurrentStatus

  ! Local variables - integers
  INTEGER :: myFile = PARAM_FILE_UNIT
  INTEGER :: iEnd, iAt, iEq, iLen, istat
  INTEGER :: iOutShift
  INTEGER :: iVisco, iLoc, iangle, mylength, nFields, i
  INTEGER :: MFILE
  INTEGER, ALLOCATABLE :: iPos(:)

  ! Local variables - characters
  CHARACTER :: letter
  CHARACTER(len=500) :: string
  CHARACTER(len=7) :: cVar
  CHARACTER(len=25) :: cPar
  CHARACTER(len=400) :: cLongString
  CHARACTER(len=8) :: cParam
  CHARACTER(len=20) :: cParam2

  ! Local variables - logicals
  LOGICAL :: bOK, bOutNMAX
  LOGICAL :: is_open

  SAVE

  inquire(unit=myFile, OPENED=is_open)
  if (.not. is_open) then
    OPEN (UNIT=myFile,FILE=TRIM(ADJUSTL(myDataFile)),action="read",iostat=istat)
    if (istat /= 0) then
      write(*,'(A)') 'ERROR: Could not open parameter data file'
      write(*,'(A,A)') '  File: ', TRIM(myDataFile)
      write(*,'(A,I0)') '  iostat = ', istat
      write(*,'(A)') 'Please check that the file exists and is readable.'
      stop 1
    end if
  end if

  bOutNMAX = .false.
  DO
  READ (UNIT=myFile,FMT='(A500)',IOSTAT=iEnd) string
  IF (iEnd == -1) EXIT
  CALL StrStuct()
  IF (bOK) THEN

    READ(string(1:iAt-1),*) cVar

    IF (TRIM(ADJUSTL(cVar)) == TRIM(ADJUSTL(cName))) THEN

      READ(string(iAt+1:iEq-1),*) cPar
      SELECT CASE (TRIM(ADJUSTL(cPar)))

      CASE ("MeshFolder")
        READ(string(iEq+1:),*) cGridFileName
      CASE ("SubMeshNumber")
        READ(string(iEq+1:),*) nSubCoarseMesh
      CASE ("ParticleFile")
        READ(string(iEq+1:),*) cFBM_File
      CASE ("ProjectFile")
        READ(string(iEq+1:),*) cProjectFile
        MMESH1 = MESH_FILE_UNIT
      CASE ("ProtocolFile")
        READ(string(iEq+1:),*) CFILE1
        if (SSE_HAS_ANGLE)then
          CFILE1=''
          iangle = int(extruder_angle)
          write(cfile1,'(a, I4.4,a)') '_data/prot.',iangle,'.txt'
        end if
        MFILE1 = PROTOCOL_FILE_UNIT
        MFILE = MFILE1
      CASE ("StartingProc")
        READ(string(iEq+1:),*) ISTART
      CASE ("Umbrella")
        READ(string(iEq+1:),*) nUmbrellaSteps
      CASE ("InitUmbrella")
        READ(string(iEq+1:),*) nInitUmbrellaSteps
      CASE ("UmbrellaStepM")
        READ(string(iEq+1:),*) nMainUmbrellaSteps
      CASE ("UmbrellaStepL")
        READ(string(iEq+1:),*) nUmbrellaStepsLvl
      CASE ("StartFile")
        READ(string(iEq+1:),*) CSTART
        MSTART = START_FILE_UNIT
      CASE ("LoadAdaptedMesh")
        bMeshAdaptation = .true.
        READ(string(iEq+1:),*) cAdaptedMeshFile
      CASE ("SolFile")
        READ(string(iEq+1:),*) CSOL
        iLen = LEN(TRIM(ADJUSTL(CSOL)))
        IF (myid < 10) THEN
          WRITE(CSOL(iLen+1:),'(A,I1)') "00",myid
        ELSEIF (myid < 100) THEN
          WRITE(CSOL(iLen+1:),'(A,I2)') "0",myid
        ELSE
          WRITE(CSOL(iLen+1:),'(I3)') myid
        END IF
        MSOL = SOL_FILE_UNIT
        ISOL = 1
      CASE ("MinMeshLevel")
        READ(string(iEq+1:),*) NLMIN
      CASE ("MaxMeshLevel")
        IF (myid /= master) THEN
          READ(string(iEq+1:),*) NLMAX
          myExport%LevelMax = NLMAX
          MaxLevelKnownToMaster = NLMAX
        else
          READ(string(iEq+1:),*) MaxLevelKnownToMaster
          myExport%LevelMax = MaxLevelKnownToMaster
        END IF
      CASE ("TimeScheme")
        cParam = " "
        READ(string(iEq+1:),*) cParam
        THETA = 0.5d0
        IF (TRIM(ADJUSTL(cParam)) == "BE") THETA = 1.0d0
        IF (TRIM(ADJUSTL(cParam)) == "FE") THETA = 0.0d0
      CASE ("TimeStep")
        READ(string(iEq+1:),*) TSTEP
      CASE ("Bench_U_mean")
        READ(string(iEq+1:),*) postParams%U_mean
        write(*,*)'umean:',postParams%U_mean
      CASE ("Bench_H")
        READ(string(iEq+1:),*) postParams%H
        write(*,*)'h:',postParams%h
      CASE ("Bench_D")
        READ(string(iEq+1:),*) postParams%D
        write(*,*)'d:',postParams%d
      CASE ("Bench_Sc_U")
        READ(string(iEq+1:),*) postParams%Sc_U
      CASE ("Bench_Sc_Mu")
        READ(string(iEq+1:),*) postParams%Sc_Mu
      CASE ("Bench_Sc_a")
        READ(string(iEq+1:),*) postParams%Sc_a
      CASE ("TimeAdaptivity")
        IADTIM = 0
        IF (read_yes_no_param(string, iEq)) IADTIM = 1
       CASE ("ProlongationDirection")
        READ(string(iEq+1:),*) ProlongationDirection
       CASE ("Tracer")
        bTracer = read_yes_no_param(string, iEq)
      CASE ("ViscoElastic")
        bViscoElastic = read_yes_no_param(string, iEq)
      CASE ("ViscoElasticBench")
        READ(string(iEq+1:),*) iVisco
        b2DViscoBench = .false.
        b3DViscoBench = .false.
        if (iVisco == 2) b2DViscoBench = .true.
        if (iVisco == 3) b3DViscoBench = .true.
      CASE ("SteadyState")
        bSteadyState = read_yes_no_param(string, iEq)
      CASE ("ReferenceFrame")
        bRefFrame = read_yes_no_param(string, iEq)
      CASE ("NoOutflow")
        bNoOutflow = read_yes_no_param(string, iEq)
      CASE ("MinTimeAdapt")
        READ(string(iEq+1:),*) DTMIN
      CASE ("MaxTimeAdapt")
        READ(string(iEq+1:),*) DTMAX
      CASE ("StartSimTime")
        READ(string(iEq+1:),*) TIMENS
      CASE ("MaxSimTime")
        READ(string(iEq+1:),*) TIMEMX
      CASE ("MaxNumStep")
        READ(string(iEq+1:),*) NITNS
      CASE ("ElemTrans")
        READ(string(iEq+1:),*) Transform%ilint
        If (Transform%ilint.lt.1.or.Transform%ilint.gt.2) Transform%ilint = 2
      CASE ("BackUpFreq")
        READ(string(iEq+1:),*) INSAV
      CASE ("BackUpNum")
        READ(string(iEq+1:),*) INSAVN
      CASE ("BoundaryCheck")
        bBoundaryCheck = read_yes_no_param(string, iEq)
      CASE ("NS_Stabilization")
        bNS_Stabilization = read_yes_no_param(string, iEq)
      CASE ("OutputFreq")
        READ(string(iEq+1:),*) DTGMV
      CASE ("MatrixRenewal")
        READ(string(iEq+1:),*) cParam2
        iLoc=INDEX (cParam2, 'M', .TRUE.)+1
        READ(cParam2(iLoc:iLoc),'(I1)') myMatrixRenewal%M
        iLoc=INDEX (cParam2, 'D', .TRUE.)+1
        READ(cParam2(iLoc:iLoc),'(I1)') myMatrixRenewal%D
        iLoc=INDEX (cParam2, 'K', .TRUE.)+1
        READ(cParam2(iLoc:iLoc),'(I1)') myMatrixRenewal%K
        iLoc=INDEX (cParam2, 'C', .TRUE.)+1
        READ(cParam2(iLoc:iLoc),'(I1)') myMatrixRenewal%C
        iLoc=INDEX (cParam2, 'S', .TRUE.)+1
        READ(cParam2(iLoc:iLoc),'(I1)') myMatrixRenewal%S
      CASE ("FlowType")
        READ(string(iEq+1:),*) cParam2
        bNonNewtonian = .true.
        IF (TRIM(ADJUSTL(cParam2)) == "Newtonian") bNonNewtonian = .false.
      CASE ("OutputLevel")
        READ(string(iEq+1:),*) cParam2
        IF (TRIM(ADJUSTL(cParam2)) == "MAX") THEN
          bOutNMAX = .true.
          iOutShift = 0
        END IF
        IF (TRIM(ADJUSTL(cParam2)) == "MAX+1") THEN
          bOutNMAX = .true.
          iOutShift = 1
        END IF
        IF (TRIM(ADJUSTL(cParam2)) == "MAX-1") THEN
          bOutNMAX = .true.
          iOutShift = -1
        END IF
        IF (.not. bOutNMAX) THEN
          READ(string(iEq+1:),*) myExport%Level
          myExport%Level = MAX(MIN(myExport%Level,myExport%LevelMax+1),1)
        END IF

      CASE ("CGALtoRealFactor")
       READ(string(iEq+1:),*)dCGALtoRealFactor

       ! --- Add new parameters here ---
      CASE ("GammaDot")
       READ(string(iEq+1:),*) GammaDot
      CASE ("AlphaRelax")
       READ(string(iEq+1:),*) AlphaRelax
      CASE ("RadParticle")
       READ(string(iEq+1:),*) RadParticle

      CASE ("aSynchComm")
        baSynch = read_yes_no_param(string, iEq)
      CASE ("CommSwitch")
        READ(string(iEq+1:),*) iCommSwitch
      CASE ("OutputFormat")
        READ(string(iEq+1:),*) myExport%Format
      CASE ("OutputFields")
        READ(string(iEq+1:),*) cLongString
        cLongString = ADJUSTL(TRIM(cLongString))
        mylength = LEN(ADJUSTL(TRIM(cLongString)))
        nFields = 0
        DO i=1,mylength
        READ(cLongString(i:i),'(A)') letter
        IF (letter == ',') nFields = nFields + 1
        END DO
        IF (ALLOCATED(myExport%Fields)) DEALLOCATE(myExport%Fields)
        IF (ALLOCATED(iPos)) DEALLOCATE(iPos)
        ALLOCATE(myExport%Fields(nFields+1),iPos(nFields+2))
        iPos(1) = 0
        nFields = 0
        DO i=1,mylength
        READ(cLongString(i:i),'(A)') letter
        IF (letter == ',') THEN
          nFields = nFields + 1
          iPos(nFields+1) = i
        END IF
        END DO
        iPos(nFields+2) = mylength+1
        DO i=1,nFields+1
        READ(cLongString(iPos(i)+1:iPos(i+1)-1),'(A)') myExport%Fields(i)
        myExport%Fields(i) = ADJUSTL(TRIM(myExport%Fields(i)))
        END DO

      END SELECT

    END IF
  END IF
  END DO

  CLOSE (myFile)

  M     = 1
  MT    = 1
  IGMV   = MaxLevelKnownToMaster
  THSTEP=TSTEP*THETA
  IF (bOutNMAX) myExport%Level = MaxLevelKnownToMaster + iOutShift



  inquire(unit=mfile1, OPENED=is_open)
  if (.not. is_open) then
    IF (myid == showid) THEN
      OPEN (UNIT=mfile1,FILE=cfile1,action="write",status="replace",iostat=istat)
      if (istat /= 0) then
        write(*,'(A)') 'ERROR: Could not open protocol file for writing'
        write(*,'(A,A)') '  File: ', TRIM(cfile1)
        write(*,'(A,I0)') '  iostat = ', istat
        write(*,'(A)') 'Please check write permissions and disk space.'
        stop 1
      end if
    end if
  end if

  IF (iCurrentStatus == 0) THEN
    IF (myid == showid) WRITE(UNIT=mterm,FMT=101) ApplicationString,VersionString
    IF (myid == showid) WRITE(UNIT=mfile,FMT=101) ApplicationString,VersionString
  END IF

  ! Printout of all loaded parameters

  IF (myid == showid) THEN
    CALL write_param_int(mfile, mterm, "CommSwitch = ", iCommSwitch)

    IF (baSynch) THEN
     CALL write_param_str(mfile, mterm, "aSynchComm = ", "YES")
    ELSE
     CALL write_param_str(mfile, mterm, "aSynchComm = ", "NO")
    END IF

    CALL write_param_int(mfile, mterm, "nSubCoarseMesh = ", nSubCoarseMesh)
    CALL write_param_str(mfile, mterm, "ParticleFile = ", cFBM_File)
    CALL write_param_str(mfile, mterm, "ProjectFile = ", TRIM(CProjectFile))
    CALL write_param_int(mfile, mterm, "StartingProc = ", ISTART)
    CALL write_param_str(mfile, mterm, "StartFile = ", CSTART)
    CALL write_param_str(mfile, mterm, "SolFile = ", CSOL)
    CALL write_param_int(mfile, mterm, "MinMeshLevel = ", NLMIN)
    CALL write_param_int(mfile, mterm, "MaxMeshLevel = ", NLMAX)
    CALL write_param_real(mfile, mterm, "TimeScheme = ", THETA)
    CALL write_param_real(mfile, mterm, "TimeStep = ", TSTEP)
    CALL write_param_int(mfile, mterm, "TimeAdaptivity = ", IADTIM)
    CALL write_param_real(mfile, mterm, "MinTimeAdapt = ", DTMIN)
    CALL write_param_real(mfile, mterm, "MaxTimeAdapt = ", DTMAX)
    CALL write_param_real(mfile, mterm, "StartSimTime = ", TIMENS)

    CALL write_param_real(mfile, mterm, "MaxSimTime = ", TIMEMX)
    CALL write_param_int(mfile, mterm, "MaxNumStep = ", NITNS)
    CALL write_param_int(mfile, mterm, "BackUpFreq = ", INSAV)
    CALL write_param_int(mfile, mterm, "BackUpNum = ", INSAVN)
    CALL write_param_real(mfile, mterm, "OutputFreq = ", DTGMV)

    WRITE(mfile,'(A,I1)') "ElemTransform = Q", Transform%ilint
    WRITE(mterm,'(A,I1)') "ElemTransform = Q", Transform%ilint

    WRITE(mfile,'(A,5(A6I1))') "Matrix Renewal scheme : ","  M = ",myMatrixRenewal%M,&
      ", D = ",myMatrixRenewal%D,", K = ",myMatrixRenewal%K,", S = ",myMatrixRenewal%S,&
      ", C = ",myMatrixRenewal%C
    WRITE(mterm,'(A,5(A6I1))') "Matrix Renewal scheme : ","  M = ",myMatrixRenewal%M,&
      ", D = ",myMatrixRenewal%D,", K = ",myMatrixRenewal%K,", S = ",myMatrixRenewal%S,&
      ", C = ",myMatrixRenewal%C

    IF (bMeshAdaptation) THEN
      CALL write_param_str(mfile, mterm, "Use initial Mesh Adaptation file: ", &
                          ADJUSTL(TRIM(cAdaptedMeshFile)))
    ELSE
      CALL write_param_str(mfile, mterm, "No Initial Mesh Adaptation", "")
    END IF

    IF (bNS_Stabilization) THEN
      CALL write_param_str(mfile, mterm, "Stabilization of Navier-Stokes is ::  ", "ON")
    ELSE
      CALL write_param_str(mfile, mterm, "Stabilization of Navier-Stokes is ::  ", "OFF")
    END IF

    CALL write_param_int(mfile, mterm, "Number of Initial Umbrella smoothening steps", &
                        nInitUmbrellaSteps)
    CALL write_param_int(mfile, mterm, "Number of Umbrella smoothening steps: ", nUmbrellaSteps)
    CALL write_param_int(mfile, mterm, "Number of Umbrella Loops: ", nMainUmbrellaSteps)

    WRITE(mfile,'(A,10I10)') "Number of Umbrella steps per levels: ", nUmbrellaStepsLvl
    WRITE(mterm,'(A,10I10)') "Number of Umbrella steps per levels: ", nUmbrellaStepsLvl

    IF (bNoOutflow) THEN
      CALL write_param_str(mfile, mterm, &
        "Matrix modification is to be performed due to the NoOuflow Condition", "")
    END IF

    IF (bTracer) THEN
      CALL write_param_str(mfile, mterm, "Tracer equation is included", "")
    END IF

    IF (bViscoElastic) THEN
      CALL write_param_str(mfile, mterm, "Visco-elastic equation is included", "")
    END IF

    IF (b2DViscoBench .or. b3DViscoBench) THEN
      WRITE(mfile,'(A,I1,A)') "Visco-elastic benchmark computation for ",iVisco,"D"
      WRITE(mterm,'(A,I1,A)') "Visco-elastic benchmark computation for ",iVisco,"D"
    END IF

    IF (bBoundaryCheck) THEN
      CALL write_param_str(mfile, mterm, "BoundaryCheck is ", "ON")
    ELSE
      CALL write_param_str(mfile, mterm, "BoundaryCheck is ", "OFF")
    END IF

    WRITE(mfile,'(A,3ES14.4)') "Newtonian FAC Benchamrk params (U,H,D) : ",postParams%U_mean,postParams%H,postParams%D
    WRITE(mterm,'(A,3ES14.4)') "Newtonian FAC Benchamrk params (U,H,D) : ",postParams%U_mean,postParams%H,postParams%D

    WRITE(mfile,'(A,3ES14.4)') "Viscoelastic Benchamrk params (Sc_U,Sc_Mu,Sc_a) : ",postParams%Sc_U,postParams%Sc_Mu,postParams%Sc_a
    WRITE(mterm,'(A,3ES14.4)') "Viscoelastic Benchamrk params (Sc_U,Sc_Mu,Sc_a) : ",postParams%Sc_U,postParams%Sc_Mu,postParams%Sc_a

    IF (bNonNewtonian) THEN
      CALL write_param_str(mfile, mterm, "FlowType = ", "non-Newtonian")
    ELSE
      CALL write_param_str(mfile, mterm, "FlowType = ", "Newtonian")
    END IF

    ! Print new parameters
    CALL write_param_real(mfile, mterm, "GammaDot = ", GammaDot)
    CALL write_param_real(mfile, mterm, "AlphaRelax = ", AlphaRelax)
    CALL write_param_real(mfile, mterm, "RadParticle = ", RadParticle)


    IF (ProlongationDirection == 0) THEN
     WRITE(mfile,'(A,D12.4)') "Mesh Prolongation is set to  = STANDARD"
     WRITE(mterm,'(A,D12.4)') "Mesh Prolongation is set to  = STANDARD"
    ELSE
     IF (ProlongationDirection == 1) THEN
      WRITE(mfile,'(A,D12.4)') "Mesh Prolongation is set to  = CYLINDRICAL in X axis"
      WRITE(mterm,'(A,D12.4)') "Mesh Prolongation is set to  = CYLINDRICAL in X axis"
     END IF
     IF (ProlongationDirection == 2) THEN
      WRITE(mfile,'(A,D12.4)') "Mesh Prolongation is set to  = CYLINDRICAL in Y axis"
      WRITE(mterm,'(A,D12.4)') "Mesh Prolongation is set to  = CYLINDRICAL in Y axis"
     END IF
     IF (ProlongationDirection == 3) THEN
      WRITE(mfile,'(A,D12.4)') "Mesh Prolongation is set to  = CYLINDRICAL in Z axis"
      WRITE(mterm,'(A,D12.4)') "Mesh Prolongation is set to  = CYLINDRICAL in Z axis"
     END IF
    END IF

    CALL write_param_real(mfile, mterm, "CGALtoRealFactor = ", dCGALtoRealFactor)


    WRITE(mfile,'(4A,I3,3A,100A)') "Exporting into ", "'",myExport%Format,"' on level: ",&
      myExport%Level," Fields: '",("["//TRIM(ADJUSTL(myExport%Fields(i)))//"]",i=1,nFields+1),"'"
    WRITE(mterm,'(4A,I3,3A,100A)') "Exporting into ", "'",myExport%Format,"' on level: ",&
      myExport%Level," Fields: '",("["//TRIM(ADJUSTL(myExport%Fields(i)))//"]",i=1,nFields+1),"'"
  END IF


  !-----------------------------------------------------------------------
  101 FORMAT(/2X,100('=')/&
    2X,"|",10X,"                                                                                        |"/&
    2X,"|",10X,"Parellel Q2/P1 FEM Fluid Dynamics code          FeatFloWer                              |"/&
    A104,/&
    A104,/&
    2X,"|",10X,"                                                                                        |"/&
    2X,"|",10X,"Developed by:                                   Otto Mierka, Raphael Münster            |"/&
    2X,"|",10X,"under the supervision of:                       Stefan Turek and Dmitri Kuzmin          |"/&
    2X,"|",10X,"additional contributions from:                  Robert Jendrny, Christoph Lohmann       |"/&
    2X,"|",10X,"                                                Tatiana Theis, Malte Schuh              |"/&
    2X,"|",10X,"                                                Michael Fast                            |"/&
    2X,"|",10X,"                                                                                        |"/&
    2X,"|",10X,"Developed at:                                                                           |"/&
    2X,"|",10X,"                                                                                        |"/&
    2X,"|",10X,"                                                                                        |"/&
    2X,"|",10X,"     #####   #    #         ###     ###   ###   ##### #   #  #    #  #   #  ###         |"/&
    2X,"|",10X,"       #     #    #         #  #   #   #  #  #    #   ## ##  #    #  ##  #  #  #        |"/&
    2X,"|",10X,"       #     #    #   ###   #   #  #   #  ###     #   # # #  #    #  # # #  #   #       |"/&
    2X,"|",10X,"       #     #    #         #  #   #   #  #  #    #   #   #  #    #  #  ##  #  #        |"/&
    2X,"|",10X,"       #      ####          ###     ###   #   #   #   #   #   ####   #   #  ###         |"/&
    2X,"|",10X,"                                                                                        |"/&
    2X,"|",10X,"                                                                                        |"/&
    2X,"|",57X,"            Chair of Mathematics III",5X,"|"/&
    2X,"|",57X,"    Applied Mathematics and Numerics",5X,"|"/&
    2X,"|",57X,"                    Vogelopthsweg 87",5X,"|"/&
    2X,"|",57X,"                      Dortmund 44225",5X,"|"/&
    2X,"|",57X,"                                    ",5X,"|"/&
    2X,"|",10X,"Based on FeatFlow (c)     ",&
    "see also: http://www.featflow.de",30X,"|"/&
    2X,"|",10X,"Correspondance:",73X,"|"/&
    2X,"|",10X,"otto.mierka@math.tu-dortmund.de, ",&
    "stefan.turek@math.tu-dortmund.de",23X,"|"/&
    2X,"|",98X,"|"/&
    2X,100('=')/)

CONTAINS

  SUBROUTINE StrStuct()
    IMPLICIT NONE
    INTEGER i,n

    n = len(string)
    iAt = 0
    iEq = 0
    DO i=1,n
    IF (string(i:i) == '@') iAt = i
    IF (string(i:i) == '=') iEq = i
    END DO

    bOk = .false.
    IF (iAt /= 0 .and. iEq /= 0) bOk = .true.

  END SUBROUTINE StrStuct

  !-----------------------------------------------------------------------
  ! Helper subroutine to write parameter to both file and terminal
  !-----------------------------------------------------------------------
  SUBROUTINE write_param_int(mf, mt, label, value)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: mf, mt, value
    CHARACTER(*), INTENT(IN) :: label

    WRITE(mf,'(A,I10)') label, value
    WRITE(mterm,'(A,I10)') label, value
  END SUBROUTINE write_param_int

  SUBROUTINE write_param_real(mf, mt, label, value)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: mf, mt
    DOUBLE PRECISION, INTENT(IN) :: value
    CHARACTER(*), INTENT(IN) :: label

    WRITE(mf,'(A,D12.4)') label, value
    WRITE(mterm,'(A,D12.4)') label, value
  END SUBROUTINE write_param_real

  SUBROUTINE write_param_str(mf, mt, label, value)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: mf, mt
    CHARACTER(*), INTENT(IN) :: label, value

    WRITE(mf,'(A,A)') label, value
    WRITE(mterm,'(A,A)') label, value
  END SUBROUTINE write_param_str

  !-----------------------------------------------------------------------
  ! Helper function to parse Yes/No string to logical value
  !-----------------------------------------------------------------------
  FUNCTION read_yes_no_param(input_string, iEq_pos) RESULT(bool_value)
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: input_string
    INTEGER, INTENT(IN) :: iEq_pos
    LOGICAL :: bool_value
    CHARACTER(len=8) :: cParam_local

    cParam_local = " "
    READ(input_string(iEq_pos+1:),*) cParam_local
    bool_value = .false.
    IF (TRIM(ADJUSTL(cParam_local)) == "Yes") bool_value = .true.
  END FUNCTION read_yes_no_param

END SUBROUTINE GDATNEW

END MODULE param_parser
