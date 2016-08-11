MODULE PLinScalar

USE QuadScalar, ONLY: QuadSc,mgDensity,Density_Primary,&
    Density_Secondary,ST_force,DiracEps
USE def_PLinScalar
USE PP3D_MPI, ONLY:myid,master,showID,COMM_Maximum,COMM_NLComplete,&
                   Comm_Summ

IMPLICIT NONE

TYPE(lScalar) PLinLS
LOGICAL bCallReInit
REAL*8 :: dFracCrit=0.1d0
REAL*8 :: dMaxSTF


CONTAINS
!
! ----------------------------------------------
!
SUBROUTINE Init_pLinScalar()
REAL*8 ttaa,ttab


 ILEV=NLMAX
 CALL SETLEV(2)

! Initialization of the LevSet field
CALL InitField(PLinLS)

PLinLS%prm%StabScheme = 1     ! 1 - upwind,   2 - Central difference
PLinLS%prm%SrfCubF   =  2     ! 1 - 1x1 CubP, 2 - 2x2CubP, 3 - 3x3CubP
PLinLS%prm%VolCubF   =  7     ! Standard 3D cubformula concepts

PLinLS%prm%defCrit   =  1d-6
PLinLS%prm%AvgGradPhiTol = 0.05d0

IF (myid.ne.0) THEN

 ! Initialization of communicators
 CALL InitCommP1(PLinLS)

 ! Initialization of mass matrix (and its inverse)
 CALL InitMassMatrix(PLinLS)

 ! Initialization of face normals
 CALL InitFNormMatrix()

 ! Initialization element midpoints
 ALLOCATE(ElemMid(3,NEL))
 CALL GetElementMids(ElemMid,DWORK(L(LCORVG)),KWORK(L(LVERT)),NEL)

 ALLOCATE (IntPhaseElem(NEL))
 ALLOCATE (dNorm(3,NVT))
 ALLOCATE (FracFieldQ0(NEL),FracFieldQ1(NVT))

END IF

END SUBROUTINE Init_pLinScalar
!
! ----------------------------------------------
!
SUBROUTINE InitCond_PLinLS()

 ! Impose initial condition on the LevSet field
 CALL InitCond(DWORK(L(LCORVG)),KWORK(L(LVERT)))

END SUBROUTINE InitCond_PLinLS
!
! ----------------------------------------------
!
SUBROUTINE TempPLinSc()

 IF (myid.ne.0) THEN
  CALL ReinitInitField(IntPhaseElem,PLinLS%val(NLMAX)%x,NEL)
 END IF

END SUBROUTINE TempPLinSc
!
! ----------------------------------------------
!
SUBROUTINE Transport_PLinLS(mfile,iNLs)
REAL*8 RESU,DEFU,DEFU0
INTEGER inl,iRI,INLComplete,iNLs,mfile
REAL*8 tstep_old


!  GOTO 1
! -------------------------------------------------
! ------------- Levelset-advection step -----------
IF (myid.eq.showID) write(MTERM,"(30('-'),A20,30('-'))") " Levelset advection "
IF (myid.eq.showID) write(MFILE,"(30('-'),A20,30('-'))") " Levelset advection "

CALL TimeSchemeSetUp("CN")

IF (myid.ne.0) THEN

 ! Initialization of convection matrix
 CALL InitConvMatrix(1,PLinLS,QuadSc%valU,QuadSc%valV,QuadSc%valW)

 thstep = tstep*(1d0-PLinLS%prm%theta)

 PLinLS%def=0d0
 ! Set up the right hand side
 CALL PLinSc_RHS(1,PLinLS,QuadSc%valU,QuadSc%valV,QuadSc%valW,thstep, 1)

 ! Store the constant right hand side
 PLinLS%rhs = PLinLS%def

 thstep = tstep*PLinLS%prm%theta

 ! Compute the defect
 CALL PLinSc_RHS(1,PLinLS,QuadSc%valU,QuadSc%valV,QuadSc%valW,thstep,-1)

 CALL LL21 (PLinLS%def,PLinLS%ndof,DEFU0)
 DEFU0=MAX(1D-16,DEFU0)

 ! Save the old solution
 PLinLS%val_old = PLinLS%val(NLMAX)%x

END IF

CALL COMM_Maximum(DEFU0)

DO INL = 1,PLinLS%prm%NLmax
 INLComplete = 0

 IF (myid.ne.0) THEN
  ! Solve the system
  CALL PLinSc_Solve(1,PLinLS)

  CALL LL21 (PLinLS%val(NLMAX)%x,PLinLS%ndof,RESU)
  RESU=MAX(1D-16,RESU)

  ! Update the new solution
  PLinLS%val(NLMAX)%x = PLinLS%val_old + PLinLS%val(NLMAX)%x
 END IF

! Explicit time stepping ?
IF  (PLinLS%prm%NLmax.EQ.1) GOTO 1

 IF (myid.ne.0) THEN
  PLinLS%def = PLinLS%rhs

  CALL PLinSc_RHS(1,PLinLS,QuadSc%valU,QuadSc%valV,QuadSc%valW,thstep,-1)

  CALL LL21 (PLinLS%def,PLinLS%ndof,DEFU)
 END IF

 CALL COMM_Maximum(DEFU)
 CALL COMM_Maximum(RESU)

 IF (myid.eq.showID) write(MTERM,'(I2,2(1XG12.4))') INL,DEFU/DEFU0,RESU
 IF (myid.eq.showID) write(MFILE,'(I2,2(1XG12.4))') INL,DEFU/DEFU0,RESU
 IF ((DEFU/DEFU0.LE.PLinLS%prm%defCrit).AND.&
     (INL.GE.PLinLS%prm%NLmin)) INLComplete = 1
 IF  (INL.EQ.PLinLS%prm%NLmax) INLComplete = 1
 CALL COMM_NLComplete(INLComplete)
 IF (INLComplete.eq.1) GOTO 1

 IF (myid.ne.0) THEN
  ! Save the old solution
  PLinLS%val_old = PLinLS%val(NLMAX)%x
 END IF

END DO

1 CONTINUE

iNLs = iNL

END SUBROUTINE Transport_PLinLS
!
! ----------------------------------------------
!
SUBROUTINE Reinitialize_PLinLS(mfile,nSteps)
REAL*8 RESU,DEFU,DEFU0
INTEGER inl,iRI,INLComplete,nSteps,iSteps,iiSteps,mfile
REAL*8 :: tstep_old,RI_factor=0.4d0
LOGICAL bRIConv,bRedo

IF (myid.eq.showID) write(MTERM,"(26('-'),A27,27('-'))")&
" Levelset reinitialization "
IF (myid.eq.showID) write(MFILE,"(26('-'),A27,27('-'))")&
" Levelset reinitialization "

IF (myid.ne.0) THEN
 CALL ReinitCloseField(IntPhaseElem,PLinLS%val(NLMAX)%x,NEL)
 CALL ReinitFarField(IntPhaseElem,ElemMid,PLinLS%val(NLMAX)%x,NEL)
!  CALL ReinitWholeField(IntPhaseElem,PLinLS%val(NLMAX)%x,NEL)
 PLinLS%prm%SlopeLimit = .FALSE.
END IF

CALL TimeSchemeSetUp("BE")

! PLinLS%val_ad = PLinLS%val(NLMAX)%x

bRIConv = .TRUE.
iiSteps = 0

DO iSteps=1,nSteps

IF (bRIConv) THEN

 CALL Init_RI()

 bRIConv = .FALSE.
 iiSteps = iiSteps + 1

END IF

IF (myid.ne.0) THEN

 thstep = RI_factor*tstep*(1d0-PLinLS%prm%theta)

 PLinLS%def=0d0
 ! Set up the right hand side
 CALL PLinSc_RHS(2,PLinLS,dNorm(1,:),dNorm(2,:),dNorm(3,:),thstep, 1)

 ! Add the source of reinitialization
 CALL RISource(PLinLS,RI_factor*tstep)

 ! Store the constant right hand side
 PLinLS%rhs = PLinLS%def

 thstep = RI_factor*tstep*PLinLS%prm%theta

 ! Compute the defect
 CALL PLinSc_RHS(2,PLinLS,dNorm(1,:),dNorm(2,:),dNorm(3,:),thstep,-1)

 CALL LL21 (PLinLS%def,PLinLS%ndof,DEFU0)

 ! Save the old solution
 PLinLS%val_old = PLinLS%val(NLMAX)%x

END IF

CALL COMM_Maximum(DEFU0)

DO INL = 1,PLinLS%prm%NLmax
 INLComplete = 0

 IF (myid.ne.0) THEN
  ! Solve the system
  CALL PLinSc_Solve(2,PLinLS)

  CALL LL21 (PLinLS%val(NLMAX)%x,PLinLS%ndof,RESU)
  RESU=MAX(1D-16,RESU)

  ! Update the new solution
  PLinLS%val(NLMAX)%x = PLinLS%val_old + PLinLS%val(NLMAX)%x

  ! Slope limiting
  IF (PLinLS%prm%SlopeLimit) CALL SlopeLimiter()
 END IF

 CALL COMM_Maximum(RESU)

 ! Explicit time stepping ?
 IF  (PLinLS%prm%NLmax.EQ.1) THEN
  IF (myid.eq.showID) write(MTERM,'(I4,G12.4)') iSteps,RESU
  IF (myid.eq.showID) write(MFILE,'(I4,G12.4)') iSteps,RESU
  GOTO 1
 END IF

 IF (myid.ne.0) THEN
  PLinLS%def = PLinLS%rhs

  CALL PLinSc_RHS(2,PLinLS,dNorm(1,:),dNorm(2,:),dNorm(3,:),thstep,-1)

  CALL LL21 (PLinLS%def,PLinLS%ndof,DEFU)
 END IF

 CALL COMM_Maximum(DEFU)

 IF (myid.eq.showID) write(MTERM,'(I2,2(1XG12.4))') INL,DEFU/DEFU0,RESU
 IF (myid.eq.showID) write(MFILE,'(I2,2(1XG12.4))') INL,DEFU/DEFU0,RESU
 IF ((DEFU/DEFU0.LE.PLinLS%prm%defCrit).AND.&
     (INL.GE.PLinLS%prm%NLmin)) INLComplete = 1
 IF  (INL.EQ.PLinLS%prm%NLmax) INLComplete = 1
 CALL COMM_NLComplete(INLComplete)
 IF (INLComplete.eq.1) GOTO 1

 IF (myid.ne.0) THEN
  ! Save the old solution
  PLinLS%val_old = PLinLS%val(NLMAX)%x
 END IF

END DO

1 CONTINUE

IF (RESU.LT.0.025d0) GOTO 2
! IF (iSteps.EQ.40) bRIConv = .TRUE.
! CALL CheckGradPhiQuality(bRedo)
! IF (iSteps.LT.0.8*nSteps.AND.iSteps.GT.40) bRIConv = .TRUE.
! IF (iiSteps.EQ.2.AND.iSteps.GT.40) GOTO 2

END DO

2 CONTINUE

nSteps = iSteps

! PLinLS%val(NLMAX)%x = PLinLS%val(NLMAX)%x - PLinLS%val_ad

END SUBROUTINE Reinitialize_PLinLS
!
! ----------------------------------------------
!
SUBROUTINE InitCond(dcorvg,kvert)
INTEGER kvert(8,*)
REAL*8 dcorvg(3,*)
REAL*8 PX,PY,PZ
INTEGER i,j,k

DO i=1,nel
 PX = 0d0
 PY = 0d0
 PZ = 0d0
 DO j=1,8
  PX = PX + 0.125d0*(dcorvg(1,kvert(j,i)))
  PY = PY + 0.125d0*(dcorvg(2,kvert(j,i)))
  PZ = PZ + 0.125d0*(dcorvg(3,kvert(j,i)))
 END DO

 k=4*(i-1)
 CALL GetInitVal(PX,PY,PZ,PLinLS%val(NLMAX)%x(k+1:k+4))

END DO

END SUBROUTINE InitCond
!
! ----------------------------------------------
!
SUBROUTINE Reinit_Interphase(mfile)
INTEGER mfile,iRe
EXTERNAL E011,E012

 IF (myid.eq.showID) write(MTERM,"(25('-'),A29,26('-'))")&
 " Interphase reinitialization "
 IF (myid.eq.showID) write(MFILE,"(25('-'),A29,26('-'))")&
 " Interphase reinitialization "
DO iRe=1,5
 IF (myid.ne.0) THEN

 ILEV=NLMAX
 CALL SETLEV(2)

 CALL PLinScP1toQ1(PLinLS)
 PLinLS%def = 0d0
 CALL GetIntQ1DistToP1(PLinLS%Q1,IntPhaseElem,&
      PLinLS%def,ElemMid,KWORK(L(LVERT)),KWORK(L(LAREA)),&
      KWORK(L(LEDGE)),VWORK(L(KLINT(ILEV))),DWORK(L(LCORVG)),&
      PLinLS%prm%VolCubF,E012,E011)
 ! Solve the system
 CALL PLinSc_Solve(3,PLinLS)

 END IF
END DO

END SUBROUTINE Reinit_Interphase


SUBROUTINE UpdateAuxVariables(mfile)
INTEGER mfile
REAL*8 dMass,dMass0,dV1,dV2
REAL*8 AvgGradPhi
LOGICAL bNEW
DATA bNEW/.TRUE./
SAVE dMass0
EXTERNAL E011,E012


IF (myid.ne.0) THEN
 ILEV=NLMAX
 CALL SETLEV(2)
 CALL IntPhaseFinder(IntPhaseElem,PLinLS%val(NLMAX)%x,ElemMid,&
      DWORK(L(LCORVG)),KWORK(L(LVERT)),KWORK(L(LAREA)),&
      KWORK(L(LADJ)),KWORK(L(LVEL)),NVEL,PLinLS%iParFace,NEL)
 CALL ReinitFarField(IntPhaseElem,ElemMid,PLinLS%val(NLMAX)%x,NEL)
END IF

CALL GetAverageGradPhi(PLinLS%val(NLMAX)%x,NEL,&
     IntPhaseElem,AvgGradPhi)

IF (myid.eq.showID) write(*,*) "Average Grad(Phi) = ",AvgGradPhi
IF (AvgGradPhi.GT.PLinLS%prm%AvgGradPhiTol) THEN
 bCallReInit=.TRUE.
ELSE
 bCallReInit=.FALSE.
END IF

IF (myid.ne.0) THEN

 ILEV=NLMAX
 CALL SETLEV(2)
 FracFieldQ0 = 0d0
 CALL GetP1Density(mgDensity(NLMAX)%x,IntPhaseElem,PLinLS%val(NLMAX)%x,&
      ElemMid,FracFieldQ0,DWORK(L(LCORVG)),KWORK(L(LVERT)),&
      VWORK(L(LVOL)),Density_Primary,Density_Secondary,dV1,dV2,NEL)

 CALL IntegrateMass(dMAss,mgDensity(NLMAX)%x,VWORK(L(KLVOL(ILEV))),NEL)

END IF

CALL Comm_Summ(dMAss)
CALL Comm_Summ(dV1)
CALL Comm_Summ(dV2)

IF (myid.eq.showID) THEN
 IF (ISTART.EQ.0.AND.bNEW) THEN
  OPEN(987,FILE='#data/mass.txt')
  dMass0 = dMass
  bNEW=.false.
 ELSE
  OPEN(987,FILE='#data/mass.txt',ACCESS='APPEND')
 END IF

 WRITE(987,'(4G13.5)') timens,1d2*dMAss/dMass0,dV1,dV2
 CLOSE(987)
END IF

IF (myid.ne.0) THEN
 CALL IntPolNormals(dNorm,PLinLS%val(NLMAX)%x,IntPhaseElem,&
      KWORK(L(LVERT)),DWORK(L(LCORVG)),VWORK(L(KLVOL(ILEV))),&
      ElemMid,NEL,NVT)

 FracFieldQ1 = 0d0
 PLinLS%aux  = 0d0
 CALL IntQ0toQ1(FracFieldQ0,FracFieldQ1,PLinLS%aux,KWORK(L(LVERT)),&
      DWORK(L(LCORVG)),VWORK(L(KLVOL(ILEV))),NEL,NVT)

END IF

CALL CorrectIntPhaseI(FracFieldQ1,IntPhaseElem,KWORK(L(LVERT)),&
     ElemMid,NEL,dFracCrit)

IF (myid.ne.0) THEN
 CALL GetSurfTensForce(PLinLS%val(NLMAX)%x,dNorm,IntPhaseElem,&
      ST_force,ElemMid,KWORK(L(LVERT)),KWORK(L(LAREA)),&
      KWORK(L(LEDGE)),VWORK(L(KLINT(ILEV))),DWORK(L(LCORVG)),&
      DiracEps,PLinLS%prm%VolCubF,E012,E011)
 CALL GetMaxSurfTensForce(ST_force,dMaxSTF)
END IF

CALL COMM_Maximum(dMaxSTF)

END SUBROUTINE UpdateAuxVariables

!
! ----------------------------------------------
!
SUBROUTINE OutputInterphase(iO)
INTEGER iO

ILEV=NLMAX
CALL SETLEV(2)
! CALL IntPhaseFinder(IntPhaseElem,PLinLS%val(NLMAX)%x,ElemMid,&
!      DWORK(L(LCORVG)),KWORK(L(LVERT)),KWORK(L(LADJ)),NEL)

CALL GetSurf(IntPhaseElem,PLinLS%val(NLMAX)%x,ElemMid,&
     DWORK(L(LCORVG)),KWORK(L(LVERT)),KWORK(L(LADJ)),iO,NEL)

END SUBROUTINE OutputInterphase
!
! ----------------------------------------------
!
SUBROUTINE CheckGradPhiQuality(bRedo)
LOGICAL bRedo
REAL*8  DGRAD

IF (myid.NE.0) THEN
 CALL CheckGradPhiQualitySub(PLinLS%val(NLMAX)%x,&
      IntPhaseElem,NEL,DGRAD)
END IF

CALL Comm_Summ(DGRAD)
IF (DGRAD.EQ.0d0) THEN
 bRedo = .TRUE.
ELSE
 bRedo = .TRUE.
END IF

END SUBROUTINE CheckGradPhiQuality
!
! ----------------------------------------------
!
SUBROUTINE Init_RI()

 ! Initialization of the normals and the convection matrix
IF (myid.ne.0) THEN
 ILEV=NLMAX
 CALL SETLEV(2)
 CALL IntPolNormals(dNorm,PLinLS%val(NLMAX)%x,IntPhaseElem,&
      KWORK(L(LVERT)),DWORK(L(LCORVG)),VWORK(L(KLVOL(ILEV))),&
      ElemMid,NEL,NVT)
 CALL InitConvMatrix(2,PLinLS,dNorm(1,:),dNorm(2,:),dNorm(3,:))

END IF

END SUBROUTINE Init_RI
!
! ----------------------------------------------
!
SUBROUTINE SlopeLimiter()

ILEV=NLMAX
CALL SETLEV(2)
CALL SlopeLimiterSub(PLinLS%val(NLMAX)%x,ElemMid,&
     DWORK(L(LCORVG)),KWORK(L(LVERT)),NEL,NVT)

END SUBROUTINE SlopeLimiter
!
! ----------------------------------------------
!
SUBROUTINE OuterRegion
INTEGER IEL,IP

DO IEL = 1,NEL
 IF (ABS(IntPhaseElem(IEL)).EQ.100) THEN
  IP=4*(IEL-1)
  IF (IntPhaseElem(IEL).GT.0) PLinLS%val(NLMAX)%x(IP+1) = 1d0
  IF (IntPhaseElem(IEL).LT.0) PLinLS%val(NLMAX)%x(IP+1) =-1d0
  PLinLS%val(NLMAX)%x(IP+2) = 0d0
  PLinLS%val(NLMAX)%x(IP+3) = 0d0
  PLinLS%val(NLMAX)%x(IP+4) = 0d0
 END IF
END DO

END SUBROUTINE OuterRegion
!
! ----------------------------------------------
!
SUBROUTINE TimeSchemeSetUp(carg)
CHARACTER*2 carg

IF     (carg.eq."BE") THEN
 PLinLS%prm%NLmin  = 1
 PLinLS%prm%NLmax  = 1
 PLinLS%prm%theta  = 0.0d0
ELSEIF (carg.eq."CN") THEN
 PLinLS%prm%NLmin  = 1
 PLinLS%prm%NLmax  = 20
 PLinLS%prm%theta  = 0.5d0
ELSEIF (carg.eq."FE") THEN
 PLinLS%prm%NLmin  = 1
 PLinLS%prm%NLmax  = 20
 PLinLS%prm%theta  = 1.0d0
ELSE
 WRITE(*,*) "Not known time stepping scheme ... ",carg
 STOP
END IF

END SUBROUTINE TimeSchemeSetUp
!
! ----------------------------------------------
!
SUBROUTINE GetMaxSurfTensForce(ST,dST)
REAL*8 ST(*),dST
INTEGER ivt

dST = 0d0
DO ivt=1,nvt
 IF (dST.lt.ST(ivt)) dST = ST(ivt)
END DO

END SUBROUTINE GetMaxSurfTensForce

END MODULE PLinScalar
