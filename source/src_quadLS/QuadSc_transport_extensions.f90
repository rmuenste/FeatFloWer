SUBROUTINE Transport_q2p1_UxyzP_ParT(mfile,inl_u,itns)
use cinterface, only: calculateDynamics,calculateFBM
use fbm, only: fbm_updateFBM

INTEGER mfile,INL,inl_u,itns
REAL*8  ResU,ResV,ResW,DefUVW,RhsUVW,DefUVWCrit
REAL*8  ResP,DefP,RhsPG,defPG,defDivU,DefPCrit
INTEGER INLComplete,I,J,IERR

INTEGER :: iOuter,iStep,nSteps,lStep,iLoop
REAL*8  :: tstep_BU,MaxInitialPressureDefect,MaxInitialPressureDefect0
character*(2) :: cEnding(0:9) = ['th','st','nd','rd','th','th','th','th','th','th']
integer       :: iEnding
logical :: bExit

TYPE tIntPolSchem
 Integer  :: Actual
 Integer  :: IntPol1,IntPol2
 Integer  :: Src1,Src2,iX1Src,iX2Src
 REAL*8   :: X1Src,X2Src,X1Dest,X2Dest
 REAL*8   :: Frac1,Frac2
END TYPE tIntPolSchem
TYPE (tIntPolSchem) IPS

TYPE tSingleSol
 REAL*8, allocatable :: U(:),V(:),W(:),P(:)
 REAL*8, allocatable :: U_aux(:),V_aux(:),W_aux(:),P_aux(:)
END TYPE tSingleSol

TYPE tSolSeq
 integer :: nOuter,nFOuter
 integer :: nSteps
 TYPE(tSingleSol), allocatable :: S(:)
END TYPE tSolSeq

TYPE(tSolSeq) :: mySolSeq

INTEGER :: LinIntPol=1,ConstIntPol=0


mySolSeq%nOuter = Properties%nTPSubSteps
mySolSeq%nFOuter = Properties%nTPFSubSteps

tstep_BU = tstep

mySolSeq%nSteps = 2**(mySolSeq%nOuter-1)

if (.not.allocated(mySolSeq%S)) THEN
 ALLOCATE(mySolSeq%S(0:mySolSeq%nSteps))
end if

DO iStep = 0,mySolSeq%nSteps
 if (.not.allocated(mySolSeq%S(iStep)%U)) ALLOCATE(mySolSeq%S(iStep)%U(QuadSc%ndof))
 if (.not.allocated(mySolSeq%S(iStep)%V)) ALLOCATE(mySolSeq%S(iStep)%V(QuadSc%ndof))
 if (.not.allocated(mySolSeq%S(iStep)%W)) ALLOCATE(mySolSeq%S(iStep)%W(QuadSc%ndof))
 if (.not.allocated(mySolSeq%S(iStep)%P)) ALLOCATE(mySolSeq%S(iStep)%P(LinSc%ndof))
 
 if (.not.allocated(mySolSeq%S(iStep)%U_aux)) ALLOCATE(mySolSeq%S(iStep)%U_aux(QuadSc%ndof))
 if (.not.allocated(mySolSeq%S(iStep)%V_aux)) ALLOCATE(mySolSeq%S(iStep)%V_aux(QuadSc%ndof))
 if (.not.allocated(mySolSeq%S(iStep)%W_aux)) ALLOCATE(mySolSeq%S(iStep)%W_aux(QuadSc%ndof))
 if (.not.allocated(mySolSeq%S(iStep)%P_aux)) ALLOCATE(mySolSeq%S(iStep)%P_aux(LinSc%ndof))
 mySolSeq%S(iStep)%U = QuadSc%ValU
 mySolSeq%S(iStep)%V = QuadSc%ValV
 mySolSeq%S(iStep)%W = QuadSc%ValW
!  mySolSeq%S(iStep)%U = 0d0
!  mySolSeq%S(iStep)%V = 0d0
!  mySolSeq%S(iStep)%W = 0d0
 mySolSeq%S(iStep)%P = LinSc%ValP(NLMAX)%x
END DO

! if (.not.(istart.eq.0.and.itns.eq.1)) then
!  mySolSeq%S(0)%U = QuadSc%ValU
!  mySolSeq%S(0)%V = QuadSc%ValV
!  mySolSeq%S(0)%W = QuadSc%ValW
! end if

!------------------------------------------------------ PREDICTION --------------------------------------------------------
MaxInitialPressureDefect = 1e-30 
iLoop =1

DO iOuter=1,mySolSeq%nOuter

 nSteps = 2**(iOuter-1)
 lStep  = mySolSeq%nSteps/nSteps
 tstep  = tstep_BU/DBLE(nSteps)
 
IF (iOuter.le.1) THEN
!  IF (iOuter.le.mySolSeq%nOuter) THEN
!  IF (iOuter.le.mySolSeq%nOuter-2) THEN
 
  CALL NonLin_BurgerStep_ParT()

  CALL PressureStep_ParT()
  
  ! Velocity update for last pressure
  CALL NonLin_BurgerStep_ParT()
  
 END IF

 IF (iOuter.ne.mySolSeq%nOuter) then
  CALL InterpolatePressure_ParT(LinIntPol)
  CALL InterpolateVelocity_ParT()
 END IF

END DO !iOuter
tstep  = tstep_BU

MaxInitialPressureDefect0 = MaxInitialPressureDefect
IF (myid.eq.1) WRITE(*,*) 'rock and roll! ',MaxInitialPressureDefect0
!------------------------------------------------------ PREDICTION --------------------------------------------------------


bEXIT = .false.
do iLoop =1 , Properties%nTPIterations

 MaxInitialPressureDefect = 0d0

 DO iOuter=mySolSeq%nFOuter,mySolSeq%nFOuter

  nSteps = 2**(iOuter-1)
  lStep  = mySolSeq%nSteps/nSteps
  tstep  = tstep_BU/DBLE(nSteps)

!   CALL Lin_BurgerStep_ParT()
   CALL NonLin_BurgerStep_ParT()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  QuadSc%ValU = mySolSeq%S(mySolSeq%nSteps)%U
  QuadSc%ValV = mySolSeq%S(mySolSeq%nSteps)%V
  QuadSc%ValW = mySolSeq%S(mySolSeq%nSteps)%W

  if (bExit) THEN
   tstep  = tstep_BU
   GOTO 100
  END IF
  
  LinSc%P_new = 1.5d0*mySolSeq%S(mySolSeq%nSteps-0)%P  - 0.5d0*mySolSeq%S(mySolSeq%nSteps-1)%P
  CALL QuadScP1ExtPoltoQ2(LinSc,QuadSc)
  CALL FAC_GetForcesParT(mfile,iLoop)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  CALL PressureStep_ParT()

 END DO !iOuter

 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! QuadSc%ValU = mySolSeq%S(mySolSeq%nSteps)%U
! QuadSc%ValV = mySolSeq%S(mySolSeq%nSteps)%V
! QuadSc%ValW = mySolSeq%S(mySolSeq%nSteps)%W
! 
! LinSc%P_new = 1.5d0*mySolSeq%S(mySolSeq%nSteps-0)%P  - 0.5d0*mySolSeq%S(mySolSeq%nSteps-1)%P
! !LinSc%P_new = mySolSeq%S(mySolSeq%nSteps-0)%P
! CALL QuadScP1ExtPoltoQ2(LinSc,QuadSc)
! 
! CALL FAC_GetForcesParT(mfile,iLoop)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  if (iLoop.eq.1) MaxInitialPressureDefect0 = MaxInitialPressureDefect

 tstep  = tstep_BU

 iEnding = mod(iLoop,10)
 if (myid.eq.1) write(mterm,'(A,I0,A,3ES12.4)') 'PressureDefectReductionIn ', iLoop,cEnding(iEnding)//' step:',MaxInitialPressureDefect0,MaxInitialPressureDefect,MaxInitialPressureDefect/MaxInitialPressureDefect0
 if (myid.eq.1) write(mfile,'(A,I0,A,3ES12.4)') 'PressureDefectReductionIn ', iLoop,cEnding(iEnding)//' step:',MaxInitialPressureDefect0,MaxInitialPressureDefect,MaxInitialPressureDefect/MaxInitialPressureDefect0
!  if (MaxInitialPressureDefect.lt.1e-7) THEN 
! if (MaxInitialPressureDefect/MaxInitialPressureDefect0.lt.Properties%DiracEps) THEN 
  if (MaxInitialPressureDefect/MaxInitialPressureDefect0.lt.Properties%DiracEps.and.MaxInitialPressureDefect.lt.1e-11) THEN 
  if (myid.eq.1) write(mterm,'(A,I0,A,3ES12.4)') 'ExitingCoarseTimeCPloopIn ', iLoop,cEnding(iEnding)//' step!|Init&FinPresDefect: ',MaxInitialPressureDefect0,MaxInitialPressureDefect,MaxInitialPressureDefect/MaxInitialPressureDefect0
  if (myid.eq.1) write(mfile,'(A,I0,A,3ES12.4)') 'ExitingCoarseTimeCPloopIn ', iLoop,cEnding(iEnding)//' step!|Init&FinPresDefect: ',MaxInitialPressureDefect0,MaxInitialPressureDefect,MaxInitialPressureDefect/MaxInitialPressureDefect0
  bEXIT = .true.
 end if
! if (bExit) GOTO 100

end do !iLoop

100 continue

if (mySolSeq%nFOuter.eq.mySolSeq%nOuter) goto 5

!!!! Extrapolation !!!!!!!
do iLoop =1 , 1 

 DO iOuter=mySolSeq%nFOuter,mySolSeq%nOuter

  nSteps = 2**(iOuter-1)
  lStep  = mySolSeq%nSteps/nSteps
  tstep  = tstep_BU/DBLE(nSteps)

  CALL NonLin_BurgerStep_ParT()

  CALL PressureStep_ParT()

 END DO !iOuter

 tstep  = tstep_BU

end do !iLoop

5 continue

QuadSc%ValU = mySolSeq%S(mySolSeq%nSteps)%U
QuadSc%ValV = mySolSeq%S(mySolSeq%nSteps)%V
QuadSc%ValW = mySolSeq%S(mySolSeq%nSteps)%W

LinSc%P_new = 1.5d0*mySolSeq%S(mySolSeq%nSteps-0)%P  - 0.5d0*mySolSeq%S(mySolSeq%nSteps-1)%P
!LinSc%P_new = mySolSeq%S(mySolSeq%nSteps-0)%P
CALL QuadScP1ExtPoltoQ2(LinSc,QuadSc)

CALL FAC_GetForces(mfile)

RETURN

 CONTAINS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE InterpolateVelocity_ParT()

  IF (myid.ne.master) THEN
  
    DO iStep = 1,nSteps
      
      IPS%Actual = (iStep-0)*lStep
      IPS%IntPol1 = (iStep-0)*lStep
      IPS%IntPol2 = (iStep-0)*lStep-int(lStep/2)
      
      IPS%Src1 = (iStep-1)*lStep
      IPS%Src2 = (iStep-0)*lStep
      
!      mySolSeq%S(IPS%IntPol1)%U = mySolSeq%S(IPS%Actual)%U
!      mySolSeq%S(IPS%IntPol1)%V = mySolSeq%S(IPS%Actual)%V
!      mySolSeq%S(IPS%IntPol1)%w = mySolSeq%S(IPS%Actual)%w
      mySolSeq%S(IPS%IntPol2)%U = 0.5d0*mySolSeq%S(IPS%Src1)%U + 0.5d0*mySolSeq%S(IPS%Src2)%U
      mySolSeq%S(IPS%IntPol2)%V = 0.5d0*mySolSeq%S(IPS%Src1)%V + 0.5d0*mySolSeq%S(IPS%Src2)%V
      mySolSeq%S(IPS%IntPol2)%W = 0.5d0*mySolSeq%S(IPS%Src1)%W + 0.5d0*mySolSeq%S(IPS%Src2)%w

      if (myid.eq.showid) WRITE(*,'(A,I0,4(A,I0))') &
      "ActualStep: ",IPS%Actual, ' IntPolSourceSteps: ',IPS%Src1,',',IPS%Src2,' UpdatedSteps: ',IPS%IntPol1,',',IPS%IntPol2
     
    END DO
  END IF
    
END SUBROUTINE InterpolateVelocity_ParT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE InterpolatePressure_ParT(IntPol)
INTEGER IntPol

  IF (myid.ne.master) THEN
  
    DO iStep = 0,nSteps
     mySolSeq%S((iStep-0)*lStep)%P_aux = mySolSeq%S((iStep-0)*lStep)%P
    END DO
  
    DO iStep = 1,nSteps
     IF (IntPol.eq.1) then
      IPS%Actual = (iStep-0)*lStep
      IPS%IntPol1 = (iStep-0)*lStep
      IPS%IntPol2 = (iStep-0)*lStep-int(lStep/2)
      IPS%Src1 = (iStep-1)*lStep
      IPS%Src2 = (iStep-0)*lStep
      IPS%iX1Src = (iStep-0)*lStep-int(lStep/2) - lStep
      IPS%iX2Src = (iStep-0)*lStep-int(lStep/2) 
      IPS%X1Src = max(-1d0/(2d0*dble(mySolSeq%nSteps)),dble(IPS%iX1Src)/dble(mySolSeq%nSteps))
      IPS%X2Src = dble(IPS%iX2Src)/dble(mySolSeq%nSteps)
      IPS%X1Dest = IPS%X2Src + dble(lStep/4d0)/dble(mySolSeq%nSteps)
      IPS%X2Dest = IPS%X2Src - dble(lStep/4d0)/dble(mySolSeq%nSteps)
     
      IPS%Frac1 = (IPS%X1Dest - IPS%X2Src )/(IPS%X1Src - IPS%X2Src)
      IPS%Frac2 = (IPS%X1Src  - IPS%X1Dest)/(IPS%X1Src - IPS%X2Src)
      mySolSeq%S(IPS%IntPol1)%P = IPS%Frac1*mySolSeq%S(IPS%Src1)%P_aux + IPS%Frac2*mySolSeq%S(IPS%Src2)%P_aux  
       
      IPS%Frac1 = (IPS%X2Dest - IPS%X2Src )/(IPS%X1Src - IPS%X2Src)
      IPS%Frac2 = (IPS%X1Src  - IPS%X2Dest)/(IPS%X1Src - IPS%X2Src)
      mySolSeq%S(IPS%IntPol2)%P = IPS%Frac1*mySolSeq%S(IPS%Src1)%P_aux + IPS%Frac2*mySolSeq%S(IPS%Src2)%P_aux  

      if (myid.eq.showid) WRITE(*,'(A,I0,4(A,I0))') &
      "ActualStep: ",IPS%Actual, ' IntPolSourceSteps: ',IPS%Src1,',',IPS%Src2,' UpdatedSteps: ',IPS%IntPol1,',',IPS%IntPol2
      if (myid.eq.showid) WRITE(*,'(2(A,I0),2(A,F6.3))') 'IntFractions:', IPS%iX1Src, ', ', IPS%iX2Src,' SrcFractions:', IPS%X1Src, ', ', IPS%X2Src
      if (myid.eq.showid) WRITE(*,'(2(A,I0),2(A,F6.3))') 'IntFractions:', IPS%IntPol1, ', ', IPS%IntPol2,' DestFractions:', IPS%X1Dest, ', ', IPS%X2Dest
     END IF

     IF (IntPol.eq.0) then
      if (myid.eq.showid) WRITE(*,*) "steps: ",(iStep-0)*lStep,(iStep-1)*lStep,(iStep-0)*lStep-int(lStep/2)
      mySolSeq%S((iStep-0)*lStep-int(lStep/2))%P = mySolSeq%S((iStep-0)*lStep)%P 
     END IF
     
    END DO
  END IF
    
END SUBROUTINE InterpolatePressure_ParT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE NonLin_BurgerStep_ParT()

DO iStep = 1,nSteps

thstep = tstep*(1d0-theta)

QuadSc%ValU = mySolSeq%S((iStep-1)*lStep)%U
QuadSc%ValV = mySolSeq%S((iStep-1)*lStep)%V
QuadSc%ValW = mySolSeq%S((iStep-1)*lStep)%W

! Set dirichlet boundary conditions on the solution
CALL Boundary_QuadScalar_Val()

LinSc%ValP(NLMAX)%x = mySolSeq%S((iStep-0)*lStep)%P

CALL OperatorRegenaration(2)

CALL OperatorRegenaration(3)

! -------------------------------------------------
! Compute the momentum equations
! -------------------------------------------------
! GOTO 1

IF (myid.ne.master) THEN

 CALL ZTIME(tttt0)

 ! Assemble the right hand side
 CALL Matdef_General_QuadScalar(QuadSc,1)

 ! Add the pressure gradient to the rhs
 CALL AddPressureGradient()
END IF

 ! Add the viscoelastic stress to the rhs
 IF(bViscoElastic)THEN
!    CALL AddViscoStress()
 END IF

IF (myid.ne.master) THEN
 ! Add the gravity force to the rhs
!  CALL AddGravForce()

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 ! Store the constant right hand side
 QuadSc%rhsU = QuadSc%defU
 QuadSc%rhsV = QuadSc%defV
 QuadSc%rhsW = QuadSc%defW

! Set dirichlet boundary conditions on the solution
 CALL Boundary_QuadScalar_Val()

END IF

thstep = tstep*theta

IF (myid.ne.master) THEN
 QuadSc%ValU = mySolSeq%S((iStep-0)*lStep)%U
 QuadSc%ValV = mySolSeq%S((iStep-0)*lStep)%V
 QuadSc%ValW = mySolSeq%S((iStep-0)*lStep)%W
END IF

! Set dirichlet boundary conditions on the solution
CALL Boundary_QuadScalar_Val()

CALL OperatorRegenaration(3)

IF (myid.ne.master) THEN

 ! Assemble the defect vector and fine level matrix
 CALL Matdef_General_QuadScalar(QuadSc,-1)

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 QuadSc%auxU = QuadSc%defU
 QuadSc%auxV = QuadSc%defV
 QuadSc%auxW = QuadSc%defW
 CALL E013Sum3(QuadSc%auxU,QuadSc%auxV,QuadSc%auxW)

 ! Save the old solution
 CALL LCP1(QuadSc%valU,QuadSc%valU_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valV,QuadSc%valV_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valW,QuadSc%valW_old,QuadSc%ndof)

 ! Compute the norm of the defect
 CALL Resdfk_General_QuadScalar(QuadSc,ResU,ResV,ResW,DefUVW,RhsUVW)

END IF

CALL COMM_Maximum(RhsUVW)
DefUVWCrit=MAX(RhsUVW*QuadSc%prm%defCrit,QuadSc%prm%MinDef)

CALL Protocol_QuadScalar(mfile,QuadSc,0,&
     ResU,ResV,ResW,DefUVW,DefUVWCrit," Momentum equation ")

CALL ZTIME(tttt1)
myStat%tDefUVW = myStat%tDefUVW + (tttt1-tttt0)


DO INL=1,QuadSc%prm%NLmax
INLComplete = 0

! ! Calling the solver
CALL Solve_General_QuadScalar(QuadSc,Boundary_QuadScalar_Val,&
Boundary_QuadScalar_Mat,Boundary_QuadScalar_Mat_9,mfile)

!!!!          Checking the quality of the result           !!!!
!!!! ----------------------------------------------------- !!!!

CALL OperatorRegenaration(3)

IF (myid.ne.master) THEN
! Restore the constant right hand side
 CALL ZTIME(tttt0)
 QuadSc%defU = QuadSc%rhsU
 QuadSc%defV = QuadSc%rhsV
 QuadSc%defW = QuadSc%rhsW
END IF

IF (myid.ne.master) THEN

 ! Assemble the defect vector and fine level matrix
 CALL Matdef_General_QuadScalar(QuadSc,-1)

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 QuadSc%auxU = QuadSc%defU
 QuadSc%auxV = QuadSc%defV
 QuadSc%auxW = QuadSc%defW
 CALL E013Sum3(QuadSc%auxU,QuadSc%auxV,QuadSc%auxW)

 ! Save the old solution
 CALL LCP1(QuadSc%valU,QuadSc%valU_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valV,QuadSc%valV_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valW,QuadSc%valW_old,QuadSc%ndof)

 ! Compute the defect
 CALL Resdfk_General_QuadScalar(QuadSc,ResU,ResV,ResW,DefUVW,RhsUVW)

END IF

! Checking convergence rates against criterions
RhsUVW=DefUVW
CALL COMM_Maximum(RhsUVW)
CALL Protocol_QuadScalar(mfile,QuadSc,INL,&
     ResU,ResV,ResW,DefUVW,RhsUVW)
IF (ISNAN(RhsUVW)) stop

IF ((DefUVW.LE.DefUVWCrit).AND.&
    (INL.GE.QuadSc%prm%NLmin)) INLComplete = 1

CALL COMM_NLComplete(INLComplete)
CALL ZTIME(tttt1)
myStat%tDefUVW = myStat%tDefUVW + (tttt1-tttt0)

IF (INLComplete.eq.1) GOTO 1

END DO

1 CONTINUE

! return
myStat%iNonLin = myStat%iNonLin + INL
inl_u = INL

IF (myid.ne.master) THEN
 mySolSeq%S((iStep-0)*lStep)%U = QuadSc%ValU
 mySolSeq%S((iStep-0)*lStep)%V = QuadSc%ValV
 mySolSeq%S((iStep-0)*lStep)%W = QuadSc%ValW
END IF

END DO ! iStep

END SUBROUTINE NonLin_BurgerStep_ParT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE Lin_BurgerStep_ParT()
integer iNonLin

! do iNonLin = 1,2

DO iStep = 0,nSteps
 mySolSeq%S(iStep)%U_aux = mySolSeq%S(iStep)%U
 mySolSeq%S(iStep)%V_aux = mySolSeq%S(iStep)%V
 mySolSeq%S(iStep)%W_aux = mySolSeq%S(iStep)%W
end do

DO iStep = 1,nSteps

thstep = tstep*(1d0-theta)

QuadSc%ValU = mySolSeq%S((iStep-1)*lStep)%U_aux
QuadSc%ValV = mySolSeq%S((iStep-1)*lStep)%V_aux
QuadSc%ValW = mySolSeq%S((iStep-1)*lStep)%W_aux

! Set dirichlet boundary conditions on the solution
CALL Boundary_QuadScalar_Val()

LinSc%ValP(NLMAX)%x = mySolSeq%S((iStep-0)*lStep)%P

CALL OperatorRegenaration(2)

CALL OperatorRegenaration(3)

! -------------------------------------------------
! Compute the momentum equations
! -------------------------------------------------
! GOTO 1

IF (myid.ne.master) THEN

 CALL ZTIME(tttt0)

 ! Assemble the right hand side
 CALL Matdef_General_QuadScalar(QuadSc,1)

 ! Add the pressure gradient to the rhs
 CALL AddPressureGradient()
END IF

 ! Add the viscoelastic stress to the rhs
 IF(bViscoElastic)THEN
!    CALL AddViscoStress()
 END IF

IF (myid.ne.master) THEN
 ! Add the gravity force to the rhs
!  CALL AddGravForce()

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 ! Store the constant right hand side
 QuadSc%rhsU = QuadSc%defU
 QuadSc%rhsV = QuadSc%defV
 QuadSc%rhsW = QuadSc%defW

! Set dirichlet boundary conditions on the solution
 CALL Boundary_QuadScalar_Val()

END IF

thstep = tstep*theta

IF (myid.ne.master) THEN
 QuadSc%ValU = mySolSeq%S((iStep-0)*lStep)%U_aux
 QuadSc%ValV = mySolSeq%S((iStep-0)*lStep)%V_aux
 QuadSc%ValW = mySolSeq%S((iStep-0)*lStep)%W_aux
END IF

! Set dirichlet boundary conditions on the solution
CALL Boundary_QuadScalar_Val()

CALL OperatorRegenaration(3)

IF (myid.ne.master) THEN

 ! Assemble the defect vector and fine level matrix
 CALL Matdef_General_QuadScalar(QuadSc,-1)

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 QuadSc%auxU = QuadSc%defU
 QuadSc%auxV = QuadSc%defV
 QuadSc%auxW = QuadSc%defW
 CALL E013Sum3(QuadSc%auxU,QuadSc%auxV,QuadSc%auxW)

 ! Save the old solution
 CALL LCP1(QuadSc%valU,QuadSc%valU_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valV,QuadSc%valV_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valW,QuadSc%valW_old,QuadSc%ndof)

 ! Compute the norm of the defect
 CALL Resdfk_General_QuadScalar(QuadSc,ResU,ResV,ResW,DefUVW,RhsUVW)

END IF

CALL COMM_Maximum(RhsUVW)
DefUVWCrit=MAX(RhsUVW*QuadSc%prm%defCrit,QuadSc%prm%MinDef)

CALL Protocol_QuadScalar(mfile,QuadSc,0,&
     ResU,ResV,ResW,DefUVW,DefUVWCrit," Momentum equation ")

CALL ZTIME(tttt1)
myStat%tDefUVW = myStat%tDefUVW + (tttt1-tttt0)

! 
! DO INL=1,QuadSc%prm%NLmax
! INLComplete = 0

! ! Calling the solver
CALL Solve_General_QuadScalar(QuadSc,Boundary_QuadScalar_Val,&
Boundary_QuadScalar_Mat,Boundary_QuadScalar_Mat_9,mfile)

!!!!          Checking the quality of the result           !!!!
!!!! ----------------------------------------------------- !!!!

! CALL OperatorRegenaration(3)
! 
IF (myid.ne.master) THEN
! Restore the constant right hand side
 CALL ZTIME(tttt0)
 QuadSc%defU = QuadSc%rhsU
 QuadSc%defV = QuadSc%rhsV
 QuadSc%defW = QuadSc%rhsW
END IF
! 
IF (myid.ne.master) THEN

 ! Assemble the defect vector and fine level matrix
 CALL Matdef_General_QuadScalar(QuadSc,-1)

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 QuadSc%auxU = QuadSc%defU
 QuadSc%auxV = QuadSc%defV
 QuadSc%auxW = QuadSc%defW
 CALL E013Sum3(QuadSc%auxU,QuadSc%auxV,QuadSc%auxW)

 ! Save the old solution
 CALL LCP1(QuadSc%valU,QuadSc%valU_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valV,QuadSc%valV_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valW,QuadSc%valW_old,QuadSc%ndof)

 ! Compute the defect
 CALL Resdfk_General_QuadScalar(QuadSc,ResU,ResV,ResW,DefUVW,RhsUVW)

END IF

! Checking convergence rates against criterions
RhsUVW=DefUVW
CALL COMM_Maximum(RhsUVW)
CALL Protocol_QuadScalar(mfile,QuadSc,INL,&
     ResU,ResV,ResW,DefUVW,RhsUVW)
! IF (ISNAN(RhsUVW)) stop
! 
! IF ((DefUVW.LE.DefUVWCrit).AND.&
!     (INL.GE.QuadSc%prm%NLmin)) INLComplete = 1
! 
! CALL COMM_NLComplete(INLComplete)
! CALL ZTIME(tttt1)
! myStat%tDefUVW = myStat%tDefUVW + (tttt1-tttt0)
! 
! IF (INLComplete.eq.1) GOTO 1
! 
! END DO
! 
! 1 CONTINUE
! 
! ! return
! myStat%iNonLin = myStat%iNonLin + INL
! inl_u = INL

IF (myid.ne.master) THEN
 mySolSeq%S((iStep-0)*lStep)%U = QuadSc%ValU
 mySolSeq%S((iStep-0)*lStep)%V = QuadSc%ValV
 mySolSeq%S((iStep-0)*lStep)%W = QuadSc%ValW
END IF

END DO ! iStep

! end do !iNonLin

END SUBROUTINE Lin_BurgerStep_ParT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE PressureStep_ParT()

DO iStep = 1,nSteps

IF (myid.ne.0) THEN

 CALL ZTIME(tttt0)
 ! Save the old solution
 ! LinSc%valP_old = LinSc%valP(NLMAX)%x
 LinSc%valP_old = mySolSeq%S((iStep-0)*lStep)%P
 LinSc%valP(NLMAX)%x = 0d0

 QuadSc%ValU = mySolSeq%S((iStep)*lStep)%U - mySolSeq%S((iStep-1)*lStep)%U 
 QuadSc%ValV = mySolSeq%S((iStep)*lStep)%V - mySolSeq%S((iStep-1)*lStep)%V
 QuadSc%ValW = mySolSeq%S((iStep)*lStep)%W - mySolSeq%S((iStep-1)*lStep)%W
! 
! Assemble the right hand side (RHS=1/k B^T U~)
 CALL Matdef_General_LinScalar(LinSc,QuadSc,PLinSc,1)

 ! Save the right hand side
 LinSc%rhsP(NLMAX)%x = LinSc%defP(NLMAX)%x

 CALL ZTIME(tttt1)
 myStat%tDefP = myStat%tDefP + (tttt1-tttt0)
END IF

! Calling the solver
CALL Solve_General_LinScalar(LinSc,PLinSc,QuadSc,Boundary_LinScalar_Mat,Boundary_LinScalar_Def,mfile)

MaxInitialPressureDefect = max(MaxInitialPressureDefect,LinSc%prm%MGprmOut%DefInitial*tstep)

CALL Protocol_LinScalar(mfile,LinSc," Pressure-Poisson equation")

2 CONTINUE

IF (myid.ne.0) THEN
 CALL ZTIME(tttt0)
 CALL Pressure_Correction()
 CALL ZTIME(tttt1)
 myStat%tCorrUVWP = myStat%tCorrUVWP + (tttt1-tttt0)
END IF

IF (myid.ne.master) THEN
 mySolSeq%S((iStep-0)*lStep)%P = LinSc%ValP(NLMAX)%x
END IF

END DO ! iStep

END SUBROUTINE PressureStep_ParT

END SUBROUTINE Transport_q2p1_UxyzP_ParT
!
! ----------------------------------------------
!
SUBROUTINE Transport_q2p1_UxyzP_fc_ext(mfile,inl_u,itns)
use cinterface, only: calculateDynamics,calculateFBM
use fbm, only: fbm_updateFBM

INTEGER mfile,INL,inl_u,itns
REAL*8  ResU,ResV,ResW,DefUVW,RhsUVW,DefUVWCrit
REAL*8  ResP,DefP,RhsPG,defPG,defDivU,DefPCrit
INTEGER INLComplete,I,J,IERR,iITER

CALL updateFBMGeometry()

thstep = tstep*(1d0-theta)

CALL OperatorRegenaration(2)

CALL OperatorRegenaration(3)

! -------------------------------------------------
! Compute the momentum equations
! -------------------------------------------------
! GOTO 1
IF (myid.ne.master) THEN

 CALL ZTIME(tttt0)

 ! Assemble the right hand side
 CALL Matdef_General_QuadScalar(QuadSc,1)

 ! Add the pressure gradient to the rhs
 CALL AddPressureGradient()
END IF

 ! Add the viscoelastic stress to the rhs
 IF(bViscoElastic)THEN
   CALL AddViscoStress()
 END IF

IF (myid.ne.master) THEN
 ! Add the gravity force to the rhs
 CALL AddGravForce()

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 ! Store the constant right hand side
 QuadSc%rhsU = QuadSc%defU
 QuadSc%rhsV = QuadSc%defV
 QuadSc%rhsW = QuadSc%defW

! Set dirichlet boundary conditions on the solution
 CALL Boundary_QuadScalar_Val()

END IF

thstep = tstep*theta

IF (myid.ne.master) THEN

 ! Assemble the defect vector and fine level matrix
 CALL Matdef_General_QuadScalar(QuadSc,-1)

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 QuadSc%auxU = QuadSc%defU
 QuadSc%auxV = QuadSc%defV
 QuadSc%auxW = QuadSc%defW
 CALL E013Sum3(QuadSc%auxU,QuadSc%auxV,QuadSc%auxW)
!  CALL E013Sum(QuadSc%auxU)
!  CALL E013Sum(QuadSc%auxV)
!  CALL E013Sum(QuadSc%auxW)

 ! Save the old solution
 CALL LCP1(QuadSc%valU,QuadSc%valU_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valV,QuadSc%valV_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valW,QuadSc%valW_old,QuadSc%ndof)

 ! Compute the norm of the defect
 CALL Resdfk_General_QuadScalar(QuadSc,ResU,ResV,ResW,DefUVW,RhsUVW)

END IF

CALL COMM_Maximum(RhsUVW)
DefUVWCrit=MAX(RhsUVW*QuadSc%prm%defCrit,QuadSc%prm%MinDef)

CALL Protocol_QuadScalar(mfile,QuadSc,0,&
     ResU,ResV,ResW,DefUVW,DefUVWCrit," Momentum equation ")

CALL ZTIME(tttt1)
myStat%tDefUVW = myStat%tDefUVW + (tttt1-tttt0)


DO INL=1,QuadSc%prm%NLmax
INLComplete = 0

! ! Calling the solver
CALL Solve_General_QuadScalar(QuadSc,Boundary_QuadScalar_Val,&
Boundary_QuadScalar_Mat,Boundary_QuadScalar_Mat_9,mfile)

!!!!          Checking the quality of the result           !!!!
!!!! ----------------------------------------------------- !!!!

CALL OperatorRegenaration(3)

IF (myid.ne.master) THEN
! Restore the constant right hand side
 CALL ZTIME(tttt0)
 QuadSc%defU = QuadSc%rhsU
 QuadSc%defV = QuadSc%rhsV
 QuadSc%defW = QuadSc%rhsW
END IF

IF (myid.ne.master) THEN

 ! Assemble the defect vector and fine level matrix
 CALL Matdef_General_QuadScalar(QuadSc,-1)

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 QuadSc%auxU = QuadSc%defU
 QuadSc%auxV = QuadSc%defV
 QuadSc%auxW = QuadSc%defW
 CALL E013Sum3(QuadSc%auxU,QuadSc%auxV,QuadSc%auxW)
! CALL E013Sum(QuadSc%auxU)
! CALL E013Sum(QuadSc%auxV)
! CALL E013Sum(QuadSc%auxW)

 ! Save the old solution
 CALL LCP1(QuadSc%valU,QuadSc%valU_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valV,QuadSc%valV_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valW,QuadSc%valW_old,QuadSc%ndof)

 ! Compute the defect
 CALL Resdfk_General_QuadScalar(QuadSc,ResU,ResV,ResW,DefUVW,RhsUVW)

END IF

! Checking convergence rates against criterions
RhsUVW=DefUVW
CALL COMM_Maximum(RhsUVW)
CALL Protocol_QuadScalar(mfile,QuadSc,INL,&
     ResU,ResV,ResW,DefUVW,RhsUVW)
IF (ISNAN(RhsUVW)) stop

IF ((DefUVW.LE.DefUVWCrit).AND.&
    (INL.GE.QuadSc%prm%NLmin)) INLComplete = 1

CALL COMM_NLComplete(INLComplete)
CALL ZTIME(tttt1)
myStat%tDefUVW = myStat%tDefUVW + (tttt1-tttt0)

IF (INLComplete.eq.1) GOTO 1
!IF (timens.lt.tstep+1d-8) GOTO 1

END DO

1 CONTINUE

! return
myStat%iNonLin = myStat%iNonLin + INL
inl_u = INL

! -------------------------------------------------
! Compute the pressure correction
! -------------------------------------------------
CALL MonitorVeloMag(QuadSc)

IF (myid.ne.0) THEN

 CALL ZTIME(tttt0)
 ! Save the old solution
 LinSc%valP_old = LinSc%valP(NLMAX)%x
 LinSc%valP(NLMAX)%x = 0d0

 ! Assemble the right hand side (RHS=1/k B^T U~)
 CALL Matdef_General_LinScalar(LinSc,QuadSc,PLinSc,1)

!  ! Assemble the right hand side (RHS:=RHS-C*Q)
!  CALL Matdef_General_LinScalar(LinSc,QuadSc,PLinSc,2)

 ! Save the right hand side
 LinSc%rhsP(NLMAX)%x = LinSc%defP(NLMAX)%x

 CALL ZTIME(tttt1)
 myStat%tDefP = myStat%tDefP + (tttt1-tttt0)
END IF

! Calling the solver
CALL Solve_General_LinScalar(LinSc,PLinSc,QuadSc,Boundary_LinScalar_Mat,Boundary_LinScalar_Def,mfile)

CALL Protocol_LinScalar(mfile,LinSc," Pressure-Poisson equation")

2 CONTINUE

IF (myid.ne.0) THEN
 CALL ZTIME(tttt0)
 !if (myid.eq.1) write(*,*) 'no correction ... '
 CALL Velocity_Correction()
 CALL Pressure_Correction()
 CALL ZTIME(tttt1)
 myStat%tCorrUVWP = myStat%tCorrUVWP + (tttt1-tttt0)
END IF

CALL QuadScP1toQ2(LinSc,QuadSc)

CALL FAC_GetForces(mfile)

CALL GetNonNewtViscosity()

IF (bNS_Stabilization) THEN
 CALL ExtractVeloGradients()
END IF

call fbm_updateFBM(Properties%Density(1),tstep,timens,&
                   Properties%Gravity,mfile,myid,&
                   QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                   LinSc%valP(NLMAX)%x,fbm_up_handler_ptr) 

!if (myid.eq.1) write(*,*) 'CommP: ',myStat%tCommP,'CommV: ',myStat%tCommV

IF (myid.ne.0) THEN
 CALL STORE_OLD_MESH(mg_mesh%level(NLMAX+1)%dcorvg)
END IF
 
CALL UmbrellaSmoother_ext(0d0,nUmbrellaSteps)
 
IF (myid.ne.0) THEN
 CALL STORE_NEW_MESH(mg_mesh%level(NLMAX+1)%dcorvg)
END IF
 
 CALL GET_MESH_VELO()
 
 ILEV=NLMAX
 CALL SETLEV(2)
 CALL SetUp_myQ2Coor( mg_mesh%level(ILEV)%dcorvg,&
                      mg_mesh%level(ILEV)%dcorag,&
                      mg_mesh%level(ILEV)%kvert,&
                      mg_mesh%level(ILEV)%karea,&
                      mg_mesh%level(ILEV)%kedge)

CALL updateFBMGeometry()

CALL MonitorVeloMag(QuadSc)

RETURN

END SUBROUTINE Transport_q2p1_UxyzP_fc_ext 
!
! ----------------------------------------------
!
SUBROUTINE Transport_q2p1_UxyzP_fc_ext_static(mfile,inl_u,itns)
use cinterface, only: calculateDynamics,calculateFBM
use fbm, only: fbm_updateFBM

INTEGER mfile,INL,inl_u,itns
REAL*8  ResU,ResV,ResW,DefUVW,RhsUVW,DefUVWCrit
REAL*8  ResP,DefP,RhsPG,defPG,defDivU,DefPCrit
INTEGER INLComplete,I,J,IERR,iITER

CALL updateFBMGeometry_Wangen()

thstep = tstep*(1d0-theta)

CALL OperatorRegenaration(2)

CALL OperatorRegenaration(3)

! -------------------------------------------------
! Compute the momentum equations
! -------------------------------------------------
! GOTO 1
IF (myid.ne.master) THEN

 CALL ZTIME(tttt0)

 ! Assemble the right hand side
 CALL Matdef_General_QuadScalar(QuadSc,1)

 ! Add the pressure gradient to the rhs
 CALL AddPressureGradient()
END IF

 ! Add the viscoelastic stress to the rhs
 IF(bViscoElastic)THEN
   CALL AddViscoStress()
 END IF

IF (myid.ne.master) THEN
 ! Add the gravity force to the rhs
 CALL AddGravForce()

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 ! Store the constant right hand side
 QuadSc%rhsU = QuadSc%defU
 QuadSc%rhsV = QuadSc%defV
 QuadSc%rhsW = QuadSc%defW

! Set dirichlet boundary conditions on the solution
 CALL Boundary_QuadScalar_Val()

END IF

thstep = tstep*theta

IF (myid.ne.master) THEN

 ! Assemble the defect vector and fine level matrix
 CALL Matdef_General_QuadScalar(QuadSc,-1)

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 QuadSc%auxU = QuadSc%defU
 QuadSc%auxV = QuadSc%defV
 QuadSc%auxW = QuadSc%defW
 CALL E013Sum3(QuadSc%auxU,QuadSc%auxV,QuadSc%auxW)
!  CALL E013Sum(QuadSc%auxU)
!  CALL E013Sum(QuadSc%auxV)
!  CALL E013Sum(QuadSc%auxW)

 ! Save the old solution
 CALL LCP1(QuadSc%valU,QuadSc%valU_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valV,QuadSc%valV_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valW,QuadSc%valW_old,QuadSc%ndof)

 ! Compute the norm of the defect
 CALL Resdfk_General_QuadScalar(QuadSc,ResU,ResV,ResW,DefUVW,RhsUVW)

END IF

CALL COMM_Maximum(RhsUVW)
DefUVWCrit=MAX(RhsUVW*QuadSc%prm%defCrit,QuadSc%prm%MinDef)

CALL Protocol_QuadScalar(mfile,QuadSc,0,&
     ResU,ResV,ResW,DefUVW,DefUVWCrit," Momentum equation ")

CALL ZTIME(tttt1)
myStat%tDefUVW = myStat%tDefUVW + (tttt1-tttt0)


DO INL=1,QuadSc%prm%NLmax
INLComplete = 0

! ! Calling the solver
CALL Solve_General_QuadScalar(QuadSc,Boundary_QuadScalar_Val,&
Boundary_QuadScalar_Mat,Boundary_QuadScalar_Mat_9,mfile)

!!!!          Checking the quality of the result           !!!!
!!!! ----------------------------------------------------- !!!!

CALL OperatorRegenaration(3)

IF (myid.ne.master) THEN
! Restore the constant right hand side
 CALL ZTIME(tttt0)
 QuadSc%defU = QuadSc%rhsU
 QuadSc%defV = QuadSc%rhsV
 QuadSc%defW = QuadSc%rhsW
END IF

IF (myid.ne.master) THEN

 ! Assemble the defect vector and fine level matrix
 CALL Matdef_General_QuadScalar(QuadSc,-1)

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 QuadSc%auxU = QuadSc%defU
 QuadSc%auxV = QuadSc%defV
 QuadSc%auxW = QuadSc%defW
 CALL E013Sum3(QuadSc%auxU,QuadSc%auxV,QuadSc%auxW)
! CALL E013Sum(QuadSc%auxU)
! CALL E013Sum(QuadSc%auxV)
! CALL E013Sum(QuadSc%auxW)

 ! Save the old solution
 CALL LCP1(QuadSc%valU,QuadSc%valU_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valV,QuadSc%valV_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valW,QuadSc%valW_old,QuadSc%ndof)

 ! Compute the defect
 CALL Resdfk_General_QuadScalar(QuadSc,ResU,ResV,ResW,DefUVW,RhsUVW)

END IF

! Checking convergence rates against criterions
RhsUVW=DefUVW
CALL COMM_Maximum(RhsUVW)
CALL Protocol_QuadScalar(mfile,QuadSc,INL,&
     ResU,ResV,ResW,DefUVW,RhsUVW)
IF (ISNAN(RhsUVW)) stop

IF ((DefUVW.LE.DefUVWCrit).AND.&
    (INL.GE.QuadSc%prm%NLmin)) INLComplete = 1

CALL COMM_NLComplete(INLComplete)
CALL ZTIME(tttt1)
myStat%tDefUVW = myStat%tDefUVW + (tttt1-tttt0)

IF (INLComplete.eq.1) GOTO 1
!IF (timens.lt.tstep+1d-8) GOTO 1

END DO

1 CONTINUE

! return
myStat%iNonLin = myStat%iNonLin + INL
inl_u = INL

! -------------------------------------------------
! Compute the pressure correction
! -------------------------------------------------

CALL MonitorVeloMag(QuadSc)

IF (myid.ne.0) THEN

 CALL ZTIME(tttt0)
 ! Save the old solution
 LinSc%valP_old = LinSc%valP(NLMAX)%x
 LinSc%valP(NLMAX)%x = 0d0

 ! Assemble the right hand side (RHS=1/k B^T U~)
 CALL Matdef_General_LinScalar(LinSc,QuadSc,PLinSc,1)

!  ! Assemble the right hand side (RHS:=RHS-C*Q)
!  CALL Matdef_General_LinScalar(LinSc,QuadSc,PLinSc,2)

 ! Save the right hand side
 LinSc%rhsP(NLMAX)%x = LinSc%defP(NLMAX)%x

 CALL ZTIME(tttt1)
 myStat%tDefP = myStat%tDefP + (tttt1-tttt0)
END IF

! Calling the solver
CALL Solve_General_LinScalar(LinSc,PLinSc,QuadSc,Boundary_LinScalar_Mat,Boundary_LinScalar_Def,mfile)

CALL Protocol_LinScalar(mfile,LinSc," Pressure-Poisson equation")

2 CONTINUE

IF (myid.ne.0) THEN
 CALL ZTIME(tttt0)
 !if (myid.eq.1) write(*,*) 'no correction ... '
 CALL Velocity_Correction()
 CALL Pressure_Correction()
 CALL ZTIME(tttt1)
 myStat%tCorrUVWP = myStat%tCorrUVWP + (tttt1-tttt0)
END IF

CALL QuadScP1toQ2(LinSc,QuadSc)

CALL FAC_GetForces(mfile)

CALL GetNonNewtViscosity()

IF (bNS_Stabilization) THEN
 CALL ExtractVeloGradients()
END IF

RETURN

END SUBROUTINE Transport_q2p1_UxyzP_fc_ext_static 
!
! ----------------------------------------------
!
SUBROUTINE Init_Q2_Structures(mfile)
implicit none
LOGICAL bExist
INTEGER I,J,ndof,mfile,LevDif
integer :: mydof
integer :: maxlevel
Real*8 :: dabl

 ILEV=NLMAX
 CALL SETLEV(2)

 ! Initialize the scalar quantity
 CALL InitializeQuadScalar(QuadSc)

 ! Initialize the scalar quantity
 CALL InitializeLinScalar(LinSc)

 ! Initialize the boundary list (QuadScBoundary)
 ALLOCATE (QuadScBoundary(mg_mesh%level(ilev)%nvt+&
                          mg_mesh%level(ilev)%net+&
                          mg_mesh%level(ilev)%nat+&
                          mg_mesh%level(ilev)%nel))

 CALL InitBoundaryList(KWORK(L(LNPR)),&
                      mg_mesh%level(ILEV)%kvert,&
                      mg_mesh%level(ILEV)%kedge,&
                      mg_mesh%level(ILEV)%karea)

 ILEV=NLMAX
 CALL SETLEV(2)

 ! Set up the Coordinate Vector
 ALLOCATE (myQ2Coor(3,mg_mesh%level(ilev)%nvt+&
                      mg_mesh%level(ilev)%net+&
                      mg_mesh%level(ilev)%nat+&
                      mg_mesh%level(ilev)%nel))

 CALL SetUp_myQ2Coor( mg_mesh%level(ILEV)%dcorvg,&
                      mg_mesh%level(ILEV)%dcorag,&
                      mg_mesh%level(ILEV)%kvert,&
                      mg_mesh%level(ILEV)%karea,&
                      mg_mesh%level(ILEV)%kedge)

 !
 !IF (myid.ne.0) CALL ParametrizeQ2Nodes(myQ2Coor)
 !

 ALLOCATE(myALE%Q2coor_old(3,&
 mg_mesh%level(ilev)%nvt+&
 mg_mesh%level(ilev)%net+&
 mg_mesh%level(ilev)%nat+&
 mg_mesh%level(ilev)%nel))

 myALE%Q2coor_old = myQ2Coor

 ALLOCATE(myALE%MeshVelo(3,&
 mg_mesh%level(ilev)%nvt+&
 mg_mesh%level(ilev)%net+&
 mg_mesh%level(ilev)%nat+&
 mg_mesh%level(ilev)%nel))

 myALE%MeshVelo = 0d0

 CALL InitBoundaryStructure(mg_mesh%level(ILEV)%kvert,&
                            mg_mesh%level(ILEV)%kedge)

 Properties%cName = "Prop"
 CALL GetPhysiclaParameters(Properties,Properties%cName,mfile)

 myPowerLawFluid(2) = 0.001d0
 myPowerLawFluid(3) = 0.75d0

 ! Initialize the arrays and the distribution of physical properties
 ALLOCATE (mgDensity(NLMIN:NLMAX))
 ALLOCATE (mgNormShearStress(NLMIN:NLMAX))
 DO ILEV=NLMIN,NLMAX

  ALLOCATE (mgDensity(ILEV)%x(mg_mesh%level(ilev)%nel))
  ALLOCATE (mgNormShearStress(ILEV)%x(mg_mesh%level(ilev)%nel))
  mgDensity(ILEV)%x          = Properties%Density(1)
  mgNormShearStress(ILEV)%x  = 0d0

 END DO

 if(myid.ne.0)then
!---------------------                          
 ALLOCATE (mgDiffCoeff(NLMIN:NLMAX+1))
 DO ILEV=NLMIN,NLMAX+1
  ALLOCATE (mgDiffCoeff(ILEV)%x(mg_mesh%level(ilev)%nel))
  mgDiffCoeff(ILEV)%x = Properties%DiffCoeff(1)
 END DO
!---------------------                          
else
 maxlevel = mg_Mesh%nlmax
 ALLOCATE (mgDiffCoeff(NLMIN:maxlevel))
 DO ILEV=NLMIN,maxlevel
  ALLOCATE (mgDiffCoeff(ILEV)%x(mg_mesh%level(ilev)%nel))
  mgDiffCoeff(ILEV)%x = Properties%DiffCoeff(1)
 END DO
end if

 ILEV = NLMAX
 ALLOCATE (Viscosity(mg_mesh%level(ilev)%nvt+&
                     mg_mesh%level(ilev)%net+&
                     mg_mesh%level(ilev)%nat+&
                     mg_mesh%level(ilev)%nel))

 Viscosity = Properties%Viscosity(1)

 mydof = mg_mesh%level(ilev)%nvt+&
         mg_mesh%level(ilev)%net+&
         mg_mesh%level(ilev)%nat+&
         mg_mesh%level(ilev)%nel

 ALLOCATE (myALE%Monitor(mydof))
 ALLOCATE (myALE%NewCoor(3,mydof))
 ALLOCATE (myALE%OldCoor(3,mydof))
 ALLOCATE (myALE%OrigCoor(3,mydof))

 myALE%Monitor   = 1d0
 myALE%MeshVelo  = 0d0

 ! Building up the E013/E013 matrix strucrures
 CALL Create_QuadMatStruct()

 ! Iteration matrix (only allocation)
 CALL Create_AMat() !(A)

 ! Building up the E012/E013 E013/E012 and matrix structures
 CALL Create_QuadLinMatStruct() 

 ! Building up the E012/E012 matrix strucrures
 CALL Create_LinMatStruct ()

 ! Pressure gradient matrix
 CALL Create_BMat() !(B,BT)

 IF (myid.EQ.ShowID) WRITE(MTERM,'(A)', advance='yes') " "

 IF (myid.ne.master) THEN
  ! Parallel E012/E013 matrix structure
  CALL Create_QuadLinParMatStruct(PLinSc) !(pB)

  ! Building up the Parallel E012/E012 matrix strucrures
  CALL Create_ParLinMatStruct ()
 END IF

! Set up the boundary condition types (knpr)
 DO ILEV=NLMIN,NLMAX
  CALL SETLEV(2)
  CALL QuadScalar_Knpr()
 END DO
 ILEV=NLMAX
 mydof = mg_mesh%level(ilev)%nvt+&
         mg_mesh%level(ilev)%net+&
         mg_mesh%level(ilev)%nat+&
         mg_mesh%level(ilev)%nel

 ALLOCATE (FictKNPR(mydof))
 FictKNPR=0
 ALLOCATE (Distance(mydof))
 Distance = 0d0

 ALLOCATE (MixerKNPR(mydof))
 MixerKNPR=0
 ALLOCATE (Distamce(mydof))
 Distamce = 0d0

 ! SEt up the knpr vector showing dofs with parallel property ...
 IF (myid.ne.0) THEN
  ALLOCATE (ParKNPR(NVT+NET+NAT+NEL))
  QuadSc%auxU = 1d0
  CALL E013Sum(QuadSc%auxU)
  DO I=1,NVT+NET+NAT+NEL
   IF (QuadSc%auxU(I).EQ.1d0) THEN
    ParKNPR(I) = 0
   ELSE
    ParKNPR(I) = 1
   END IF
  END DO
 END IF

 IF (myid.eq.showID) THEN
  INQUIRE (FILE="_data/BenchValues.txt", EXIST=bExist)
  IF (ISTART.EQ.0.OR.(.NOT.bExist)) THEN
   OPEN(666,FILE="_data/BenchValues.txt")
   WRITE(666,'(4A16)') "Time","Drag","Lift","ZForce"
  ELSE
   OPEN(666,FILE="_data/BenchValues.txt",ACCESS='APPEND')
  END IF
 END IF

 CALL InitializeProlRest(QuadSc,LinSc)

 CALL OperatorRegenaration(1)

END SUBROUTINE Init_Q2_Structures
!
! ----------------------------------------------
!
SUBROUTINE InitCond_Disp_QuadScalar()

 ILEV=NLMAX
 CALL SETLEV(2)

IF (myid.ne.0) THEN

 ! Set initial conditions
 CALL QuadScalar_InitCond()

 IF (myFBM%nParticles.GT.0) THEN
  CALL updateFBMGeometry()
 END IF

 ! Set dirichlet boundary conditions on the solution
 CALL Boundary_QuadScalar_Val()

 ! Set initial conditions
 CALL LinScalar_InitCond(mg_mesh%level(ilev)%dcorvg,&
                         mg_mesh%level(ilev)%kvert)

END IF

end subroutine InitCond_Disp_QuadScalar
!
! ----------------------------------------------
!
SUBROUTINE InitCond_GenLinSc_Q1_QuadScalar()

 ILEV=NLMAX
 CALL SETLEV(2)

IF (myid.ne.0) THEN

 ! Set initial conditions
 CALL QuadScalar_InitCond()

 IF (myFBM%nParticles.GT.0) THEN
  CALL updateFBMGeometry()
 END IF

 ! Set dirichlet boundary conditions on the solution
 CALL Boundary_QuadScalar_Val()

 ! Set initial conditions
 CALL LinScalar_InitCond(mg_mesh%level(ilev)%dcorvg,&
                         mg_mesh%level(ilev)%kvert)

END IF

end subroutine InitCond_GenLinSc_Q1_QuadScalar
!
! ----------------------------------------------
!
subroutine Transport_q2p1_UxyzP_sse(mfile, inl_u, itns)

use, intrinsic :: ieee_arithmetic

INTEGER mfile,INL,inl_u,itns
REAL*8  ResU,ResV,ResW,DefUVW,RhsUVW,DefUVWCrit
REAL*8  ResP,DefP,RhsPG,defPG,defDivU,DefPCrit
INTEGER INLComplete,I,J,IERR,iITER

if (.not.allocated(MGSteps%n)) THEN 
 allocate(MGSteps%n(MGSteps%m))
 allocate(MGSteps%r(MGSteps%m))
 MGSteps%i = 0
 MGSteps%n = 0
 MGSteps%r = 0d0
END IF

thstep = tstep*(1d0-theta)

ILEV = NLMAX
CALL SETLEV(2)
CALL OperatorRegenaration(2)

!CALL Setup_PeriodicPressureRHS(LinSc,PLinSc,QuadSc%knprU,QuadSc%knprV,QuadSc%knprW)
ILEV = NLMAX
CALL SETLEV(2)
CALL Setup_PeriodicVelocityRHS()

CALL OperatorRegenaration(3)

! -------------------------------------------------
! Compute the momentum equations
! -------------------------------------------------
IF (myid.ne.master) THEN

 CALL ZTIME(tttt0)

 ! Assemble the right hand side
 CALL Matdef_General_QuadScalar(QuadSc,1)

 ! Add the pressure gradient
 CALL AddPeriodicPressureGradient()

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 ! Store the constant right hand side
 QuadSc%rhsU = QuadSc%defU
 QuadSc%rhsV = QuadSc%defV
 QuadSc%rhsW = QuadSc%defW

! Set dirichlet boundary conditions on the solution
 CALL Boundary_QuadScalar_Val()

END IF

  
thstep = tstep*theta

IF (myid.ne.master) THEN

 ! Assemble the defect vector and fine level matrix
 CALL Matdef_General_QuadScalar(QuadSc,-1)

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 QuadSc%auxU = QuadSc%defU
 QuadSc%auxV = QuadSc%defV
 QuadSc%auxW = QuadSc%defW
 CALL E013Sum(QuadSc%auxU)
 CALL E013Sum(QuadSc%auxV)
 CALL E013Sum(QuadSc%auxW)

 ! Save the old solution
 CALL LCP1(QuadSc%valU,QuadSc%valU_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valV,QuadSc%valV_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valW,QuadSc%valW_old,QuadSc%ndof)

 ! Compute the norm of the defect
 CALL Resdfk_General_QuadScalar(QuadSc,ResU,ResV,ResW,DefUVW,RhsUVW)

END IF

CALL COMM_Maximum(RhsUVW)
DefUVWCrit=MAX(RhsUVW*QuadSc%prm%defCrit,QuadSc%prm%MinDef)

CALL Protocol_QuadScalar(mfile,QuadSc,0,&
     ResU,ResV,ResW,DefUVW,DefUVWCrit," Momentum equation ")

CALL ZTIME(tttt1)
myStat%tDefUVW = myStat%tDefUVW + (tttt1-tttt0)

DO INL=1,QuadSc%prm%NLmax
  INLComplete = 0

  ! ! Calling the solver
  CALL Solve_General_QuadScalar(QuadSc,Boundary_QuadScalar_Val,&
  Boundary_QuadScalar_Mat,Boundary_QuadScalar_Mat_9,mfile)

  MGSteps%j = MOD(MGSteps%i,MGSteps%m) + 1
  MGSteps%n = 0
  MGSteps%r(MGSteps%j) = QuadSc%prm%MGprmOut(1)%RhoMG1
  MGSteps%n(MGSteps%j) = QuadSc%prm%MGprmIn%nSmootherSteps
  MGSteps%i = MGSteps%i + 1
  MGSteps%daux = 0d0
!   if (myid.eq.1) WRITE(*,*) 'MG data:', MGSteps%i,MGSteps%m
  if (MGSteps%i.ge.MGSteps%m) THEN
   DO i = 1,MGSteps%m
    MGSteps%daux = MGSteps%daux + MGSteps%r(i)
   END DO
   MGSteps%daux = MGSteps%daux/dble(MGSteps%m)
!    if (myid.eq.1) WRITE(*,*) 'AVG MG rate:', MGSteps%daux,MGSteps%r
   IF (MGSteps%daux.lt.3d-3) THEN
    QuadSc%prm%MGprmIn%nSmootherSteps = MAX(MIN(INT(DBLE(QuadSc%prm%MGprmIn%nSmootherSteps)*0.75d0),QuadSc%prm%MGprmIn%nSmootherSteps-1),4)
    if (myid.eq.1) WRITE(*,*) 'MG smoothening step reduction to:', QuadSc%prm%MGprmIn%nSmootherSteps
   END IF
   IF (MGSteps%daux.gt.3d-1) THEN
    QuadSc%prm%MGprmIn%nSmootherSteps = MIN(INT(DBLE(QuadSc%prm%MGprmIn%nSmootherSteps)*1.25d0),32)
    if (myid.eq.1) WRITE(*,*) 'MG smoothening step increasement to:', QuadSc%prm%MGprmIn%nSmootherSteps
   END IF
  END IF

  !!!!!          Checking the quality of the result           !!!!
  !!!!! ----------------------------------------------------- !!!!

  CALL OperatorRegenaration(3)

  IF (myid.ne.master) THEN
  ! Restore the constant right hand side
   CALL ZTIME(tttt0)
   QuadSc%defU = QuadSc%rhsU
   QuadSc%defV = QuadSc%rhsV
   QuadSc%defW = QuadSc%rhsW
  END IF

  IF (myid.ne.master) THEN

   ! Assemble the defect vector and fine level matrix
   CALL Matdef_General_QuadScalar(QuadSc,-1)

   ! Set dirichlet boundary conditions on the defect
   CALL Boundary_QuadScalar_Def()

   QuadSc%auxU = QuadSc%defU
   QuadSc%auxV = QuadSc%defV
   QuadSc%auxW = QuadSc%defW
   CALL E013Sum(QuadSc%auxU)
   CALL E013Sum(QuadSc%auxV)
   CALL E013Sum(QuadSc%auxW)

   ! Save the old solution
   CALL LCP1(QuadSc%valU,QuadSc%valU_old,QuadSc%ndof)
   CALL LCP1(QuadSc%valV,QuadSc%valV_old,QuadSc%ndof)
   CALL LCP1(QuadSc%valW,QuadSc%valW_old,QuadSc%ndof)

   ! Compute the defect
   CALL Resdfk_General_QuadScalar(QuadSc,ResU,ResV,ResW,DefUVW,RhsUVW)

  END IF

  ! Checking convergence rates against criterions
  RhsUVW=DefUVW
  CALL COMM_Maximum(RhsUVW)
  CALL Protocol_QuadScalar(mfile,QuadSc,INL,&
       ResU,ResV,ResW,DefUVW,RhsUVW)
  IF (ISNAN(RhsUVW)) stop

  IF ((DefUVW.LE.DefUVWCrit).AND.&
     (INL.GE.QuadSc%prm%NLmin)) INLComplete = 1

  CALL COMM_NLComplete(INLComplete)
  CALL ZTIME(tttt1)
  myStat%tDefUVW = myStat%tDefUVW + (tttt1-tttt0)

  IF (INLComplete.eq.1) exit 

END DO

myStat%iNonLin = myStat%iNonLin + INL
inl_u = INL

! -------------------------------------------------
! Compute the pressure correction
! -------------------------------------------------
IF (myid.ne.0) THEN

 CALL ZTIME(tttt0)
 ! Save the old solution
 LinSc%valP_old = LinSc%valP(NLMAX)%x
 LinSc%valP(NLMAX)%x = 0d0

 ! Assemble the right hand side (RHS=1/k B^T U~)
 CALL Matdef_General_LinScalar(LinSc,QuadSc,PLinSc,1)

 ! Save the right hand side
 LinSc%rhsP(NLMAX)%x = LinSc%defP(NLMAX)%x

 CALL ZTIME(tttt1)
 myStat%tDefP = myStat%tDefP + (tttt1-tttt0)

END IF

! Calling the solver
CALL Solve_General_LinScalar(LinSc,PLinSc,QuadSc,Boundary_LinScalar_Mat,Boundary_LinScalar_Def,mfile)

CALL Protocol_LinScalar(mfile,LinSc," Pressure-Poisson equation")

2 CONTINUE

IF (myid.ne.0) THEN
 CALL ZTIME(tttt0)
 CALL Velocity_Correction()
 CALL Pressure_Correction()
 CALL ZTIME(tttt1)
 myStat%tCorrUVWP = myStat%tCorrUVWP + (tttt1-tttt0)
END IF

CALL CorrectPressure_WRT_KNPRP()

IF (ieee_is_finite(myProcess%dPress)) THEN
 CALL QuadScP1toQ2Periodic(LinSc,QuadSc)
ELSE
 CALL QuadScP1toQ2(LinSc,QuadSc)
END IF

CALL GetNonNewtViscosity_sse()

IF (.not.bKTPRelease) then
 CALL Calculate_Torque(mfile)
END IF

CALL LL21(LinSc%valP(NLMAX)%x,LinSc%ndof,defP)
call COMM_SUMM(defp)
if (ieee_is_nan(defP)) DivergedSolution = .true.
if (.not.ieee_is_finite(defP)) DivergedSolution = .true.

ILEV = NLMAX
CALL SETLEV(2)
CALL Comm_Solution(QuadSc%ValU,QuadSc%ValV,QuadSc%ValW,QuadSc%auxU,QuadSc%ndof)

end subroutine Transport_q2p1_UxyzP_sse
!
! ----------------------------------------------
!
SUBROUTINE TemporalFieldInterpolator(iL,iS)
USE Sigma_User, Only : myTransientSolution
INTEGER iL,iS,iL1,iS1,iL2,iS2
INTEGER i,nnn
REAL*8 dS,t_bu,MeshVelo(3)
integer nLengthV,nLengthE,LevDif
REAL*8 , ALLOCATABLE :: SendVect(:,:,:)

iL1 = iL
iL2 = iL+1
if (iL2.eq.myProcess%nTimeLevels) iL2 = 0

ILEV = NLMAX

dS = DBLE(iS)/DBLE(myTransientSolution%nTimeSubStep)
if (myid.eq.1) write(*,*) 'Time fraction: ',dS,LinSc%prm%MGprmIn%MedLev
IF (myid.ne.0) then
 DO i=1,QuadSc%ndof
 
  QuadSc%ValU(i)                  = (1d0-dS)*myTransientSolution%Velo(1,iL1)%x(i) + (dS)*myTransientSolution%Velo(1,iL2)%x(i)
  QuadSc%ValV(i)                  = (1d0-dS)*myTransientSolution%Velo(2,iL1)%x(i) + (dS)*myTransientSolution%Velo(2,iL2)%x(i)
  QuadSc%ValW(i)                  = (1d0-dS)*myTransientSolution%Velo(3,iL1)%x(i) + (dS)*myTransientSolution%Velo(3,iL2)%x(i)

  mg_mesh%level(ILEV)%dcorvg(1,i) = (1d0-dS)*myTransientSolution%Coor(1,iL1)%x(i) + (dS)*myTransientSolution%Coor(1,iL2)%x(i)
  mg_mesh%level(ILEV)%dcorvg(2,i) = (1d0-dS)*myTransientSolution%Coor(2,iL1)%x(i) + (dS)*myTransientSolution%Coor(2,iL2)%x(i)
  mg_mesh%level(ILEV)%dcorvg(3,i) = (1d0-dS)*myTransientSolution%Coor(3,iL1)%x(i) + (dS)*myTransientSolution%Coor(3,iL2)%x(i)
 
  Screw(i) = (1d0-dS)*myTransientSolution%Dist(iL1)%x(i) + (dS)*myTransientSolution%Dist(iL2)%x(i)
  
  IF (myProcess%SegmentThermoPhysProps) THEN
   mySegmentIndicator(2,i) = DBLE(NINT((1d0-dS)*myTransientSolution%iSeg(iL1)%x(i) + (dS)*myTransientSolution%iSeg(iL2)%x(i)))
  END IF
  
 END DO
END IF

IF (myid.EQ.0) THEN
  CALL CreateDumpStructures(0)
ELSE
  LevDif = LinSc%prm%MGprmIn%MedLev - NLMAX
  CALL CreateDumpStructures(LevDif)
END IF
 
ILEV = LinSc%prm%MGprmIn%MedLev

nLengthV = (2**(ILEV-1)+1)**3
nLengthE = mg_mesh%level(NLMIN)%nel

ALLOCATE(SendVect(3,nLengthV,nLengthE))

CALL SendNodeValuesToCoarse(SendVect,mg_mesh%level(NLMAX)%dcorvg,&
                            mg_mesh%level(ILEV)%kvert,&
                            nLengthV,&
                            nLengthE,&
                            mg_mesh%level(ILEV)%nel,&
                            mg_mesh%level(ILEV)%nvt)
DEALLOCATE(SendVect)

CALL Create_MRhoMat()
IF (myid.EQ.ShowID) WRITE(MTERM,'(A)', advance='yes') " "

CALL GetNonNewtViscosity_sse()

!!!!!!!!!!!!!!!!Correction of the Velocities due to the ALE components !!!!!!!!!!!!!!!!!!!

IF (myid.ne.0) then
 DO i=1,QuadSc%ndof
 
  MeshVelo(1)     = (myTransientSolution%Coor(1,iL2)%x(i)-myTransientSolution%Coor(1,iL1)%x(i))/(DBLE(myTransientSolution%nTimeSubStep)*tstep)
  MeshVelo(2)     = (myTransientSolution%Coor(2,iL2)%x(i)-myTransientSolution%Coor(2,iL1)%x(i))/(DBLE(myTransientSolution%nTimeSubStep)*tstep)
  MeshVelo(3)     = (myTransientSolution%Coor(3,iL2)%x(i)-myTransientSolution%Coor(3,iL1)%x(i))/(DBLE(myTransientSolution%nTimeSubStep)*tstep)
  
  QuadSc%ValU(i)  = QuadSc%ValU(i) - MeshVelo(1)
  QuadSc%ValV(i)  = QuadSc%ValV(i) - MeshVelo(2)
  QuadSc%ValW(i)  = QuadSc%ValW(i) - MeshVelo(3)
  
 END DO
END IF
!!!!!!!!!!!!!!!!Correction of the Velocities due to the ALE components !!!!!!!!!!!!!!!!!!!

END SUBROUTINE TemporalFieldInterpolator


SUBROUTINE CorrectPressure_WRT_KNPRP()
integer iel

if (myid.ne.0) then
 ilev = nlmax
 do iel = 1,mg_mesh%level(ilev)%nel
  if (LinSc%knprP(ilev)%x(iel).eq.1) then
   LinSc%valP(NLMAX)%x(4*(iel-1)+1) = 0d0
   LinSc%valP(NLMAX)%x(4*(iel-1)+2) = 0d0
   LinSc%valP(NLMAX)%x(4*(iel-1)+3) = 0d0
   LinSc%valP(NLMAX)%x(4*(iel-1)+4) = 0d0
  end if
 end do
end if

END SUBROUTINE CorrectPressure_WRT_KNPRP