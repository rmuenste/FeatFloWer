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
 integer :: ndofP,ndofU
 integer :: nOuter,nFOuter
 integer :: nSteps
 TYPE(tSingleSol), allocatable :: S(:)
END TYPE tSolSeq

TYPE(tSolSeq),TARGET :: Sol2,Sol3
TYPE(tSolSeq), POINTER :: mySolSeq

INTEGER :: LinIntPol=1,ConstIntPol=0,iOX,ndofP,ndofU,nLinBurgers
CHARACTER :: PredictionScheme*(2)='PP'
REAL*8, allocatable :: PressureDefect(:)

tstep_BU = tstep

CALL InitSolutions(Sol3,0)
mySolSeq => Sol3

CALL InitSolutions(Sol2,1)
mySolSeq => Sol2
if (myid.ne.master) then
 nlmax = nlmax - 1
 QuadSc%ndof = knvt(nlmax)+knel(nlmax)+knet(nlmax)+knat(nlmax)
 PLinSc%ndof = 4*knel(nlmax)
 LinSc%ndof = 4*knel(nlmax)
end if

!------------------------------------------------------ PREDICTION --------------------------------------------------------
MaxInitialPressureDefect = 1e-30 
iLoop =1

! DO iOuter=1,mySolSeq%nOuter
DO iOuter=mySolSeq%nOuter,mySolSeq%nOuter

 nSteps = 2**(iOuter-1)
 lStep  = mySolSeq%nSteps/nSteps
 tstep  = tstep_BU/DBLE(nSteps)
 
!  IF (iOuter.le.1) THEN
 IF (iOuter.le.mySolSeq%nOuter) THEN
!  IF (iOuter.le.mySolSeq%nOuter-2) THEN
 
  IF (PredictionScheme.eq.'CP') THEN
   CALL NonLin_SequentialBurgerStep_ParT()

   CALL PressureStepCP_ParT()
  
   ! Velocity update for last pressure
   CALL NonLin_SequentialBurgerStep_ParT()
  END IF
  
  IF (PredictionScheme.eq.'PP') THEN

   DO iStep = 1,nSteps
    CALL NonLin_BurgerStep_SingleT()

    CALL PressureStep_SingleT()
    
    IF (myid.ne.master) THEN
     mySolSeq%S(iStep)%P(1:LinSc%ndof)  = LinSc%ValP(NLMAX)%x(1:LinSc%ndof)
     mySolSeq%S(iStep)%U(1:QuadSc%ndof) = QuadSc%ValU(1:QuadSc%ndof)
     mySolSeq%S(iStep)%V(1:QuadSc%ndof) = QuadSc%ValV(1:QuadSc%ndof)
     mySolSeq%S(iStep)%W(1:QuadSc%ndof) = QuadSc%ValW(1:QuadSc%ndof)
    end if
    
   END DO
  
  END IF
  
 END IF
  
 IF (iOuter.ne.mySolSeq%nOuter) then
  CALL InterpolatePressure_ParT(LinIntPol)
  CALL InterpolateVelocity_ParT(LinIntPol)
 END IF

END DO !iOuter
tstep  = tstep_BU


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IF (myid.ne.master) THEN
!  QuadSc%ValU(1:QuadSc%ndof) = mySolSeq%S(mySolSeq%nSteps)%U(1:QuadSc%ndof)
!  QuadSc%ValV(1:QuadSc%ndof) = mySolSeq%S(mySolSeq%nSteps)%V(1:QuadSc%ndof)
!  QuadSc%ValW(1:QuadSc%ndof) = mySolSeq%S(mySolSeq%nSteps)%W(1:QuadSc%ndof)
! end if
! 
! if (bExit) THEN
!  tstep  = tstep_BU
!  GOTO 100
! END IF
! 
! IF (myid.ne.master) THEN
!  LinSc%P_new = 1.5d0*mySolSeq%S(mySolSeq%nSteps-0)%P(1:LinSc%ndof)  - 0.5d0*mySolSeq%S(mySolSeq%nSteps-1)%P(1:LinSc%ndof)
! end if
! CALL QuadScP1ExtPoltoQ2(LinSc,QuadSc)
! CALL FAC_GetForcesParT(mfile,iLoop)
! 
! return

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MaxInitialPressureDefect0 = MaxInitialPressureDefect
IF (myid.eq.1) WRITE(*,*) 'rock and roll! ',MaxInitialPressureDefect0
! 
! DO iStep = 0,mySolSeq%nSteps
! if (myid.ne.master) then
!  QuadSc%ValU = Sol3%S(iStep)%U(1:Sol3%ndofU)
!  QuadSc%ValV = Sol3%S(iStep)%V(1:Sol3%ndofU)
!  QuadSc%ValW = Sol3%S(iStep)%W(1:Sol3%ndofU)
!  LinSc%ValP(NLMAX)%x = Sol3%S(iStep)%P(1:Sol3%ndofP)
!  
!  LinSc%P_new = LinSc%ValP(NLMAX)%x
!  CALL QuadScP1ExtPoltoQ2(LinSc,QuadSc)
! end if
! CALL Output_Profiles(iStep)
! end do
! pause

if (myid.ne.master) then
 nlmax = nlmax + 1
 QuadSc%ndof = knvt(nlmax)+knel(nlmax)+knet(nlmax)+knat(nlmax)
 PLinSc%ndof = 4*knel(nlmax)
 LinSc%ndof = 4*knel(nlmax)
end if

DO iStep = 1,mySolSeq%nSteps

 if (myid.ne.master) then
  QuadSc%ValU = Sol2%S(iStep)%U(1:Sol2%ndofU)
  QuadSc%ValV = Sol2%S(iStep)%V(1:Sol2%ndofU)
  QuadSc%ValW = Sol2%S(iStep)%W(1:Sol2%ndofU)
  LinSc%ValP(NLMAX-1)%x = mySolSeq%S(iStep)%P(1:mySolSeq%ndofP)
  
  ndofP = 4*KNEL(NLMAX-1)
  LinSc%valP(NLMAX)%x(1:ndofP) = LinSc%valP(NLMAX-1)%x(1:ndofP)
  
  CALL ProlongateSolution()
  
  Sol3%S(iStep)%U = QuadSc%ValU(1:Sol3%ndofU)
  Sol3%S(iStep)%V = QuadSc%ValV(1:Sol3%ndofU)
  Sol3%S(iStep)%W = QuadSc%ValW(1:Sol3%ndofU)
  Sol3%S(iStep)%P = LinSc%ValP(NLMAX)%x(1:Sol3%ndofP)
  
!   LinSc%P_new = LinSc%ValP(NLMAX)%x
!   CALL QuadScP1ExtPoltoQ2(LinSc,QuadSc)
 end if
  
!  CALL Output_Profiles(iStep)
  
END DO
!------------------------------------------------------ PREDICTION --------------------------------------------------------
! pause
mySolSeq => Sol3

bEXIT = .false.
do iLoop =1 , Properties%nTPIterations

 MaxInitialPressureDefect = 0d0

 DO iOuter=mySolSeq%nFOuter,mySolSeq%nFOuter

  nSteps = 2**(iOuter-1)
  lStep  = mySolSeq%nSteps/nSteps
  tstep  = tstep_BU/DBLE(nSteps)

!   nLinBurgers = 1
!   do iOX=1,nLinBurgers
!    IF (myid.eq.1) WRITE(*,'(A,I0,A,i0)') "NumberOfLinBurgers: ",iOX," / ",nLinBurgers 
!    CALL Lin_BurgerStep_ParT()
!   end do
  CALL NonLin_SimultanBurgerStep_ParT()

!  CALL NonLin_SequentialBurgerStep_ParT()
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (myid.ne.master) THEN
   QuadSc%ValU(1:QuadSc%ndof) = mySolSeq%S(mySolSeq%nSteps)%U(1:QuadSc%ndof)
   QuadSc%ValV(1:QuadSc%ndof) = mySolSeq%S(mySolSeq%nSteps)%V(1:QuadSc%ndof)
   QuadSc%ValW(1:QuadSc%ndof) = mySolSeq%S(mySolSeq%nSteps)%W(1:QuadSc%ndof)
  end if

  if (bExit) THEN
   tstep  = tstep_BU
   GOTO 100
  END IF
  
  IF (myid.ne.master) THEN
   LinSc%P_new = 1.5d0*mySolSeq%S(mySolSeq%nSteps-0)%P(1:LinSc%ndof)  - 0.5d0*mySolSeq%S(mySolSeq%nSteps-1)%P(1:LinSc%ndof)
  end if
  CALL QuadScP1ExtPoltoQ2(LinSc,QuadSc)
  CALL FAC_GetForcesParT(mfile,iLoop)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  CALL PressureStepCP_ParT()

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
 if (myid.eq.1) write(mterm,'(A,I0,A,256ES12.4)') 'PressureDefects ', iLoop,cEnding(iEnding)//' : ',PressureDefect(:)
 if (myid.eq.1) write(mfile,'(A,I0,A,256ES12.4)') 'PressureDefects ', iLoop,cEnding(iEnding)//' : ',PressureDefect(:)
!  if (MaxInitialPressureDefect.lt.1e-7) THEN 
! if (MaxInitialPressureDefect/MaxInitialPressureDefect0.lt.Properties%DiracEps) THEN 
 if (MaxInitialPressureDefect/MaxInitialPressureDefect0.lt.Properties%DiracEps.and.MaxInitialPressureDefect.lt.1e-10) THEN 
  if (myid.eq.1) write(mterm,'(A,I0,A,3ES12.4)') 'ExitingCoarseTimeCPloopIn ', iLoop,cEnding(iEnding)//' step!|Init&FinPresDefect: ',MaxInitialPressureDefect0,MaxInitialPressureDefect,MaxInitialPressureDefect/MaxInitialPressureDefect0
  if (myid.eq.1) write(mfile,'(A,I0,A,3ES12.4)') 'ExitingCoarseTimeCPloopIn ', iLoop,cEnding(iEnding)//' step!|Init&FinPresDefect: ',MaxInitialPressureDefect0,MaxInitialPressureDefect,MaxInitialPressureDefect/MaxInitialPressureDefect0
  bEXIT = .true.
 end if
! if (bExit) GOTO 100

end do !iLoop

100 continue

! DO iStep = 0,mySolSeq%nSteps
! if (myid.ne.master) then
!  QuadSc%ValU = Sol3%S(iStep)%U(1:Sol3%ndofU)
!  QuadSc%ValV = Sol3%S(iStep)%V(1:Sol3%ndofU)
!  QuadSc%ValW = Sol3%S(iStep)%W(1:Sol3%ndofU)
!  LinSc%ValP(NLMAX)%x = Sol3%S(iStep)%P(1:Sol3%ndofP)
!  
!  LinSc%P_new = LinSc%ValP(NLMAX)%x
!  CALL QuadScP1ExtPoltoQ2(LinSc,QuadSc)
! end if
! CALL Output_Profiles(iStep)
! end do
! pause

if (mySolSeq%nFOuter.eq.mySolSeq%nOuter) goto 5

!!!! Extrapolation !!!!!!!
do iLoop =1 , 1 

 DO iOuter=mySolSeq%nFOuter,mySolSeq%nOuter

  nSteps = 2**(iOuter-1)
  lStep  = mySolSeq%nSteps/nSteps
  tstep  = tstep_BU/DBLE(nSteps)

  CALL NonLin_SequentialBurgerStep_ParT()

  CALL PressureStepCP_ParT()

 END DO !iOuter

 tstep  = tstep_BU

end do !iLoop

5 continue

IF (myid.ne.master) THEN
 QuadSc%ValU = mySolSeq%S(mySolSeq%nSteps)%U(1:QuadSc%ndof)
 QuadSc%ValV = mySolSeq%S(mySolSeq%nSteps)%V(1:QuadSc%ndof)
 QuadSc%ValW = mySolSeq%S(mySolSeq%nSteps)%W(1:QuadSc%ndof)

 LinSc%P_new = 1.5d0*mySolSeq%S(mySolSeq%nSteps-0)%P(1:LinSc%ndof)  - 0.5d0*mySolSeq%S(mySolSeq%nSteps-1)%P(1:LinSc%ndof)
end if
!LinSc%P_new = mySolSeq%S(mySolSeq%nSteps-0)%P
CALL QuadScP1ExtPoltoQ2(LinSc,QuadSc)

CALL FAC_GetForces(mfile)

RETURN

 CONTAINS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE InitSolutions(Sol,iCase)
integer, intent(in) :: iCase
TYPE(tSolSeq) :: Sol

integer iC,iS,iE,indice,nU,nP
real*8 daux

ILEV = NLMAX-iCase
nP = 4*KNEL(ILEV)
nU = KNEL(ILEV)+KNET(ILEV)+KNAT(ILEV)+KNVT(ILEV)

Sol%ndofP = nP
Sol%ndofU = nU
Sol%nOuter  = Properties%nTPSubSteps
Sol%nFOuter = Properties%nTPFSubSteps
Sol%nSteps  = 2**(Sol%nOuter-1)

if (.not.allocated(PressureDefect)) allocate(PressureDefect(Sol%nSteps))

IF (myid.ne.master) then
 if (.not.allocated(Sol%S)) THEN
  ALLOCATE(Sol%S(0:Sol%nSteps))
 end if

 DO iStep = 0,Sol%nSteps
  if (.not.allocated(Sol%S(iStep)%U)) ALLOCATE(Sol%S(iStep)%U(nU))
  if (.not.allocated(Sol%S(iStep)%V)) ALLOCATE(Sol%S(iStep)%V(nU))
  if (.not.allocated(Sol%S(iStep)%W)) ALLOCATE(Sol%S(iStep)%W(nU))
  if (.not.allocated(Sol%S(iStep)%P)) ALLOCATE(Sol%S(iStep)%P(nP))
 
  if (.not.allocated(Sol%S(iStep)%U_aux)) ALLOCATE(Sol%S(iStep)%U_aux(nU))
  if (.not.allocated(Sol%S(iStep)%V_aux)) ALLOCATE(Sol%S(iStep)%V_aux(nU))
  if (.not.allocated(Sol%S(iStep)%W_aux)) ALLOCATE(Sol%S(iStep)%W_aux(nU))
  if (.not.allocated(Sol%S(iStep)%P_aux)) ALLOCATE(Sol%S(iStep)%P_aux(nP))
 
  Sol%S(iStep)%U = QuadSc%ValU(1:nU)
  Sol%S(iStep)%V = QuadSc%ValV(1:nU)
  Sol%S(iStep)%W = QuadSc%ValW(1:nU)
  
  if (icase.eq.0) then
   Sol%S(iStep)%P = LinSc%ValP(NLMAX)%x(1:nP)
  else
   do iC = 1,4
    do iE=1,nP/4
     daux = LinSc%ValP(NLMAX)%x(4*(iE-1)+(iC-1)+1)
     do iS = 1,7
      indice = nP + 4*7*(iE-1) + 4*(iS-1) + (iC-1) + 1 
      daux = daux + LinSc%ValP(NLMAX)%x(indice)
     end do
    Sol%S(iStep)%P(4*(iE-1)+(iC-1)+1) = 0.125d0*daux
    end do
   end do
  end if
 END DO

 IF (itns.eq.1.and.istart.eq.0) then
  Sol%S(0)%U = 0d0
  Sol%S(0)%V = 0d0
  Sol%S(0)%W = 0d0
 end if
END IF

END SUBROUTINE InitSolutions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE InterpolateVelocity_ParT(IntPol)
INTEGER IntPol

  IF (myid.ne.master) THEN
  
    DO iStep = 1,nSteps
      
     IPS%Actual = (iStep-0)*lStep
     IPS%IntPol1 = (iStep-0)*lStep
     IPS%IntPol2 = (iStep-0)*lStep-int(lStep/2)
     
     IPS%Src1 = (iStep-1)*lStep
     IPS%Src2 = (iStep-0)*lStep
     
     IF (IntPol.eq.1) then
      
      if (myid.eq.showid) WRITE(*,'(A,I0,4(A,I0))') &
      "ActualStep: ",IPS%Actual, ' IntPolSourceSteps: ',IPS%Src1,',',IPS%Src2,' UpdatedSteps: ',IPS%IntPol1,',',IPS%IntPol2
      mySolSeq%S(IPS%IntPol2)%U = 0.5d0*mySolSeq%S(IPS%Src1)%U(1:QuadSc%ndof) + 0.5d0*mySolSeq%S(IPS%Src2)%U(1:QuadSc%ndof)
      mySolSeq%S(IPS%IntPol2)%V = 0.5d0*mySolSeq%S(IPS%Src1)%V(1:QuadSc%ndof) + 0.5d0*mySolSeq%S(IPS%Src2)%V(1:QuadSc%ndof)
      mySolSeq%S(IPS%IntPol2)%W = 0.5d0*mySolSeq%S(IPS%Src1)%W(1:QuadSc%ndof) + 0.5d0*mySolSeq%S(IPS%Src2)%w(1:QuadSc%ndof)

     END IF
     
     IF (IntPol.eq.0) then
      if (myid.eq.showid) WRITE(*,'(A,I0,A,i0)') "steps: from: ",IPS%Src2,' to: ', IPS%IntPol2
      mySolSeq%S(IPS%IntPol2)%U = mySolSeq%S(IPS%Src2)%U(1:QuadSc%ndof)
      mySolSeq%S(IPS%IntPol2)%V = mySolSeq%S(IPS%Src2)%V(1:QuadSc%ndof)
      mySolSeq%S(IPS%IntPol2)%W = mySolSeq%S(IPS%Src2)%W(1:QuadSc%ndof)
     END IF
     
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
     mySolSeq%S((iStep-0)*lStep)%P_aux = mySolSeq%S((iStep-0)*lStep)%P(1:LinSc%ndof)
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
      mySolSeq%S(IPS%IntPol1)%P = IPS%Frac1*mySolSeq%S(IPS%Src1)%P_aux(1:LinSc%ndof) + IPS%Frac2*mySolSeq%S(IPS%Src2)%P_aux(1:LinSc%ndof)
       
      IPS%Frac1 = (IPS%X2Dest - IPS%X2Src )/(IPS%X1Src - IPS%X2Src)
      IPS%Frac2 = (IPS%X1Src  - IPS%X2Dest)/(IPS%X1Src - IPS%X2Src)
      mySolSeq%S(IPS%IntPol2)%P = IPS%Frac1*mySolSeq%S(IPS%Src1)%P_aux(1:LinSc%ndof) + IPS%Frac2*mySolSeq%S(IPS%Src2)%P_aux(1:LinSc%ndof)

      if (myid.eq.showid) WRITE(*,'(A,I0,4(A,I0))') &
      "ActualStep: ",IPS%Actual, ' IntPolSourceSteps: ',IPS%Src1,',',IPS%Src2,' UpdatedSteps: ',IPS%IntPol1,',',IPS%IntPol2
      if (myid.eq.showid) WRITE(*,'(2(A,I0),2(A,F6.3))') 'IntFractions:', IPS%iX1Src, ', ', IPS%iX2Src,' SrcFractions:', IPS%X1Src, ', ', IPS%X2Src
      if (myid.eq.showid) WRITE(*,'(2(A,I0),2(A,F6.3))') 'IntFractions:', IPS%IntPol1, ', ', IPS%IntPol2,' DestFractions:', IPS%X1Dest, ', ', IPS%X2Dest
     END IF

     IF (IntPol.eq.0) then
      if (myid.eq.showid) WRITE(*,*) "steps: ",(iStep-0)*lStep,(iStep-1)*lStep,(iStep-0)*lStep-int(lStep/2)
      mySolSeq%S((iStep-0)*lStep-int(lStep/2))%P = mySolSeq%S((iStep-0)*lStep)%P(1:LinSc%ndof)
     END IF
     
    END DO
  END IF
    
END SUBROUTINE InterpolatePressure_ParT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE NonLin_BurgerStep_SingleT()

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

IF (myid.ne.master) THEN

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

END SUBROUTINE NonLin_BurgerStep_SingleT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE NonLin_SimultanBurgerStep_ParT()
REAL*8 u_rel_max(6),u_def0(3),u_def0_max(3),u_def0_max_Val
integer :: nLinBurgers=128, iLinBurgers

do iLinBurgers=1,nLinBurgers

u_rel_max = -1d0
u_def0_max = -1d0

IF (myid.ne.master) THEN
 DO iStep = 0,nSteps
  mySolSeq%S((iStep-0)*lStep)%U_aux = mySolSeq%S((iStep-0)*lStep)%U(1:QuadSc%ndof)
  mySolSeq%S((iStep-0)*lStep)%V_aux = mySolSeq%S((iStep-0)*lStep)%V(1:QuadSc%ndof)
  mySolSeq%S((iStep-0)*lStep)%W_aux = mySolSeq%S((iStep-0)*lStep)%W(1:QuadSc%ndof)
 end do
 DO iStep = 0,nSteps
  mySolSeq%S((iStep-0)*lStep)%U(1:QuadSc%ndof) = 0d0
  mySolSeq%S((iStep-0)*lStep)%V(1:QuadSc%ndof) = 0d0
  mySolSeq%S((iStep-0)*lStep)%W(1:QuadSc%ndof) = 0d0
 end do
end if

DO iStep = 1,nSteps

thstep = tstep*(1d0-theta)

IF (myid.ne.master) THEN
 QuadSc%ValU = mySolSeq%S((iStep-1)*lStep)%U_aux(1:QuadSc%ndof)
 QuadSc%ValV = mySolSeq%S((iStep-1)*lStep)%V_aux(1:QuadSc%ndof)
 QuadSc%ValW = mySolSeq%S((iStep-1)*lStep)%W_aux(1:QuadSc%ndof)
 LinSc%ValP(NLMAX)%x = mySolSeq%S((iStep-0)*lStep)%P(1:LinSc%ndof)
end if

! Set dirichlet boundary conditions on the solution
CALL Boundary_QuadScalar_Val()

! write(*,*) size(LinSc%ValP(NLMAX)%x),size(mySolSeq%S(0)%P),size(QuadSc%ValU),size(mySolSeq%S(0)%U)
! pause
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
 QuadSc%ValU = mySolSeq%S((iStep-0)*lStep)%U_aux(1:QuadSc%ndof)
 QuadSc%ValV = mySolSeq%S((iStep-0)*lStep)%V_aux(1:QuadSc%ndof)
 QuadSc%ValW = mySolSeq%S((iStep-0)*lStep)%W_aux(1:QuadSc%ndof)
END IF

! Set dirichlet boundary conditions on the solution
CALL Boundary_QuadScalar_Val()

CALL OperatorRegenaration(3)

IF (myid.ne.master) THEN

 ! Assemble the defect vector and fine level matrix
 CALL Matdef_General_QuadScalar(QuadSc,-1)
 
 ! Assemble the defect vector and fine level matrix
 IF (iStep.gt.1) then
  CALL AddTimeNewtonAccelerator(QuadSc,&
  mySolSeq%S((iStep-1)*lStep)%U,&
  mySolSeq%S((iStep-1)*lStep)%V,&
  mySolSeq%S((iStep-1)*lStep)%W)
 end if

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 QuadSc%auxU = QuadSc%defU
 QuadSc%auxV = QuadSc%defV
 QuadSc%auxW = QuadSc%defW
 CALL E013Sum3(QuadSc%auxU,QuadSc%auxV,QuadSc%auxW)

 ! Save the old solution
!  CALL LCP1(QuadSc%valU,QuadSc%valU_old,QuadSc%ndof)
!  CALL LCP1(QuadSc%valV,QuadSc%valV_old,QuadSc%ndof)
!  CALL LCP1(QuadSc%valW,QuadSc%valW_old,QuadSc%ndof)

 ! Compute the norm of the defect
 CALL Resdfk_General_QuadScalar(QuadSc,ResU,ResV,ResW,DefUVW,RhsUVW)
 
 CALL MeasureDefectNorms(QuadSc,u_def0)

END IF

CALL COMM_MaximumN(u_def0,3)
u_def0_max = max(u_def0_max,u_def0)

CALL COMM_Maximum(RhsUVW)
DefUVWCrit=MAX(RhsUVW*QuadSc%prm%defCrit,QuadSc%prm%MinDef)

CALL Protocol_QuadScalar(mfile,QuadSc,0,&
     ResU,ResV,ResW,DefUVW,DefUVWCrit," Momentum equation ")

CALL ZTIME(tttt1)
myStat%tDefUVW = myStat%tDefUVW + (tttt1-tttt0)

INL=1

QuadSc%valU_old = 0d0
QuadSc%valV_old = 0d0
QuadSc%valW_old = 0d0

! ! Calling the solver
CALL Solve_General_QuadScalar(QuadSc,Boundary_QuadScalar_Val,&
Boundary_QuadScalar_Mat,Boundary_QuadScalar_Mat_9,mfile)

!!!!          Checking the quality of the result           !!!!
!!!! ----------------------------------------------------- !!!!

u_rel_max = max(u_rel_max,QuadSc%prm%MGprmOut(1)%u_rel)

IF (myid.ne.master) THEN
 mySolSeq%S((iStep-0)*lStep)%U = QuadSc%ValU(1:QuadSc%ndof)
 mySolSeq%S((iStep-0)*lStep)%V = QuadSc%ValV(1:QuadSc%ndof)
 mySolSeq%S((iStep-0)*lStep)%W = QuadSc%ValW(1:QuadSc%ndof)
END IF

END DO ! iStep

if (myid.eq.1) write(MTERM,"(104('X'))") 
if (myid.eq.1) write(mFILE,"(104('X'))") 
if (myid.eq.1) write(MTERM,'(2(A,I4),3ES12.4)') 'RelativeChangesOfU:', iLoop,'|',iLinBurgers,u_rel_max(4)/u_rel_max(1),u_rel_max(5)/u_rel_max(2),u_rel_max(6)/u_rel_max(3)
if (myid.eq.1) write(MFILE,'(2(A,I4),3ES12.4)') 'RelativeChangesOfU:', iLoop,'|',iLinBurgers,u_rel_max(4)/u_rel_max(1),u_rel_max(5)/u_rel_max(2),u_rel_max(6)/u_rel_max(3)
if (myid.eq.1) write(MTERM,'(2(A,I4),3ES12.4)') 'MaxInitDefectsOfU: ', iLoop,'|',iLinBurgers,u_def0_max(1),u_def0_max(2),u_def0_max(3)
if (myid.eq.1) write(MFILE,'(2(A,I4),3ES12.4)') 'MaxInitDefectsOfU: ', iLoop,'|',iLinBurgers,u_def0_max(1),u_def0_max(2),u_def0_max(3)
if (myid.eq.1) write(MTERM,"(104('X'))") 
if (myid.eq.1) write(mFILE,"(104('X'))") 

! Update the solution
DO iStep = 0,nSteps
 IF (myid.ne.0) THEN
  QuadSc%ValU(1:QuadSc%ndof) = mySolSeq%S((iStep-0)*lStep)%U_aux(1:QuadSc%ndof) + mySolSeq%S((iStep-0)*lStep)%U(1:QuadSc%ndof)
  QuadSc%ValV(1:QuadSc%ndof) = mySolSeq%S((iStep-0)*lStep)%V_aux(1:QuadSc%ndof) + mySolSeq%S((iStep-0)*lStep)%V(1:QuadSc%ndof)
  QuadSc%ValW(1:QuadSc%ndof) = mySolSeq%S((iStep-0)*lStep)%W_aux(1:QuadSc%ndof) + mySolSeq%S((iStep-0)*lStep)%W(1:QuadSc%ndof)
  
  ! Set dirichlet boundary conditions on the solution
  CALL Boundary_QuadScalar_Val()
  
  mySolSeq%S((iStep-0)*lStep)%U(1:QuadSc%ndof) = QuadSc%ValU(1:QuadSc%ndof)
  mySolSeq%S((iStep-0)*lStep)%V(1:QuadSc%ndof) = QuadSc%ValV(1:QuadSc%ndof)
  mySolSeq%S((iStep-0)*lStep)%W(1:QuadSc%ndof) = QuadSc%ValW(1:QuadSc%ndof)
 END IF
END DO

u_def0_max_Val=max(u_def0_max(1),u_def0_max(2),u_def0_max(3))

IF (u_def0_max_Val.lt.1e-15) exit

END DO ! iLinBurgers

END SUBROUTINE NonLin_SimultanBurgerStep_ParT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE NonLin_SequentialBurgerStep_ParT()
REAL*8 u_rel_max(6)

u_rel_max = -1d0

DO iStep = 1,nSteps

thstep = tstep*(1d0-theta)

IF (myid.ne.master) THEN
 QuadSc%ValU = mySolSeq%S((iStep-1)*lStep)%U(1:QuadSc%ndof)
 QuadSc%ValV = mySolSeq%S((iStep-1)*lStep)%V(1:QuadSc%ndof)
 QuadSc%ValW = mySolSeq%S((iStep-1)*lStep)%W(1:QuadSc%ndof)
 LinSc%ValP(NLMAX)%x = mySolSeq%S((iStep-0)*lStep)%P(1:LinSc%ndof)
end if

! Set dirichlet boundary conditions on the solution
CALL Boundary_QuadScalar_Val()

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

! if (myid.eq.1) write(*,*) 'A0', size(LinSc%ValP(NLMAX)%x),size(mySolSeq%S(iStep)%P),size(QuadSc%ValU),size(mySolSeq%S(iStep)%U)
IF (myid.ne.master) THEN
 QuadSc%ValU = mySolSeq%S((iStep-0)*lStep)%U(1:QuadSc%ndof)
 QuadSc%ValV = mySolSeq%S((iStep-0)*lStep)%V(1:QuadSc%ndof)
 QuadSc%ValW = mySolSeq%S((iStep-0)*lStep)%W(1:QuadSc%ndof)
END IF
! if (myid.eq.1) write(*,*) 'A1', size(LinSc%ValP(NLMAX)%x),size(mySolSeq%S(iStep)%P),size(QuadSc%ValU),size(mySolSeq%S(iStep)%U)

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

! if (myid.eq.1) write(*,*) 'Z0', size(LinSc%ValP(NLMAX)%x),size(mySolSeq%S(iStep)%P),size(QuadSc%ValU),size(mySolSeq%S(iStep)%U)
! ! Calling the solver
CALL Solve_General_QuadScalar(QuadSc,Boundary_QuadScalar_Val,&
Boundary_QuadScalar_Mat,Boundary_QuadScalar_Mat_9,mfile)

u_rel_max = max(u_rel_max,QuadSc%prm%MGprmOut(1)%u_rel)

! if (myid.eq.1) write(*,*) 'Z1', size(LinSc%ValP(NLMAX)%x),size(mySolSeq%S(iStep)%P),size(QuadSc%ValU),size(mySolSeq%S(iStep)%U)

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

! if (myid.eq.1) write(*,*) 'X0', size(LinSc%ValP(NLMAX)%x),size(mySolSeq%S(iStep)%P),size(QuadSc%ValU),size(mySolSeq%S(iStep)%U)
IF (myid.ne.master) THEN
 mySolSeq%S((iStep-0)*lStep)%U = QuadSc%ValU(1:QuadSc%ndof)
 mySolSeq%S((iStep-0)*lStep)%V = QuadSc%ValV(1:QuadSc%ndof)
 mySolSeq%S((iStep-0)*lStep)%W = QuadSc%ValW(1:QuadSc%ndof)
END IF
! if (myid.eq.1) write(*,*) 'X1', size(LinSc%ValP(NLMAX)%x),size(mySolSeq%S(iStep)%P),size(QuadSc%ValU),size(mySolSeq%S(iStep)%U)

END DO ! iStep

if (myid.eq.1) write(MTERM,'(A,I4,3ES12.4)') 'RelativeChangesOfU:', iLoop,u_rel_max(4)/u_rel_max(1),u_rel_max(5)/u_rel_max(2),u_rel_max(6)/u_rel_max(3)
if (myid.eq.1) write(MFILE,'(A,I4,3ES12.4)') 'RelativeChangesOfU:', iLoop,u_rel_max(4)/u_rel_max(1),u_rel_max(5)/u_rel_max(2),u_rel_max(6)/u_rel_max(3)

END SUBROUTINE NonLin_SequentialBurgerStep_ParT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE Lin_BurgerStep_ParT()
integer iNonLin

! do iNonLin = 1,2

IF (myid.ne.master) THEN
 DO iStep = 0,nSteps
  mySolSeq%S(iStep)%U_aux = mySolSeq%S(iStep)%U(1:QuadSc%ndof)
  mySolSeq%S(iStep)%V_aux = mySolSeq%S(iStep)%V(1:QuadSc%ndof)
  mySolSeq%S(iStep)%W_aux = mySolSeq%S(iStep)%W(1:QuadSc%ndof)
 end do
end if

DO iStep = 1,nSteps

thstep = tstep*(1d0-theta)

IF (myid.ne.master) THEN
 QuadSc%ValU = mySolSeq%S((iStep-1)*lStep)%U_aux(1:QuadSc%ndof)
 QuadSc%ValV = mySolSeq%S((iStep-1)*lStep)%V_aux(1:QuadSc%ndof)
 QuadSc%ValW = mySolSeq%S((iStep-1)*lStep)%W_aux(1:QuadSc%ndof)
 LinSc%ValP(NLMAX)%x = mySolSeq%S((iStep-0)*lStep)%P(1:LinSc%ndof)
end if

! Set dirichlet boundary conditions on the solution
CALL Boundary_QuadScalar_Val()

! write(*,*) size(LinSc%ValP(NLMAX)%x),size(mySolSeq%S(0)%P),size(QuadSc%ValU),size(mySolSeq%S(0)%U)
! pause
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
 QuadSc%ValU = mySolSeq%S((iStep-0)*lStep)%U_aux(1:QuadSc%ndof)
 QuadSc%ValV = mySolSeq%S((iStep-0)*lStep)%V_aux(1:QuadSc%ndof)
 QuadSc%ValW = mySolSeq%S((iStep-0)*lStep)%W_aux(1:QuadSc%ndof)
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

INL=1
! 
! DO INL=1,2!QuadSc%prm%NLmax
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
 mySolSeq%S((iStep-0)*lStep)%U = QuadSc%ValU(1:QuadSc%ndof)
 mySolSeq%S((iStep-0)*lStep)%V = QuadSc%ValV(1:QuadSc%ndof)
 mySolSeq%S((iStep-0)*lStep)%W = QuadSc%ValW(1:QuadSc%ndof)
END IF

END DO ! iStep

! end do !iNonLin

END SUBROUTINE Lin_BurgerStep_ParT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE PressureStep_SingleT()

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

MaxInitialPressureDefect = max(MaxInitialPressureDefect,LinSc%prm%MGprmOut%DefInitial*tstep)

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

END SUBROUTINE PressureStep_SingleT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE PressureStepCP_ParT()

DO iStep = 1,nSteps

IF (myid.ne.master) THEN

 CALL ZTIME(tttt0)
 ! Save the old solution
 ! LinSc%valP_old = LinSc%valP(NLMAX)%x
 LinSc%valP_old = mySolSeq%S((iStep-0)*lStep)%P(1:LinSc%ndof)
 LinSc%valP(NLMAX)%x = 0d0

 QuadSc%ValU = mySolSeq%S((iStep)*lStep)%U(1:QuadSc%ndof) - mySolSeq%S((iStep-1)*lStep)%U(1:QuadSc%ndof) 
 QuadSc%ValV = mySolSeq%S((iStep)*lStep)%V(1:QuadSc%ndof) - mySolSeq%S((iStep-1)*lStep)%V(1:QuadSc%ndof)
 QuadSc%ValW = mySolSeq%S((iStep)*lStep)%W(1:QuadSc%ndof) - mySolSeq%S((iStep-1)*lStep)%W(1:QuadSc%ndof)
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
PressureDefect(iStep) = LinSc%prm%MGprmOut%DefInitial*tstep

CALL Protocol_LinScalar(mfile,LinSc," Pressure-Poisson equation")

2 CONTINUE

IF (myid.ne.0) THEN
 CALL ZTIME(tttt0)
 CALL Pressure_Correction()
 CALL ZTIME(tttt1)
 myStat%tCorrUVWP = myStat%tCorrUVWP + (tttt1-tttt0)
END IF
 
IF (myid.ne.master) THEN
 mySolSeq%S((iStep-0)*lStep)%P = LinSc%ValP(NLMAX)%x(1:LinSc%ndof)
END IF

END DO ! iStep

END SUBROUTINE PressureStepCP_ParT

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

old_TSTEP = tstep
tstep = xTimeStepMultiplier*tstep

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
    QuadSc%prm%MGprmIn%nSmootherSteps = MIN(INT(DBLE(QuadSc%prm%MGprmIn%nSmootherSteps)*1.25d0),MaxSmootheningSteps)
    if (myid.eq.1) WRITE(*,*) 'MG smoothening step increasement to:', QuadSc%prm%MGprmIn%nSmootherSteps
   END IF
   
   !!!!!!!!!!!!!!!!!!!!!!!!! Artificial TimeStep Control !!!!!!!!!!!!!!!!!!!!!!!!!
   IF (QuadSc%prm%MGprmIn%nSmootherSteps.ge.MaxSmootheningSteps) THEN
    IF (TimeStepIncrease.ge.TimeStepIncreaseCrit) then
     xTimeStepMultiplier = MAX(0.066667d0,0.2d0*xTimeStepMultiplier)
     TimeStepIncrease = 0
     if (myid.eq.1) WRITE(MTERM,*) 'Artificial TimeStep Reduction:', xTimeStepMultiplier
     if (myid.eq.1) WRITE(MFILE,*) 'Artificial TimeStep Reduction:', xTimeStepMultiplier
     QuadSc%prm%MGprmIn%nSmootherSteps = INT(DBLE(QuadSc%prm%MGprmIn%nSmootherSteps)*0.8d0)
     if (myid.eq.1) WRITE(*,*) 'MG smoothening step reduction to:', QuadSc%prm%MGprmIn%nSmootherSteps
     if (myid.eq.1) WRITE(*,*) 'MG smoothening step reduction to:', QuadSc%prm%MGprmIn%nSmootherSteps
    else
     TimeStepIncrease = TimeStepIncrease + 1
    END IF
   END IF
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
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

if (bMultiMat) then
 CALL GetAlphaNonNewtViscosity_sse()
else
 CALL GetNonNewtViscosity_sse()
end if

IF (.not.bKTPRelease) then
 CALL Calculate_Torque(mfile)
END IF

if ((.not.bKTPRelease).and.bMultiMat) then
 CALL GetPressureAtInflows(mfile)
end if

! checking convergence of the pressure (worst case scenario)
CALL LL21(LinSc%valP(NLMAX)%x,LinSc%ndof,defP)
call COMM_SUMM(defp)
if (ieee_is_nan(defP)) DivergedSolution = .true.
if (.not.ieee_is_finite(defP)) DivergedSolution = .true.

! checking convergence of the pressure (no improvement of the defect)
if (LinSc%prm%MGprmOut%DefInitial.lt.LinSc%prm%MGprmOut%DefFinal) then
 DivergedSolution = .true.
end if

ILEV = NLMAX
CALL SETLEV(2)
CALL Comm_Solution(QuadSc%ValU,QuadSc%ValV,QuadSc%ValW,QuadSc%auxU,QuadSc%ndof)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ConvergenceCheck !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (bMultiMat) then


end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ConvergenceCheck !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

tstep = old_TSTEP 

end subroutine Transport_q2p1_UxyzP_sse
!
! ----------------------------------------------
!
subroutine Transport_q2p1_UxyzP_sse_PF(mfile, inl_u, itns)

use, intrinsic :: ieee_arithmetic

INTEGER mfile,INL,inl_u,itns
REAL*8  ResU,ResV,ResW,DefUVW,RhsUVW,DefUVWCrit
REAL*8  ResP,DefP,RhsPG,defPG,defDivU,DefPCrit
INTEGER INLComplete,I,J,IERR,iITER

old_TSTEP = tstep
tstep = xTimeStepMultiplier*tstep

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

 ! Add the gravity force to the rhs
 CALL AddGravForce()
 
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
    QuadSc%prm%MGprmIn%nSmootherSteps = MIN(INT(DBLE(QuadSc%prm%MGprmIn%nSmootherSteps)*1.25d0),MaxSmootheningSteps)
    if (myid.eq.1) WRITE(*,*) 'MG smoothening step increasement to:', QuadSc%prm%MGprmIn%nSmootherSteps
   END IF
  END IF

  !!!!!!!!!!!!!!!!!!!!!!!!!! Artificial TimeStep Control !!!!!!!!!!!!!!!!!!!!!!!!!
  IF (QuadSc%prm%MGprmIn%nSmootherSteps.ge.MaxSmootheningSteps) THEN
   IF (TimeStepIncrease.ge.TimeStepIncreaseCrit) then
    xTimeStepMultiplier = MAX(0.066667d0,0.2d0*xTimeStepMultiplier)
    TimeStepIncrease = 0
    QuadSc%prm%MGprmIn%nSmootherSteps = INT(DBLE(QuadSc%prm%MGprmIn%nSmootherSteps)*0.8d0)
    if (myid.eq.1) WRITE(*,*) 'MG smoothening step reduction to:', QuadSc%prm%MGprmIn%nSmootherSteps
    if (myid.eq.1) WRITE(*,*) 'MG smoothening step reduction to:', QuadSc%prm%MGprmIn%nSmootherSteps
    if (myid.eq.1) WRITE(MTERM,*) 'Artificial TimeStep Reduction:', xTimeStepMultiplier
    if (myid.eq.1) WRITE(MFILE,*) 'Artificial TimeStep Reduction:', xTimeStepMultiplier
   else
    TimeStepIncrease = TimeStepIncrease + 1
   END IF
  ELSE
   TimeStepIncrease = 0
  END IF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
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

if (bMultiMat) then
 CALL GetAlphaNonNewtViscosity_sse()
else
 CALL GetNonNewtViscosity_sse()
end if

IF (.not.bKTPRelease) then
 CALL Calculate_Torque(mfile)
END IF

if ((.not.bKTPRelease).and.bMultiMat) then
 CALL GetPressureAtInflows(mfile)
end if

CALL LL21(LinSc%valP(NLMAX)%x,LinSc%ndof,defP)
call COMM_SUMM(defp)
if (ieee_is_nan(defP)) DivergedSolution = .true.
if (.not.ieee_is_finite(defP)) DivergedSolution = .true.

ILEV = NLMAX
CALL SETLEV(2)
CALL Comm_Solution(QuadSc%ValU,QuadSc%ValV,QuadSc%ValW,QuadSc%auxU,QuadSc%ndof)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ConvergenceCheck !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (bMultiMat) then


end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ConvergenceCheck !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

tstep = old_TSTEP 

end subroutine Transport_q2p1_UxyzP_sse_PF
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

SUBROUTINE GetPressureAtInflows(mfile)
integer mfile
!---------------------------------
INTEGER NeighA(4,6),iInflow
REAL*8 dPress,dArea,dAN(3),stdNorm
DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4,ivt5,iMem
INTEGER :: nMem = 20
REAL*8, allocatable :: dIntPress(:),dIntArea(:)

if (.not.allocated(dIntPress)) allocate(dIntPress(myProcess%nOfInflows))
if (.not.allocated(dIntArea))  allocate(dIntArea(myProcess%nOfInflows))

dIntPress = 0d0
dIntArea = 0d0

DO iInflow=1,myProcess%nOfInflows
 k=1
 DO i=1,mg_mesh%level(nlmax)%nel
  DO j=1,6
   IF (k.eq.mg_mesh%level(nlmax)%karea(j,i)) THEN
    IF (ABS(myBoundary%iInflow(mg_mesh%level(nlmax)%nvt+mg_mesh%level(nlmax)%net+k)).eq.iInflow) THEN
     ivt1 = mg_mesh%level(nlmax)%kvert(NeighA(1,j),i)
     ivt2 = mg_mesh%level(nlmax)%kvert(NeighA(2,j),i)
     ivt3 = mg_mesh%level(nlmax)%kvert(NeighA(3,j),i)
     ivt4 = mg_mesh%level(nlmax)%kvert(NeighA(4,j),i)
     ivt5 = mg_mesh%level(nlmax)%nvt + mg_mesh%level(nlmax)%net + mg_mesh%level(nlmax)%nat + i
     
     if (shell(ivt5).gt.0d0) THEN  ! Only faces which belong to inner (no-fbm) elements 
      CALL GET_NormalArea(mg_mesh%level(nlmax)%dcorvg(1:3,ivt1),&
                          mg_mesh%level(nlmax)%dcorvg(1:3,ivt2),&
                          mg_mesh%level(nlmax)%dcorvg(1:3,ivt3),&
                          mg_mesh%level(nlmax)%dcorvg(1:3,ivt4),&
                          mg_mesh%level(nlmax)%dcorvg(1:3,ivt5),dAN)
      dPress = 0.25d0*(LinSc%Q2(ivt1) + LinSc%Q2(ivt2) + LinSc%Q2(ivt3) + LinSc%Q2(ivt4))
      dArea = sqrt(dAN(1)**2d0 + dAN(2)**2d0 + dAN(3)**2d0)
!       write(*,*) dPress, dArea
      dIntArea(iInflow)  = dIntArea(iInflow)  + dArea
      dIntPress(iInflow) = dIntPress(iInflow) + dArea*dPress
     end if
    END IF
    k = k + 1
   END IF
  END DO
 END DO
END DO 
 
CALL COMM_SUMMn(dIntArea,myProcess%nOfInflows)
CALL COMM_SUMMn(dIntPress,myProcess%nOfInflows)

mySetup%bPressureConvergence = .true.
DO iInflow=1,myProcess%nOfInflows

 if (.not.allocated(myProcess%myInflow(iInflow)%PressureEvolution)) then 
  allocate(myProcess%myInflow(iInflow)%PressureEvolution(0:nMem))
  myProcess%myInflow(iInflow)%PressureEvolution = 0d0
 end if
 DO iMem = nMem-1,1,-1
  myProcess%myInflow(iInflow)%PressureEvolution(iMem+1) = myProcess%myInflow(iInflow)%PressureEvolution(iMem)
 end do
 myProcess%myInflow(iInflow)%PressureEvolution(1) = 1e-6*dIntPress(iInflow)/dIntArea(iInflow)

 ! TimeAveragedPressureComputation
 myProcess%myInflow(iInflow)%PressureEvolution(0) = 0d0
 DO iMem = 1,nMem
  myProcess%myInflow(iInflow)%PressureEvolution(0) = myProcess%myInflow(iInflow)%PressureEvolution(0) + myProcess%myInflow(iInflow)%PressureEvolution(iMem)
 END DO
 myProcess%myInflow(iInflow)%PressureEvolution(0) = myProcess%myInflow(iInflow)%PressureEvolution(0) / DBLE(nMem)

 ! Standard deviation for the last nMem steps
 stdNorm = 0d0
 DO iMem = 1,nMem
  stdNorm = stdNorm + (myProcess%myInflow(iInflow)%PressureEvolution(iMem)-myProcess%myInflow(iInflow)%PressureEvolution(0))**2d0
 END DO
 stdNorm = stdNorm / DBLE(nMem)
 stdNorm = stdNorm**0.5d0
 stdNorm = 1d2*stdNorm/myProcess%myInflow(iInflow)%PressureEvolution(0)
 
 myProcess%myInflow(iInflow)%PressureEvolution(0) = stdNorm !! percentage error for the last nMem steps
 
!  if (myid.eq.1) write(mfile,'(A,100ES12.4)') 'test :',myProcess%myInflow(iInflow)%PressureEvolution(0:nMem)
!  if (myid.eq.1) write(mterm,'(A,100ES12.4)') 'test :',myProcess%myInflow(iInflow)%PressureEvolution(0:nMem)
 
 if (itns.gt.nMem.and.istart.eq.1) then
  if (myProcess%myInflow(iInflow)%PressureEvolution(0).gt.mySetup%PressureConvergenceTolerance) mySetup%bPressureConvergence = .false.
  if (myid.eq.1) write(mfile,'(A,I0,A,4ES12.4)') 'Pressure[bar]AtInflowArea[cm2]stdErr[%]',iInflow,':',1e-6*dIntPress(iInflow)/dIntArea(iInflow),dIntArea(iInflow),myProcess%myInflow(iInflow)%PressureEvolution(0)
  if (myid.eq.1) write(mterm,'(A,I0,A,4ES12.4)') 'Pressure[bar]AtInflowArea[cm2]stdErr[%]',iInflow,':',1e-6*dIntPress(iInflow)/dIntArea(iInflow),dIntArea(iInflow),myProcess%myInflow(iInflow)%PressureEvolution(0)
 else
  mySetup%bPressureConvergence = .false.
  if (myid.eq.1) write(mfile,'(A,I0,A,4ES12.4)') 'Pressure[bar]AtInflowArea[cm2]stdErr[%]',iInflow,':',1e-6*dIntPress(iInflow)/dIntArea(iInflow),dIntArea(iInflow),1d0
  if (myid.eq.1) write(mterm,'(A,I0,A,4ES12.4)') 'Pressure[bar]AtInflowArea[cm2]stdErr[%]',iInflow,':',1e-6*dIntPress(iInflow)/dIntArea(iInflow),dIntArea(iInflow),1d0
 end if
END DO

mySetup%bPressureConvergence = mySetup%bPressureConvergence.and.(istart.eq.1)

END SUBROUTINE GetPressureAtInflows
