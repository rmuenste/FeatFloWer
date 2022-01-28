MODULE mg_LinScalar

USE PP3D_MPI, ONLY : myid,showID,&
              COMM_Maximum,COMM_SUMM,COMM_NLComplete
USE var_QuadScalar
USE UMFPackSolver, ONLY : myUmfPack_Solve,myUmfPack_Free

IMPLICIT NONE


INTEGER :: iMGL
! -****--------------*-*********-*----------
INTEGER IterCycle,mgLev,CoarseIter
REAL*8 DefNorm
REAL*8 time0,time1


TYPE tJCB
 LOGICAL :: bActivated = .FALSE.
 REAL*8, DIMENSION(:), ALLOCATABLE :: d1
END TYPE tJCB

TYPE (tJCB) :: myJCB

CONTAINS
!
! ----------------------------------------------
!
SUBROUTINE MG_Solver(mfile,mterm)
INTEGER mfile,mterm
INTEGER INLComplete,i
REAL*8 DefI1,DefI2,DefImpr,DDD,AccCoarseIter

! IF (myMG%MaxLev.EQ.myMG%MinLev) THEN
! 
!  RETURN
! END IF

CALL mgInit()

IF (myid.ne.0) CALL mgProlRestInit()

myMG%DefInitial = DefNorm
!  IF (myid.eq.1) WRITE(*,*) "Initial ",DefNorm
IF (DefNorm.LT.1d-32*MyMG%Criterion2) GOTO 88

DefI1 = 0d0
DefI2 = 0d0
AccCoarseIter = 0

DO IterCycle=1,MyMG%MaxIterCycle

 INLComplete = 0
 CALL mg_cycle()
!  pause

 ! Evaluation of new defect
 CALL mgUpdateDefect(myMG%MaxLev,.TRUE.)

 AccCoarseIter = AccCoarseIter + CoarseIter
 DefImpr = DefNorm/myMG%DefInitial
 IF (DefNorm.LE.myMG%DefInitial*MyMG%Criterion1.AND.&
     DefNorm.LE.MyMG%Criterion2.AND.&
     IterCycle.GE.MyMG%MinIterCycle) INLComplete = 1
 IF (IterCycle.GE.MyMG%MaxIterCycle) INLComplete = 1
 CALL COMM_NLComplete(INLComplete)
!  IF (MyMG%cVariable.EQ."Pressure".AND.myid.eq.showid) THEN
!   write(mfile,'(I4,2G12.4,A3,I5,A3,9I4)') &
!   IterCycle,DefNorm,DefNorm/myMG%DefInitial,&
!   " | ",CoarseIter," | ",MyMG%nSmootherSteps
!   write(mterm,'(I4,2G12.4,A3,I5,A3,9I4)') &
!   IterCycle,DefNorm,DefNorm/myMG%DefInitial,&
!   " | ",CoarseIter," | ",MyMG%nSmootherSteps
!  END IF
 IF (INLComplete.eq.1) GOTO 88
 DEFI2 = DEFI1
 DEFI1 = DefNorm

!       pause
END DO

88 CONTINUE

! IF (MyMG%cVariable.EQ."Pressure".AND.myid.eq.0.AND.myMatrixRenewal%C.GE.2) THEN
!   CALL myUmfPack_Free()
! END IF

!  IF (myid.eq.1) WRITE(*,*) "Final ",DefNorm

myMG%DefFinal = DefNorm
myMG%RhoMG1 = (myMG%DefFinal/myMG%DefInitial)**(1d0/DBLE(MAX(IterCycle,1)))
myMG%RhoMG2 = (myMG%DefFinal/DEFI2)**(0.5d0)
MyMG%UsedIterCycle = IterCycle
MyMG%nIterCoarse = INT(AccCoarseIter/MAX(IterCycle,1))

END SUBROUTINE MG_Solver
!
! ----------------------------------------------
!
SUBROUTINE mg_cycle()

IF (MyMG%CycleType.EQ."W")  CALL mg_W_cycle()
IF (MyMG%CycleType.EQ."V")  CALL mg_V_cycle()
IF (MyMG%CycleType.EQ."F")  CALL mg_F_cycle()

END SUBROUTINE mg_cycle
!
! ----------------------------------------------
!
SUBROUTINE mg_W_cycle()
INTEGER imgLev

 CALL mg_down(myMG%MaxLev)

 CALL mg_W_subcycle(myMG%MaxLev-1)

 CALL mg_up(myMG%MaxLev)

END SUBROUTINE mg_W_cycle
!
! ----------------------------------------------
!
RECURSIVE SUBROUTINE mg_W_subcycle(imgLev)
INTEGER imgLev

 IF (imgLev.NE.2) CALL mg_W_subcycle(imgLev-1)

 CALL mg_up(imgLev)
 CALL mg_down(imgLev)

 IF (imgLev.NE.2) CALL mg_W_subcycle(imgLev-1)

END SUBROUTINE mg_W_subcycle
!
! ----------------------------------------------
!
SUBROUTINE mg_F_cycle()
INTEGER imgLev

 CALL mg_down(myMG%MaxLev)

 DO imgLev = myMG%MedLev+1,myMG%MaxLev-1
  CALL mg_up(imgLev)
  CALL mg_down(imgLev)
 END DO

 CALL mg_up(myMG%MaxLev)

END SUBROUTINE mg_F_cycle
!
! ----------------------------------------------
!
SUBROUTINE mg_V_cycle()

 CALL mg_down(myMG%MaxLev)

 CALL mg_up(myMG%MaxLev)

END SUBROUTINE mg_V_cycle
!
! ----------------------------------------------
!
SUBROUTINE mg_down(nimgLev)
INTEGER imgLev,nimgLev

! Presmoothing + restriction
 DO imgLev = nimgLev,myMG%MedLev+1,-1

!  IF (myid.eq.1) write(*,*) "D",imgLev
  mgLev = imgLev
  CALL mgSmoother()                                          ! takes B as RHS and smoothes X

  CALL mgUpdateDefect(imgLev,.FALSE.)                        ! puts B into D and gets the new defect

!   IF (MyMG%cVariable.EQ."Pressure".and.mglev.eq.myMG%MaxLev) THEN
!     write(*,*) "asdasdasd - 4 -- 44",myMG%MaxLev,mgLev
!     CALL outputsol(myMG%D(mgLev)%x,myQ2coor,KWORK(L(KLVERT(mgLev))),KNEL(mgLev),KNVT(mgLev),0)
!   END IF

  mgLev = imgLev - 1
  CALL mgRestriction()                                       ! brings D down by 1 level and stores it as B

  CALL mgGuessSolution()                                     ! choose zero initial vector
 END DO

!  CALL mgSmoother ! Corse Grid Solver
 CALL mgCoarseGridSolver()                                   ! computes linear system with X=(A^-1)B

!   pause
END SUBROUTINE mg_down
!
! ----------------------------------------------
!
SUBROUTINE mg_up(nimgLev)
INTEGER imgLev,nimgLev

 ! Postsmoothing + prolongation
 DO imgLev = myMG%MedLev+1, nimgLev

!  IF (myid.eq.1) write(*,*) "U",imgLev
  mgLev = imgLev

!   IF (MyMG%cVariable.EQ."Pressure".and.mglev.eq.myMG%MaxLev) THEN
!     write(*,*) "asdasdasd - 2 -- 22",myMG%MaxLev,mgLev-1
!     CALL outputsol(myMG%X(mgLev-1)%x,myQ2coor,KWORK(L(KLVERT(mgLev-1))),KNEL(mgLev-1),KNVT(mgLev-1),3)
!   END IF

  CALL mgProlongation()                                      ! brings X up by one level and stores it as AUX

!   IF (MyMG%cVariable.EQ."Pressure".and.mglev.eq.myMG%MaxLev) THEN
!     write(*,*) "asdasdasd - 3 -- 33",myMG%MaxLev,mgLev
!     CALL outputsol(myMG%AUX(mgLev)%x,myQ2coor,KWORK(L(KLVERT(mgLev))),KNEL(mgLev),KNVT(mgLev),2)
!   END IF

  CALL mgUpdateSolution()                                    ! updates solution X = X_old + AUX(update)

!   IF (MyMG%cVariable.EQ."Pressure".and.mglev.eq.myMG%MaxLev) THEN
!     write(*,*) "asdasdasd - 4 -- 44",myMG%MaxLev,mgLev
!     CALL outputsol(myMG%X(mgLev)%x,myQ2coor,KWORK(L(KLVERT(mgLev))),KNEL(mgLev),KNVT(mgLev),0)
!   END IF

  CALL mgSmoother()                                          ! takes B as RHS And smoothes the solution X further

!   IF (MyMG%cVariable.EQ."Pressure".and.mglev.eq.myMG%MaxLev) THEN
!     write(*,*) "asdasdasd - 5 -- 55",myMG%MaxLev,mgLev
!     CALL outputsol(myMG%X(mgLev)%x,myQ2coor,KWORK(L(KLVERT(mgLev))),KNEL(mgLev),KNVT(mgLev),1)
!     pause
!   END IF

 END DO
 
END SUBROUTINE mg_up
!
! ----------------------------------------------
!
SUBROUTINE mgGuessSolution()

IF (mgLev.EQ.myMG%MedLev) RETURN

IF (myid.ne.0) THEN
  myMG%X(mgLev)%x = 0d0
END IF

END SUBROUTINE mgGuessSolution
!
! ----------------------------------------------
!
SUBROUTINE mgUpdateSolution()
INTEGER ndof

IF (myid.ne.0) THEN
 ndof = SIZE(myMG%X(mgLev)%x) 
 CALL LLC1(myMG%AUX(mgLev)%x,myMG%X(mgLev)%x,ndof,1D0,1D0)      ! myMG%X(mgLev)%x = myMG%X(mgLev)%x + myMG%AUX(mgLev)%x
END IF

END SUBROUTINE mgUpdateSolution
!
! ----------------------------------------------
!
SUBROUTINE mgInit()

mgLev = myMG%MaxLev

CALL JCB_Activation()

CALL mgUpdateDefect(mgLev,.TRUE.)

IF (myid.ne.0) THEN
 myMG%AUX(mgLev)%x = 0d0
 myMG%X(mgLev)%x = 0d0
END IF

IterCycle = 0

END SUBROUTINE mgInit
!
! ----------------------------------------------
!
SUBROUTINE mgDefectNorm(II)
INTEGER II,ndof,iFLD,i,j,NN

IF (myid.ne.0) THEN
 mgLev = NLMAX
 ndof  = SIZE(myMG%X(mgLev)%x)
 
 IF (MyMG%cVariable.EQ."Displac") THEN
  CALL LL21(myMG%AUX(II)%x,ndof,DefNorm)
 END IF
 
 IF (MyMG%cVariable.EQ."Temper") THEN !done
  DO iFld=1,myMG%nOfFields
   NN = KNVT(mgLev)
   DO I=1,NN
    J=(iFld-1)*NN + I
    if (myMG%knpr(J).ne.0) myMG%AUX(II)%x(J) = 0d0
   END DO
  END DO

  CALL LL21(myMG%AUX(II)%x,ndof,DefNorm)
!   write(*,*) DefNorm,myid
 END IF
 
END IF

CALL COMM_Maximum(DefNorm)

END SUBROUTINE mgDefectNorm
!
! ----------------------------------------------
!
SUBROUTINE mgUpdateDefect(imgLev,bDef)
USE PP3D_MPI, only: E011SUM
INTEGER i,j,ndof,imgLev,neq,iFld
REAL*8  daux
LOGICAL bDef

IF (myid.ne.0) THEN
 mgLev = imgLev
 ndof  = SIZE(myMG%X(mgLev)%x)
 CALL LCP1(myMG%B(mgLev)%x,myMG%D(mgLev)%x,ndof)
 IF (MyMG%cVariable.EQ."Displac") THEN
  neq = KNVT(mgLev)
  CALL E011_DAX(myMG%X(mgLev)%x,myMG%D(mgLev)%x,-1D0,1D0)
 END IF
 IF (MyMG%cVariable.EQ."Temper") THEN !done
  neq = KNVT(mgLev)
  CALL E011_GenLinSc_Q1_DAX(myMG%X(mgLev)%x,myMG%D(mgLev)%x,-1D0,1D0,MyMG%nOfFields)
 END IF
END IF

IF (bDef) THEN 
 IF (myid.ne.0) THEN
  CALL LCP1(myMG%D(mgLev)%x,myMG%AUX(mgLev)%x,ndof)
  
  IF (MyMG%cVariable.EQ."Displac") THEN
   ILEV = mgLev
   NDOF = KNVT(mgLev)
   CALL E011Sum(myMG%AUX(mgLev)%x(0*NDOF+1))
   CALL E011Sum(myMG%AUX(mgLev)%x(1*NDOF+1))
   CALL E011Sum(myMG%AUX(mgLev)%x(2*NDOF+1))
  END IF
  
  IF (MyMG%cVariable.EQ."Temper") THEN !done
   ILEV = mgLev
   NDOF = KNVT(mgLev)
   do iFld=1,MyMG%nOfFields
    CALL E011Sum(myMG%AUX(mgLev)%x((iFld-1)*NDOF+1))
   END DO
  END IF
  
 END IF
 CALL mgDefectNorm(mgLev)
 IF (myid.ne.0) THEN
  myMG%AUX(mgLev)%x = 0d0
 END IF
END IF

END SUBROUTINE mgUpdateDefect
!
! ----------------------------------------------
!
SUBROUTINE mgRestriction()
INTEGER I,NDOF

IF (myid.ne.0) THEN
 IF (MyMG%cVariable.EQ."Displac") THEN
  CALL ZTIME(time0)
  CALL E011_Restriction(myMG%D(mgLev+1)%x,myMG%B(mgLev)%x,mg_E011Rest(mgLev)%a,&
       mg_E011RestM(mgLev)%LdA,mg_E011RestM(mgLev)%ColA,myMG%KNPRU,myMG%KNPRV,myMG%KNPRW)
  CALL ZTIME(time1)
  myStat%tRestUVW = myStat%tRestUVW + (time1-time0)
 END IF
 
 IF (MyMG%cVariable.EQ."Temper") THEN !done
  CALL ZTIME(time0)
  CALL E011_GenLinSc_Q1_Restriction(myMG%D(mgLev+1)%x,myMG%B(mgLev)%x,mg_E011Rest(mgLev)%a,&
       mg_E011RestM(mgLev)%LdA,mg_E011RestM(mgLev)%ColA,myMG%KNPR,myMG%nOfFields)
  CALL ZTIME(time1)
  myStat%tRestUVW = myStat%tRestUVW + (time1-time0)
 END IF
END IF

END SUBROUTINE mgRestriction
!
! ----------------------------------------------
!
SUBROUTINE mgProlongation()
INTEGER I,NDOF

IF (myid.ne.0) THEN
 IF (MyMG%cVariable.EQ."Displac") THEN
  CALL ZTIME(time0)
  CALL E011_Prolongation(myMG%AUX(mgLev)%x,myMG%X(mgLev-1)%x,mg_E011Prol(mgLev-1)%a,&
       mg_E011ProlM(mgLev-1)%LdA,mg_E011ProlM(mgLev-1)%ColA,myMG%KNPRU,myMG%KNPRV,myMG%KNPRW)
  CALL ZTIME(time1)
  myStat%tProlUVW = myStat%tProlUVW + (time1-time0)
 END IF
 
 IF (MyMG%cVariable.EQ."Temper") THEN !done
  CALL ZTIME(time0)
  CALL E011_GenLinSc_Q1_Prolongation(myMG%AUX(mgLev)%x,myMG%X(mgLev-1)%x,mg_E011Prol(mgLev-1)%a,&
       mg_E011ProlM(mgLev-1)%LdA,mg_E011ProlM(mgLev-1)%ColA,myMG%KNPR,myMG%nOfFields)
  CALL ZTIME(time1)
  myStat%tProlUVW = myStat%tProlUVW + (time1-time0)
 END IF
END IF

END SUBROUTINE mgProlongation
!
! ----------------------------------------------
!
SUBROUTINE mgCoarseGridSolver()
INTEGER ndof,neq,nnSteps
INTEGER ITE,I,j
REAL*8 def,def0

 mgLev = myMG%MedLev
 ILEV = myMG%MedLev
!   write(*,*) myMG%nIterCoarse,MyMG%cVariable

 IF (MyMG%cVariable.EQ."Displac") THEN
  CALL ZTIME(time0)
  IF (myid.ne.0) THEN
   CoarseIter  = myMG%nIterCoarse
   ILEV = mgLev
   neq = KNVT(mgLev) !+ KNAT(mgLev) + KNET(mgLev) + KNEL(mgLev)
   myMG%X(mgLev)%x = 0d0
   CALL SSORSolverYYY(myMG%A11(mgLev)%a,myMG%A22(mgLev)%a,myMG%A33(mgLev)%a,&
        myMG%L(mgLEV)%ColA,myMG%L(mgLEV)%LdA,&
        myMG%X(mgLev)%x,myMG%B(mgLev)%x,myJCB%d1,ParKNPR,neq,0,myMG%RLX,def0)
  END IF
  def0 = MAX(1d-30,def0)
  CALL COMM_maximum(def0)

  nnSteps = myMG%nIterCoarse
  ITE = 0
  DO I=1,10
   ITE = ITE + 1
   IF (myid.ne.0) THEN
     CALL E011_SSORSolver(myMG%A11(mgLev)%a,myMG%A22(mgLev)%a,myMG%A33(mgLev)%a,&
          myMG%L(mgLEV)%ColA,myMG%L(mgLEV)%LdA,&
          myMG%X(mgLev)%x,myMG%B(mgLev)%x,myJCB%d1,ParKNPR,neq,nnSteps,myMG%RLX,def)
   END IF
   def = MAX(1d-33,def)
   CALL COMM_maximum(def)
   
!    if (myid.eq.showid) write(*,*) "def/def0",def/def0
   IF (def/def0.LT.MyMG%DefImprCoarse) GOTO 1
  END DO
1 CONTINUE
  CoarseIter = ITE*nnSteps
  CALL ZTIME(time1)
  myStat%tSolvUVW = myStat%tSolvUVW + (time1-time0)
 END IF
 
 IF (MyMG%cVariable.EQ."Temper") THEN
  CALL ZTIME(time0)
  IF (myid.ne.0) THEN
   CoarseIter  = myMG%nIterCoarse
   ILEV = mgLev
   neq = KNVT(mgLev) !+ KNAT(mgLev) + KNET(mgLev) + KNEL(mgLev)
   myMG%X(mgLev)%x = 0d0
   CALL SSORSolver_GenLinSc_Q1(myMG%AXX,myMG%L(mgLEV)%ColA,myMG%L(mgLEV)%LdA,&
        myMG%X(mgLev)%x,myMG%B(mgLev)%x,myJCB%d1,ParKNPR,neq,0,myMG%nOfFields,myMG%RLX,def0)
  END IF
  def0 = MAX(1d-30,def0)
  CALL COMM_maximum(def0)

  nnSteps = myMG%nIterCoarse
  ITE = 0
  DO I=1,10
   ITE = ITE + 1
   IF (myid.ne.0) THEN
    CALL SSORSolver_GenLinSc_Q1(myMG%AXX,myMG%L(mgLEV)%ColA,myMG%L(mgLEV)%LdA,&
             myMG%X(mgLev)%x,myMG%B(mgLev)%x,myJCB%d1,ParKNPR,neq,nnSteps,myMG%nOfFields,myMG%RLX,def)
   END IF
   def = MAX(1d-33,def)
   CALL COMM_maximum(def)
   
!    if (myid.eq.showid) write(*,*) "def/def0",def/def0
   IF (def/def0.LT.MyMG%DefImprCoarse) GOTO 2
  END DO
2 CONTINUE
  CoarseIter = ITE*nnSteps
  CALL ZTIME(time1)
  myStat%tSolvUVW = myStat%tSolvUVW + (time1-time0)
 END IF

 !  IF (myid.eq.0) write(*,*) myMG%X(mgLev)%x
END SUBROUTINE mgCoarseGridSolver
!
! ----------------------------------------------
!
SUBROUTINE mgSmoother()
INTEGER Iter,i,j,ndof,neq
REAL*8 daux

 IF (MyMG%cVariable.EQ."Displac") THEN
  ILEV = mgLev
  if (myid.ne.0) then
   Iter  = myMG%nSmootherSteps
   ndof  = SIZE(myMG%X(mgLev)%x)

   CALL ZTIME(time0)
   neq = KNVT(mgLev) !+ KNAT(mgLev) + KNET(mgLev) + KNEL(mgLev)
   CALL SSORSolverXXX(myMG%A11(mgLev)%a,myMG%A22(mgLev)%a,myMG%A33(mgLev)%a,&
             myMG%L(mgLEV)%ColA,myMG%L(mgLEV)%LdA,&
             myMG%X(mgLev)%x,myMG%B(mgLev)%x,myJCB%d1,ParKNPR,neq,Iter,myMG%RLX)
   CALL ZTIME(time1)
   myStat%tSmthUVW = myStat%tSmthUVW + (time1-time0)
  end if
 END IF

 IF (MyMG%cVariable.EQ."Temper") THEN !done
  ILEV = mgLev
  if (myid.ne.0) then
   Iter  = myMG%nSmootherSteps
   ndof  = SIZE(myMG%X(mgLev)%x)

   CALL ZTIME(time0)
   neq = KNVT(mgLev) !+ KNAT(mgLev) + KNET(mgLev) + KNEL(mgLev)
!       SUBROUTINE SSORSmoother_GenLinSc_Q1(AXX,KCOL,KLD,
!      *           DX,DB,DD,KNPR,NEQ,NIT,NFLD,OMEGA)
   CALL SSORSmoother_GenLinSc_Q1(myMG%AXX,myMG%L(mgLEV)%ColA,myMG%L(mgLEV)%LdA,&
             myMG%X(mgLev)%x,myMG%B(mgLev)%x,myJCB%d1,ParKNPR,neq,Iter,myMG%nOfFields,myMG%RLX)
   CALL ZTIME(time1)
   myStat%tSmthUVW = myStat%tSmthUVW + (time1-time0)
  end if
 END IF
 
END SUBROUTINE mgSmoother
!
! ----------------------------------------------
!
SUBROUTINE E011_GenLinSc_Q1_Prolongation(D2,D1,A,KLD,KCOL,KNPR,NFLD)
IMPLICIT NONE
REAL*8 A(*)
INTEGER KLD(*),KCOL(*),KNPR(*)
REAL*8 D1(*),D2(*)
INTEGER I,J,ICOL,IFLD,NFLD
INTEGER MEQ,iMEQ,iNEQ
INTEGER LEQ,iLEQ

IF (myid.eq.0) RETURN
MEQ = KNVT(mgLev  )!+KNAT(mgLev  )+KNET(mgLev  )+KNEL(mgLev  )
LEQ = KNVT(mgLev-1)!+KNAT(mgLev-1)+KNET(mgLev-1)+KNEL(mgLev-1)

DO iFld=1,nFLD

 iMEQ = (iFld-1)*MEQ
 iNEQ = (iFld-1)*myMG%nOfSubsystemEqs
 iLEQ = (iFld-1)*LEQ

 DO I=1,MEQ
  D2(iMEQ+I) = 0d0
!   IF (KNPR(iNEQ+I).EQ.0) THEN
   DO J=KLD(I),KLD(I+1)-1
    ICOL = KCOL(J)
    D2(iMEQ+I) = D2(iMEQ+I) + A(J)*D1(iLEQ+ICOL)
   END DO
!   END IF
 END DO

END DO
 
END SUBROUTINE
!
! ----------------------------------------------
!
SUBROUTINE E011_Prolongation(D2,D1,A,KLD,KCOL,KNPRU,KNPRV,KNPRW)
IMPLICIT NONE
REAL*8 A(*)
INTEGER KLD(*),KCOL(*),KNPRU(*),KNPRV(*),KNPRW(*)
REAL*8 D1(*),D2(*)
INTEGER I,J,ICOL
INTEGER MEQ,MEQ1,MEQ2,MEQ3
INTEGER LEQ,LEQ1,LEQ2,LEQ3

IF (myid.eq.0) RETURN
MEQ = KNVT(mgLev  )!+KNAT(mgLev  )+KNET(mgLev  )+KNEL(mgLev  )
LEQ = KNVT(mgLev-1)!+KNAT(mgLev-1)+KNET(mgLev-1)+KNEL(mgLev-1)

MEQ1 = 0
MEQ2 = MEQ
MEQ3 = 2*MEQ
LEQ1 = 0
LEQ2 = LEQ
LEQ3 = 2*LEQ

DO I=1,MEQ
 D2(MEQ1+I) = 0d0
 D2(MEQ2+I) = 0d0
 D2(MEQ3+I) = 0d0
 IF (KNPRU(I).NE.1) THEN
  DO J=KLD(I),KLD(I+1)-1
   ICOL = KCOL(J)
   D2(MEQ1+I) = D2(MEQ1+I) + A(J)*D1(LEQ1+ICOL)
  END DO
 END IF
 IF (KNPRV(I).NE.1) THEN
  DO J=KLD(I),KLD(I+1)-1
   ICOL = KCOL(J)
   D2(MEQ2+I) = D2(MEQ2+I) + A(J)*D1(LEQ2+ICOL)
  END DO
 END IF
 IF (KNPRW(I).NE.1) THEN
  DO J=KLD(I),KLD(I+1)-1
   ICOL = KCOL(J)
   D2(MEQ3+I) = D2(MEQ3+I) + A(J)*D1(LEQ3+ICOL)
  END DO
 END IF
END DO

END SUBROUTINE
!
! ----------------------------------------------
!
SUBROUTINE E011_GenLinSc_Q1_Restriction(D2,D1,A,KLD,KCOL,KNPR,nFLD)
IMPLICIT NONE
REAL*8 A(*)
INTEGER N,KLD(*),KCOL(*),KNPR(*)
INTEGER nFLD,iFLD
REAL*8 D1(*),D2(*)
INTEGER I,J,ICOL
INTEGER MEQ,iMEQ,iNEQ
INTEGER LEQ,iLEQ

IF (myid.eq.0) RETURN
MEQ = KNVT(mgLev  )!+KNAT(mgLev  )+KNET(mgLev  )+KNEL(mgLev  )
LEQ = KNVT(mgLev+1)!+KNAT(mgLev+1)+KNET(mgLev+1)+KNEL(mgLev+1)


DO iFld=1,nFLD

 iMEQ = (iFld-1)*MEQ
 iNEQ = (iFld-1)*myMG%nOfSubsystemEqs
 iLEQ = (iFld-1)*LEQ

 DO I=1,MEQ
  D1(iMEQ+I) = 0d0
!   IF (KNPR(iNEQ+I).eq.0) THEN
   DO J=KLD(I),KLD(I+1)-1
    ICOL = KCOL(J)
    D1(iMEQ+I) = D1(iMEQ+I) + A(J)*D2(iLEQ+ICOL)
   END DO
!   END IF
 END DO

END DO

ILEV = mgLev
! if (myid.eq.1) WRITE(*,*) ilev
! CALL E013Sum(D1)

END SUBROUTINE
!
! ----------------------------------------------
!
SUBROUTINE E011_Restriction(D2,D1,A,KLD,KCOL,KNPRU,KNPRV,KNPRW)
IMPLICIT NONE
REAL*8 A(*)
INTEGER N,KLD(*),KCOL(*),KNPRU(*),KNPRV(*),KNPRW(*)
REAL*8 D1(*),D2(*)
INTEGER I,J,ICOL
INTEGER MEQ,MEQ1,MEQ2,MEQ3
INTEGER LEQ,LEQ1,LEQ2,LEQ3

IF (myid.eq.0) RETURN
MEQ = KNVT(mgLev  )!+KNAT(mgLev  )+KNET(mgLev  )+KNEL(mgLev  )
LEQ = KNVT(mgLev+1)!+KNAT(mgLev+1)+KNET(mgLev+1)+KNEL(mgLev+1)

MEQ1 = 0
MEQ2 = MEQ
MEQ3 = 2*MEQ
LEQ1 = 0
LEQ2 = LEQ
LEQ3 = 2*LEQ

DO I=1,MEQ
 D1(MEQ1+I) = 0d0
 D1(MEQ2+I) = 0d0
 D1(MEQ3+I) = 0d0
 IF (KNPRU(I).NE.1) THEN
  DO J=KLD(I),KLD(I+1)-1
   ICOL = KCOL(J)
   D1(MEQ1+I) = D1(MEQ1+I) + A(J)*D2(LEQ1+ICOL)
  END DO
 END IF
 IF (KNPRV(I).NE.1) THEN
  DO J=KLD(I),KLD(I+1)-1
   ICOL = KCOL(J)
   D1(MEQ2+I) = D1(MEQ2+I) + A(J)*D2(LEQ2+ICOL)
  END DO
 END IF
 IF (KNPRW(I).NE.1) THEN
  DO J=KLD(I),KLD(I+1)-1
   ICOL = KCOL(J)
   D1(MEQ3+I) = D1(MEQ3+I) + A(J)*D2(LEQ3+ICOL)
  END DO
 END IF
END DO

ILEV = mgLev
! if (myid.eq.1) WRITE(*,*) ilev
! CALL E013Sum(D1)

END SUBROUTINE
!
! ----------------------------------------------
!
SUBROUTINE E011_DAX(DX,DD,A1,A2)
IMPLICIT NONE
REAL*8  DX(*),DD(*),A1,A2
INTEGER MEQ,MEQ1,MEQ2,MEQ3
INTEGER IEQ,IA,ICOL

ILEV = mgLev

MEQ = KNVT(mgLev  ) !+KNAT(mgLev  )+KNET(mgLev  )+KNEL(mgLev  )
MEQ1= 0
MEQ2= MEQ
MEQ3= 2*MEQ

DO IEQ=1,MEQ
 DD(MEQ1+IEQ) = A2*DD(MEQ1+IEQ)
 DD(MEQ2+IEQ) = A2*DD(MEQ2+IEQ)
 DD(MEQ3+IEQ) = A2*DD(MEQ3+IEQ)
 DO IA=myMG%L(mgLEV)%LdA(IEQ),myMG%L(mgLEV)%LdA(IEQ+1)-1
  ICOL = myMG%L(mgLEV)%ColA(IA)
  DD(MEQ1+IEQ) = DD(MEQ1+IEQ) + A1*myMG%A11(mgLEV)%a(IA)*DX(MEQ1+ICOL)
  DD(MEQ2+IEQ) = DD(MEQ2+IEQ) + A1*myMG%A22(mgLEV)%a(IA)*DX(MEQ2+ICOL)
  DD(MEQ3+IEQ) = DD(MEQ3+IEQ) + A1*myMG%A33(mgLEV)%a(IA)*DX(MEQ3+ICOL)
 END DO
END DO

END SUBROUTINE E011_DAX
!
! ----------------------------------------------
!
SUBROUTINE E011_GenLinSc_Q1_DAX(DX,DD,A1,A2,NFLD)
IMPLICIT NONE
REAL*8  DX(*),DD(*),A1,A2
INTEGER MEQ,iMEQ,NFLD
INTEGER IEQ,IA,ICOL,IFLD

ILEV = mgLev

MEQ = KNVT(mgLev  ) !+KNAT(mgLev  )+KNET(mgLev  )+KNEL(mgLev  )

DO IFLD=1,NFLD
 iMEQ= (IFLD-1)*MEQ
 DO IEQ=1,MEQ
  DD(iMEQ+IEQ) = A2*DD(iMEQ+IEQ)
  DO IA=myMG%L(mgLEV)%LdA(IEQ),myMG%L(mgLEV)%LdA(IEQ+1)-1
   ICOL = myMG%L(mgLEV)%ColA(IA)
   DD(iMEQ+IEQ) = DD(iMEQ+IEQ) + A1*myMG%AXX(IFLD)%fld(mgLEV)%a(IA)*DX(iMEQ+ICOL)
  END DO
 END DO
END DO

END SUBROUTINE E011_GenLinSc_Q1_DAX
!
! ----------------------------------------------
!
SUBROUTINE mgProlRestInit

!  IF (myMG%bProlRest) RETURN
  IF (MyMG%cVariable.EQ."Displac") THEN
   DO mgLev = myMG%MedLev+1, myMG%MaxLev
    ILEV = mgLev
    CALL InitE011ProlMat(mg_E011ProlM(mgLev-1)%na,&
         mg_E011Prol(mgLev-1)%a,mg_E011ProlM(mgLev-1)%LdA,mg_E011ProlM(mgLev-1)%ColA,&
         mg_E011Rest(mgLev-1)%a,mg_E011RestM(mgLev-1)%LdA,mg_E011RestM(mgLev-1)%ColA,&
         mg_mesh%level(mgLev-1)%kadj,mg_mesh%level(mgLev-1)%kvert,&
         mg_mesh%level(mgLev-1)%kedge,&
         mg_mesh%level(mgLev-1)%karea,mg_mesh%level(mgLev-1)%nvt,&
         mg_mesh%level(mgLev-1)%net,mg_mesh%level(mgLev-1)%nat,&
         mg_mesh%level(mgLev-1)%nel,&
         mg_mesh%level(mgLev)%kadj,mg_mesh%level(mgLev)%kvert,&
         mg_mesh%level(mgLev)%kedge,&
         mg_mesh%level(mgLev)%karea,mg_mesh%level(mgLev)%nvt,&
         mg_mesh%level(mgLev)%net,mg_mesh%level(mgLev)%nat,&
         mg_mesh%level(mgLev)%nel)
   END DO
  END IF

  IF (MyMG%cVariable.EQ."Temper") THEN !done
   DO mgLev = myMG%MedLev+1, myMG%MaxLev
    ILEV = mgLev
    CALL InitE011ProlMat(mg_E011ProlM(mgLev-1)%na,&
         mg_E011Prol(mgLev-1)%a,mg_E011ProlM(mgLev-1)%LdA,mg_E011ProlM(mgLev-1)%ColA,&
         mg_E011Rest(mgLev-1)%a,mg_E011RestM(mgLev-1)%LdA,mg_E011RestM(mgLev-1)%ColA,&
         mg_mesh%level(mgLev-1)%kadj,mg_mesh%level(mgLev-1)%kvert,&
         mg_mesh%level(mgLev-1)%kedge,&
         mg_mesh%level(mgLev-1)%karea,mg_mesh%level(mgLev-1)%nvt,&
         mg_mesh%level(mgLev-1)%net,mg_mesh%level(mgLev-1)%nat,&
         mg_mesh%level(mgLev-1)%nel,&
         mg_mesh%level(mgLev)%kadj,mg_mesh%level(mgLev)%kvert,&
         mg_mesh%level(mgLev)%kedge,&
         mg_mesh%level(mgLev)%karea,mg_mesh%level(mgLev)%nvt,&
         mg_mesh%level(mgLev)%net,mg_mesh%level(mgLev)%nat,&
         mg_mesh%level(mgLev)%nel)
   END DO
  END IF

  MyMG%bProlRest=.TRUE.

END SUBROUTINE mgProlRestInit
!
! ----------------------------------------------
!
SUBROUTINE InitE011ProlMat(NA,PA,PLDA,PCOLA,RA,RLDA,RCOLA,&
   KADJ1,KVERT1,KEDGE1,KAREA1,NVT1,NET1,NAT1,NEL1,&
   KADJ2,KVERT2,KEDGE2,KAREA2,NVT2,NET2,NAT2,NEL2)

IMPLICIT NONE
REAL*8 PA(*),RA(*)
INTEGER NA,PLDA(*),PCOLA(*),RLDA(*),RCOLA(*)
INTEGER, ALLOCATABLE, DIMENSION(:) :: RLDB,PLDB
REAL*8 A(27,8)

INTEGER KVERT2(8,*),KADJ2(6,*),KEDGE2(12,*),KAREA2(6,*),NVT2,NET2,NAT2,NEL2
INTEGER KVERT1(8,*),KADJ1(6,*),KEDGE1(12,*),KAREA1(6,*),NVT1,NET1,NAT1,NEL1
INTEGER IEL,JEL(8),I,J,II,JJ,KK,IND_I(27),IND_J(8),II_POS,JJ_POS
REAL*8  DAUX


CALL InitE011ElemProlMat(A)


PA(1:NA) = 0d0
RA(1:NA) = 0d0
PCOLA(1:NA) = 0
RCOLA(1:NA) = 0
PLDA(1:NVT2+1) = 0
RLDA(1:NVT1+1) = 0

DO IEL=1,NEL1
 JEL(1)  = IEL
 JEL(2)  = KADJ2(3,JEL(1))
 JEL(3)  = KADJ2(3,JEL(2))
 JEL(4)  = KADJ2(3,JEL(3))
 JEL(5)  = KADJ2(6,JEL(1))
 JEL(6)  = KADJ2(6,JEL(2))
 JEL(7)  = KADJ2(6,JEL(3))
 JEL(8)  = KADJ2(6,JEL(4))

 IND_I( 1) = KVERT2( 1,JEL(1))
 IND_I( 2) = KVERT2( 2,JEL(1))
 IND_I( 3) = KVERT2( 1,JEL(2))
 IND_I( 4) = KVERT2( 4,JEL(1))
 IND_I( 5) = KVERT2( 3,JEL(1))
 IND_I( 6) = KVERT2( 2,JEL(2))
 IND_I( 7) = KVERT2( 1,JEL(4))
 IND_I( 8) = KVERT2( 4,JEL(4))
 IND_I( 9) = KVERT2( 1,JEL(3))

 IND_I(10) = KVERT2( 5,JEL(1))
 IND_I(11) = KVERT2( 6,JEL(1))
 IND_I(12) = KVERT2( 5,JEL(2))
 IND_I(13) = KVERT2( 8,JEL(1))
 IND_I(14) = KVERT2( 7,JEL(1))
 IND_I(15) = KVERT2( 6,JEL(2))
 IND_I(16) = KVERT2( 5,JEL(4))
 IND_I(17) = KVERT2( 8,JEL(4))
 IND_I(18) = KVERT2( 5,JEL(3))

 IND_I(19) = KVERT2( 1,JEL(5))
 IND_I(20) = KVERT2( 2,JEL(5))
 IND_I(21) = KVERT2( 1,JEL(6))
 IND_I(22) = KVERT2( 4,JEL(5))
 IND_I(23) = KVERT2( 3,JEL(5))
 IND_I(24) = KVERT2( 2,JEL(6))
 IND_I(25) = KVERT2( 1,JEL(8))
 IND_I(26) = KVERT2( 4,JEL(8))
 IND_I(27) = KVERT2( 1,JEL(7))


 IND_J( 1) = KVERT1( 1,IEL)
 IND_J( 2) = KVERT1( 2,IEL)
 IND_J( 3) = KVERT1( 3,IEL)
 IND_J( 4) = KVERT1( 4,IEL)
 IND_J( 5) = KVERT1( 5,IEL)
 IND_J( 6) = KVERT1( 6,IEL)
 IND_J( 7) = KVERT1( 7,IEL)
 IND_J( 8) = KVERT1( 8,IEL)

 DO I=1,27
   II = IND_I(I)
   IF (PLDA(II+1).EQ.0) THEN
    KK = 0
    DO J=1,8
     JJ = IND_J(J)
     IF (ABS(A(I,J)).GT.1d-16) THEN
      RLDA(JJ+1) = RLDA(JJ+1) + 1
      KK = KK + 1
     END IF
    END DO
    PLDA(II+1) = KK
   END IF
 END DO

END DO !IEL

PLDA(1) = 1
RLDA(1) = 1
DO I=2,NVT1+1
 RLDA(I) = RLDA(I) + RLDA(I-1)
END DO
DO I=2,NVT2+1
 PLDA(I) = PLDA(I) + PLDA(I-1)
END DO

ALLOCATE (RLDB(NVT1+1),PLDB(NVT2+1))
RLDB = 0
PLDB = 0

DO IEL=1,NEL1
 JEL(1)  = IEL
 JEL(2)  = KADJ2(3,JEL(1))
 JEL(3)  = KADJ2(3,JEL(2))
 JEL(4)  = KADJ2(3,JEL(3))
 JEL(5)  = KADJ2(6,JEL(1))
 JEL(6)  = KADJ2(6,JEL(2))
 JEL(7)  = KADJ2(6,JEL(3))
 JEL(8)  = KADJ2(6,JEL(4))


 IND_I( 1) = KVERT2( 1,JEL(1))
 IND_I( 2) = KVERT2( 2,JEL(1))
 IND_I( 3) = KVERT2( 1,JEL(2))
 IND_I( 4) = KVERT2( 4,JEL(1))
 IND_I( 5) = KVERT2( 3,JEL(1))
 IND_I( 6) = KVERT2( 2,JEL(2))
 IND_I( 7) = KVERT2( 1,JEL(4))
 IND_I( 8) = KVERT2( 4,JEL(4))
 IND_I( 9) = KVERT2( 1,JEL(3))

 IND_I(10) = KVERT2( 5,JEL(1))
 IND_I(11) = KVERT2( 6,JEL(1))
 IND_I(12) = KVERT2( 5,JEL(2))
 IND_I(13) = KVERT2( 8,JEL(1))
 IND_I(14) = KVERT2( 7,JEL(1))
 IND_I(15) = KVERT2( 6,JEL(2))
 IND_I(16) = KVERT2( 5,JEL(4))
 IND_I(17) = KVERT2( 8,JEL(4))
 IND_I(18) = KVERT2( 5,JEL(3))

 IND_I(19) = KVERT2( 1,JEL(5))
 IND_I(20) = KVERT2( 2,JEL(5))
 IND_I(21) = KVERT2( 1,JEL(6))
 IND_I(22) = KVERT2( 4,JEL(5))
 IND_I(23) = KVERT2( 3,JEL(5))
 IND_I(24) = KVERT2( 2,JEL(6))
 IND_I(25) = KVERT2( 1,JEL(8))
 IND_I(26) = KVERT2( 4,JEL(8))
 IND_I(27) = KVERT2( 1,JEL(7))


 IND_J( 1) = KVERT1( 1,IEL)
 IND_J( 2) = KVERT1( 2,IEL)
 IND_J( 3) = KVERT1( 3,IEL)
 IND_J( 4) = KVERT1( 4,IEL)
 IND_J( 5) = KVERT1( 5,IEL)
 IND_J( 6) = KVERT1( 6,IEL)
 IND_J( 7) = KVERT1( 7,IEL)
 IND_J( 8) = KVERT1( 8,IEL)

 DO I=1,27
   II = IND_I(I)
   IF (PLDB(II+1).EQ.0) THEN
    DAUX = 1d0
    KK = 0
    DO J=1,8
     JJ = IND_J(J)
     IF (ABS(A(I,J)).GT.1d-16) THEN
      JJ_POS = RLDA(JJ) + RLDB(JJ+1)
      II_POS = PLDA(II) + KK
      RCOLA(JJ_POS) = II
      PCOLA(II_POS) = JJ
      RA(JJ_POS) = DAUX*A(I,J)
      PA(II_POS) = A(I,J)
      RLDB(JJ+1) = RLDB(JJ+1) + 1
      KK = KK + 1
     END IF
    END DO
    PLDB(II+1) = KK
   END IF
 END DO

END DO !IEL

DEALLOCATE (RLDB,PLDB)

END SUBROUTINE InitE011ProlMat
!
! ----------------------------------------------
!
SUBROUTINE InitE011ElemProlMat(A)
REAL*8 A(27,8)
REAL*8 F(8),PI,PJ,PK
INTEGER I,J,K,N,IEL,II

! OPEN(789,FILE=outfile)

A(:,:) = 0d0

II = 0
DO I=-1,1
 PI = DBLE(I)/1d0
 DO J=-1,1
  PJ = DBLE(J)/1d0
  DO K=-1,1
   PK = DBLE(K)/1d0
   CALL E011_Lin(PK,PJ,PI,F)
   II = II + 1
   A(II,:) = F(:)
!   WRITE(789,'(I3,27D12.4)') II,F(:)
  END DO
 END DO
END DO

! CLOSE(789)
! stop

END SUBROUTINE InitE011ElemProlMat
!
! ----------------------------------------------
!
SUBROUTINE E011_Lin(X1,X2,X3,D)
REAL*8 D(8),X1,X2,X3
REAL*8 :: Q8=0.125D0

D(1)=Q8*(1D0-X1)*(1D0-X2)*(1D0-X3)
D(2)=Q8*(1D0+X1)*(1D0-X2)*(1D0-X3)
D(3)=Q8*(1D0+X1)*(1D0+X2)*(1D0-X3)
D(4)=Q8*(1D0-X1)*(1D0+X2)*(1D0-X3)
D(5)=Q8*(1D0-X1)*(1D0-X2)*(1D0+X3)
D(6)=Q8*(1D0+X1)*(1D0-X2)*(1D0+X3)
D(7)=Q8*(1D0+X1)*(1D0+X2)*(1D0+X3)
D(8)=Q8*(1D0-X1)*(1D0+X2)*(1D0+X3)

END SUBROUTINE E011_Lin
!
! ----------------------------------------------
!
SUBROUTINE JCB_Activation()
INTEGER ndof,ndof_orig,NN

IF (myid.ne.0) THEN 

 ndof  = SIZE(myMG%X(MyMG%MaxLev)%x)

 IF (MyMG%cVariable.EQ."Disp") THEN ! done
  IF (.NOT.myJCB%bActivated) THEN

    ALLOCATE(myJCB%d1(3*ndof))

   myJCB%bActivated = .TRUE.
  ELSE
   ndof_orig  = SIZE(myJCB%d1)
   IF (ndof_orig.LT.ndof.AND.myid.ne.0) THEN
    DEALLOCATE(myJCB%d1)
    ALLOCATE(myJCB%d1(3*ndof))
   END IF
  END IF
 END IF

 IF (MyMG%cVariable.EQ."Temper") THEN ! done
  NN = KNVT(mgLev)
  IF (.NOT.myJCB%bActivated) THEN

    ALLOCATE(myJCB%d1((8**(1-MyMG%MaxDifLev))*MyMG%nOfFields*NN))

   myJCB%bActivated = .TRUE.
  ELSE
   ndof_orig  = SIZE(myJCB%d1)
   IF (ndof_orig.LT.ndof.AND.myid.ne.0) THEN
    DEALLOCATE(myJCB%d1)
    ALLOCATE(myJCB%d1((8**(1-MyMG%MaxDifLev))*MyMG%nOfFields*NN))
   END IF
  END IF
 END IF

END IF

END SUBROUTINE JCB_Activation

END MODULE mg_LinScalar
