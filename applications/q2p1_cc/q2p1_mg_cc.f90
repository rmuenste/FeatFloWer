MODULE mg_cc

USE PP3D_MPI, ONLY:E011Sum,E011DMat,myid,showID,MGE013,&
                   COMM_Maximum,COMM_SUMM,COMM_NLComplete,myMPI_Barrier
USE var_QuadScalar
USE var_QuadScalar_newton
USE UMFPackSolver_CC, ONLY : myUmfPack_CCSolve,myUmfPack_CCSolveMaster,&
                   myUmfPack_CCSolveLocalMat    

USE UMFPackSolver, only : myUmfPack_Solve
USE MumpsSolver, only : MUMPS_solver_Distributed


IMPLICIT NONE

TYPE tCG
 LOGICAL :: bActivated = .FALSE.
 REAL*8, DIMENSION(:), ALLOCATABLE :: d1,d2,d3,d4,d5,d6
END TYPE tCG

TYPE tJCB
 LOGICAL :: bActivated = .FALSE.
 REAL*8, DIMENSION(:), ALLOCATABLE :: d1
END TYPE tJCB

TYPE (tCG) :: myCG
TYPE (tJCB) :: myJCB

INTEGER :: iMGL
! -****--------------*-*********-*----------
INTEGER IterCycle,mgLev,CoarseIter
REAL*8 mgDefNorm
REAL*8 time0,time1

CONTAINS
!
! ----------------------------------------------
!
SUBROUTINE MG_Solver_CC(mfile,mterm,i)
implicit none
INTEGER mfile,mterm
INTEGER INLComplete
REAL*8 DefI1,DefI2,DefImpr,DDD,AccCoarseIter,DDNorm,i

!IF (MyMg%MaxLev.EQ.MyMg%MinLev) THEN
!
! RETURN
!END IF


CALL mgInit_cc()
CALL mgComputeDefect_cc(MyMg%MaxLev,.TRUE.)

DefI1 = 0d0
DefI2 = 0d0
AccCoarseIter = 0
i=0d0

DO IterCycle=1,MyMg%MaxIterCycle
 i=i+1
 INLComplete = 0

 CALL mg_cycle_cc()

!  CALL ImposeBC(MyMg%MaxLev)
 
 ! Evaluation of new defect
 CALL mgComputeDefect_cc(MyMg%MaxLev,.TRUE.)

 AccCoarseIter = AccCoarseIter + CoarseIter
 DefImpr = mgDefNorm/MyMg%DefInitial

 IF (mgDefNorm.LE.MyMg%DefInitial*MyMg%Criterion1.AND.&
     mgDefNorm.LE.MyMg%Criterion2.AND.&
     IterCycle.GE.MyMg%MinIterCycle) INLComplete = 1
 IF (IterCycle.GE.MyMg%MaxIterCycle) INLComplete = 1
 CALL COMM_NLComplete(INLComplete)
 IF (MyMg%cVariable.EQ."Pressure".AND.myid.eq.showid) THEN
  write(mfile,'(I4,2G12.4,A3,I5,A3,9I4)') &
  IterCycle,mgDefNorm,mgDefNorm/MyMg%DefInitial,&
  " | ",CoarseIter," | ",MyMg%nSmootherSteps
  write(mterm,'(I4,2G12.4,A3,I5,A3,9I4)') &
  IterCycle,mgDefNorm,mgDefNorm/MyMg%DefInitial,&
  " | ",CoarseIter," | ",MyMg%nSmootherSteps
 END IF
 IF (INLComplete.eq.1) GOTO 88
 DEFI2 = DEFI1
 DEFI1 = mgDefNorm

!       pause
END DO

88 CONTINUE

MyMg%DefFinal = mgDefNorm
MyMg%RhoMG1 = (MyMg%DefFinal/MyMg%DefInitial)**(1d0/DBLE(MAX(IterCycle,1)))
MyMg%RhoMG2 = (MyMg%DefFinal/DEFI2)**(0.5d0)
MyMg%UsedIterCycle = IterCycle
MyMg%nIterCoarse = INT(AccCoarseIter/MAX(IterCycle,1))

END SUBROUTINE MG_Solver_cc
!
! ----------------------------------------------
!
SUBROUTINE mg_cycle_cc()

IF (MyMG%CycleType.EQ."W")  CALL mg_W_cycle_cc()
IF (MyMG%CycleType.EQ."V")  CALL mg_V_cycle_cc()
IF (MyMG%CycleType.EQ."F")  CALL mg_F_cycle_cc()
IF (MyMG%CycleType.EQ."N")  CALL mgCoarseGridSolverShouldBe_cc()

END SUBROUTINE mg_cycle_cc
!
! ----------------------------------------------
!
SUBROUTINE mg_W_cycle_cc()
INTEGER imgLev

 CALL mg_down_cc(myMG%MaxLev)

 CALL mg_W_subcycle_cc(myMG%MaxLev-1)

 CALL mg_up_cc(myMG%MaxLev)

END SUBROUTINE mg_W_cycle_cc
!
! ----------------------------------------------
!
RECURSIVE SUBROUTINE mg_W_subcycle_cc(imgLev)
INTEGER imgLev

 IF (imgLev.NE.2) CALL mg_W_subcycle_cc(imgLev-1)

 CALL mg_up_cc(imgLev)
 CALL mg_down_cc(imgLev)
!  CALL mg_up(myMG%MaxLev)
!  CALL mg_down(myMG%MaxLev)

 IF (imgLev.NE.2) CALL mg_W_subcycle_cc(imgLev-1)

END SUBROUTINE mg_W_subcycle_cc
!
! ----------------------------------------------
!
SUBROUTINE mg_F_cycle_cc()
INTEGER imgLev

 CALL mg_down_cc(myMG%MaxLev)

 DO imgLev = myMG%MinLev+1,myMG%MaxLev-1
  CALL mg_up_cc(imgLev)
  CALL mg_down_cc(imgLev)
 END DO

 CALL mg_up_cc(myMG%MaxLev)

END SUBROUTINE mg_F_cycle_cc
!
! ----------------------------------------------
!
SUBROUTINE mg_V_cycle_cc()

 CALL mg_down_cc(myMG%MaxLev)

 CALL mg_up_cc(myMG%MaxLev)

END SUBROUTINE mg_V_cycle_cc
!
! ----------------------------------------------
!
SUBROUTINE mg_down_cc(nimgLev)
INTEGER imgLev,nimgLev

! Presmoothing + restriction
 DO imgLev = nimgLev,myMG%MinLev+1,-1

  mgLev = imgLev
  CALL mgSmoother_cc()                                          ! takes B as RHS and smoothes X

  CALL mgComputeDefect_cc(imgLev,.FALSE.)                        ! puts B into D and gets the new defect
!
  mgLev = imgLev - 1
  CALL mgRestriction_cc()                                       ! brings D down by 1 level and stores it as B
!
  CALL mgGuessSolution_cc()                          ! choose zero initial vector

 END DO

 CALL mgCoarseGridSolverShouldBe_cc()                ! computes linear system with X=(A^-1)B

END SUBROUTINE mg_down_cc
!
! ----------------------------------------------
!
SUBROUTINE mg_up_cc(nimgLev)
INTEGER imgLev,nimgLev

 ! Postsmoothing + prolongation
 DO imgLev = myMG%MinLev+1, nimgLev

  mgLev = imgLev

  CALL mgProlongation_cc()                                      ! brings X up by one level and stores it as AUX

  CALL mgUpdateSolution_cc()                                    ! updates solution X = X_old + AUX(update)

  CALL mgSmoother_cc()                                          ! takes B as RHS And smoothes the solution X further

 END DO

END SUBROUTINE mg_up_cc
!
! ----------------------------------------------
!
SUBROUTINE mgGuessSolution_cc()

IF (myid.ne.0) THEN
  myMG%dX_u(mgLev)%x = 0d0
  myMG%dX_p(mgLev)%x = 0d0
END IF

END SUBROUTINE mgGuessSolution_cc
!
! ----------------------------------------------
!
SUBROUTINE mgUpdateSolution_cc()
INTEGER ndof

IF (myid.ne.0) THEN
 myMG%dX_u(mgLev)%x = myMG%dX_u(mgLev)%x + myMG%A_u(mgLev)%x
 myMG%dX_p(mgLev)%x = myMG%dX_p(mgLev)%x + myMG%A_p(mgLev)%x
END IF

END SUBROUTINE mgUpdateSolution_cc
!
! ----------------------------------------------
!
SUBROUTINE mgInit_cc()
INTEGER ndof_u

mgLev = myMG%MaxLev

! CALL mgComputeDefect(mgLev,.TRUE.)

IF (myid.ne.0) THEN

 ILEV = mgLev
 CALL SETLEV(2)
!  CALL E013UVWSUM(MyMG%B_u(mgLev)%x)

 myMG%A_u(mgLev)%x  = 0d0
!  myMG%X_u(mgLev)%x  = 0d0
 myMG%dX_u(mgLev)%x = 0d0
 myMG%A_p(mgLev)%x  = 0d0
!  myMG%X_p(mgLev)%x  = 0d0
 myMG%dX_p(mgLev)%x = 0d0
 ndof_u  = KNEL(mgLev) + KNVT(mgLev) + KNET(mgLev) + KNAT(mgLev)

END IF

IterCycle = 0

END SUBROUTINE mgInit_cc
!
! ----------------------------------------------
!
SUBROUTINE mgDefectNorm_cc(II)
INTEGER II,ndof

IF (myid.ne.0) THEN
 mgLev = NLMAX
 ndof  = SIZE(myMG%X(mgLev)%x)
 CALL LL21(myMG%AUX(II)%x,ndof,mgDefNorm)
END IF

CALL COMM_Maximum(mgDefNorm)

END SUBROUTINE mgDefectNorm_cc
!
! ----------------------------------------------
!
SUBROUTINE mgComputeDefect_cc(imgLev,bDef)
INTEGER i,j,ndof_u,ndof_p,imgLev,neq,jVelo,jPres
REAL*8  daux_u(3),daux_v(3),dt,dP,daux_p,DefNorm(4)
LOGICAL bDef

IF (myid.ne.0) THEN
 mgLev = imgLev
 ndof_u  = KNEL(mgLev) + KNVT(mgLev) + KNET(mgLev) + KNAT(mgLev)
 ndof_p  = 4*KNEL(mgLev)
 dt = TSTEP

 DO i=1,ndof_u

  daux_u(1) = MyMG%B_u(mgLev)%x(         i)
  daux_u(2) = MyMG%B_u(mgLev)%x(  ndof_u+i)
  daux_u(3) = MyMG%B_u(mgLev)%x(2*ndof_u+i)

  DO j=MyMG%Lq(mgLev)%LdA(i),MyMG%Lq(mgLev)%LdA(i+1)-1
   jVelo = MyMG%Lq(mgLev)%ColA(j)
   daux_v(1) = MyMG%dX_u(mgLev)%x(         jVelo)
   daux_v(2) = MyMG%dX_u(mgLev)%x(  ndof_u+jVelo)
   daux_v(3) = MyMG%dX_u(mgLev)%x(2*ndof_u+jVelo)
   daux_u(1) = daux_u(1) - (MyMG%A11(mgLev)%a(j)*daux_v(1) + MyMG%A12(mgLev)%a(j)*daux_v(2) + MyMG%A13(mgLev)%a(j)*daux_v(3))
   daux_u(2) = daux_u(2) - (MyMG%A21(mgLev)%a(j)*daux_v(1) + MyMG%A22(mgLev)%a(j)*daux_v(2) + MyMG%A23(mgLev)%a(j)*daux_v(3))
   daux_u(3) = daux_u(3) - (MyMG%A31(mgLev)%a(j)*daux_v(1) + MyMG%A32(mgLev)%a(j)*daux_v(2) + MyMG%A33(mgLev)%a(j)*daux_v(3))
  END DO

  DO j=MyMG%Lql(mgLev)%LdA(i),MyMG%Lql(mgLev)%LdA(i+1)-1
   jPres = MyMG%Lql(mgLev)%ColA(j)
   dP    = MyMG%dX_p(mgLev)%x(jPres)
   daux_u(1) = daux_u(1) + dt*MyMG%BX(mgLev)%a(j)*dP
   daux_u(2) = daux_u(2) + dt*MyMG%BY(mgLev)%a(j)*dP
   daux_u(3) = daux_u(3) + dt*MyMG%BZ(mgLev)%a(j)*dP
  END DO

  IF (myMG%KNPRU(i).EQ.1) daux_u(1) = 0d0
  IF (myMG%KNPRV(i).EQ.1) daux_u(2) = 0d0
  IF (myMG%KNPRW(i).EQ.1) daux_u(3) = 0d0
  MyMG%D_u(mgLev)%x(         i) = daux_u(1)
  MyMG%D_u(mgLev)%x(  ndof_u+i) = daux_u(2)
  MyMG%D_u(mgLev)%x(2*ndof_u+i) = daux_u(3)

 END DO

 CALL E013UVWSUM(MyMG%D_u(mgLev)%x)

 DO i=1,ndof_p
  daux_p = MyMG%B_p(mgLev)%x(i)
  DO j=MyMG%Llq(mgLev)%LdA(i),MyMG%Llq(mgLev)%LdA(i+1)-1
   jVelo = MyMG%Llq(mgLev)%ColA(j)
   daux_v(1) = MyMG%dX_u(mgLev)%x(         jVelo)
   daux_v(2) = MyMG%dX_u(mgLev)%x(  ndof_u+jVelo)
   daux_v(3) = MyMG%dX_u(mgLev)%x(2*ndof_u+jVelo)
   daux_p = daux_p - MyMG%BTX(mgLev)%a(j)*daux_v(1) - MyMG%BTY(mgLev)%a(j)*daux_v(2) - MyMG%BTZ(mgLev)%a(j)*daux_v(3)
  END DO
  MyMG%D_p(mgLev)%x(         i) = daux_p
 END DO
END IF

IF (bDef) THEN 
 IF (myid.ne.0) THEN
  CALL LL21(myMG%D_u(mgLev)%x(         1),ndof_u,DefNorm(1))
  CALL LL21(myMG%D_u(mgLev)%x(  ndof_u+1),ndof_u,DefNorm(2))
  CALL LL21(myMG%D_u(mgLev)%x(2*ndof_u+1),ndof_u,DefNorm(3))
  CALL LL21(myMG%D_p(mgLev)%x(         1),ndof_p,DefNorm(4))
 END IF
 CALL COMM_Maximum(DefNorm(1))
 CALL COMM_Maximum(DefNorm(2))
 CALL COMM_Maximum(DefNorm(3))
 CALL COMM_Maximum(DefNorm(4))

 IF (myid.eq.showid) WRITE(*,'(A,I3,4ES12.3)') "linear Defect",mgLev,DefNorm
 IF (IterCycle.eq.0) THEN
  myMG%DefInitial = MAX(DefNorm(1),DefNorm(2),DefNorm(3),DefNorm(4))
 ELSE
  mgDefNorm = MAX(DefNorm(1),DefNorm(2),DefNorm(3),DefNorm(4))
 END IF 
 
END IF

END SUBROUTINE mgComputeDefect_cc
!
! ----------------------------------------------
!
SUBROUTINE mgUpdateDefect_cc()
INTEGER i,j,ndof_u,ndof_p,neq,jVelo,jPres
REAL*8  daux_u(3),daux_v(3),dt,dP,daux_p,DefNorm(4)

IF (myid.ne.0) THEN
 ndof_u  = KNEL(mgLev) + KNVT(mgLev) + KNET(mgLev) + KNAT(mgLev)
 ndof_p  = 4*KNEL(mgLev)
 dt = TSTEP

 DO i=1,ndof_u

  daux_u(1) = 0d0
  daux_u(2) = 0d0
  daux_u(3) = 0d0

  DO j=MyMG%Lq(mgLev)%LdA(i),MyMG%Lq(mgLev)%LdA(i+1)-1
   jVelo = MyMG%Lq(mgLev)%ColA(j)
   daux_v(1) = MyMG%dX_u(mgLev)%x(         jVelo)
   daux_v(2) = MyMG%dX_u(mgLev)%x(  ndof_u+jVelo)
   daux_v(3) = MyMG%dX_u(mgLev)%x(2*ndof_u+jVelo)
   daux_u(1) = daux_u(1) + (MyMG%A11(mgLev)%a(j)*daux_v(1) + MyMG%A12(mgLev)%a(j)*daux_v(2) + MyMG%A13(mgLev)%a(j)*daux_v(3))
   daux_u(2) = daux_u(2) + (MyMG%A21(mgLev)%a(j)*daux_v(1) + MyMG%A22(mgLev)%a(j)*daux_v(2) + MyMG%A23(mgLev)%a(j)*daux_v(3))
   daux_u(3) = daux_u(3) + (MyMG%A31(mgLev)%a(j)*daux_v(1) + MyMG%A32(mgLev)%a(j)*daux_v(2) + MyMG%A33(mgLev)%a(j)*daux_v(3))
  END DO

  DO j=MyMG%Lql(mgLev)%LdA(i),MyMG%Lql(mgLev)%LdA(i+1)-1
   jPres = MyMG%Lql(mgLev)%ColA(j)
   dP    = MyMG%dX_p(mgLev)%x(jPres)
   daux_u(1) = daux_u(1) - dt*MyMG%BX(mgLev)%a(j)*dP
   daux_u(2) = daux_u(2) - dt*MyMG%BY(mgLev)%a(j)*dP
   daux_u(3) = daux_u(3) - dt*MyMG%BZ(mgLev)%a(j)*dP
  END DO

  IF (myMG%KNPRU(i).EQ.1) daux_u(1) = 0d0
  IF (myMG%KNPRV(i).EQ.1) daux_u(2) = 0d0
  IF (myMG%KNPRW(i).EQ.1) daux_u(3) = 0d0
  MyMG%A_u(mgLev)%x(         i) = daux_u(1)
  MyMG%A_u(mgLev)%x(  ndof_u+i) = daux_u(2)
  MyMG%A_u(mgLev)%x(2*ndof_u+i) = daux_u(3)

 END DO

 CALL E013UVWSUM(MyMG%A_u(mgLev)%x)

 DO i=1,ndof_p
  daux_p = 0d0
  DO j=MyMG%Llq(mgLev)%LdA(i),MyMG%Llq(mgLev)%LdA(i+1)-1
   jVelo = MyMG%Llq(mgLev)%ColA(j)
   daux_v(1) = MyMG%dX_u(mgLev)%x(         jVelo)
   daux_v(2) = MyMG%dX_u(mgLev)%x(  ndof_u+jVelo)
   daux_v(3) = MyMG%dX_u(mgLev)%x(2*ndof_u+jVelo)
   daux_p = daux_p + MyMG%BTX(mgLev)%a(j)*daux_v(1) &
                   + MyMG%BTY(mgLev)%a(j)*daux_v(2) &   
                   + MyMG%BTZ(mgLev)%a(j)*daux_v(3)    
  END DO
  MyMG%A_p(mgLev)%x(         i) = daux_p
 END DO

!  MyMG%D_u(mgLev)%x = MyMG%B_u(mgLev)%x
!  CALL E013UVWSUM(MyMG%D_u(mgLev)%x)
!  MyMG%D_p(mgLev)%x = MyMG%B_p(mgLev)%x

 DO i=1,ndof_u
  MyMG%D_u(mgLev)%x(         i) = MyMG%D_u(mgLev)%x(         i) - MyMG%A_u(mgLev)%x(         i)
  MyMG%D_u(mgLev)%x(  ndof_u+i) = MyMG%D_u(mgLev)%x(  ndof_u+i) - MyMG%A_u(mgLev)%x(  ndof_u+i)
  MyMG%D_u(mgLev)%x(2*ndof_u+i) = MyMG%D_u(mgLev)%x(2*ndof_u+i) - MyMG%A_u(mgLev)%x(2*ndof_u+i)
  MyMG%X_u(mgLev)%x(         i) = MyMG%X_u(mgLev)%x(         i) + MyMG%dX_u(mgLev)%x(         i)
  MyMG%X_u(mgLev)%x(  ndof_u+i) = MyMG%X_u(mgLev)%x(  ndof_u+i) + MyMG%dX_u(mgLev)%x(  ndof_u+i)
  MyMG%X_u(mgLev)%x(2*ndof_u+i) = MyMG%X_u(mgLev)%x(2*ndof_u+i) + MyMG%dX_u(mgLev)%x(2*ndof_u+i)
 END DO

 DO i=1,ndof_p
  MyMG%D_p(mgLev)%x(         i) = MyMG%D_p(mgLev)%x(         i) - MyMG%A_p(mgLev)%x(         i)
  MyMG%X_p(mgLev)%x(         i) = MyMG%X_p(mgLev)%x(         i) + MyMG%dX_p(mgLev)%x(         i)
 END DO

END IF

IF (myid.ne.0) THEN
 CALL LL21(myMG%D_u(mgLev)%x(         1),ndof_u,DefNorm(1))
 CALL LL21(myMG%D_u(mgLev)%x(  ndof_u+1),ndof_u,DefNorm(2))
 CALL LL21(myMG%D_u(mgLev)%x(2*ndof_u+1),ndof_u,DefNorm(3))
 CALL LL21(myMG%D_p(mgLev)%x(         1),ndof_p,DefNorm(4))
END IF
CALL COMM_Maximum(DefNorm(1))
CALL COMM_Maximum(DefNorm(2))
CALL COMM_Maximum(DefNorm(3))
CALL COMM_Maximum(DefNorm(4))

IF (myid.eq.showid) WRITE(*,'(A,I3,4ES12.3)') "DefNorm",mgLev,DefNorm

END SUBROUTINE mgUpdateDefect_cc
!
! ----------------------------------------------
!
SUBROUTINE mgccRestriction_cc(DU1,DU2,DP1,DP2,KNPRU,KNPRV,KNPRW)
REAL*8  DU1(*),DU2(*),DP1(*),DP2(*)
INTEGER KNPRU(*),KNPRV(*),KNPRW(*)

IF (myid.ne.0) THEN
  CALL E012_Restriction(DP2,DP1,mg_E012Prol(mgLev)%a,KNEL(mgLev))
  CALL E013_Restriction(DU2,DU1,mg_E013Rest(mgLev)%a,&
       mg_E013RestM(mgLev)%LdA,mg_E013RestM(mgLev)%ColA,KNPRU,KNPRV,KNPRW)
END IF

END SUBROUTINE mgccRestriction_cc
!
! ----------------------------------------------
!
SUBROUTINE mgRestriction_cc()
INTEGER I,NDOF

IF (myid.ne.0) THEN

 CALL ZTIME(time0)
 CALL E012_Restriction(myMG%D_p(mgLev+1)%x,myMG%B_p(mgLev)%x,mg_E012Prol(mgLev)%a,KNEL(mgLev))
 CALL ZTIME(time1)
 myStat%tRestP = myStat%tRestP + (time1-time0)

 CALL ZTIME(time0)
 CALL E013_Restriction(myMG%D_u(mgLev+1)%x,myMG%B_u(mgLev)%x,mg_E013Rest(mgLev)%a,&
      mg_E013RestM(mgLev)%LdA,mg_E013RestM(mgLev)%ColA,myMG%KNPRU,myMG%KNPRV,myMG%KNPRW)
 CALL ZTIME(time1)
 myStat%tRestUVW = myStat%tRestUVW + (time1-time0)

END IF

END SUBROUTINE mgRestriction_cc
!
! ----------------------------------------------
!
SUBROUTINE mgccProlongation_cc(DU1,DU2,DP1,DP2,KNPRU,KNPRV,KNPRW)
INTEGER I,NDOF
REAL*8  DU1(*),DU2(*),DP1(*),DP2(*)
INTEGER KNPRU(*),KNPRV(*),KNPRW(*)

IF (myid.ne.0) THEN
  CALL E012_Prolongation(DP2,DP1,mg_E012Prol(mgLev-1)%a,KNEL(mgLev-1))
  CALL E013_Prolongation(DU2,DU1,mg_E013Prol(mgLev-1)%a,&
       mg_E013ProlM(mgLev-1)%LdA,mg_E013ProlM(mgLev-1)%ColA,KNPRU,KNPRV,KNPRW)
END IF

END SUBROUTINE mgccProlongation_cc
!
! ----------------------------------------------
!
SUBROUTINE mgProlongation_cc()
INTEGER I,NDOF
REAL*8, ALLOCATABLE :: aux_X1(:),aux_X2(:)
EXTERNAL E013

IF (myid.ne.0) THEN
  CALL ZTIME(time0)
  CALL E012_Prolongation(myMG%A_p(mgLev)%x,myMG%dX_p(mgLev-1)%x,mg_E012Prol(mgLev-1)%a,KNEL(mgLev-1))
  CALL ZTIME(time1)
  myStat%tProlP = myStat%tProlP + (time1-time0)

  CALL ZTIME(time0)
  CALL E013_Prolongation(myMG%A_u(mgLev)%x,myMG%dX_u(mgLev-1)%x,mg_E013Prol(mgLev-1)%a,&
       mg_E013ProlM(mgLev-1)%LdA,mg_E013ProlM(mgLev-1)%ColA,myMG%KNPRU,myMG%KNPRV,myMG%KNPRW)
  CALL ZTIME(time1)
  myStat%tProlUVW = myStat%tProlUVW + (time1-time0)
END IF

END SUBROUTINE mgProlongation_cc
!
! ----------------------------------------------
!
SUBROUTINE mgCoarseGridSolverShouldBe_cc()
INTEGER Iter,i,j,ndof_u,ndof_p,neq
REAL*8 daux

CALL ZTIME(time0)
ILEV = mgLev
CALL SETLEV(2)

ndof_u  = KNEL(mgLev) + KNVT(mgLev) + KNET(mgLev) + KNAT(mgLev)
ndof_p  = 4*KNEL(mgLev)

CALL MasterVanka()

CALL ZTIME(time1)
myStat%tSolvP = myStat%tSolvP + (time1-time0)

END SUBROUTINE mgCoarseGridSolverShouldBe_cc
!
! ----------------------------------------------
!
 SUBROUTINE mgCoarseGridSolver_cc()
 INTEGER Iter,i,j,ndof_u,ndof_p,neq
 REAL*8 daux
 
  IF (myid.ne.0) THEN
   CALL ZTIME(time0)
 
   ILEV = mgLev
   CALL SETLEV(2)
 
   ndof_u  = KNEL(mgLev) + KNVT(mgLev) + KNET(mgLev) + KNAT(mgLev)
   ndof_p  = 4*KNEL(mgLev)
 
   myMG%dX_u(mgLev)%x = 0d0
   myMG%dX_p(mgLev)%x = 0d0
 
   ! the pressure rhs is constant .... 
   MyMG%D_p(mgLev)%x = MyMG%B_p(mgLev)%x
 
   DO iter = 1,MyMG%nIterCoarse
 
    ! the velocity rhs needs to be updated due to the parallelization ...
    CALL mgUpdateDefectForVanka()

    CALL VANKA(mg_mesh%level(ilev)%kvert,&
               mg_mesh%level(ilev)%karea,&
               mg_mesh%level(ilev)%kedge,&
               ndof_u,ndof_p,&
               myMG%dX_u(mgLev)%x,&
               myMG%dX_p(mgLev)%x,&
               myMG%D_u(mgLev)%x,&
               myMG%D_p(mgLev)%x,&
               myMG%KNPRU,&
               myMG%KNPRV,&
               myMG%KNPRW)

    END DO
 
   CALL ZTIME(time1)
   myStat%tSmthP = myStat%tSmthP + (time1-time0)
 
  END IF
 
 END SUBROUTINE mgCoarseGridSolver_cc
!
! ----------------------------------------------
!
SUBROUTINE mgSmoother_cc()
INTEGER Iter,i,j,ndof_u,ndof_p,neq
REAL*8 daux

! RETURN
 IF (myid.ne.0) THEN
  CALL ZTIME(time0)

  ILEV = mgLev
  CALL SETLEV(2)

  ndof_u  = KNEL(mgLev) + KNVT(mgLev) + KNET(mgLev) + KNAT(mgLev)
  ndof_p  = 4*KNEL(mgLev)

  ! the pressure rhs is constant .... 
  MyMG%D_p(mgLev)%x = MyMG%B_p(mgLev)%x

!  IF (myid.eq.1) WRITE(*,*) "operating on Level: ", mgLev
  DO iter = 1,MyMG%nSmootherSteps

!    ! the velocity rhs needs to be updated due to the parallelization ...
   CALL mgUpdateDefectForVanka()

   IF (myMG%VANKA.eq.1) THEN
   CALL VANKA_new(mg_mesh%level(ilev)%kvert,&
                  mg_mesh%level(ilev)%karea,&
                  mg_mesh%level(ilev)%kedge,&
                  ndof_u,ndof_p,&
                  myMG%dX_u(mgLev)%x,&
                  myMG%dX_p(mgLev)%x,&
                  myMG%D_u(mgLev)%x,&
                  myMG%D_p(mgLev)%x,&
                  myMG%KNPRU,&
                  myMG%KNPRV,&
                  myMG%KNPRW)

   ELSE
   CALL VANKA(mg_mesh%level(ilev)%kvert,&
              mg_mesh%level(ilev)%karea,&
              mg_mesh%level(ilev)%kedge,&
              ndof_u,ndof_p,&
              myMG%dX_u(mgLev)%x,&
              myMG%dX_p(mgLev)%x,&
              myMG%D_u(mgLev)%x,&
              myMG%D_p(mgLev)%x,&
              myMG%KNPRU,&
              myMG%KNPRV,&
              myMG%KNPRW)
   END IF

 ! ----------------------------------------------------------------------------------------------------- !
  myMG%A_u(mgLev)%x = 1d0
  CALL E013SUM(myMG%A_u(mgLev)%x)
  CALL E013UVWSUM(myMG%dX_u(mgLev)%x)

  DO i=1,3*ndof_u
   daux = myMG%A_u(mgLev)%x(mod(i-1,ndof_u)+1)
   myMG%dX_u(mgLev)%x(         i) = myMG%dX_u(mgLev)%x(         i)/daux
  END DO
! ----------------------------------------------------------------------------------------------------- !

  END DO

  CALL ZTIME(time1)
  myStat%tSmthP = myStat%tSmthP + (time1-time0)

 END IF

END SUBROUTINE mgSmoother_cc
!
! ----------------------------------------------
!
SUBROUTINE CG_Activation()
INTEGER ndof,ndof_orig

 ndof  = SIZE(myMG%X(MyMG%MaxLev)%x)
 IF (myid.eq.0) ndof  = SIZE(myMG%X(MyMG%MinLev)%x) 

 IF (.NOT.myCG%bActivated) THEN
!   IF (myid.ne.0) THEN 
   ALLOCATE(myCG%d1(ndof))
   ALLOCATE(myCG%d2(ndof))
   ALLOCATE(myCG%d3(ndof))
   ALLOCATE(myCG%d4(ndof))
   ALLOCATE(myCG%d5(ndof))
!   END IF
  myCG%bActivated = .TRUE.
 ELSE
  ndof_orig  = SIZE(myCG%d1)
  IF (ndof_orig.LT.ndof) THEN
   DEALLOCATE(myCG%d1,myCG%d2,myCG%d3,myCG%d4,myCG%d5)
   ALLOCATE(myCG%d1(ndof))
   ALLOCATE(myCG%d2(ndof))
   ALLOCATE(myCG%d3(ndof))
   ALLOCATE(myCG%d4(ndof))
   ALLOCATE(myCG%d5(ndof))
  END IF
 END IF

END SUBROUTINE CG_Activation
!
! ----------------------------------------------
!
SUBROUTINE JCB_Activation()
INTEGER ndof,ndof_orig

 ndof  = SIZE(myMG%X(MyMG%MaxLev)%x)

 IF (.NOT.myJCB%bActivated) THEN
  IF (myid.ne.0) THEN 
   ALLOCATE(myJCB%d1(ndof))
  END IF
  myJCB%bActivated = .TRUE.
 ELSE
  ndof_orig  = SIZE(myJCB%d1)
  IF (ndof_orig.LT.ndof.AND.myid.ne.0) THEN
   DEALLOCATE(myJCB%d1)
   ALLOCATE(myJCB%d1(ndof))
  END IF
 END IF

END SUBROUTINE JCB_Activation
!
! ----------------------------------------------
!
SUBROUTINE mgProlRestInit

 IF (myMG%bProlRest) RETURN

 IF (myid.ne.0) THEN
  IF (MyMG%cVariable.EQ."Pressure") THEN
   DO mgLev = myMG%MinLev+1, myMG%MaxLev
    ILEV = mgLev
    CALL InitE012ProlMat(mg_E012Prol(mgLev-1)%a,KWORK(L(KLADJ(mgLev))),&
         KWORK(L(KLVERT(mgLev))),DWORK(L(KLCVG(mgLev))),KNEL(mgLev-1))
   END DO
  END IF

  IF (MyMG%cVariable.EQ."Velocity") THEN
   DO mgLev = myMG%MinLev+1, myMG%MaxLev
    ILEV = mgLev
    CALL InitE013ProlMat(mg_E013ProlM(mgLev-1)%na,&
         mg_E013Prol(mgLev-1)%a,mg_E013ProlM(mgLev-1)%LdA,mg_E013ProlM(mgLev-1)%ColA,&
         mg_E013Rest(mgLev-1)%a,mg_E013RestM(mgLev-1)%LdA,mg_E013RestM(mgLev-1)%ColA,&
         KWORK(L(KLADJ(mgLev-1))),KWORK(L(KLVERT(mgLev-1))),KWORK(L(KLEDGE(mgLev-1))),&
         KWORK(L(KLAREA(mgLev-1))),KNVT(mgLev-1),KNET(mgLev-1),KNAT(mgLev-1),KNEL(mgLev-1),&
         KWORK(L(KLADJ(mgLev))),KWORK(L(KLVERT(mgLev))),KWORK(L(KLEDGE(mgLev))),&
         KWORK(L(KLAREA(mgLev))),KNVT(mgLev),KNET(mgLev),KNAT(mgLev),KNEL(mgLev))
   END DO
  END IF

 ELSE
 IF (MyMG%cVariable.EQ."Pressure") THEN
  DO mgLev = myMG%MinLev+1, myMG%MinLev
   ILEV = mgLev
   CALL InitE012ProlMat(mg_E012Prol(mgLev-1)%a,KWORK(L(KLADJ(mgLev))),&
        KWORK(L(KLVERT(mgLev))),DWORK(L(KLCVG(mgLev))),KNEL(mgLev-1))
  END DO
!    write(*,*) mg_E012Prol(mgLev-1)%a
!    write(*,*) 'asdasdas'
!    pause
 END IF
 
 END IF

 MyMG%bProlRest=.TRUE.

END SUBROUTINE mgProlRestInit
!
! ============================= E012 ELEMENT =================================
!
SUBROUTINE E012_SSOR(DX,DXP,DB,NEQ,nIter)
INTEGER NEQ,nIter
REAL*8 DX(*),DXP(*),DB(*)
REAL*8 drlx
INTEGER i

IF (myid.ne.0) THEN
 drlx=0.5d0
 ILEV = mgLev

 DO i=1,nIter
  CALL GetParPressure(DX,DXP)
  CALL E012_SOR_Smoother(DB,myMG%A(mgLev)%a,myMG%AP(mgLev)%a,myMG%L(mgLEV)%ColA,&
      myMG%LP(mgLEV)%ColA,myMG%L(mgLEV)%LdA,myMG%LP(mgLEV)%LdA,&
      DX,DXP,neq,drlx)
  CALL GetParPressure(DX,DXP)
  CALL E012_ROS_Smoother(DB,myMG%A(mgLev)%a,myMG%AP(mgLev)%a,myMG%L(mgLEV)%ColA,&
      myMG%LP(mgLEV)%ColA,myMG%L(mgLEV)%LdA,myMG%LP(mgLEV)%LdA,&
      DX,DXP,neq,drlx)
 END DO
END IF

END SUBROUTINE E012_SSOR
!
! ----------------------------------------------
!
SUBROUTINE E012_SOR(DX,DXP,DB,NEQ,nIter)
INTEGER NEQ,nIter
REAL*8 DX(*),DXP(*),DB(*)
REAL*8 drlx
INTEGER i

IF (myid.ne.0) THEN
 drlx=0.9d0
 ILEV = mgLev

 DO i=1,2*nIter
  CALL GetParPressure(DX,DXP)
  CALL E012_SOR_Smoother(DB,myMG%A(mgLev)%a,myMG%AP(mgLev)%a,myMG%L(mgLEV)%ColA,&
      myMG%LP(mgLEV)%ColA,myMG%L(mgLEV)%LdA,myMG%LP(mgLEV)%LdA,&
      DX,DXP,neq,drlx)
 END DO
END IF

END SUBROUTINE E012_SOR
!
! ----------------------------------------------
!
SUBROUTINE crs_E012_SOR(DX,DB,NEQ,nIter)
INTEGER NEQ,nIter
REAL*8 DX(*),DB(*)
REAL*8 drlx
INTEGER i

drlx=1.25d0
ILEV = mgLev

DO i=1,2*nIter
 CALL crs_E012_SOR_Smoother(DB,myMG%A(mgLev)%a,myMG%L(mgLEV)%ColA,&
     myMG%L(mgLEV)%LdA,DX,neq,drlx)
END DO

END SUBROUTINE crs_E012_SOR
!
! ----------------------------------------------
!
SUBROUTINE E012_DCG(DX,DXP,NEQ)
INTEGER NEQ
REAL*8 DX(*),DXP(*)
REAL*8 drlx

drlx=1.0d0
ILEV = mgLev
CALL GetParPressure(DX,DXP)
CALL PARID117(myMG%A(mgLev)%a,myMG%AP(mgLev)%a,myMG%L(mgLEV)%ColA,&
     myMG%LP(mgLEV)%ColA,myMG%L(mgLEV)%LdA,myMG%LP(mgLEV)%LdA,&
     dx,dxp,neq,drlx)

END SUBROUTINE E012_DCG
!
! ----------------------------------------------
!
SUBROUTINE E012_DCG_Master(DX,NEQ)
INTEGER NEQ
REAL*8 DX(*)
REAL*8 drlx

drlx=0.5d0
ILEV = mgLev
CALL PARID117_Master(myMG%A(mgLev)%a,myMG%L(mgLEV)%ColA,&
     myMG%L(mgLEV)%LdA,dx,neq,drlx)

END SUBROUTINE E012_DCG_Master
!
! ----------------------------------------------
!
SUBROUTINE E012_DAX(DX,DXP,DA,NEQ,A1,A2)
IMPLICIT NONE
REAL*8  DX(*),DXP(*),DA(*),A1,A2
INTEGER NEQ
INTEGER IEQ,IA,ICOL

ILEV = mgLev
CALL GetParPressure(DX,DXP)

DO IEQ=1,NEQ
 DA(IEQ) = A2*DA(IEQ)
END DO
DO IEQ=1,NEQ
 DO IA=myMG%L(mgLEV)%LdA(IEQ),myMG%L(mgLEV)%LdA(IEQ+1)-1
  ICOL = myMG%L(mgLEV)%ColA(IA)
  DA(IEQ) = DA(IEQ) + A1*myMG%A(mgLEV)%a(IA)*DX(ICOL)
 END DO
 DO IA=myMG%LP(mgLEV)%LdA(IEQ),myMG%LP(mgLEV)%LdA(IEQ+1)-1
  ICOL = myMG%LP(mgLEV)%ColA(IA)
  DA(IEQ) = DA(IEQ) + A1*myMG%AP(mgLEV)%a(IA)*DXP(ICOL)
 END DO
END DO

END SUBROUTINE E012_DAX
!
! ----------------------------------------------
!
SUBROUTINE E012_DAX_Master(DX,DA,NEQ,A1,A2)
IMPLICIT NONE
REAL*8  DX(*),DA(*),A1,A2
INTEGER NEQ
INTEGER IEQ,IA,ICOL

ILEV = mgLev

DO IEQ=1,NEQ
 DA(IEQ) = A2*DA(IEQ)
END DO
DO IEQ=1,NEQ
 DO IA=myMG%L(mgLEV)%LdA(IEQ),myMG%L(mgLEV)%LdA(IEQ+1)-1
  ICOL = myMG%L(mgLEV)%ColA(IA)
  DA(IEQ) = DA(IEQ) + A1*myMG%A(mgLEV)%a(IA)*DX(ICOL)
 END DO
END DO

END SUBROUTINE E012_DAX_Master
!
! ----------------------------------------------
!
SUBROUTINE InitE012ProlMat(A,KADJ2,KVERT2,DCORVG2,N)
REAL*8 A(8,4,4,*)
INTEGER KVERT2(8,*),KADJ2(6,*),N
REAL*8 DCORVG2(3,*)
REAL*8 DAUX(5),DX,DY,DZ,DV,DX0,DY0,DZ0,DXP,DYP,DZP
INTEGER IEL,JEL(8),I,J,K

DO IEL=1,N

 A(:,:,:,IEL) = 0d0

 DX0 = DCORVG2(1,KVERT2(7,IEL))
 DY0 = DCORVG2(2,KVERT2(7,IEL))
 DZ0 = DCORVG2(3,KVERT2(7,IEL))

 JEL(1)  = IEL
 JEL(2)  = KADJ2(3,JEL(1))
 JEL(3)  = KADJ2(3,JEL(2))
 JEL(4)  = KADJ2(3,JEL(3))
 JEL(5)  = KADJ2(6,JEL(1))
 JEL(6)  = KADJ2(6,JEL(2))
 JEL(7)  = KADJ2(6,JEL(3))
 JEL(8)  = KADJ2(6,JEL(4))

 DO J=1,8
  DXP = 0d0
  DYP = 0d0
  DZP = 0d0
  DO K=1,8
   DXP = DXP + 0.125D0*DCORVG2(1,KVERT2(K,JEL(J)))
   DYP = DYP + 0.125D0*DCORVG2(2,KVERT2(K,JEL(J)))
   DZP = DZP + 0.125D0*DCORVG2(3,KVERT2(K,JEL(J)))
  END DO

  A(J,1,1,IEL) = 1d0
  A(J,1,2,IEL) = DXP-DX0
  A(J,1,3,IEL) = DYP-DY0
  A(J,1,4,IEL) = DZP-DZ0

  A(J,2,2,IEL) = 1d0

  A(J,3,3,IEL) = 1d0

  A(J,4,4,IEL) = 1d0

 END DO

END DO

END SUBROUTINE InitE012ProlMat
!
! ----------------------------------------------
!
SUBROUTINE E012_Prolongation(D2,D1,A,N)
IMPLICIT NONE
REAL*8 A(8,4,4,*)
INTEGER N
REAL*8 D1(4,*),D2(4,*)
INTEGER IEL,JEL(8),I,J,K,KEL

DO IEL=1,N
 JEL(1)  = IEL
 JEL(2)  = N + 7*(IEL-1) + 1 !KADJ2(3,JEL(1))
 JEL(3)  = JEL(2) + 1        !N + 7*(IEL-1) + 2 !KADJ2(3,JEL(2))
 JEL(4)  = JEL(3) + 1        !N + 7*(IEL-1) + 3 !KADJ2(3,JEL(3))
 JEL(5)  = JEL(4) + 1        !N + 7*(IEL-1) + 4 !KADJ2(6,JEL(1))
 JEL(6)  = JEL(5) + 1        !N + 7*(IEL-1) + 5 !KADJ2(6,JEL(2))
 JEL(7)  = JEL(6) + 1        !N + 7*(IEL-1) + 6 !KADJ2(6,JEL(3))
 JEL(8)  = JEL(7) + 1        !N + 7*(IEL-1) + 7 !KADJ2(6,JEL(4))

 DO I=1,4
  DO J=1,8
   KEL = JEL(J)
   D2(I,KEL) = 0d0
   DO K=1,4
    D2(I,KEL) =  D2(I,KEL) + D1(K,IEL)*A(J,I,K,IEL)
   END DO
  END DO
 END DO
END DO

END SUBROUTINE E012_Prolongation
!
! ----------------------------------------------
!
SUBROUTINE E012_Restriction(D2,D1,A,N)
IMPLICIT NONE
REAL*8 A(8,4,4,*)
INTEGER N
REAL*8 D1(4,*),D2(4,*)
INTEGER IEL,JEL(8),I,J,K,KEL

DO IEL=1,N
 JEL(1)  = IEL
 JEL(2)  = N + 7*(IEL-1) + 1 !KADJ2(3,JEL(1))
 JEL(3)  = JEL(2) + 1        !N + 7*(IEL-1) + 2 !KADJ2(3,JEL(2))
 JEL(4)  = JEL(3) + 1        !N + 7*(IEL-1) + 3 !KADJ2(3,JEL(3))
 JEL(5)  = JEL(4) + 1        !N + 7*(IEL-1) + 4 !KADJ2(6,JEL(1))
 JEL(6)  = JEL(5) + 1        !N + 7*(IEL-1) + 5 !KADJ2(6,JEL(2))
 JEL(7)  = JEL(6) + 1        !N + 7*(IEL-1) + 6 !KADJ2(6,JEL(3))
 JEL(8)  = JEL(7) + 1        !N + 7*(IEL-1) + 7 !KADJ2(6,JEL(4))

 DO I=1,4
  D1(I,IEL) = 0d0
  DO J=1,8
   KEL = JEL(J)
   DO K=1,4
    D1(I,IEL) =  D1(I,IEL) + D2(K,KEL)*A(J,I,K,IEL)
   END DO
  END DO
 END DO
END DO

END SUBROUTINE E012_Restriction
!
! ============================================================================
!
SUBROUTINE E013_Prolongation(D2,D1,A,KLD,KCOL,KNPRU,KNPRV,KNPRW)
IMPLICIT NONE
REAL*8 A(*)
INTEGER KLD(*),KCOL(*),KNPRU(*),KNPRV(*),KNPRW(*)
REAL*8 D1(*),D2(*)
INTEGER I,J,ICOL
INTEGER MEQ,MEQ1,MEQ2,MEQ3
INTEGER LEQ,LEQ1,LEQ2,LEQ3

IF (myid.eq.0) RETURN
MEQ = KNVT(mgLev  )+KNAT(mgLev  )+KNET(mgLev  )+KNEL(mgLev  )
LEQ = KNVT(mgLev-1)+KNAT(mgLev-1)+KNET(mgLev-1)+KNEL(mgLev-1)

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
SUBROUTINE E013_ProlongationX(D2,D1,A,KLD,KCOL)
IMPLICIT NONE
REAL*8 A(*)
INTEGER KLD(*),KCOL(*)
REAL*8 D1(*),D2(*)
INTEGER I,J,ICOL
INTEGER MEQ
INTEGER LEQ

IF (myid.eq.0) RETURN
MEQ = KNVT(mgLev  )+KNAT(mgLev  )+KNET(mgLev  )+KNEL(mgLev  )
LEQ = KNVT(mgLev-1)+KNAT(mgLev-1)+KNET(mgLev-1)+KNEL(mgLev-1)

DO I=1,MEQ
 D2(I) = 0d0
 DO J=KLD(I),KLD(I+1)-1
  ICOL = KCOL(J)
  D2(I) = D2(I) + A(J)*D1(ICOL)
 END DO
END DO

END SUBROUTINE
!
! ----------------------------------------------
!
SUBROUTINE E013_Restriction(D2,D1,A,KLD,KCOL,KNPRU,KNPRV,KNPRW)
IMPLICIT NONE
REAL*8 A(*)
INTEGER N,KLD(*),KCOL(*),KNPRU(*),KNPRV(*),KNPRW(*)
REAL*8 D1(*),D2(*)
INTEGER I,J,ICOL
INTEGER MEQ,MEQ1,MEQ2,MEQ3
INTEGER LEQ,LEQ1,LEQ2,LEQ3

IF (myid.eq.0) RETURN
MEQ = KNVT(mgLev  )+KNAT(mgLev  )+KNET(mgLev  )+KNEL(mgLev  )
LEQ = KNVT(mgLev+1)+KNAT(mgLev+1)+KNET(mgLev+1)+KNEL(mgLev+1)

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
SUBROUTINE InitE013ProlMat(NA,PA,PLDA,PCOLA,RA,RLDA,RCOLA,&
   KADJ1,KVERT1,KEDGE1,KAREA1,NVT1,NET1,NAT1,NEL1,&
   KADJ2,KVERT2,KEDGE2,KAREA2,NVT2,NET2,NAT2,NEL2)

IMPLICIT NONE
REAL*8 PA(*),RA(*)
INTEGER NA,PLDA(*),PCOLA(*),RLDA(*),RCOLA(*)
INTEGER, ALLOCATABLE, DIMENSION(:) :: RLDB,PLDB
REAL*8 A(125,27)
INTEGER KVERT2(8,*),KADJ2(6,*),KEDGE2(12,*),KAREA2(6,*),NVT2,NET2,NAT2,NEL2
INTEGER KVERT1(8,*),KADJ1(6,*),KEDGE1(12,*),KAREA1(6,*),NVT1,NET1,NAT1,NEL1
INTEGER IEL,JEL(8),I,J,II,JJ,KK,IND_I(125),IND_J(27),II_POS,JJ_POS
INTEGER N00,NV1,NA1,NE1,NN1,NV2,NA2,NE2,NN2
REAL*8  DAUX

CALL InitE013ElemProlMat(A)

N00 = 0

NV1 = NVT1 + N00
NE1 = NET1 + NV1
NA1 = NAT1 + NE1
NN1 = NEL1 + NA1

NV2 = NVT2 + N00
NE2 = NET2 + NV2
NA2 = NAT2 + NE2
NN2 = NEL2 + NA2

PA(1:NA) = 0d0
RA(1:NA) = 0d0
PCOLA(1:NA) = 0
RCOLA(1:NA) = 0
PLDA(1:NN2+1) = 0
RLDA(1:NN1+1) = 0

DO IEL=1,NEL1
 JEL(1)  = IEL
 JEL(2)  = KADJ2(3,JEL(1))
 JEL(3)  = KADJ2(3,JEL(2))
 JEL(4)  = KADJ2(3,JEL(3))
 JEL(5)  = KADJ2(6,JEL(1))
 JEL(6)  = KADJ2(6,JEL(2))
 JEL(7)  = KADJ2(6,JEL(3))
 JEL(8)  = KADJ2(6,JEL(4))

 IND_I(  1) = KVERT2( 1,JEL(1)) + N00
 IND_I(  2) = KEDGE2( 1,JEL(1)) + NV2
 IND_I(  3) = KVERT2( 2,JEL(1)) + N00
 IND_I(  4) = KEDGE2( 4,JEL(2)) + NV2
 IND_I(  5) = KVERT2( 1,JEL(2)) + N00
 IND_I(  6) = KEDGE2( 4,JEL(1)) + NV2
 IND_I(  7) = KAREA2( 1,JEL(1)) + NE2
 IND_I(  8) = KEDGE2( 2,JEL(1)) + NV2
 IND_I(  9) = KAREA2( 1,JEL(2)) + NE2
 IND_I( 10) = KEDGE2( 1,JEL(2)) + NV2
 IND_I( 11) = KVERT2( 4,JEL(1)) + N00
 IND_I( 12) = KEDGE2( 3,JEL(1)) + NV2
 IND_I( 13) = KVERT2( 3,JEL(1)) + N00
 IND_I( 14) = KEDGE2( 2,JEL(2)) + NV2
 IND_I( 15) = KVERT2( 2,JEL(2)) + N00
 IND_I( 16) = KEDGE2( 1,JEL(4)) + NV2
 IND_I( 17) = KAREA2( 1,JEL(4)) + NE2
 IND_I( 18) = KEDGE2( 3,JEL(4)) + NV2
 IND_I( 19) = KAREA2( 1,JEL(3)) + NE2
 IND_I( 20) = KEDGE2( 4,JEL(3)) + NV2
 IND_I( 21) = KVERT2( 1,JEL(4)) + N00
 IND_I( 22) = KEDGE2( 4,JEL(4)) + NV2
 IND_I( 23) = KVERT2( 4,JEL(4)) + N00
 IND_I( 24) = KEDGE2( 1,JEL(3)) + NV2
 IND_I( 25) = KVERT2( 1,JEL(3)) + N00

 IND_I( 26) = KEDGE2( 5,JEL(1)) + NV2
 IND_I( 27) = KAREA2( 2,JEL(1)) + NE2
 IND_I( 28) = KEDGE2( 6,JEL(1)) + NV2
 IND_I( 29) = KAREA2( 5,JEL(2)) + NE2
 IND_I( 30) = KEDGE2( 5,JEL(2)) + NV2
 IND_I( 31) = KAREA2( 5,JEL(1)) + NE2
 IND_I( 32) =           JEL(1)  + NA2
 IND_I( 33) = KAREA2( 3,JEL(1)) + NE2
 IND_I( 34) =           JEL(2)  + NA2
 IND_I( 35) = KAREA2( 2,JEL(2)) + NE2
 IND_I( 36) = KEDGE2( 8,JEL(1)) + NV2
 IND_I( 37) = KAREA2( 4,JEL(1)) + NE2
 IND_I( 38) = KEDGE2( 7,JEL(1)) + NV2
 IND_I( 39) = KAREA2( 3,JEL(2)) + NE2
 IND_I( 40) = KEDGE2( 6,JEL(2)) + NV2
 IND_I( 41) = KAREA2( 2,JEL(4)) + NE2
 IND_I( 42) =           JEL(4)  + NA2
 IND_I( 43) = KAREA2( 4,JEL(4)) + NE2
 IND_I( 44) =           JEL(3)  + NA2
 IND_I( 45) = KAREA2( 5,JEL(3)) + NE2
 IND_I( 46) = KEDGE2( 5,JEL(4)) + NV2
 IND_I( 47) = KAREA2( 5,JEL(4)) + NE2
 IND_I( 48) = KEDGE2( 8,JEL(4)) + NV2
 IND_I( 49) = KAREA2( 2,JEL(3)) + NE2
 IND_I( 50) = KEDGE2( 5,JEL(3)) + NV2

 IND_I( 51) = KVERT2( 5,JEL(1)) + N00
 IND_I( 52) = KEDGE2( 9,JEL(1)) + NV2
 IND_I( 53) = KVERT2( 6,JEL(1)) + N00
 IND_I( 54) = KEDGE2(12,JEL(2)) + NV2
 IND_I( 55) = KVERT2( 5,JEL(2)) + N00
 IND_I( 56) = KEDGE2(12,JEL(1)) + NV2
 IND_I( 57) = KAREA2( 6,JEL(1)) + NE2
 IND_I( 58) = KEDGE2(10,JEL(1)) + NV2
 IND_I( 59) = KAREA2( 6,JEL(2)) + NE2
 IND_I( 60) = KEDGE2( 9,JEL(2)) + NV2
 IND_I( 61) = KVERT2( 8,JEL(1)) + N00
 IND_I( 62) = KEDGE2(11,JEL(1)) + NV2
 IND_I( 63) = KVERT2( 7,JEL(1)) + N00
 IND_I( 64) = KEDGE2(10,JEL(2)) + NV2
 IND_I( 65) = KVERT2( 6,JEL(2)) + N00
 IND_I( 66) = KEDGE2( 9,JEL(4)) + NV2
 IND_I( 67) = KAREA2( 6,JEL(4)) + NE2
 IND_I( 68) = KEDGE2(11,JEL(4)) + NV2
 IND_I( 69) = KAREA2( 6,JEL(3)) + NE2
 IND_I( 70) = KEDGE2(12,JEL(3)) + NV2
 IND_I( 71) = KVERT2( 5,JEL(4)) + N00
 IND_I( 72) = KEDGE2(12,JEL(4)) + NV2
 IND_I( 73) = KVERT2( 8,JEL(4)) + N00
 IND_I( 74) = KEDGE2( 9,JEL(3)) + NV2
 IND_I( 75) = KVERT2( 5,JEL(3)) + N00

 IND_I( 76) = KEDGE2( 5,JEL(5)) + NV2
 IND_I( 77) = KAREA2( 2,JEL(5)) + NE2
 IND_I( 78) = KEDGE2( 6,JEL(5)) + NV2
 IND_I( 79) = KAREA2( 5,JEL(6)) + NE2
 IND_I( 80) = KEDGE2( 5,JEL(6)) + NV2
 IND_I( 81) = KAREA2( 5,JEL(5)) + NE2
 IND_I( 82) =           JEL(5)  + NA2
 IND_I( 83) = KAREA2( 3,JEL(5)) + NE2
 IND_I( 84) =           JEL(6)  + NA2
 IND_I( 85) = KAREA2( 2,JEL(6)) + NE2
 IND_I( 86) = KEDGE2( 8,JEL(5)) + NV2
 IND_I( 87) = KAREA2( 4,JEL(5)) + NE2
 IND_I( 88) = KEDGE2( 7,JEL(5)) + NV2
 IND_I( 89) = KAREA2( 3,JEL(6)) + NE2
 IND_I( 90) = KEDGE2( 6,JEL(6)) + NV2
 IND_I( 91) = KAREA2( 2,JEL(8)) + NE2
 IND_I( 92) =           JEL(8)  + NA2
 IND_I( 93) = KAREA2( 4,JEL(8)) + NE2
 IND_I( 94) =           JEL(7)  + NA2
 IND_I( 95) = KAREA2( 5,JEL(7)) + NE2
 IND_I( 96) = KEDGE2( 5,JEL(8)) + NV2
 IND_I( 97) = KAREA2( 5,JEL(8)) + NE2
 IND_I( 98) = KEDGE2( 8,JEL(8)) + NV2
 IND_I( 99) = KAREA2( 2,JEL(7)) + NE2
 IND_I(100) = KEDGE2( 5,JEL(7)) + NV2

 IND_I(101) = KVERT2( 1,JEL(5)) + N00
 IND_I(102) = KEDGE2( 1,JEL(5)) + NV2
 IND_I(103) = KVERT2( 2,JEL(5)) + N00
 IND_I(104) = KEDGE2( 4,JEL(6)) + NV2
 IND_I(105) = KVERT2( 1,JEL(6)) + N00
 IND_I(106) = KEDGE2( 4,JEL(5)) + NV2
 IND_I(107) = KAREA2( 1,JEL(5)) + NE2
 IND_I(108) = KEDGE2( 2,JEL(5)) + NV2
 IND_I(109) = KAREA2( 1,JEL(6)) + NE2
 IND_I(110) = KEDGE2( 1,JEL(6)) + NV2
 IND_I(111) = KVERT2( 4,JEL(5)) + N00
 IND_I(112) = KEDGE2( 3,JEL(5)) + NV2
 IND_I(113) = KVERT2( 3,JEL(5)) + N00
 IND_I(114) = KEDGE2( 2,JEL(6)) + NV2
 IND_I(115) = KVERT2( 2,JEL(6)) + N00
 IND_I(116) = KEDGE2( 1,JEL(8)) + NV2
 IND_I(117) = KAREA2( 1,JEL(8)) + NE2
 IND_I(118) = KEDGE2( 3,JEL(8)) + NV2
 IND_I(119) = KAREA2( 1,JEL(7)) + NE2
 IND_I(120) = KEDGE2( 4,JEL(7)) + NV2
 IND_I(121) = KVERT2( 1,JEL(8)) + N00
 IND_I(122) = KEDGE2( 4,JEL(8)) + NV2
 IND_I(123) = KVERT2( 4,JEL(8)) + N00
 IND_I(124) = KEDGE2( 1,JEL(7)) + NV2
 IND_I(125) = KVERT2( 1,JEL(7)) + N00


 IND_J( 1) = KVERT1( 1,IEL) + N00
 IND_J( 2) = KVERT1( 2,IEL) + N00
 IND_J( 3) = KVERT1( 3,IEL) + N00
 IND_J( 4) = KVERT1( 4,IEL) + N00
 IND_J( 5) = KVERT1( 5,IEL) + N00
 IND_J( 6) = KVERT1( 6,IEL) + N00
 IND_J( 7) = KVERT1( 7,IEL) + N00
 IND_J( 8) = KVERT1( 8,IEL) + N00
 IND_J( 9) = KEDGE1( 1,IEL) + NV1
 IND_J(10) = KEDGE1( 2,IEL) + NV1
 IND_J(11) = KEDGE1( 3,IEL) + NV1
 IND_J(12) = KEDGE1( 4,IEL) + NV1
 IND_J(13) = KEDGE1( 5,IEL) + NV1
 IND_J(14) = KEDGE1( 6,IEL) + NV1
 IND_J(15) = KEDGE1( 7,IEL) + NV1
 IND_J(16) = KEDGE1( 8,IEL) + NV1
 IND_J(17) = KEDGE1( 9,IEL) + NV1
 IND_J(18) = KEDGE1(10,IEL) + NV1
 IND_J(19) = KEDGE1(11,IEL) + NV1
 IND_J(20) = KEDGE1(12,IEL) + NV1
 IND_J(21) = KAREA1( 1,IEL) + NE1
 IND_J(22) = KAREA1( 2,IEL) + NE1
 IND_J(23) = KAREA1( 3,IEL) + NE1
 IND_J(24) = KAREA1( 4,IEL) + NE1
 IND_J(25) = KAREA1( 5,IEL) + NE1
 IND_J(26) = KAREA1( 6,IEL) + NE1
 IND_J(27) =           IEL  + NA1

 DO I=1,125
   II = IND_I(I)
   IF (PLDA(II+1).EQ.0) THEN
    KK = 0
    DO J=1,27
     JJ = IND_J(J)
     IF (ABS(A(I,J)).GT.1d-5) THEN
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
DO I=2,NN1+1
 RLDA(I) = RLDA(I) + RLDA(I-1)
END DO
DO I=2,NN2+1
 PLDA(I) = PLDA(I) + PLDA(I-1)
END DO

ALLOCATE (RLDB(NN1+1),PLDB(NN2+1))
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

 IND_I(  1) = KVERT2( 1,JEL(1)) + N00
 IND_I(  2) = KEDGE2( 1,JEL(1)) + NV2
 IND_I(  3) = KVERT2( 2,JEL(1)) + N00
 IND_I(  4) = KEDGE2( 4,JEL(2)) + NV2
 IND_I(  5) = KVERT2( 1,JEL(2)) + N00
 IND_I(  6) = KEDGE2( 4,JEL(1)) + NV2
 IND_I(  7) = KAREA2( 1,JEL(1)) + NE2
 IND_I(  8) = KEDGE2( 2,JEL(1)) + NV2
 IND_I(  9) = KAREA2( 1,JEL(2)) + NE2
 IND_I( 10) = KEDGE2( 1,JEL(2)) + NV2
 IND_I( 11) = KVERT2( 4,JEL(1)) + N00
 IND_I( 12) = KEDGE2( 3,JEL(1)) + NV2
 IND_I( 13) = KVERT2( 3,JEL(1)) + N00
 IND_I( 14) = KEDGE2( 2,JEL(2)) + NV2
 IND_I( 15) = KVERT2( 2,JEL(2)) + N00
 IND_I( 16) = KEDGE2( 1,JEL(4)) + NV2
 IND_I( 17) = KAREA2( 1,JEL(4)) + NE2
 IND_I( 18) = KEDGE2( 3,JEL(4)) + NV2
 IND_I( 19) = KAREA2( 1,JEL(3)) + NE2
 IND_I( 20) = KEDGE2( 4,JEL(3)) + NV2
 IND_I( 21) = KVERT2( 1,JEL(4)) + N00
 IND_I( 22) = KEDGE2( 4,JEL(4)) + NV2
 IND_I( 23) = KVERT2( 4,JEL(4)) + N00
 IND_I( 24) = KEDGE2( 1,JEL(3)) + NV2
 IND_I( 25) = KVERT2( 1,JEL(3)) + N00

 IND_I( 26) = KEDGE2( 5,JEL(1)) + NV2
 IND_I( 27) = KAREA2( 2,JEL(1)) + NE2
 IND_I( 28) = KEDGE2( 6,JEL(1)) + NV2
 IND_I( 29) = KAREA2( 5,JEL(2)) + NE2
 IND_I( 30) = KEDGE2( 5,JEL(2)) + NV2
 IND_I( 31) = KAREA2( 5,JEL(1)) + NE2
 IND_I( 32) =           JEL(1)  + NA2
 IND_I( 33) = KAREA2( 3,JEL(1)) + NE2
 IND_I( 34) =           JEL(2)  + NA2
 IND_I( 35) = KAREA2( 2,JEL(2)) + NE2
 IND_I( 36) = KEDGE2( 8,JEL(1)) + NV2
 IND_I( 37) = KAREA2( 4,JEL(1)) + NE2
 IND_I( 38) = KEDGE2( 7,JEL(1)) + NV2
 IND_I( 39) = KAREA2( 3,JEL(2)) + NE2
 IND_I( 40) = KEDGE2( 6,JEL(2)) + NV2
 IND_I( 41) = KAREA2( 2,JEL(4)) + NE2
 IND_I( 42) =           JEL(4)  + NA2
 IND_I( 43) = KAREA2( 4,JEL(4)) + NE2
 IND_I( 44) =           JEL(3)  + NA2
 IND_I( 45) = KAREA2( 5,JEL(3)) + NE2
 IND_I( 46) = KEDGE2( 5,JEL(4)) + NV2
 IND_I( 47) = KAREA2( 5,JEL(4)) + NE2
 IND_I( 48) = KEDGE2( 8,JEL(4)) + NV2
 IND_I( 49) = KAREA2( 2,JEL(3)) + NE2
 IND_I( 50) = KEDGE2( 5,JEL(3)) + NV2

 IND_I( 51) = KVERT2( 5,JEL(1)) + N00
 IND_I( 52) = KEDGE2( 9,JEL(1)) + NV2
 IND_I( 53) = KVERT2( 6,JEL(1)) + N00
 IND_I( 54) = KEDGE2(12,JEL(2)) + NV2
 IND_I( 55) = KVERT2( 5,JEL(2)) + N00
 IND_I( 56) = KEDGE2(12,JEL(1)) + NV2
 IND_I( 57) = KAREA2( 6,JEL(1)) + NE2
 IND_I( 58) = KEDGE2(10,JEL(1)) + NV2
 IND_I( 59) = KAREA2( 6,JEL(2)) + NE2
 IND_I( 60) = KEDGE2( 9,JEL(2)) + NV2
 IND_I( 61) = KVERT2( 8,JEL(1)) + N00
 IND_I( 62) = KEDGE2(11,JEL(1)) + NV2
 IND_I( 63) = KVERT2( 7,JEL(1)) + N00
 IND_I( 64) = KEDGE2(10,JEL(2)) + NV2
 IND_I( 65) = KVERT2( 6,JEL(2)) + N00
 IND_I( 66) = KEDGE2( 9,JEL(4)) + NV2
 IND_I( 67) = KAREA2( 6,JEL(4)) + NE2
 IND_I( 68) = KEDGE2(11,JEL(4)) + NV2
 IND_I( 69) = KAREA2( 6,JEL(3)) + NE2
 IND_I( 70) = KEDGE2(12,JEL(3)) + NV2
 IND_I( 71) = KVERT2( 5,JEL(4)) + N00
 IND_I( 72) = KEDGE2(12,JEL(4)) + NV2
 IND_I( 73) = KVERT2( 8,JEL(4)) + N00
 IND_I( 74) = KEDGE2( 9,JEL(3)) + NV2
 IND_I( 75) = KVERT2( 5,JEL(3)) + N00

 IND_I( 76) = KEDGE2( 5,JEL(5)) + NV2
 IND_I( 77) = KAREA2( 2,JEL(5)) + NE2
 IND_I( 78) = KEDGE2( 6,JEL(5)) + NV2
 IND_I( 79) = KAREA2( 5,JEL(6)) + NE2
 IND_I( 80) = KEDGE2( 5,JEL(6)) + NV2
 IND_I( 81) = KAREA2( 5,JEL(5)) + NE2
 IND_I( 82) =           JEL(5)  + NA2
 IND_I( 83) = KAREA2( 3,JEL(5)) + NE2
 IND_I( 84) =           JEL(6)  + NA2
 IND_I( 85) = KAREA2( 2,JEL(6)) + NE2
 IND_I( 86) = KEDGE2( 8,JEL(5)) + NV2
 IND_I( 87) = KAREA2( 4,JEL(5)) + NE2
 IND_I( 88) = KEDGE2( 7,JEL(5)) + NV2
 IND_I( 89) = KAREA2( 3,JEL(6)) + NE2
 IND_I( 90) = KEDGE2( 6,JEL(6)) + NV2
 IND_I( 91) = KAREA2( 2,JEL(8)) + NE2
 IND_I( 92) =           JEL(8)  + NA2
 IND_I( 93) = KAREA2( 4,JEL(8)) + NE2
 IND_I( 94) =           JEL(7)  + NA2
 IND_I( 95) = KAREA2( 5,JEL(7)) + NE2
 IND_I( 96) = KEDGE2( 5,JEL(8)) + NV2
 IND_I( 97) = KAREA2( 5,JEL(8)) + NE2
 IND_I( 98) = KEDGE2( 8,JEL(8)) + NV2
 IND_I( 99) = KAREA2( 2,JEL(7)) + NE2
 IND_I(100) = KEDGE2( 5,JEL(7)) + NV2

 IND_I(101) = KVERT2( 1,JEL(5)) + N00
 IND_I(102) = KEDGE2( 1,JEL(5)) + NV2
 IND_I(103) = KVERT2( 2,JEL(5)) + N00
 IND_I(104) = KEDGE2( 4,JEL(6)) + NV2
 IND_I(105) = KVERT2( 1,JEL(6)) + N00
 IND_I(106) = KEDGE2( 4,JEL(5)) + NV2
 IND_I(107) = KAREA2( 1,JEL(5)) + NE2
 IND_I(108) = KEDGE2( 2,JEL(5)) + NV2
 IND_I(109) = KAREA2( 1,JEL(6)) + NE2
 IND_I(110) = KEDGE2( 1,JEL(6)) + NV2
 IND_I(111) = KVERT2( 4,JEL(5)) + N00
 IND_I(112) = KEDGE2( 3,JEL(5)) + NV2
 IND_I(113) = KVERT2( 3,JEL(5)) + N00
 IND_I(114) = KEDGE2( 2,JEL(6)) + NV2
 IND_I(115) = KVERT2( 2,JEL(6)) + N00
 IND_I(116) = KEDGE2( 1,JEL(8)) + NV2
 IND_I(117) = KAREA2( 1,JEL(8)) + NE2
 IND_I(118) = KEDGE2( 3,JEL(8)) + NV2
 IND_I(119) = KAREA2( 1,JEL(7)) + NE2
 IND_I(120) = KEDGE2( 4,JEL(7)) + NV2
 IND_I(121) = KVERT2( 1,JEL(8)) + N00
 IND_I(122) = KEDGE2( 4,JEL(8)) + NV2
 IND_I(123) = KVERT2( 4,JEL(8)) + N00
 IND_I(124) = KEDGE2( 1,JEL(7)) + NV2
 IND_I(125) = KVERT2( 1,JEL(7)) + N00


 IND_J( 1) = KVERT1( 1,IEL) + N00
 IND_J( 2) = KVERT1( 2,IEL) + N00
 IND_J( 3) = KVERT1( 3,IEL) + N00
 IND_J( 4) = KVERT1( 4,IEL) + N00
 IND_J( 5) = KVERT1( 5,IEL) + N00
 IND_J( 6) = KVERT1( 6,IEL) + N00
 IND_J( 7) = KVERT1( 7,IEL) + N00
 IND_J( 8) = KVERT1( 8,IEL) + N00
 IND_J( 9) = KEDGE1( 1,IEL) + NV1
 IND_J(10) = KEDGE1( 2,IEL) + NV1
 IND_J(11) = KEDGE1( 3,IEL) + NV1
 IND_J(12) = KEDGE1( 4,IEL) + NV1
 IND_J(13) = KEDGE1( 5,IEL) + NV1
 IND_J(14) = KEDGE1( 6,IEL) + NV1
 IND_J(15) = KEDGE1( 7,IEL) + NV1
 IND_J(16) = KEDGE1( 8,IEL) + NV1
 IND_J(17) = KEDGE1( 9,IEL) + NV1
 IND_J(18) = KEDGE1(10,IEL) + NV1
 IND_J(19) = KEDGE1(11,IEL) + NV1
 IND_J(20) = KEDGE1(12,IEL) + NV1
 IND_J(21) = KAREA1( 1,IEL) + NE1
 IND_J(22) = KAREA1( 2,IEL) + NE1
 IND_J(23) = KAREA1( 3,IEL) + NE1
 IND_J(24) = KAREA1( 4,IEL) + NE1
 IND_J(25) = KAREA1( 5,IEL) + NE1
 IND_J(26) = KAREA1( 6,IEL) + NE1
 IND_J(27) =           IEL  + NA1

 DO I=1,125
   II = IND_I(I)
   IF (PLDB(II+1).EQ.0) THEN
    DAUX = 1d0
    KK = 0
    DO J=1,27
     JJ = IND_J(J)
     IF (ABS(A(I,J)).GT.1d-5) THEN
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

!  WRITE(outfile(8:9),'(A1,I1)') "P",myid
!  OPEN(FILE=outfile,UNIT=987)
!   DO I=1,nn2
!     WRITE(987,'(2I8,A3,500D12.4)') I,PLDA(I)," | ",(PA(J),J=PLDA(I),PLDA(I+1)-1)
!   END DO
!  CLOSE(987)
! 
!  WRITE(outfile(8:9),'(A1,I1)') "R",myid
!  OPEN(FILE=outfile,UNIT=987)
!   DO I=1,nn1
!     WRITE(987,'(2I8,A3,500D12.4)') I,RLDA(I)," | ",(RA(J),J=RLDA(I),RLDA(I+1)-1)
!   END DO
!  CLOSE(987)

END SUBROUTINE InitE013ProlMat
!
! ----------------------------------------------
!
SUBROUTINE InitE013ElemProlMat(A)
REAL*8 A(125,27)
REAL*8 F(27),PI,PJ,PK
INTEGER I,J,K,N,IEL,II

! OPEN(789,FILE=outfile)

A(:,:) = 0d0

II = 0
DO I=-2,2
 PI = DBLE(I)/2d0
 DO J=-2,2
  PJ = DBLE(J)/2d0
  DO K=-2,2
   PK = DBLE(K)/2d0
   CALL E013_Lin(PK,PJ,PI,F)
   II = II + 1
   A(II,:) = F(:)
!   WRITE(789,'(I3,27D12.4)') II,F(:)
  END DO
 END DO
END DO

! CLOSE(789)
! stop

END SUBROUTINE InitE013ElemProlMat
!
! ----------------------------------------------
!
SUBROUTINE E013_Lin(X1,X2,X3,D)
REAL*8 D(27),X1,X2,X3
REAL*8 :: Q8=0.125D0

D( 1)=-Q8*X1*(1D0-X1)*X2*(1D0-X2)*X3*(1D0-X3)
D( 2)= Q8*X1*(1D0+X1)*X2*(1D0-X2)*X3*(1D0-X3)
D( 3)=-Q8*X1*(1D0+X1)*X2*(1D0+X2)*X3*(1D0-X3)
D( 4)= Q8*X1*(1D0-X1)*X2*(1D0+X2)*X3*(1D0-X3)
D( 5)= Q8*X1*(1D0-X1)*X2*(1D0-X2)*X3*(1D0+X3)
D( 6)=-Q8*X1*(1D0+X1)*X2*(1D0-X2)*X3*(1D0+X3)
D( 7)= Q8*X1*(1D0+X1)*X2*(1D0+X2)*X3*(1D0+X3)
D( 8)=-Q8*X1*(1D0-X1)*X2*(1D0+X2)*X3*(1D0+X3)
D( 9)= Q8*(2D0-2D0*X1**2D0)*X2*(1D0-X2)*X3*(1D0-X3)
D(10)=-Q8*X1*(1D0+X1)*(2D0-2D0*X2**2D0)*X3*(1D0-X3)
D(11)=-Q8*(2D0-2D0*X1**2D0)*X2*(1D0+X2)*X3*(1D0-X3)
D(12)= Q8*X1*(1D0-X1)*(2D0-2D0*X2**2D0)*X3*(1D0-X3)
D(13)= Q8*X1*(1D0-X1)*X2*(1D0-X2)*(2D0-2D0*X3**2D0)
D(14)=-Q8*X1*(1D0+X1)*X2*(1D0-X2)*(2D0-2D0*X3**2D0)
D(15)= Q8*X1*(1D0+X1)*X2*(1D0+X2)*(2D0-2D0*X3**2D0)
D(16)=-Q8*X1*(1D0-X1)*X2*(1D0+X2)*(2D0-2D0*X3**2D0)
D(17)=-Q8*(2D0-2D0*X1**2D0)*X2*(1D0-X2)*X3*(1D0+X3)
D(18)= Q8*X1*(1D0+X1)*(2D0-2D0*X2**2D0)*X3*(1D0+X3)
D(19)= Q8*(2D0-2D0*X1**2D0)*X2*(1D0+X2)*X3*(1D0+X3)
D(20)=-Q8*X1*(1D0-X1)*(2D0-2D0*X2**2D0)*X3*(1D0+X3)
D(21)=-Q8*(2D0-2D0*X1**2D0)*(2D0-2D0*X2**2D0)*X3*(1D0-X3)
D(22)=-Q8*(2D0-2D0*X1**2D0)*X2*(1D0-X2)*(2D0-2D0*X3**2D0)
D(23)= Q8*X1*(1D0+X1)*(2D0-2D0*X2**2D0)*(2D0-2D0*X3**2D0)
D(24)= Q8*(2D0-2D0*X1**2D0)*X2*(1D0+X2)*(2D0-2D0*X3**2D0)
D(25)=-Q8*X1*(1D0-X1)*(2D0-2D0*X2**2D0)*(2D0-2D0*X3**2D0)
D(26)= Q8*(2D0-2D0*X1**2D0)*(2D0-2D0*X2**2D0)*X3*(1D0+X3)
D(27)= Q8*(2D0-2D0*X1**2D0)*(2D0-2D0*X2**2D0)*(2D0-2D0*X3**2D0)

END SUBROUTINE E013_Lin
!
! ----------------------------------------------
!
! SUBROUTINE outputsol(x,dcoor,kvert,NoOfElem,NoOfVert,iInd)
! REAL*8 X(*),dcoor(3,*)
! INTEGER kvert(8,*),NoOfElem,NoOfVert
! INTEGER I,J,iOutUnit,iInd
! CHARACTER*12 cf
! 
! ! RETURN
! iOutUnit = 442
! 
! WRITE(cf,'(A,I1.1,A,I2.2,A)') "gmv_",iInd,'_',myid,".gmv"
! OPEN (UNIT=iOutUnit,FILE=cf,buffered="yes")
! 
! WRITE(iOutUnit,'(A)')'gmvinput ascii'
! WRITE(iOutUnit,*)'nodes ',NoOfVert
! 
! DO i=1,NoOfVert
!  WRITE(iOutUnit,1200) REAL(dcoor(1,i))
! END DO
! DO i=1,NoOfVert
!  WRITE(iOutUnit,1200) REAL(dcoor(2,i))
! END DO
! DO i=1,NoOfVert
!  WRITE(iOutUnit,1200) REAL(dcoor(3,i))
! END DO
! 
! WRITE(iOutUnit,*)'cells ',NoOfElem
! DO i=1,NoOfElem
!  WRITE(iOutUnit,*)'hex 8'
!  WRITE(iOutUnit,1300) (kvert(j,i),j=1,8)
! END DO
! 
! WRITE(iOutUnit,*)  'variable'
! WRITE(iOutUnit,*)  'pressure'
! DO i=1,NoOfElem
!  j = 4*(i-1) + 1
!  WRITE(iOutUnit,1000) REAL(x(j))
! END DO
! 
! WRITE(iOutUnit,*)  'endvars'
! WRITE(iOutUnit,*)  'probtime',timens
! 
! WRITE(iOutUnit,*)  'endgmv'
! 
! CLOSE  (iOutUnit)
! 
! ! pause
! 
! 1000  FORMAT(E12.5)
! 1200  FORMAT(E12.5)
! 1300  FORMAT(8I8)
! 
! END SUBROUTINE outputsol
!
! ----------------------------------------------
!
       SUBROUTINE mgCoarseGridSolver_P()
       INTEGER Iter,i,j,ndof
       REAL*8 daux
       
         IF (myMG%MedLev.EQ.1) CALL E012DISTR_L1(myMG%B(mgLev)%x,KNEL(mgLev))
         IF (myMG%MedLev.EQ.2) CALL E012DISTR_L2(myMG%B(mgLev)%x,KNEL(mgLev))
         IF (myMG%MedLev.EQ.3) CALL E012DISTR_L3(myMG%B(mgLev)%x,KNEL(mgLev))
       
         IF (myid.eq.0) THEN
       
          myMG%X(mgLev)%x = 0d0
       
          IF (myMG%MinLev.EQ.myMG%MedLev) THEN
           CALL myUmfPack_Solve(myMG%X(mgLev)%x,myMG%B(mgLev)%x,UMF_CMat,UMF_lMat,1)
           CoarseIter = 1
       !    CALL outputsol(myMG%X(mgLev)%x,myQ2coor,KWORK(L(KLVERT(mgLev))),KNEL(mgLev),KNVT(mgLev))
       
       !    CALL E012_BiCGStabSolverMaster(myMG%X(mgLev)%x,myMG%B(mgLev)%x,&
       !         ndof,CoarseIter,E012_DAX_Master,E012_DCG_Master,.true.,&
       !         myCG%d1,myCG%d2,myCG%d3,myCG%d4,myCG%d5,myMG%DefImprCoarse)
          ELSE
       
           CALL crs_cycle()
       !    CALL crs_oneStep()
        
         END IF
        END IF
       
        IF (myMG%MedLev.EQ.1) CALL E012GATHR_L1(myMG%X(mgLev)%x,KNEL(mgLev))
        IF (myMG%MedLev.EQ.2) CALL E012GATHR_L2(myMG%X(mgLev)%x,KNEL(mgLev))
        IF (myMG%MedLev.EQ.3) CALL E012GATHR_L3(myMG%X(mgLev)%x,KNEL(mgLev))
       
       !  CALL outputsol(myMG%X(mgLev)%x,myQ2coor,KWORK(L(KLVERT(mgLev))),KNEL(mgLev),KNVT(mgLev))
       
        CALL E013SendK(0,showid,CoarseIter)
       
       END SUBROUTINE mgCoarseGridSolver_P
!
! ----------------------------------------------
!
SUBROUTINE crs_cycle()

IF (MyMG%CycleType.EQ."W")  CALL crs_W_cycle()
IF (MyMG%CycleType.EQ."V")  CALL crs_V_cycle()
IF (MyMG%CycleType.EQ."F")  CALL crs_F_cycle()

END SUBROUTINE crs_cycle
!
! ----------------------------------------------
!
SUBROUTINE crs_W_cycle()
INTEGER imgLev

 CALL crs_down(myMG%MedLev)

 CALL crs_W_subcycle_cc(myMG%MedLev-1)

 CALL crs_up(myMG%MedLev)

END SUBROUTINE crs_W_cycle
!
! ----------------------------------------------
!
RECURSIVE SUBROUTINE crs_W_subcycle_cc(imgLev)
INTEGER imgLev

 IF (imgLev.NE.2) CALL crs_W_subcycle_cc(imgLev-1)

 CALL crs_up(imgLev)
 CALL crs_down(imgLev)

 IF (imgLev.NE.2) CALL crs_W_subcycle_cc(imgLev-1)

END SUBROUTINE crs_W_subcycle_cc
!
! ----------------------------------------------
!
SUBROUTINE crs_F_cycle()
INTEGER imgLev

 CALL crs_down(myMG%MedLev)

 DO imgLev = myMG%MedLev+1,myMG%MedLev-1
  CALL crs_up(imgLev)
  CALL crs_down(imgLev)
 END DO

 CALL crs_up(myMG%MedLev)

END SUBROUTINE crs_F_cycle
!
! ----------------------------------------------
!
SUBROUTINE crs_V_cycle()

 CALL crs_down(myMG%MedLev)

 CALL crs_up(myMG%MedLev)

END SUBROUTINE crs_V_cycle
!
! ----------------------------------------------
!
SUBROUTINE crs_down(nimgLev)
INTEGER imgLev,nimgLev

! Presmoothing + restriction
 DO imgLev = nimgLev,myMG%MinLev+1,-1

!  IF (myid.eq.1) write(*,*) "D",imgLev
  mgLev = imgLev
  CALL crsSmoother()                                          ! takes B as RHS and smoothes X
  CALL crsUpdateDefect(imgLev,.FALSE.)                        ! puts B into D and gets the new defect

  mgLev = imgLev - 1
  CALL crsRestriction()                                       ! brings D down by 1 level and stores it as B

  CALL crsGuessSolution()                                     ! choose zero initial vector
 END DO

 ! Corse Grid Solver
 CALL crsCoarseGridSolver()                                   ! computes linear system with X=(A^-1)B

END SUBROUTINE crs_down
!
! ----------------------------------------------
!
SUBROUTINE crs_up(nimgLev)
INTEGER imgLev,nimgLev

 ! Postsmoothing + prolongation
 DO imgLev = myMG%MinLev+1, nimgLev

!  IF (myid.eq.1) write(*,*) "U",imgLev
  mgLev = imgLev
  CALL crsProlongation()                                      ! brings X up by one level and stores it as AUX

  CALL crsUpdateSolution()                                    ! updates solution X = X_old + AUX(update)

  CALL crsSmoother()                                          ! takes B as RHS And smoothes the solution X further

 END DO

END SUBROUTINE crs_up
!
! ----------------------------------------------
!
SUBROUTINE crsSmoother()
INTEGER Iter,i,j,ndof,neq
REAL*8 daux

 Iter  = myMG%nSmootherSteps
 ndof  = SIZE(myMG%X(mgLev)%x)

 CALL ZTIME(time0)
 CALL crs_E012_SOR(myMG%X(mgLev)%x,myMG%B(mgLev)%x,ndof,Iter)
!  CALL E012_SOR(myMG%X(mgLev)%x,myMG%XP,myMG%B(mgLev)%x,ndof,Iter)
 CALL ZTIME(time1)
 myStat%tSmthP = myStat%tSmthP + (time1-time0)

END SUBROUTINE crsSmoother
!
! ----------------------------------------------
!
SUBROUTINE crsUpdateDefect(imgLev,bDef)
INTEGER i,j,ndof,imgLev,neq
REAL*8  daux
LOGICAL bDef

 mgLev = imgLev
 ndof  = SIZE(myMG%X(mgLev)%x)
 CALL LCP1(myMG%B(mgLev)%x,myMG%D(mgLev)%x,ndof) !myMG%D(mgLev)%x = myMG%B(mgLev)%x
! CALL E012_DAX(myMG%X(mgLev)%x,myMG%XP,myMG%D(mgLev)%x,ndof,-1D0,1D0)
 CALL E012_DAX_Master(myMG%X(mgLev)%x,myMG%D(mgLev)%x,ndof,-1D0,1D0)

 IF (bDef) THEN 
  CALL LCP1(myMG%D(mgLev)%x,myMG%AUX(mgLev)%x,ndof) !myMG%AUX(mgLev)%x = myMG%D(mgLev)%x
  CALL mgDefectNorm_cc(mgLev)
  IF (myid.ne.0) THEN
   myMG%AUX(mgLev)%x = 0d0
  END IF
 END IF

END SUBROUTINE crsUpdateDefect
!
! ----------------------------------------------
!
SUBROUTINE crsRestriction()
INTEGER I,NDOF

 CALL ZTIME(time0)
 CALL E012_Restriction(myMG%D(mgLev+1)%x,myMG%B(mgLev)%x,mg_E012Prol(mgLev)%a,KNEL(mgLev))
 CALL ZTIME(time1)
 myStat%tRestP = myStat%tRestP + (time1-time0)

END SUBROUTINE crsRestriction
!
! ----------------------------------------------
!
SUBROUTINE crsProlongation()
INTEGER I,NDOF

 CALL ZTIME(time0)
 CALL E012_Prolongation(myMG%AUX(mgLev)%x,myMG%X(mgLev-1)%x,mg_E012Prol(mgLev-1)%a,KNEL(mgLev-1))
 CALL ZTIME(time1)
 myStat%tProlP = myStat%tProlP + (time1-time0)

END SUBROUTINE crsProlongation
!
! ----------------------------------------------
!
SUBROUTINE crsGuessSolution()

 IF (mgLev.EQ.myMG%MinLev) RETURN
 myMG%X(mgLev)%x = 0d0

END SUBROUTINE crsGuessSolution
!
! ----------------------------------------------
!
SUBROUTINE crsUpdateSolution()
INTEGER ndof

 ndof = SIZE(myMG%X(mgLev)%x) 
 CALL LLC1(myMG%AUX(mgLev)%x,myMG%X(mgLev)%x,ndof,1D0,1D0)      ! myMG%X(mgLev)%x = myMG%X(mgLev)%x + myMG%AUX(mgLev)%x

END SUBROUTINE crsUpdateSolution
!
! ----------------------------------------------
!
SUBROUTINE crsCoarseGridSolver()
INTEGER ndof,neq
INTEGER ITE,I,j
REAL*8 def,def0

 mgLev = myMG%MinLev
 ILEV = myMG%MinLev
 CoarseIter  = myMG%nIterCoarse
 ndof  = SIZE(myMG%X(mgLev)%x)

 CALL ZTIME(time0)
  CALL myUmfPack_Solve(myMG%X(mgLev)%x,myMG%B(mgLev)%x,UMF_CMat,UMF_lMat,1)
  CoarseIter = 1

!  CALL E012_BiCGStabSolverMaster(myMG%X(mgLev)%x,myMG%B(mgLev)%x,&
!       ndof,CoarseIter,E012_DAX_Master,E012_DCG_Master,.true.,&
!       myCG%d1,myCG%d2,myCG%d3,myCG%d4,myCG%d5,myMG%DefImprCoarse)

 CALL ZTIME(time1)
 myStat%tSolvP = myStat%tSolvP + (time1-time0)

END SUBROUTINE crsCoarseGridSolver
!
! ----------------------------------------------
!

SUBROUTINE crs_oneStep()
INTEGER Iter,i,j,ndof
REAL*8 daux

 Iter  = myMG%nSmootherSteps
 ndof  = SIZE(myMG%X(mgLev)%x)

 CALL ZTIME(time0)
 CALL crs_E012_SOR(myMG%X(mgLev)%x,myMG%B(mgLev)%x,ndof,Iter)
 CALL ZTIME(time1)
 myStat%tSmthP = myStat%tSmthP + (time1-time0)

 CALL LCP1(myMG%B(mgLev)%x,myMG%D(mgLev)%x,ndof)
 CALL E012_DAX_Master(myMG%X(mgLev)%x,myMG%D(mgLev)%x,ndof,-1D0,1D0)

 mgLev = mgLev - 1
 ndof  = SIZE(myMG%X(mgLev)%x)
 
 CALL ZTIME(time0)
 CALL E012_Restriction(myMG%D(mgLev+1)%x,myMG%B(mgLev)%x,mg_E012Prol(mgLev)%a,KNEL(mgLev))
 CALL ZTIME(time1)
 myStat%tRestP = myStat%tRestP + (time1-time0)

 myMG%X(mgLev)%x = 0d0

!  CALL ZTIME(time0)
!  CALL crs_E012_SOR(myMG%X(mgLev)%x,myMG%B(mgLev)%x,ndof,250*Iter)
! ! CALL crss_E012_SOR(myMG%X(mgLev)%x,myMG%B(mgLev)%x,ndof,100*Iter)
!  CALL ZTIME(time1)

!  write(*,*) myMG%B(mgLev)%x
!  write(*,*) " - - -  - -" 
!  write(*,*) myMG%X(mgLev)%x
!  write(*,*) " - - -  - -" 
!  CALL outputsol(myMG%X(mgLev)%x,myQ2coor,KWORK(L(KLVERT(mgLev))),KNEL(mgLev),KNVT(mgLev))
!  pause

 CALL myUmfPack_Solve(myMG%X(mgLev)%x,myMG%B(mgLev)%x,UMF_CMat,UMF_lMat,1)

 mgLev = mgLev + 1
 ndof  = SIZE(myMG%X(mgLev)%x)

 CALL ZTIME(time0)
 CALL E012_Prolongation(myMG%AUX(mgLev)%x,myMG%X(mgLev-1)%x,mg_E012Prol(mgLev-1)%a,KNEL(mgLev-1))
 CALL ZTIME(time1)
 myStat%tProlP = myStat%tProlP + (time1-time0)

 CALL LLC1(myMG%AUX(mgLev)%x,myMG%X(mgLev)%x,ndof,1D0,1D0)

 CALL ZTIME(time0)
 CALL crs_E012_SOR(myMG%X(mgLev)%x,myMG%B(mgLev)%x,ndof,Iter)
 CALL ZTIME(time1)
 myStat%tSmthP = myStat%tSmthP + (time1-time0)

END SUBROUTINE crs_oneStep
!
! ----------------------------------------------
!
SUBROUTINE mgUpdateDefectForVanka()
INTEGER i,j,ndof_u,ndof_p,neq,jVelo,jPres
REAL*8  daux_u(3),daux_v(3),dt,dP,daux_p,DefNorm(4)

IF (myid.ne.0) THEN
 ndof_u  = KNEL(mgLev) + KNVT(mgLev) + KNET(mgLev) + KNAT(mgLev)
 ndof_p  = 4*KNEL(mgLev)
 dt = TSTEP

 DO i=1,ndof_u

  j=MyMG%Lq(mgLev)%LdA(i)
  jVelo = MyMG%Lq(mgLev)%ColA(j)
  daux_v(1) = MyMG%dX_u(mgLev)%x(         jVelo)
  daux_v(2) = MyMG%dX_u(mgLev)%x(  ndof_u+jVelo)
  daux_v(3) = MyMG%dX_u(mgLev)%x(2*ndof_u+jVelo)
  daux_u(1) = (MyMG%A12(mgLev)%a(j)*daux_v(2) + MyMG%A13(mgLev)%a(j)*daux_v(3))
  daux_u(2) = (MyMG%A21(mgLev)%a(j)*daux_v(1) + MyMG%A23(mgLev)%a(j)*daux_v(3))
  daux_u(3) = (MyMG%A31(mgLev)%a(j)*daux_v(1) + MyMG%A32(mgLev)%a(j)*daux_v(2))

  DO j=MyMG%Lq(mgLev)%LdA(i)+1,MyMG%Lq(mgLev)%LdA(i+1)-1
   jVelo = MyMG%Lq(mgLev)%ColA(j)
   daux_v(1) = MyMG%dX_u(mgLev)%x(         jVelo)
   daux_v(2) = MyMG%dX_u(mgLev)%x(  ndof_u+jVelo)
   daux_v(3) = MyMG%dX_u(mgLev)%x(2*ndof_u+jVelo)
   daux_u(1) = daux_u(1) + (MyMG%A11(mgLev)%a(j)*daux_v(1) + MyMG%A12(mgLev)%a(j)*daux_v(2) + MyMG%A13(mgLev)%a(j)*daux_v(3))
   daux_u(2) = daux_u(2) + (MyMG%A21(mgLev)%a(j)*daux_v(1) + MyMG%A22(mgLev)%a(j)*daux_v(2) + MyMG%A23(mgLev)%a(j)*daux_v(3))
   daux_u(3) = daux_u(3) + (MyMG%A31(mgLev)%a(j)*daux_v(1) + MyMG%A32(mgLev)%a(j)*daux_v(2) + MyMG%A33(mgLev)%a(j)*daux_v(3))
  END DO

  DO j=MyMG%Lql(mgLev)%LdA(i),MyMG%Lql(mgLev)%LdA(i+1)-1
   jPres = MyMG%Lql(mgLev)%ColA(j)
   dP    = MyMG%dX_p(mgLev)%x(jPres)
   daux_u(1) = daux_u(1) - dt*MyMG%BX(mgLev)%a(j)*dP
   daux_u(2) = daux_u(2) - dt*MyMG%BY(mgLev)%a(j)*dP
   daux_u(3) = daux_u(3) - dt*MyMG%BZ(mgLev)%a(j)*dP
  END DO

  IF (myMG%KNPRU(i).EQ.1) daux_u(1) = 0d0
  IF (myMG%KNPRV(i).EQ.1) daux_u(2) = 0d0
  IF (myMG%KNPRW(i).EQ.1) daux_u(3) = 0d0
  MyMG%A_u(mgLev)%x(         i) = daux_u(1)
  MyMG%A_u(mgLev)%x(  ndof_u+i) = daux_u(2)
  MyMG%A_u(mgLev)%x(2*ndof_u+i) = daux_u(3)

 END DO

 MyMG%D_u(mgLev)%x = MyMG%B_u(mgLev)%x
 CALL E013UVWSUM(MyMG%D_u(mgLev)%x)
 MyMG%D_u(mgLev)%x = MyMG%D_u(mgLev)%x + MyMG%A_u(mgLev)%x
 CALL E013UVWSUM(MyMG%A_u(mgLev)%x)
 MyMG%D_u(mgLev)%x = MyMG%D_u(mgLev)%x - MyMG%A_u(mgLev)%x

END IF

END SUBROUTINE mgUpdateDefectForVanka
!
! ----------------------------------------------
!
SUBROUTINE VANKA_new(KVERT,KAREA,KEDGE,ndof_u,ndof_p,&
           VALU,VALP,RHSU,RHSP,KNU,KNV,KNW)

REAL*8 VALU(*),VALP(*),RHSU(*),RHSP(*)
INTEGER ndof_u,ndof_p
INTEGER KVERT(8,*),KAREA(6,*),KEDGE(12,*),KNU(*),KNV(*),KNW(*)
INTEGER KDFG1(729),KDFL1(729),KDFG2(256),KDFL2(256)
REAL*8  DEF(2443),SOL(2443)
INTEGER I,J,IL1,IL2,IL3,IG1,IG2,IG3,IL,IG,IEL,IDFL,JDFL
INTEGER jVelo,jPres,nEq,jel
REAL*8 dU,dV,dW,dP,dt,dNorm(2)

dt = TSTEP

IG1 = 0
IG2 = ndof_u
IG3 = 2*ndof_u

CALL E013_PUTMAT(MyMG%A11(mgLev)%a,MyMG%A22(mgLev)%a,MyMG%A33(mgLev)%a,myMg%Lq(ILEV)%LdA,myMg%Lq(ILEV)%nu)

DO IEL = 1,KNEL(1)

 DO i=1,my_mg_CCPiece(ILEV)%nu_qq
  KDFL1(i) = my_mg_CCPiece(ILEV)%E(IEL)%pairV(1,i)
  KDFG1(i) = my_mg_CCPiece(ILEV)%E(IEL)%pairV(2,i)
 END DO
 
 DO i=1,my_mg_CCPiece(ILEV)%nu_lq/4
  KDFL2(4*(i-1)+1) = 4*(i-1)+1
  KDFL2(4*(i-1)+2) = 4*(i-1)+2
  KDFL2(4*(i-1)+3) = 4*(i-1)+3
  KDFL2(4*(i-1)+4) = 4*(i-1)+4

  jel = my_mg_CCPiece(ILEV)%E(IEL)%pairE(2,i)
  KDFG2(4*(i-1)+1) = 4*(jel-1)+1
  KDFG2(4*(i-1)+2) = 4*(jel-1)+2
  KDFG2(4*(i-1)+3) = 4*(jel-1)+3
  KDFG2(4*(i-1)+4) = 4*(jel-1)+4
 END DO

 DO IDFL=1,my_mg_CCPiece(ILEV)%nu_qq

  IL1 = 0*my_mg_CCPiece(ILEV)%nu_qq + KDFL1(IDFL)
  IL2 = 1*my_mg_CCPiece(ILEV)%nu_qq + KDFL1(IDFL)
  IL3 = 2*my_mg_CCPiece(ILEV)%nu_qq + KDFL1(IDFL)
  IG  = KDFG1(IDFL)

  DEF(IL1) =  RHSU(IG1+IG)
  DEF(IL2) =  RHSU(IG2+IG)
  DEF(IL3) =  RHSU(IG3+IG)

  DO JDFL = myMg%Lq(ILEV)%LdA(IG),myMg%Lq(ILEV)%LdA(IG+1)-1
   jVelo = myMg%Lq(ILEV)%ColA(JDFL)
   dU = VALU(IG1+jVelo)
   dV = VALU(IG2+jVelo)
   dW = VALU(IG3+jVelo)
   DEF(IL1) = DEF(IL1) - (MyMG%A11(mgLev)%a(JDFL)*dU + MyMG%A12(mgLev)%a(JDFL)*dV + MyMG%A13(mgLev)%a(JDFL)*dW)
   DEF(IL2) = DEF(IL2) - (MyMG%A21(mgLev)%a(JDFL)*dU + MyMG%A22(mgLev)%a(JDFL)*dV + MyMG%A23(mgLev)%a(JDFL)*dW)
   DEF(IL3) = DEF(IL3) - (MyMG%A31(mgLev)%a(JDFL)*dU + MyMG%A32(mgLev)%a(JDFL)*dV + MyMG%A33(mgLev)%a(JDFL)*dW)
  END DO

  DO JDFL = myMg%Lql(ILEV)%LdA(IG),myMg%Lql(ILEV)%LdA(IG+1)-1
   jPres = myMg%Lql(ILEV)%ColA(JDFL)
   dP    = VALP(jPres)
   DEF(IL1) =  DEF(IL1) + dt*myMg%BX(ILEV)%a(JDFL)*dP
   DEF(IL2) =  DEF(IL2) + dt*myMg%BY(ILEV)%a(JDFL)*dP
   DEF(IL3) =  DEF(IL3) + dt*myMg%BZ(ILEV)%a(JDFL)*dP
  END DO

  IF (KNU(IG).EQ.1) DEF(IL1) = 0d0
  IF (KNV(IG).EQ.1) DEF(IL2) = 0d0
  IF (KNW(IG).EQ.1) DEF(IL3) = 0d0

 END DO

 DO IDFL=1,my_mg_CCPiece(ILEV)%nu_lq
  IL = 3*my_mg_CCPiece(ILEV)%nu_qq + KDFL2(IDFL)
  IG = KDFG2(IDFL)

  DEF(IL) = RHSP(IG)

  DO JDFL = myMg%Llq(ILEV)%LdA(IG),myMg%Llq(ILEV)%LdA(IG+1)-1
   jVelo = myMg%Llq(ILEV)%ColA(JDFL)
   dU = VALU(IG1+jVelo)
   dV = VALU(IG2+jVelo)
   dW = VALU(IG3+jVelo)
   DEF(IL) =  DEF(IL) - myMg%BTX(ILEV)%a(JDFL)*dU - myMg%BTY(ILEV)%a(JDFL)*dV - myMg%BTZ(ILEV)%a(JDFL)*dW
  END DO

 END DO

 nEq = 3*my_mg_CCPiece(ILEV)%nu_qq + my_mg_CCPiece(ILEV)%nu_lq
 SOL=0d0
 CALL myUmfPack_CCSolveLocalMat(SOL,DEF,my_mg_CCPiece(ILEV)%E(iel)%A,my_mg_CCPiece(ILEV)%MPatchCopy,&
                                my_mg_CCPiece(ILEV)%E(iel)%sym,my_mg_CCPiece(ILEV)%E(iel)%num,nEq)
                                     
 DO IDFL=1,my_mg_CCPiece(ILEV)%nu_qq
  IL1 = KDFL1(IDFL) + 0*my_mg_CCPiece(ILEV)%nu_qq
  IL2 = KDFL1(IDFL) + 1*my_mg_CCPiece(ILEV)%nu_qq
  IL3 = KDFL1(IDFL) + 2*my_mg_CCPiece(ILEV)%nu_qq
  IG  = KDFG1(IDFL)
  IF (KNU(IG).NE.1) VALU(IG1+IG) = VALU(IG1+IG) + MyMG%RLX*SOL(IL1)
  IF (KNV(IG).NE.1) VALU(IG2+IG) = VALU(IG2+IG) + MyMG%RLX*SOL(IL2)
  IF (KNW(IG).NE.1) VALU(IG3+IG) = VALU(IG3+IG) + MyMG%RLX*SOL(IL3)
 END DO

 DO IDFL=1,my_mg_CCPiece(ILEV)%nu_lq
  IL = 3*my_mg_CCPiece(ILEV)%nu_qq + KDFL2(IDFL)
  IG = KDFG2(IDFL)
  VALP(IG) = VALP(IG) + SOL(IL)!*MyMG%RLX
 END DO       

END DO

CALL E013_GETMAT(MyMG%A11(mgLev)%a,MyMG%A22(mgLev)%a,MyMG%A33(mgLev)%a,myMg%Lq(ILEV)%LdA,myMg%Lq(ILEV)%nu)

END SUBROUTINE VANKA_new
!
! ----------------------------------------------
!
SUBROUTINE VANKA(KVERT,KAREA,KEDGE,ndof_u,ndof_p,&
           VALU,VALP,RHSU,RHSP,KNU,KNV,KNW)

REAL*8 VALU(*),VALP(*),RHSU(*),RHSP(*)
INTEGER ndof_u,ndof_p
INTEGER KVERT(8,*),KAREA(6,*),KEDGE(12,*),KNU(*),KNV(*),KNW(*)
INTEGER KDFG1(27),KDFL1(27),KDFG2(4),KDFL2(4)
REAL*8  DEF(85),SOL(85)
INTEGER LDA(86),COLA(85,85)
INTEGER I,J,IL1,IL2,IL3,IG1,IG2,IG3,IL,IG,IEL,IDFL,JDFL
INTEGER jVelo,jPres
REAL*8 dU,dV,dW,dP,dt,dNorm(2)

dt = TSTEP

IG1 = 0
IG2 = ndof_u
IG3 = 2*ndof_u

DO i = 1, 86
 LdA(i) =  ((i-1)*85 + 1) - 1
END DO
DO i = 1, 85
 DO j = 1, 85
  ColA (i,j) = (i)-1
 END DO
END DO

!########################################################
DO IEL = 1,KNEL(ILEV)

 CALL NDFGL(IEL,1,13,KVERT,KEDGE,KAREA,KDFG1,KDFL1)
 KDFL2 = [1,2,3,4]
 KDFG2 = [4*iel-3,4*iel-2,4*iel-1,4*iel]

 DO IDFL=1,27
  IL1 = KDFL1(IDFL)
  IL2 = KDFL1(IDFL)+27
  IL3 = KDFL1(IDFL)+54
  IG  = KDFG1(IDFL)

  DEF(IL1) =  RHSU(IG1+IG)
  DEF(IL2) =  RHSU(IG2+IG)
  DEF(IL3) =  RHSU(IG3+IG)

  DO JDFL = myMg%Lq(ILEV)%LdA(IG),myMg%Lq(ILEV)%LdA(IG+1)-1
   jVelo = myMg%Lq(ILEV)%ColA(JDFL)
   dU = VALU(IG1+jVelo)
   dV = VALU(IG2+jVelo)
   dW = VALU(IG3+jVelo)
   IF (jVelo.eq.IG) THEN

    DEF(IL1) = DEF(IL1) - (  MGE013(ILEV)%UE11(IG)*dU + MyMG%A12(mgLev)%a(JDFL)*dV + MyMG%A13(mgLev)%a(JDFL)*dW)

    DEF(IL2) = DEF(IL2) - (MyMG%A21(mgLev)%a(JDFL)*dU +   MGE013(ILEV)%UE22(IG)*dV + MyMG%A23(mgLev)%a(JDFL)*dW)

    DEF(IL3) = DEF(IL3) - (MyMG%A31(mgLev)%a(JDFL)*dU + MyMG%A32(mgLev)%a(JDFL)*dV +   MGE013(ILEV)%UE33(IG)*dW)
   ELSE
    DEF(IL1) = DEF(IL1) - (MyMG%A11(mgLev)%a(JDFL)*dU + MyMG%A12(mgLev)%a(JDFL)*dV + MyMG%A13(mgLev)%a(JDFL)*dW)
    DEF(IL2) = DEF(IL2) - (MyMG%A21(mgLev)%a(JDFL)*dU + MyMG%A22(mgLev)%a(JDFL)*dV + MyMG%A23(mgLev)%a(JDFL)*dW)
    DEF(IL3) = DEF(IL3) - (MyMG%A31(mgLev)%a(JDFL)*dU + MyMG%A32(mgLev)%a(JDFL)*dV + MyMG%A33(mgLev)%a(JDFL)*dW)

!     DEF(IL1) = DEF(IL1) - myMg%A11(ILEV)%a(JDFL)*dU
!     DEF(IL2) = DEF(IL2) - myMg%A22(ILEV)%a(JDFL)*dV
!     DEF(IL3) = DEF(IL3) - myMg%A33(ILEV)%a(JDFL)*dW

   END IF
  END DO

  DO JDFL = myMg%Lql(ILEV)%LdA(IG),myMg%Lql(ILEV)%LdA(IG+1)-1
   jPres = myMg%Lql(ILEV)%ColA(JDFL)
   dP    = VALP(jPres)
   DEF(IL1) =  DEF(IL1) + dt*myMg%BX(ILEV)%a(JDFL)*dP
   DEF(IL2) =  DEF(IL2) + dt*myMg%BY(ILEV)%a(JDFL)*dP
   DEF(IL3) =  DEF(IL3) + dt*myMg%BZ(ILEV)%a(JDFL)*dP
  END DO

  IF (KNU(IG).EQ.1) DEF(IL1) = 0d0
  IF (KNV(IG).EQ.1) DEF(IL2) = 0d0
  IF (KNW(IG).EQ.1) DEF(IL3) = 0d0

 END DO ! idfl 1,27

 DO IDFL=1,4
  IL = 81 + IDFL
  IG = KDFG2(IDFL)

  DEF(IL) = RHSP(IG)

  DO JDFL = myMg%Llq(ILEV)%LdA(IG),myMg%Llq(ILEV)%LdA(IG+1)-1
   jVelo = myMg%Llq(ILEV)%ColA(JDFL)
   dU = VALU(IG1+jVelo)
   dV = VALU(IG2+jVelo)
   dW = VALU(IG3+jVelo)
   DEF(IL) =  DEF(IL) - myMg%BTX(ILEV)%a(JDFL)*dU - myMg%BTY(ILEV)%a(JDFL)*dV - myMg%BTZ(ILEV)%a(JDFL)*dW
  END DO

 END DO ! idfl 1,4
!########################################################

 CALL ll21(def,85,dNorm(1))

 SOL=0d0
 CALL myUmfPack_CCSolve(SOL,DEF,CC_EMat(ILEV)%E(IEL)%a,LdA,ColA,&
      CC_EMat(ILEV)%E(IEL)%H(1),CC_EMat(ILEV)%E(IEL)%H(2),85)

 DO IDFL=1,27
  IL1 = KDFL1(IDFL)
  IL2 = KDFL1(IDFL)+27
  IL3 = KDFL1(IDFL)+54
  IG  = KDFG1(IDFL)
  IF (KNU(IG).NE.1) VALU(IG1+IG) = VALU(IG1+IG) + MyMG%RLX*SOL(IL1)
  IF (KNV(IG).NE.1) VALU(IG2+IG) = VALU(IG2+IG) + MyMG%RLX*SOL(IL2)
  IF (KNW(IG).NE.1) VALU(IG3+IG) = VALU(IG3+IG) + MyMG%RLX*SOL(IL3)
 END DO
 DO IDFL=1,4
  VALP(4*(iel-1)+IDFL) = VALP(4*(iel-1)+IDFL) + SOL(81+IDFL)!*MyMG%RLX
 END DO       
 

END DO

END SUBROUTINE VANKA
!
! ----------------------------------------------
!
SUBROUTINE MasterVanka()
INTEGER iiITER,Iter,i,j,ndof_u,ndof_p,neq,jVelo,jPres
REAL*8  daux_u(3),daux_v(3),dt,dP,daux_p
REAL*8 daux
CHARACTER*20 cFile
REAL*8, ALLOCATABLE :: rhs_cc(:),sol_cc(:)

ILEV = mgLev
CALL SETLEV(2)

ndof_u  = KNEL(mgLev) + KNVT(mgLev) + KNET(mgLev) + KNAT(mgLev)
ndof_p  = 4*KNEL(mgLev)
dt = TSTEP

ALLOCATE(rhs_cc(3*ndof_u+ndof_p),sol_cc(3*ndof_u+ndof_p)) 

neq = 3*ndof_u + ndof_p

IF (myid.ne.0) THEN

 DO i=1,ndof_u

  daux_u(1) = MyMG%B_u(mgLev)%x(         i)
  daux_u(2) = MyMG%B_u(mgLev)%x(  ndof_u+i)
  daux_u(3) = MyMG%B_u(mgLev)%x(2*ndof_u+i)

  DO j=MyMG%Lq(mgLev)%LdA(i),MyMG%Lq(mgLev)%LdA(i+1)-1
   jVelo = MyMG%Lq(mgLev)%ColA(j)
   daux_v(1) = MyMG%dX_u(mgLev)%x(         jVelo)
   daux_v(2) = MyMG%dX_u(mgLev)%x(  ndof_u+jVelo)
   daux_v(3) = MyMG%dX_u(mgLev)%x(2*ndof_u+jVelo)
   daux_u(1) = daux_u(1) - (MyMG%A11(mgLev)%a(j)*daux_v(1) + MyMG%A12(mgLev)%a(j)*daux_v(2) + MyMG%A13(mgLev)%a(j)*daux_v(3))
   daux_u(2) = daux_u(2) - (MyMG%A21(mgLev)%a(j)*daux_v(1) + MyMG%A22(mgLev)%a(j)*daux_v(2) + MyMG%A23(mgLev)%a(j)*daux_v(3))
   daux_u(3) = daux_u(3) - (MyMG%A31(mgLev)%a(j)*daux_v(1) + MyMG%A32(mgLev)%a(j)*daux_v(2) + MyMG%A33(mgLev)%a(j)*daux_v(3))
  END DO

  DO j=MyMG%Lql(mgLev)%LdA(i),MyMG%Lql(mgLev)%LdA(i+1)-1
   jPres = MyMG%Lql(mgLev)%ColA(j)
   dP    = MyMG%dX_p(mgLev)%x(jPres)
   daux_u(1) = daux_u(1) + dt*MyMG%BX(mgLev)%a(j)*dP
   daux_u(2) = daux_u(2) + dt*MyMG%BY(mgLev)%a(j)*dP
   daux_u(3) = daux_u(3) + dt*MyMG%BZ(mgLev)%a(j)*dP
  END DO

  IF (myMG%KNPRU(i).EQ.1) daux_u(1) = 0d0
  IF (myMG%KNPRV(i).EQ.1) daux_u(2) = 0d0
  IF (myMG%KNPRW(i).EQ.1) daux_u(3) = 0d0
  MyMG%D_u(mgLev)%x(         i) = daux_u(1)
  MyMG%D_u(mgLev)%x(  ndof_u+i) = daux_u(2)
  MyMG%D_u(mgLev)%x(2*ndof_u+i) = daux_u(3)

 END DO

 DO i=1,ndof_p
  daux_p = MyMG%B_p(mgLev)%x(i)
  DO j=MyMG%Llq(mgLev)%LdA(i),MyMG%Llq(mgLev)%LdA(i+1)-1
   jVelo = MyMG%Llq(mgLev)%ColA(j)
   daux_v(1) = MyMG%dX_u(mgLev)%x(         jVelo)
   daux_v(2) = MyMG%dX_u(mgLev)%x(  ndof_u+jVelo)
   daux_v(3) = MyMG%dX_u(mgLev)%x(2*ndof_u+jVelo)
   daux_p = daux_p - MyMG%BTX(mgLev)%a(j)*daux_v(1) - MyMG%BTY(mgLev)%a(j)*daux_v(2) - MyMG%BTZ(mgLev)%a(j)*daux_v(3)
  END DO
  MyMG%D_p(mgLev)%x(         i) = daux_p
 END DO

 DO i=1,3*ndof_u
  rhs_cc(i) = myMG%D_u(mgLev)%x(i)
 END DO
 DO i=1,ndof_p
  rhs_cc(3*ndof_u + i) = -myMG%D_p(mgLev)%x(i)
 END DO
 CALL ll21(rhs_cc,neq,daux)
!  WRITE(*,*) myid,'ll21 norm:',daux
END IF

CALL COMM_cc_def(rhs_cc,neq)

!  WRITE(Cfile,'(A,I2.2,A)')'input_',myid,'.txt'
!  OPEN(FILE=cFile,unit=336)
!  IF (myid.eq.0) THEN
!   write(336,*) myCrsMat%nu
!   write(336,*) myCrsMat%na
!   DO i =1,myCrsMat%na
!    write(336,*) myCrsMat%Row(i),myCrsMat%Col(i)
!   END DO
!   do i=1,myCrsMat%nu
!    write(336,*) rhs_cc(i)
!   enddo
!  ELSE
!   write(336,*) myCrsMat%na
!   DO i =1,myCrsMat%na
!    write(336,*) myCrsMat%Row(i),myCrsMat%Col(i),myCrsMat%A(i)
!   END DO
!  END IF
!  CLOSE(336)
!  WRITE(*,*) "myid#",myid,' is fine!'
!  pause
 
 IF (myid.eq.0) myCrsMat%D = rhs_cc
!  CALL MUMPS_solver_Central()
 CALL MUMPS_solver_Distributed()
 IF (myid.eq.0) sol_cc = myCrsMat%D
!  IF (myid.eq.0) CALL myUmfPack_CCSolveMaster(sol_cc,rhs_cc,CC_crs_AMat,CC_crs_lMat%LdA,CC_crs_lMat%ColA,CC_H(1),CC_H(2),CC_crs_lMat%nu)
 

!  WRITE(Cfile,'(A,I2.2,A)')'output_',myid,'.txt'
!  OPEN(FILE=ADJUSTL(TRIM(cFile)),unit=336)
!  do i=1,neq
!   write(336,*) sol_cc(i)
!  enddo
!  CLOSE(336)
!  WRITE(*,*) "myid#",myid,' is fine!'

 !  CALL outputsol(sol_cc(1),sol_cc(3*ndof_u+1),myq2coor,KWORK(L(KLVERT(mgLev))),KNEL(mgLev),KNVT(mgLev),ndof_u,0)
!  pause
 !  OPEN (file='sol.txt',unit=199)
!  WRITE(199,*) ((sol_cc(0*ndof_u+i),sol_cc(1*ndof_u+i),sol_cc(2*ndof_u+i)),i=1,ndof_u)
!  WRITE(199,*) 
!  WRITE(199,*) (sol_cc(3*ndof_u+i),i=1,4*ndof_p)
!  close(199)
! END IF

CALL COMM_cc_sol(sol_cc,neq)

! WRITE(*,*) "myid#",myid,' is very fine!'
! pause

IF (myid.ne.0) THEN
 DO i=1,3*ndof_u
  myMG%dX_u(mgLev)%x(i) = myMG%dX_u(mgLev)%x(i) + sol_cc(i)
 END DO
!  DO i=1,ndof_u
!   IF (myMG%KNPRU(i).NE.1) myMG%dX_u(mgLev)%x(0*ndof_u+i) = myMG%dX_u(mgLev)%x(0*ndof_u+i) + sol_cc(0*ndof_u+i)
!   IF (myMG%KNPRV(i).NE.1) myMG%dX_u(mgLev)%x(1*ndof_u+i) = myMG%dX_u(mgLev)%x(1*ndof_u+i) + sol_cc(1*ndof_u+i)
!   IF (myMG%KNPRW(i).NE.1) myMG%dX_u(mgLev)%x(2*ndof_u+i) = myMG%dX_u(mgLev)%x(2*ndof_u+i) + sol_cc(2*ndof_u+i)
!  END DO
 DO i=1,ndof_p
  myMG%dX_p(mgLev)%x(i) = myMG%dX_p(mgLev)%x(i) + sol_cc(3*ndof_u + i)
 END DO
END IF

DEALLOCATE(rhs_cc,sol_cc)


END SUBROUTINE MasterVanka
!
! ----------------------------------------------
!
SUBROUTINE outputsol(u,x,dcoor,kvert,NoOfElem,NoOfVert,ndof_u,iInd)
REAL*8 X(*),dcoor(3,*),u(*)
INTEGER kvert(8,*),NoOfElem,NoOfVert
INTEGER I,J,iOutUnit,iInd,ndof_u
CHARACTER*12 cf

! RETURN
iOutUnit = 442

WRITE(cf,'(A,I1.1,A,I2.2,A)') "gmv_",iInd,'_',myid,".gmv"
OPEN (UNIT=iOutUnit,FILE=cf,buffered="yes")

WRITE(iOutUnit,'(A)')'gmvinput ascii'
WRITE(iOutUnit,*)'nodes ',NoOfVert

DO i=1,NoOfVert
 WRITE(iOutUnit,1200) REAL(dcoor(1,i))
END DO
DO i=1,NoOfVert
 WRITE(iOutUnit,1200) REAL(dcoor(2,i))
END DO
DO i=1,NoOfVert
 WRITE(iOutUnit,1200) REAL(dcoor(3,i))
END DO

WRITE(iOutUnit,*)'cells ',NoOfElem
DO i=1,NoOfElem
 WRITE(iOutUnit,*)'hex 8'
 WRITE(iOutUnit,1300) (kvert(j,i),j=1,8)
END DO

WRITE(iOutUnit,*)  'velocity 1'
DO i=1,NoOfVert
 WRITE(iOutUnit,1000) REAL(U(0*ndof_u+i))
END DO
DO i=1,NoOfVert
 WRITE(iOutUnit,1000) REAL(u(1*ndof_u+i))
END DO
DO i=1,NoOfVert
 WRITE(iOutUnit,1000) REAL(u(2*ndof_u+i))
END DO

WRITE(iOutUnit,*)  'variable'
WRITE(iOutUnit,*)  'pressure 0'
DO i=1,NoOfElem
 j = 4*(i-1) + 1
 WRITE(iOutUnit,1000) REAL(x(j))
END DO

WRITE(iOutUnit,*)  'endvars'
WRITE(iOutUnit,*)  'probtime',timens

WRITE(iOutUnit,*)  'endgmv'

CLOSE  (iOutUnit)

! pause

1000  FORMAT(E12.5)
1200  FORMAT(E12.5)
1300  FORMAT(8I8)

END SUBROUTINE outputsol

END MODULE mg_cc
