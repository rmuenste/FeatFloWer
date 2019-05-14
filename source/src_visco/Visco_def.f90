MODULE def_ViscoScalar

USE PP3D_MPI, ONLY:myid,showID,MGE013,COMM_Maximum,COMM_SUMM,COMM_NLComplete,&
    subnodes,SENDI_myMPI,SENDD_myMPI,RECVI_myMPI,RECVD_myMPI,myMPI_Barrier
USE var_QuadScalar
USE mg_QuadScalar, ONLY : MG_Solver,mgProlRestInit,mgProlongation,myMG,mgLev

IMPLICIT NONE

CONTAINS
!
! ----------------------------------------------
!
SUBROUTINE GetRHS_Visco(Scalar)
TYPE(TViscoScalar) Scalar
INTEGER I,J,ndof
REAL*8 daux

EXTERNAL E013

ILEV = NLMAX
CALL SETLEV(2)

 VisMat_11 => mg_VisMat_11(ILEV)%a
 VisMat_22 => mg_VisMat_22(ILEV)%a
 VisMat_33 => mg_VisMat_33(ILEV)%a
 VisMat_12 => mg_VisMat_12(ILEV)%a
 VisMat_13 => mg_VisMat_13(ILEV)%a
 VisMat_23 => mg_VisMat_23(ILEV)%a
 !ConstDMat => mg_ConstDMat(ILEV)%a

MlMat     => mg_MlMat(ILEV)%a
KMat      => mg_KMat(ILEV)%a
qMat      => mg_qMat(ILEV)

! Include the convection and mass matrix into the RHS

ndof = Scalar%ndof

Scalar%def(ILEV)%x((0*ndof + 1):(1*ndof)) = MlMat*Scalar%val11
Scalar%def(ILEV)%x((1*ndof + 1):(2*ndof)) = MlMat*Scalar%val22
Scalar%def(ILEV)%x((2*ndof + 1):(3*ndof)) = MlMat*Scalar%val33
Scalar%def(ILEV)%x((3*ndof + 1):(4*ndof)) = MlMat*Scalar%val12
Scalar%def(ILEV)%x((4*ndof + 1):(5*ndof)) = MlMat*Scalar%val13
Scalar%def(ILEV)%x((5*ndof + 1):(6*ndof)) = MlMat*Scalar%val23

CALL LAX17(KMat,qMat%ColA,qMat%LdA,qMat%nu,&
Scalar%val11,Scalar%def(ILEV)%x(0*ndof + 1),-thstep,1d0)
CALL LAX17(KMat,qMat%ColA,qMat%LdA,qMat%nu,&
Scalar%val22,Scalar%def(ILEV)%x(1*ndof + 1),-thstep,1d0)
CALL LAX17(KMat,qMat%ColA,qMat%LdA,qMat%nu,&
Scalar%val33,Scalar%def(ILEV)%x(2*ndof + 1),-thstep,1d0)
CALL LAX17(KMat,qMat%ColA,qMat%LdA,qMat%nu,&
Scalar%val12,Scalar%def(ILEV)%x(3*ndof + 1),-thstep,1d0)
CALL LAX17(KMat,qMat%ColA,qMat%LdA,qMat%nu,&
Scalar%val13,Scalar%def(ILEV)%x(4*ndof + 1),-thstep,1d0)
CALL LAX17(KMat,qMat%ColA,qMat%LdA,qMat%nu,&
Scalar%val23,Scalar%def(ILEV)%x(5*ndof + 1),-thstep,1d0)

! CALL LAX17(ConstDMat,qMat%ColA,qMat%LdA,qMat%nu,&
! Scalar%val11,Scalar%def(ILEV)%x(0*ndof + 1),-tstep*Properties%ViscoAlphaExp,1d0)
! CALL LAX17(ConstDMat,qMat%ColA,qMat%LdA,qMat%nu,&
! Scalar%val22,Scalar%def(ILEV)%x(1*ndof + 1),-tstep*Properties%ViscoAlphaExp,1d0)
! CALL LAX17(ConstDMat,qMat%ColA,qMat%LdA,qMat%nu,&
! Scalar%val33,Scalar%def(ILEV)%x(2*ndof + 1),-tstep*Properties%ViscoAlphaExp,1d0)
! CALL LAX17(ConstDMat,qMat%ColA,qMat%LdA,qMat%nu,&
! Scalar%val12,Scalar%def(ILEV)%x(3*ndof + 1),-tstep*Properties%ViscoAlphaExp,1d0)
! CALL LAX17(ConstDMat,qMat%ColA,qMat%LdA,qMat%nu,&
! Scalar%val13,Scalar%def(ILEV)%x(4*ndof + 1),-tstep*Properties%ViscoAlphaExp,1d0)
! CALL LAX17(ConstDMat,qMat%ColA,qMat%LdA,qMat%nu,&
! Scalar%val23,Scalar%def(ILEV)%x(5*ndof + 1),-tstep*Properties%ViscoAlphaExp,1d0)

CALL   DivGradStress(Scalar%val11,&
       Scalar%grad11%x,Scalar%grad11%y,Scalar%grad11%z,&
       Scalar%def(ILEV)%x(0*ndof + 1),&
       mg_mesh%level(ILEV)%kvert,&
       mg_mesh%level(ILEV)%karea,&
       mg_mesh%level(ILEV)%kedge,&
       mg_mesh%level(ILEV)%dcorvg,&
       E013,Properties%ViscoAlphaExp,1d0)

CALL   DivGradStress(Scalar%val22,&
       Scalar%grad22%x,Scalar%grad22%y,Scalar%grad22%z,&
       Scalar%def(ILEV)%x(1*ndof + 1),&
       mg_mesh%level(ILEV)%kvert,&
       mg_mesh%level(ILEV)%karea,&
       mg_mesh%level(ILEV)%kedge,&
       mg_mesh%level(ILEV)%dcorvg,&
       E013,Properties%ViscoAlphaExp,1d0)

CALL   DivGradStress(Scalar%val33,&
       Scalar%grad33%x,Scalar%grad33%y,Scalar%grad33%z,&
       Scalar%def(ILEV)%x(2*ndof + 1),&
       mg_mesh%level(ILEV)%kvert,&
       mg_mesh%level(ILEV)%karea,&
       mg_mesh%level(ILEV)%kedge,&
       mg_mesh%level(ILEV)%dcorvg,&
       E013,Properties%ViscoAlphaExp,1d0)

CALL   DivGradStress(Scalar%val12,&
       Scalar%grad12%x,Scalar%grad12%y,Scalar%grad12%z,&
       Scalar%def(ILEV)%x(3*ndof + 1),&
       mg_mesh%level(ILEV)%kvert,&
       mg_mesh%level(ILEV)%karea,&
       mg_mesh%level(ILEV)%kedge,&
       mg_mesh%level(ILEV)%dcorvg,&
       E013,Properties%ViscoAlphaExp,1d0)

CALL   DivGradStress(Scalar%val13,&
       Scalar%grad13%x,Scalar%grad13%y,Scalar%grad13%z,&
       Scalar%def(ILEV)%x(4*ndof + 1),&
       mg_mesh%level(ILEV)%kvert,&
       mg_mesh%level(ILEV)%karea,&
       mg_mesh%level(ILEV)%kedge,&
       mg_mesh%level(ILEV)%dcorvg,&
       E013,Properties%ViscoAlphaExp,1d0)
       
CALL   DivGradStress(Scalar%val23,&
       Scalar%grad23%x,Scalar%grad23%y,Scalar%grad23%z,&
       Scalar%def(ILEV)%x(5*ndof + 1),&
       mg_mesh%level(ILEV)%kvert,&
       mg_mesh%level(ILEV)%karea,&
       mg_mesh%level(ILEV)%kedge,&
       mg_mesh%level(ILEV)%dcorvg,&
       E013,Properties%ViscoAlphaExp,1d0)

Scalar%rhs0 = Scalar%def(ILEV)%x

END SUBROUTINE GetRHS_Visco
!
! ----------------------------------------------
!
SUBROUTINE GetMat_Visco(Scalar)
TYPE(TViscoScalar) Scalar
INTEGER I,J,ndof
REAL*8 daux

DO ILEV = NLMIN,NLMAX

 CALL SETLEV(2)

 VisMat_11 => mg_VisMat_11(ILEV)%a
 VisMat_22 => mg_VisMat_22(ILEV)%a
 VisMat_33 => mg_VisMat_33(ILEV)%a
 VisMat_12 => mg_VisMat_12(ILEV)%a
 VisMat_13 => mg_VisMat_13(ILEV)%a
 VisMat_23 => mg_VisMat_23(ILEV)%a
 !ConstDMat => mg_ConstDMat(ILEV)%a
 hDMat => mg_hDMat(ILEV)%a
 MlMat     => mg_MlMat(ILEV)%a
 KMat      => mg_KMat(ILEV)%a
 qMat      => mg_qMat(ILEV)

 DO I=1,qMat%nu
  J = qMat%LdA(I)
  daux = MlMat(I) + thstep*KMat(J) + tstep*Properties%ViscoAlphaImp*hDMat(J)
  VisMat_11(J) = daux
  VisMat_22(J) = daux
  VisMat_33(J) = daux
  VisMat_12(J) = daux
  VisMat_13(J) = daux
  VisMat_23(J) = daux
  DO J=qMat%LdA(I)+1,qMat%LdA(I+1)-1
   daux = thstep*KMat(J) + tstep*Properties%ViscoAlphaImp*hDMat(J)
   VisMat_11(J) = daux
   VisMat_22(J) = daux
   VisMat_33(J) = daux
   VisMat_12(J) = daux
   VisMat_13(J) = daux
   VisMat_23(J) = daux
  END DO
 END DO
END DO

END SUBROUTINE GetMat_Visco
!
! ----------------------------------------------
!
SUBROUTINE MatDef_Visco(ViscoSc,VeloSc,dLambda)
TYPE(TViscoScalar) ViscoSc
TYPE(TQuadScalar) VeloSc
INTEGER I,J,ndof
REAL*8 daux,dLambda
EXTERNAL E013

IF (myid.eq.0) RETURN

DO ILEV = NLMIN,NLMAX

 CALL SETLEV(2)

 VisMat_11 => mg_VisMat_11(ILEV)%a
 VisMat_22 => mg_VisMat_22(ILEV)%a
 VisMat_33 => mg_VisMat_33(ILEV)%a
 VisMat_12 => mg_VisMat_12(ILEV)%a
 VisMat_13 => mg_VisMat_13(ILEV)%a
 VisMat_23 => mg_VisMat_23(ILEV)%a
 MlMat     => mg_MlMat(ILEV)%a
 KMat      => mg_KMat(ILEV)%a
 qMat      => mg_qMat(ILEV)

!  VisMat_11 = 0d0
!  VisMat_22 = 0d0
!  VisMat_33 = 0d0
!  VisMat_12 = 0d0
!  VisMat_13 = 0d0
!  VisMat_23 = 0d0
!  CALL ViscoDriver_Mat(VeloSc%ValU,VeloSc%ValV,VeloSc%ValW,&
!       ViscoSc%val11,ViscoSc%val22,ViscoSc%val33,ViscoSc%val12,ViscoSc%val13,ViscoSc%val23,&
!       VisMat_11,VisMat_22,VisMat_33,VisMat_12,VisMat_13,VisMat_23,&
!       qMat%na,qMat%ColA,qMat%LdA,&
!       KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),DWORK(L(LCORVG)),&
!       E013,dLambda,tstep)
END DO

!IF (myid.eq.8) write(*,*)  VisMat_11

ILEV = NLMAX
CALL SETLEV(2)

ndof = ViscoSc%ndof
ViscoSc%def(ILEV)%x = ViscoSc%rhs0

CALL LAX17(VisMat_11,qMat%ColA,qMat%LdA,qMat%nu,&
ViscoSc%val11,ViscoSc%def(ILEV)%x(0*ndof + 1),-1d0,1d0)
CALL LAX17(VisMat_22,qMat%ColA,qMat%LdA,qMat%nu,&
ViscoSc%val22,ViscoSc%def(ILEV)%x(1*ndof + 1),-1d0,1d0)
CALL LAX17(VisMat_33,qMat%ColA,qMat%LdA,qMat%nu,&
ViscoSc%val33,ViscoSc%def(ILEV)%x(2*ndof + 1),-1d0,1d0)
CALL LAX17(VisMat_12,qMat%ColA,qMat%LdA,qMat%nu,&
ViscoSc%val12,ViscoSc%def(ILEV)%x(3*ndof + 1),-1d0,1d0)
CALL LAX17(VisMat_13,qMat%ColA,qMat%LdA,qMat%nu,&
ViscoSc%val13,ViscoSc%def(ILEV)%x(4*ndof + 1),-1d0,1d0)
CALL LAX17(VisMat_23,qMat%ColA,qMat%LdA,qMat%nu,&
ViscoSc%val23,ViscoSc%def(ILEV)%x(5*ndof + 1),-1d0,1d0)


CALL ViscoDriver_Def(VeloSc%ValU,VeloSc%ValV,VeloSc%ValW,myQ2coor,&
     ViscoSc%val11,ViscoSc%val22,ViscoSc%val33,ViscoSc%val12,ViscoSc%val13,ViscoSc%val23,&
     ViscoSc%def(ILEV)%x(0*ndof + 1),ViscoSc%def(ILEV)%x(1*ndof + 1),ViscoSc%def(ILEV)%x(2*ndof + 1),&
     ViscoSc%def(ILEV)%x(3*ndof + 1),ViscoSc%def(ILEV)%x(4*ndof + 1),ViscoSc%def(ILEV)%x(5*ndof + 1),&
     mg_mesh%level(ILEV)%kvert,&
     mg_mesh%level(ILEV)%karea,&
     mg_mesh%level(ILEV)%kedge,&
     mg_mesh%level(ILEV)%dcorvg,&
     E013,dLambda,tstep)


END SUBROUTINE MatDef_Visco
!
! ----------------------------------------------
!
SUBROUTINE CopyOldSolution(ViscoSc)
TYPE(TViscoScalar) ViscoSc
INTEGER I,J,ndof

ndof = ViscoSc%ndof
ViscoSc%ValOld((0*ndof + 1):(1*ndof)) = ViscoSc%val11
ViscoSc%ValOld((1*ndof + 1):(2*ndof)) = ViscoSc%val22
ViscoSc%ValOld((2*ndof + 1):(3*ndof)) = ViscoSc%val33
ViscoSc%ValOld((3*ndof + 1):(4*ndof)) = ViscoSc%val12
ViscoSc%ValOld((4*ndof + 1):(5*ndof)) = ViscoSc%val13
ViscoSc%ValOld((5*ndof + 1):(6*ndof)) = ViscoSc%val23

END SUBROUTINE CopyOldSolution
!
! ----------------------------------------------
!
SUBROUTINE UpdateSolution(ViscoSc,Bndry_Val)
TYPE(TViscoScalar) ViscoSc
EXTERNAL Bndry_Val
INTEGER I,J,ndof

ndof = ViscoSc%ndof
ViscoSc%val11 = ViscoSc%val11 + ViscoSc%ValOld((0*ndof + 1):(1*ndof))
ViscoSc%val22 = ViscoSc%val22 + ViscoSc%ValOld((1*ndof + 1):(2*ndof))
ViscoSc%val33 = ViscoSc%val33 + ViscoSc%ValOld((2*ndof + 1):(3*ndof))
ViscoSc%val12 = ViscoSc%val12 + ViscoSc%ValOld((3*ndof + 1):(4*ndof))
ViscoSc%val13 = ViscoSc%val13 + ViscoSc%ValOld((4*ndof + 1):(5*ndof))
ViscoSc%val23 = ViscoSc%val23 + ViscoSc%ValOld((5*ndof + 1):(6*ndof))

CALL Bndry_Val()

END SUBROUTINE UpdateSolution
!
! ----------------------------------------------
!
SUBROUTINE Solve_Visco(ViscoSc,Bndry_Mat,def0,def1)
TYPE(TViscoScalar) ViscoSc
EXTERNAL Bndry_Mat
INTEGER I,J,ndof,nIter
REAL*8 dRelax,def0(6),def1(6)

IF (myid.eq.0) GOTO 1

ILEV = NLMAX
CALL SETLEV(2)

ndof = ViscoSc%ndof

VisMat_11 => mg_VisMat_11(ILEV)%a
VisMat_22 => mg_VisMat_22(ILEV)%a
VisMat_33 => mg_VisMat_33(ILEV)%a
VisMat_12 => mg_VisMat_12(ILEV)%a
VisMat_13 => mg_VisMat_13(ILEV)%a
VisMat_23 => mg_VisMat_23(ILEV)%a
qMat      => mg_qMat(ILEV)

CALL Bndry_Mat(VisMat_11,VisMat_22,VisMat_33,VisMat_12,VisMat_13,VisMat_23,qMat%LdA,qMat%nu)

CALL E013UVWMAT(VisMat_11,VisMat_22,VisMat_33,qMat%LdA,qMat%nu)
ViscoSc%diag((0*ndof + 1):(1*ndof)) = MGE013(ILEV)%UE11
ViscoSc%diag((1*ndof + 1):(2*ndof)) = MGE013(ILEV)%UE22
ViscoSc%diag((2*ndof + 1):(3*ndof)) = MGE013(ILEV)%UE33
CALL E013UVWMAT(VisMat_12,VisMat_13,VisMat_23,qMat%LdA,qMat%nu)
ViscoSc%diag((3*ndof + 1):(4*ndof)) = MGE013(ILEV)%UE11
ViscoSc%diag((4*ndof + 1):(5*ndof)) = MGE013(ILEV)%UE22
ViscoSc%diag((5*ndof + 1):(6*ndof)) = MGE013(ILEV)%UE33

! CALL LL21(ViscoSc%diag,6*ndof,DEF1)
! write(*,'(A,3G16.8)') 'norm of the diagonal entries ', def1
! ViscoSc%aux(ILEV)%x = 1d0/ViscoSc%diag
! CALL LL21(ViscoSc%aux(ILEV)%x,6*ndof,DEF1)
! write(*,'(A,3G16.8)') 'norm of the 1/diagonal entries ', def1

dRelax = 0.7d0
nIter = 4
ViscoSc%sol(ILEV)%x = 0d0
CALL Visco_SSOR_Solver(VisMat_11,VisMat_22,VisMat_33,VisMat_12,VisMat_13,VisMat_23,qMat%ColA,qMat%LdA,&
                       ViscoSc%sol(ILEV)%x,ViscoSc%def(ILEV)%x,ViscoSc%aux(ILEV)%x,ViscoSc%diag,&
                       ParKNPR,ndof,nIter,dRelax,def0,def1)

1 CONTINUE

! CALL LL21(ViscoSc%sol(ILEV)%x,6*ndof,DEF1)
! write(*,'(A,3G16.8)') 'norm of the solution vector ', def1

CALL COMM_Maximum(def0(1))
CALL COMM_Maximum(def0(2))
CALL COMM_Maximum(def0(3))
CALL COMM_Maximum(def0(4))
CALL COMM_Maximum(def0(5))
CALL COMM_Maximum(def0(6))
CALL COMM_Maximum(def1(1))
CALL COMM_Maximum(def1(2))
CALL COMM_Maximum(def1(3))
CALL COMM_Maximum(def1(4))
CALL COMM_Maximum(def1(5))
CALL COMM_Maximum(def1(6))

IF (myid.ne.0) THEN
 ViscoSc%val11 = ViscoSc%sol(ILEV)%x((0*ndof + 1):(1*ndof))
 ViscoSc%val22 = ViscoSc%sol(ILEV)%x((1*ndof + 1):(2*ndof))
 ViscoSc%val33 = ViscoSc%sol(ILEV)%x((2*ndof + 1):(3*ndof))
 ViscoSc%val12 = ViscoSc%sol(ILEV)%x((3*ndof + 1):(4*ndof))
 ViscoSc%val13 = ViscoSc%sol(ILEV)%x((4*ndof + 1):(5*ndof))
 ViscoSc%val23 = ViscoSc%sol(ILEV)%x((5*ndof + 1):(6*ndof))
END IF

! ViscoSc%val13 = 0d0
! ViscoSc%val23 = 0d0
! ViscoSc%val33 = 1d0

END SUBROUTINE Solve_Visco
!
! ----------------------------------------------
!
SUBROUTINE GetDef_Visco(ViscoSc,Bndry_Mat,def)
TYPE(TViscoScalar) ViscoSc
EXTERNAL Bndry_Mat
INTEGER I,J,ndof
REAL*8 def(6)

IF (myid.eq.0) GOTO 1

ILEV = NLMAX
CALL SETLEV(2)

ndof = ViscoSc%ndof

VisMat_11 => mg_VisMat_11(ILEV)%a
VisMat_22 => mg_VisMat_22(ILEV)%a
VisMat_33 => mg_VisMat_33(ILEV)%a
VisMat_12 => mg_VisMat_12(ILEV)%a
VisMat_13 => mg_VisMat_13(ILEV)%a
VisMat_23 => mg_VisMat_23(ILEV)%a
qMat      => mg_qMat(ILEV)

CALL Bndry_Mat(VisMat_11,VisMat_22,VisMat_33,VisMat_12,VisMat_13,VisMat_23,qMat%LdA,qMat%nu)

CALL Visco_CompDefect(VisMat_11,VisMat_22,VisMat_33,VisMat_12,VisMat_13,VisMat_23,qMat%ColA,qMat%LdA,&
                      ViscoSc%sol(ILEV)%x,ViscoSc%def(ILEV)%x,ViscoSc%aux(ILEV)%x,ViscoSc%diag,&
                      ParKNPR,ndof,def)

1 CONTINUE

CALL COMM_Maximum(def(1))
CALL COMM_Maximum(def(2))
CALL COMM_Maximum(def(3))
CALL COMM_Maximum(def(4))
CALL COMM_Maximum(def(5))
CALL COMM_Maximum(def(6))

END SUBROUTINE GetDef_Visco
!
!----------------------------------------------
!
SUBROUTINE Protocol_Visco(iT,def,mfile)
INTEGER I,J,iT,mfile
REAL*8 def(3,6)

IF (myid.eq.1) THEN
 IF (iT.EQ.1) THEN
  write(MTERM,'(104("-"))')
  write(MFILE,'(104("-"))')
  write(MTERM,'(A3,6(A9,("        ")))') "-","S11 ","S22 ","S33 ","S12 ","S13 ","S23 "
  write(MFILE,'(A3,6(A9,("        ")))') "-","S11 ","S22 ","S33 ","S12 ","S13 ","S23 "
  write(MTERM,'(104("-"))')
  write(MFILE,'(104("-"))')
!  write(MTERM,'(A3,("    "),6(ES8.1,A9))') "0",((def(1,i)," "),i=1,6)
  write(MTERM,'(A3,("    "),6(ES8.1,("         ")))') "0",((def(1,i)),i=1,6)
  write(MFILE,'(A3,("    "),6(ES8.1,("         ")))') "0",((def(1,i)),i=1,6)
!  write(MFILE,'(A3,("    "),6(ES8.1,A9))') "0",((def(1,i)," "),i=1,6)
 END IF
! write(MTERM,'(I3,6(ES8.1,ES8.1,A1))') iT,((def(2,i),def(3,i)," "),i=1,6)
! write(MFILE,'(I3,6(ES8.1,ES8.1,A1))') iT,((def(2,i),def(3,i)," "),i=1,6)
 write(MTERM,'(I3,6(ES8.1,ES8.1,(" ")))') iT,((def(2:3,i)),i=1,6)
 write(MFILE,'(I3,6(ES8.1,ES8.1,(" ")))') iT,((def(2:3,i)),i=1,6)
END IF

END SUBROUTINE Protocol_Visco
!
! ----------------------------------------------
!
SUBROUTINE GatherValues(dV,nV)
INTEGER nV,pV
REAL*8 dV(*)
INTEGER iLoc,pID,i

! RETURN
IF (myid.ne.0) THEN

 CALL SENDI_myMPI(nV,0)
 IF (nV.gt.0) THEN
  CALL SENDD_myMPI(dV,2*nV,0)
 END IF

ELSE

 iLoc = 0
 DO pID=1,subnodes
  CALL RECVI_myMPI(pV,pID)
  IF (pV.gt.0) THEN
!    WRITE(*,*) 'receiving',pID
   CALL RECVD_myMPI(dV(iLoc+1),2*pV,pID)
!    WRITE(*,*) 'receiving',pID,"OK"
   iLoc = iLoc + 2*pV
  END IF
 END DO

 CALL SortOutputValues(dV,iLoc/2)

 OPEN(FILE='_data/Result.txt',UNIT=621)
 DO i=1,iLoc/2
  WRITE(621,*) dV(2*(i-1)+1),dV(2*(i-1)+2)
 END DO
 CLOSE(621)


END IF

END SUBROUTINE GatherValues
!
! ----------------------------------------------
!
SUBROUTINE SortOutputValues(DW,N)
INTEGER N
REAL*8 DW(2,*)
INTEGER I,J
REAL*8 DW1,DW2

 DO I=2,N
  DO J=N,I,-1
   IF (DW(1,J).LT.DW(1,J-1)) THEN
    DW1       = DW(1,J)
    DW2       = DW(2,J)
    DW(1,J)   = DW(1,J-1)
    DW(2,J)   = DW(2,J-1)
    DW(1,J-1) = DW1
    DW(2,J-1) = DW2
   END IF
  END DO
 END DO

END SUBROUTINE SortOutputValues
!
! ----------------------------------------------
!
SUBROUTINE ProlongateViscoSolutionSub(qScalar,ViscoSc,Bndry_Val)
TYPE(TViscoScalar) ViscoSc
TYPE(TQuadScalar), INTENT(INOUT), TARGET :: qScalar
INTEGER J,K,ndof
REAL*8 dnorm
EXTERNAL Bndry_Val

IF (myid.ne.0) THEN

 ILEV=NLMAX
 CALL SETLEV(2)

 J = NLMAX
 K = NLMAX-1
 mgLev = J

 IF(myid.eq.showid) WRITE(*,*) "Prolongation of the viscoelastic stress solution to a higher level"

 MyMG%cVariable = "ViscoEla"
 MyMG%KNPRU => qScalar%knprU(J)%x
 MyMG%KNPRV => qScalar%knprV(J)%x
 MyMG%KNPRW => qScalar%knprW(J)%x

 ndof = KNVT(K) + KNAT(K) + KNET(K) + KNEL(K)
 qScalar%sol(K)%x(0*ndof+1:1*ndof) = ViscoSc%Val11(1:ndof)
 qScalar%sol(K)%x(1*ndof+1:2*ndof) = ViscoSc%Val22(1:ndof)
 qScalar%sol(K)%x(2*ndof+1:3*ndof) = ViscoSc%Val33(1:ndof)
 MyMG%X    => qScalar%sol
 MyMG%AUX  => qScalar%aux

 CALL mgProlongation()

 ndof = KNVT(J) + KNAT(J) + KNET(J) + KNEL(J)
 ViscoSc%Val11 = qScalar%aux(J)%x(0*ndof+1:1*ndof)
 ViscoSc%Val22 = qScalar%aux(J)%x(1*ndof+1:2*ndof)
 ViscoSc%Val33 = qScalar%aux(J)%x(2*ndof+1:3*ndof)

 ndof = KNVT(K) + KNAT(K) + KNET(K) + KNEL(K)
 qScalar%sol(K)%x(0*ndof+1:1*ndof) = ViscoSc%Val12(1:ndof)
 qScalar%sol(K)%x(1*ndof+1:2*ndof) = ViscoSc%Val13(1:ndof)
 qScalar%sol(K)%x(2*ndof+1:3*ndof) = ViscoSc%Val23(1:ndof)
 MyMG%X    => qScalar%sol
 MyMG%AUX  => qScalar%aux

 CALL mgProlongation()

 ndof = KNVT(J) + KNAT(J) + KNET(J) + KNEL(J)
 ViscoSc%Val12 = qScalar%aux(J)%x(0*ndof+1:1*ndof)
 ViscoSc%Val13 = qScalar%aux(J)%x(1*ndof+1:2*ndof)
 ViscoSc%Val23 = qScalar%aux(J)%x(2*ndof+1:3*ndof)

 ! Set dirichlet boundary conditions on the solution
 CALL Bndry_Val()

END IF

END SUBROUTINE ProlongateViscoSolutionSub
!
! ----------------------------------------------
!
END MODULE
