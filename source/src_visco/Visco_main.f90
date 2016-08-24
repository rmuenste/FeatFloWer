MODULE ViscoScalar

USE QuadScalar
USE def_ViscoScalar

IMPLICIT NONE

INTEGER :: ninl=3
CONTAINS
!
! ----------------------------------------------
!
SUBROUTINE Transport_ViscoScalar(mfile)
INTEGER mfile
INTEGER inl
REAL*8 def(3,6)

CALL Create_KMat(QuadSc)
IF (myid.EQ.ShowID) WRITE(MTERM,'(A)', advance='yes') " "

thstep = 0.5d0*tstep
CALL GetRHS_Visco(ViscoSc)

CALL Boundary_ViscoScalar_Def(ViscoSc%def(NLMAX)%x)

CALL GetMat_Visco(ViscoSc)

CALL MatDef_Visco(ViscoSc,QuadSc,Properties%ViscoLambda)

CALL Boundary_ViscoScalar_Def(ViscoSc%def(NLMAX)%x)

! !IF (myid.eq.3) WRITE(*,*) 'the value : ', ViscoSc%def(NLMAX)%x(1*ViscoSc%ndof + 113), ViscoSc%def(NLMAX)%x(1*ViscoSc%ndof + 3553)
! IF (myid.eq.22) WRITE(*,*) 'the value : ', ViscoSc%def(NLMAX)%x(1*ViscoSc%ndof + 1), ViscoSc%def(NLMAX)%x(1*ViscoSc%ndof + 693)
! ! pause

CALL CopyOldSolution(ViscoSc)


DO inl=1,ninl

 CALL Solve_Visco(ViscoSc,Boundary_ViscoScalar_Mat,def(1,:),def(2,:))

 CALL UpdateSolution(ViscoSc,Boundary_ViscoScalar_Val)

 CALL GetMat_Visco(ViscoSc)

 CALL MatDef_Visco(ViscoSc,QuadSc,Properties%ViscoLambda)

 CALL Boundary_ViscoScalar_Def(ViscoSc%def(NLMAX)%x)

 CALL CopyOldSolution(ViscoSc)

 CALL GetDef_Visco(ViscoSc,Boundary_ViscoScalar_Mat,def(3,:))

 CALL Protocol_Visco(inl,def,mfile)

END DO

! ViscoSc%Val33 = 0d0
! ViscoSc%Val13 = 0d0
! ViscoSc%Val23 = 0d0

!CALL EqualizeSolutions()

 CALL Output_BenchQuantity()

END SUBROUTINE Transport_ViscoScalar
!
! ----------------------------------------------
!
SUBROUTINE Init_ViscoScalar_Stuctures(mfile)
INTEGER mfile
INTEGER ndof

ViscoSc%ndof = QuadSc%ndof
ViscoSc%na   = QuadSc%na

ndof = ViscoSc%ndof

ALLOCATE(ViscoSc%Val11(ndof))
ALLOCATE(ViscoSc%Val22(ndof))
ALLOCATE(ViscoSc%Val33(ndof))
ALLOCATE(ViscoSc%Val12(ndof))
ALLOCATE(ViscoSc%Val13(ndof))
ALLOCATE(ViscoSc%Val23(ndof))
ViscoSc%Val11 = 1d0
ViscoSc%Val22 = 1d0
ViscoSc%Val33 = 1d0
ViscoSc%Val12 = 0d0
ViscoSc%Val13 = 0d0
ViscoSc%Val23 = 0d0
ALLOCATE(ViscoSc%rhs0(6*ndof))
ALLOCATE(ViscoSc%ValOld(6*ndof))
ALLOCATE(ViscoSc%diag(6*ndof))
ALLOCATE(ViscoSc%knpr(6*ndof))
ViscoSc%knpr = 0

ALLOCATE(ViscoSc%def(NLMIN:NLMAX))
ALLOCATE(ViscoSc%rhs(NLMIN:NLMAX))
ALLOCATE(ViscoSc%aux(NLMIN:NLMAX))
ALLOCATE(ViscoSc%sol(NLMIN:NLMAX))
DO ILEV=NLMIN,NLMAX
 NDOF = KNVT(ILEV)+KNET(ILEV)+KNAT(ILEV)+KNEL(ILEV)
 ALLOCATE(ViscoSc%def(ILEV)%x(6*NDOF))
 ALLOCATE(ViscoSc%rhs(ILEV)%x(6*NDOF))
 ALLOCATE(ViscoSc%aux(ILEV)%x(6*NDOF))
 ALLOCATE(ViscoSc%sol(ILEV)%x(6*NDOF))
END DO


IF (.not.ALLOCATED(mg_VisMat_11))    ALLOCATE(mg_VisMat_11(NLMIN:NLMAX))
IF (.not.ALLOCATED(mg_VisMat_22))    ALLOCATE(mg_VisMat_22(NLMIN:NLMAX))
IF (.not.ALLOCATED(mg_VisMat_33))    ALLOCATE(mg_VisMat_33(NLMIN:NLMAX))
IF (.not.ALLOCATED(mg_VisMat_12))    ALLOCATE(mg_VisMat_12(NLMIN:NLMAX))
IF (.not.ALLOCATED(mg_VisMat_13))    ALLOCATE(mg_VisMat_13(NLMIN:NLMAX))
IF (.not.ALLOCATED(mg_VisMat_23))    ALLOCATE(mg_VisMat_23(NLMIN:NLMAX))

DO ILEV=NLMIN,NLMAX
 NDOF = KNVT(ILEV)+KNET(ILEV)+KNAT(ILEV)+KNEL(ILEV)

 IF (.not.ALLOCATED(mg_VisMat_11(ILEV)%a)) THEN
  ALLOCATE(mg_VisMat_11(ILEV)%a(qMat%na))
 END IF

 IF (.not.ALLOCATED(mg_VisMat_22(ILEV)%a)) THEN
  ALLOCATE(mg_VisMat_22(ILEV)%a(qMat%na))
 END IF

 IF (.not.ALLOCATED(mg_VisMat_33(ILEV)%a)) THEN
  ALLOCATE(mg_VisMat_33(ILEV)%a(qMat%na))
 END IF

 IF (.not.ALLOCATED(mg_VisMat_12(ILEV)%a)) THEN
  ALLOCATE(mg_VisMat_12(ILEV)%a(qMat%na))
 END IF

 IF (.not.ALLOCATED(mg_VisMat_13(ILEV)%a)) THEN
  ALLOCATE(mg_VisMat_13(ILEV)%a(qMat%na))
 END IF

 IF (.not.ALLOCATED(mg_VisMat_23(ILEV)%a)) THEN
  ALLOCATE(mg_VisMat_23(ILEV)%a(qMat%na))
 END IF

END DO

CALL Create_MMat()

END SUBROUTINE Init_ViscoScalar_Stuctures
!
! ----------------------------------------------
!
SUBROUTINE IniProf_ViscoScalar()
REAL*8 PX,PY,daux
REAL*8 tau(6),psi(6)
INTEGER i

! ViscoSc%Val11 = 0d0
! ViscoSc%Val22 = 0d0
! ViscoSc%Val33 = 0d0
! ViscoSc%Val12 = 0d0
! ViscoSc%Val13 = 0d0
! ViscoSc%Val23 = 0d0
! 
! return
DO i=1,ViscoSc%ndof
 PX = myQ2Coor(1,i)
 PY = myQ2Coor(2,i)
 daux = -Properties%ViscoLambda*0.75*PY
! tau = [1d0 + 2d0*daux*daux,1d0,1d0,daux,0d0,0d0]
 tau = [1d0,1d0,1d0,0d0,0d0,0d0]
 CALL ConvertTauToPsi(tau,psi)

 ViscoSc%Val11(i) = psi(1)
 ViscoSc%Val22(i) = psi(2)
 ViscoSc%Val33(i) = psi(3)
 ViscoSc%Val12(i) = psi(4)
 ViscoSc%Val13(i) = psi(5)
 ViscoSc%Val23(i) = psi(6)
END DO
! ViscoSc%Val11 = 1d0
! ViscoSc%Val22 = 1d0
! ViscoSc%Val33 = 1d0
! ViscoSc%Val12 = 0d0
! ViscoSc%Val13 = 0d0
! ViscoSc%Val23 = 0d0

END SUBROUTINE IniProf_ViscoScalar
!
! ----------------------------------------------
!
SUBROUTINE Boundary_ViscoScalar_Def(def)
REAL*8 def(*)
REAL*8 PX,PZ
INTEGER i

! return
DO i=1,ViscoSc%ndof
 PX = myQ2Coor(1,i)
 PZ = myQ2Coor(3,i)
 IF (PZ.lt.-11.99d0) THEN
  def(0*ViscoSc%ndof+i) = 0d0
  def(1*ViscoSc%ndof+i) = 0d0
  def(2*ViscoSc%ndof+i) = 0d0
  def(3*ViscoSc%ndof+i) = 0d0
  def(4*ViscoSc%ndof+i) = 0d0
  def(5*ViscoSc%ndof+i) = 0d0
 END IF
END DO

END SUBROUTINE Boundary_ViscoScalar_Def
!
! ----------------------------------------------
!
SUBROUTINE Boundary_ViscoScalar_Val
REAL*8 PX,PY,PZ,daux
REAL*8 tau(6),psi(6)
INTEGER i

! return
DO i=1,ViscoSc%ndof
 PX = myQ2Coor(1,i)
 PY = myQ2Coor(2,i)
 PZ = myQ2Coor(3,i)
 IF (PZ.lt.-11.99d0) THEN
  daux = -Properties%ViscoLambda*0.75d0*PY
!   tau = [1d0 + 2d0*daux*daux,1d0,1d0,daux,0d0,0d0]
  tau = [1d0,1d0,1d0,0d0,0d0,0d0]
  CALL ConvertTauToPsi(tau,psi)

  ViscoSc%Val11(i) = psi(1)
  ViscoSc%Val22(i) = psi(2)
  ViscoSc%Val33(i) = psi(3)
  ViscoSc%Val12(i) = psi(4)
  ViscoSc%Val13(i) = psi(5)
  ViscoSc%Val23(i) = psi(6)
!   ViscoSc%Val11(i) = 1d0 + 2d0*daux*daux
!   ViscoSc%Val22(i) = 1d0
!   ViscoSc%Val33(i) = 1d0
!   ViscoSc%Val12(i) = daux
!   ViscoSc%Val13(i) = 0d0
!   ViscoSc%Val23(i) = 0d0
 END IF
END DO

END SUBROUTINE Boundary_ViscoScalar_Val
!
! ----------------------------------------------
!
SUBROUTINE Boundary_ViscoScalar_Mat(DA11,DA22,DA33,DA12,DA13,DA23,KLD,NDOF)
REAL*8  DA11(*),DA22(*),DA33(*),DA12(*),DA13(*),DA23(*)
INTEGER KLD(*),ICOL,I,NDOF
REAL*8  PX,PZ

! return
DO I=1,NDOF
 PX = myQ2Coor(1,i)
 PZ = myQ2Coor(3,i)
 IF (PZ.LT.-11.99d0) THEN
   ICOL = KLD(I)
   DO ICOL=KLD(I)+1,KLD(I+1)-1
    DA11(ICOL) = 0d0
    DA22(ICOL) = 0d0
    DA33(ICOL) = 0d0
    DA12(ICOL) = 0d0
    DA13(ICOL) = 0d0
    DA23(ICOL) = 0d0
   END DO
  END IF
END DO

END SUBROUTINE Boundary_ViscoScalar_Mat
!
! ----------------------------------------------
!
SUBROUTINE Output_BenchQuantity()
REAL*8 PX,PY,PZ,daux(2)
INTEGER i,j,k,l
REAL*8, ALLOCATABLE :: dVals(:,:)
REAL*8 :: PI=3.141592654d0,dAlpha,psi(6),tau(6),dSigma11

j=0
k=0
DO i=1,KNVT(NLMAX+1)!ViscoSc%ndof
 PZ = myQ2Coor(1,i)  
 PY = myQ2Coor(2,i)  
 PX = myQ2Coor(3,i)  
 IF (myBoundary%bWall(i)) THEN
  IF (PY.LT.1.25d0.AND.ABS(PZ).LT.1d-4) THEN
   j = j + 1
  END IF
 END IF
 IF (PY.LT.1d-4.AND.PX.GT.1d0.AND.PX.LT.5d0.AND.ABS(PZ).LT.1d-4) THEN
  k = k + 1
 END IF
END DO

daux(1) = DBLE(j)
daux(2) = DBLE(k)
CALL Comm_SummN(daux,2)
l = INT(daux(1) + daux(2))
ALLOCATE(dVals(2,l))

l=0
DO i=1,KNVT(NLMAX+1)!ViscoSc%ndof
 PZ = myQ2Coor(1,i)  
 PY = myQ2Coor(2,i)  
 PX = myQ2Coor(3,i)  
 IF (myBoundary%bWall(i)) THEN
  IF (PY.LT.1.25d0.AND.ABS(PZ).LT.1d-4) THEN
   l = l + 1
   IF (PX.LT.0d0) THEN
    dAlpha = 0.0d0 - ATAN(PY/PX)/PI
!     WRITE(*,*) "-----", l,PX,dAlpha
   ELSE
    dAlpha = 1.0d0 - ATAN(PY/PX)/PI
!     WRITE(*,*) "+++++", l,PX,dAlpha
   END IF

   psi = [ViscoSc%Val11(i),ViscoSc%Val22(i),ViscoSc%Val33(i),&
          ViscoSc%Val12(i),ViscoSc%Val13(i),ViscoSc%Val23(i)]
   CALL ConvertPsiToTau(psi,tau)
   dSigma11 = (tau(3) - 1d0)*Properties%Viscosity(2)/Properties%ViscoLambda
   
   dVals(1,l) = dAlpha*PI
   dVals(2,l) = dSigma11
  END IF
 END IF
 IF (PY.LT.1d-4.AND.PX.GT.1d0.AND.PX.LT.5d0.AND.ABS(PZ).LT.1d-4) THEN
  l = l + 1

   psi = [ViscoSc%Val11(i),ViscoSc%Val22(i),ViscoSc%Val33(i),&
          ViscoSc%Val12(i),ViscoSc%Val13(i),ViscoSc%Val23(i)]
   CALL ConvertPsiToTau(psi,tau)
   dSigma11 = (tau(3) - 1d0)*Properties%Viscosity(2)/Properties%ViscoLambda

  dVals(1,l) = 1d0*PI + (PX-1d0)
  dVals(2,l) = dSigma11
 END IF
END DO

j = SIZE(dVals(1,:))
! WRITE(*,*) "Number of sample points: ",myid,j,l
CALL GatherValues(dVals,l)

! IF (myid.eq.1) WRITE(*,*) INT(daux)

DEALLOCATE(dVals)

END SUBROUTINE Output_BenchQuantity
!
! ----------------------------------------------
!
SUBROUTINE EqualizeSolutions()
INTEGER ndof,i

 ndof = ViscoSc%ndof

 ViscoSc%rhs0 = 1d0
 CALL E013Sum(ViscoSc%rhs0)

 CALL E013Sum(ViscoSc%val11)
 CALL E013Sum(ViscoSc%val22)
 CALL E013Sum(ViscoSc%val33)
 CALL E013Sum(ViscoSc%val12)
 CALL E013Sum(ViscoSc%val13)
 CALL E013Sum(ViscoSc%val23)

 DO i=1,ndof
  ViscoSc%val11(i)  = ViscoSc%val11(i)/ViscoSc%rhs0(i)
  ViscoSc%val22(i)  = ViscoSc%val22(i)/ViscoSc%rhs0(i)
  ViscoSc%val33(i)  = ViscoSc%val33(i)/ViscoSc%rhs0(i)
  ViscoSc%val12(i)  = ViscoSc%val12(i)/ViscoSc%rhs0(i)
  ViscoSc%val13(i)  = ViscoSc%val13(i)/ViscoSc%rhs0(i)
  ViscoSc%val23(i)  = ViscoSc%val23(i)/ViscoSc%rhs0(i)
 END DO

END SUBROUTINE EqualizeSolutions
!
! ----------------------------------------------
!
SUBROUTINE ProlongateViscoSolution()

 CALL ProlongateViscoSolutionSub(QuadSc,ViscoSc,Boundary_ViscoScalar_Val)

END SUBROUTINE ProlongateViscoSolution
!
! ----------------------------------------------
!
END MODULE ViscoScalar
