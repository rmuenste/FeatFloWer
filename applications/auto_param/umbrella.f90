SUBROUTINE EdgeRunner(f,x,y,z,w,v,dcorvg,kvert,kedge,nel,nvt,net,nProjStep)
USE MeshProcDef, ONLY : myUmbrella_Relax
USE Parametrization, ONLY: ParametrizeBndr

! USE Sigma_User,  ONLY : Shell_dist,mySigma
IMPLICIT NONE

REAL*8 f(*),x(*),y(*),z(*),w(*),v(*),dcorvg(3,*)
INTEGER kedge(2,*),kvert(8,*),nel,nvt,net,nProjStep
INTEGER i,j,k,ivt1,ivt2,iProjStep,iaux,iel
REAL*8 WeightE,P1(3),P2(3),daux2,daux1,PX,PY,PZ,dScale1,dScale2
REAL*8 DIST,dIII
REAL*8 :: dCrit1,dCrit2 
REAL*8 dFactor,dKernel,dPower
REAL*4, ALLOCATABLE :: myVol(:)

!Get the monitorfunction

DO k=nvt+1,nvt+net
 v(k) = 1d0
END DO

ALLOCATE(myVol(nel))

 myVol = 0e0
 CALL  SETARE(myVol,nel,kvert,dcorvg)

 f(1:nvt) = 0d0
 w(1:nvt) = 0d0

 DO iel=1,nel
  DO i=1,8
   j = kvert(i,iel)
   f(j) = f(j) + abs(myVol(iel))
   w(j) = w(j) + 1d0
  END DO
 END DO

 DO i=1,nvt
  f(i) = f(i)/w(i)
 END DO

 DO i=1,nvt
   PX = dcorvg(1,i)
   PY = dcorvg(2,i)
   PZ = dcorvg(3,i)

   dFactor = 1d0
   f(i) = dFactor*f(i)

 END DO

 ! real umbrella
 
 DO iProjStep=1,nProjStep

 x(1:nvt) = 0d0
 y(1:nvt) = 0d0
 z(1:nvt) = 0d0
 w(1:nvt) = 0d0
! 
 k=1
 DO i=1,net

   ivt1 = kedge(1,i)
   ivt2 = kedge(2,i)
   P1(:) = dcorvg(:,ivt1)
   P2(:) = dcorvg(:,ivt2)

    daux1 = ABS(f(ivt1))
    daux2 = ABS(f(ivt2))
    WeightE = 1d0/(v(nvt + k))

    x(ivt1) = x(ivt1) + WeightE*P2(1)*daux2
    y(ivt1) = y(ivt1) + WeightE*P2(2)*daux2
    z(ivt1) = z(ivt1) + WeightE*P2(3)*daux2
    w(ivt1) = w(ivt1) + WeightE*daux2

    x(ivt2) = x(ivt2) + WeightE*P1(1)*daux1
    y(ivt2) = y(ivt2) + WeightE*P1(2)*daux1
    z(ivt2) = z(ivt2) + WeightE*P1(3)*daux1
    w(ivt2) = w(ivt2) + WeightE*daux1

    k = k + 1
 END DO

 DO i=1,nvt
   PX = x(i)/w(i)
   PY = y(i)/w(i)
   PZ = z(i)/w(i)
   dcorvg(1,i) = MAX(0d0,(1d0-myUmbrella_Relax))*dcorvg(1,i) + myUmbrella_Relax*PX
   dcorvg(2,i) = MAX(0d0,(1d0-myUmbrella_Relax))*dcorvg(2,i) + myUmbrella_Relax*PY
   dcorvg(3,i) = MAX(0d0,(1d0-myUmbrella_Relax))*dcorvg(3,i) + myUmbrella_Relax*PZ
 END DO

 CALL ParametrizeBndr()

END DO

DEALLOCATE(myVol)

END SUBROUTINE EdgeRunner

SUBROUTINE SETARE(AVOL,NEL,KVERT,DCORVG)
IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
REAL*4 AVOL
PARAMETER (NNVE=8)
PARAMETER (A1=1D0/6D0)
DIMENSION  AVOL(*),KVERT(NNVE,*),DCORVG(3,*)

DO  IEL=1,NEL

I1=KVERT(1,IEL)
I2=KVERT(2,IEL)
I3=KVERT(3,IEL)
I4=KVERT(4,IEL)
I5=KVERT(5,IEL)
I6=KVERT(6,IEL)
I7=KVERT(7,IEL)
I8=KVERT(8,IEL)

X1=DCORVG(1,I1)
X2=DCORVG(1,I2)
X3=DCORVG(1,I3)
X4=DCORVG(1,I4)
X5=DCORVG(1,I5)
X6=DCORVG(1,I6)
X7=DCORVG(1,I7)
X8=DCORVG(1,I8)

Y1=DCORVG(2,I1)
Y2=DCORVG(2,I2)
Y3=DCORVG(2,I3)
Y4=DCORVG(2,I4)
Y5=DCORVG(2,I5)
Y6=DCORVG(2,I6)
Y7=DCORVG(2,I7)
Y8=DCORVG(2,I8)

Z1=DCORVG(3,I1)
Z2=DCORVG(3,I2)
Z3=DCORVG(3,I3)
Z4=DCORVG(3,I4)
Z5=DCORVG(3,I5)
Z6=DCORVG(3,I6)
Z7=DCORVG(3,I7)
Z8=DCORVG(3,I8)

AAA=A1*((DABS((X4-X1)*(Y4-Y3)*(Z4-Z8)+(Y4-Y1)*  &
       (Z4-Z3)*(X4-X8)+(Z4-Z1)*(X4-X3)*(Y4-Y8)- &
       (X4-X8)*(Y4-Y3)*(Z4-Z1)-(Y4-Y8)*(Z4-Z3)* &
       (X4-X1)-(Z4-Z8)*(X4-X3)*(Y4-Y1)))+       &
       (DABS((X2-X3)*(Y2-Y1)*(Z2-Z6)+(Y2-Y3)*   &
       (Z2-Z1)*(X2-X6)+(Z2-Z3)*(X2-X1)*(Y2-Y6)- &
       (X2-X6)*(Y2-Y1)*(Z2-Z3)-(Y2-Y6)*(Z2-Z1)* &
       (X2-X3)-(Z2-Z6)*(X2-X1)*(Y2-Y3)))+       &
       (DABS((X5-X8)*(Y5-Y6)*(Z5-Z1)+(Y5-Y8)*   &
       (Z5-Z6)*(X5-X1)+(Z5-Z8)*(X5-X6)*(Y5-Y1)- &
       (X5-X1)*(Y5-Y6)*(Z5-Z8)-(Y5-Y1)*(Z5-Z6)* &
       (X5-X8)-(Z5-Z1)*(X5-X6)*(Y5-Y8)))+       &
       (DABS((X7-X6)*(Y7-Y8)*(Z7-Z3)+(Y7-Y6)*   &
       (Z7-Z8)*(X7-X3)+(Z7-Z6)*(X7-X8)*(Y7-Y3)- &
       (X7-X3)*(Y7-Y8)*(Z7-Z6)-(Y7-Y3)*(Z7-Z8)* &
       (X7-X6)-(Z7-Z3)*(X7-X8)*(Y7-Y6)))+       &
       (DABS((X1-X3)*(Y1-Y8)*(Z1-Z6)+(Y1-Y3)*   &
       (Z1-Z8)*(X1-X6)+(Z1-Z3)*(X1-X8)*(Y1-Y6)- &
       (X1-X6)*(Y1-Y8)*(Z1-Z3)-(Y1-Y6)*(Z1-Z8)* &
       (X1-X3)-(Z1-Z6)*(X1-X8)*(Y1-Y3))))
 AVOL(IEL)=REAL(AAA)
END DO

END SUBROUTINE SETARE                                            
