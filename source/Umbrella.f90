SUBROUTINE UmbrellaSmoother(myTime,nSteps)
USE var_QuadScalar
USE QuadScalar, ONLY : QuadSc,LinSc
USE PP3D_MPI, ONLY: myid,coarse
REAL*8 myTime
INTEGEr nSteps,nUsedSteps
CHARACTER*60 :: cFile
REAL*8 , ALLOCATABLE :: SendVect(:,:,:)
REAL*8 , ALLOCATABLE :: a1(:),a2(:),a3(:),a4(:),a5(:),a6(:)

IF (nSteps.EQ.0) RETURN

IF (myid.eq.0) GOTO 1

ndof = KNEL(NLMAX) + KNVT(NLMAX) + KNET(NLMAX) + KNAT(NLMAX)
ALLOCATE(a1(ndof))
ALLOCATE(a2(ndof))
ALLOCATE(a3(ndof))
ALLOCATE(a4(ndof))
ALLOCATE(a5(ndof))
ALLOCATE(a6(ndof))
a1 = 0d0
a2 = 0d0
a3 = 0d0
a4 = 0d0
a5 = 0d0
a6 = 0d0

DO ILEV = NLMIN+1,NLMAX

 CALL SETLEV(2)

 nUsedSteps = max(1,nSteps/(4**(ILEV-(NLMIN+1))))
 CALL EdgeRunner(a1,a2,a3,a4,a5,a6,&
      DWORK(L(LCORVG)),KWORK(L(LVERT)),KWORK(L(LEDGE)),NEL,NVT,NET,nUsedSteps,myTime)
 CALL ProlongateCoordinates(DWORK(L(LCORVG)),KWORK(L(LAREA)),KWORK(L(LVERT)),&
      KWORK(L(LEDGE)),nel,nvt,net,nat)

END DO

DEALLOCATE(a1,a2,a3,a4,a5,a6)

1 CONTINUE


!  ILEV = LinSc%prm%MGprmIn%MedLev
!  CALL SETLEV(2)
!  nLengthV = (2**(ILEV-1)+1)**3
!  nLengthE = KNEL(NLMIN)
!  ALLOCATE(SendVect(3,nLengthV,nLengthE))
!  CALL SendNodeValuesToCoarse(SendVect,DWORK(L(LCORVG)),KWORK(L(LVERT)),nLengthV,nLengthE,NEL,NVT)
!  DEALLOCATE(SendVect)

ILEV = NLMIN
CALL SETLEV(2)

CALL ExchangeNodeValuesOnCoarseLevel(DWORK(L(LCORVG)),KWORK(L(LVERT)),NVT,NEL)


! IF (myid.ne.0) THEN
!  cFile="vertLink_00.txt"
!  WRITE(cFile(10:11),'(I2.2)') myid
!  OPEN (FILE=ADJUSTL(TRIM(cFile)),UNIT=555)
!  WRITE(555,*) coarse%myVertLink(1:NVT)
!  CLOSE(555)
! END IF

END
!
! --------------------------------------------------------------
!
SUBROUTINE SendNodeValuesToCoarse(sv,dcoor,kvert,nV,nE,nel,nvt)
 USE QuadScalar, ONLY : myDump
 USE PP3D_MPI

 INTEGER kvert(8,*)
 INTEGER nV,nE,pnE,nEmin,nvt,pN,pID
 REAL*8  sv(3,nV,nE),dcoor(3,*)
 INTEGER iAUX(nvt)

 IF (myid.NE.0) THEN

  DO iel = 1,nE
   DO ivt=1,nV
    jvt = myDump%Vertices(iel,ivt)
    SV(:,ivt,iel) = dcoor(:,jvt)
   END DO
  END DO

 END IF

 IF (myid.NE.0) THEN
   pN = nE*nV
   CALL SENDI_myMPI(pN,0)
   CALL SENDI_myMPI(nE,0)
   CALL SENDD_myMPI(SV,3*pN,0)

!   OPEN (UNIT=555,FILE='mesha.gmv')
! 
!   WRITE(555,'(A)')'gmvinput ascii'
!   WRITE(555,*)'nodes ',nvt
! 
!   DO i=1,nvt
!    WRITE(555,'(E13.5)') REAL(dcoor(1,i))
!   END DO
!   DO i=1,nvt
!    WRITE(555,'(E13.5)') REAL(dcoor(2,i))
!   END DO
!   DO i=1,nvt
!    WRITE(555,'(E13.5)') REAL(dcoor(3,i))
!   END DO
! 
!   WRITE(555,*)'cells ',nel
!   DO i=1,nel
!    WRITE(555,*)'hex 8'
!    WRITE(555,'(8I8)') (kvert(j,i),j=1,8)
!   END DO
! 
!   WRITE(555,*)  'endgmv'
! 
!   CLOSE  (555)
! 
!   end if
 ELSE

  SV = 0d0
  iAux = 0
  DO ivt=1,nvt
   dcoor(:,ivt) = 0d0
  END DO

   DO pID=1,subnodes
!   pID = 1
   CALL RECVI_myMPI(pN ,pID)
   CALL RECVI_myMPI(pnE ,pID)
   CALL RECVD_myMPI(sv,3*pN,pID)

   DO iel=1,pnE
    jel = coarse%pELEMLINK(pID,iel)
    DO ivt = 1,nV
     jvt = myDump%Vertices(jel,ivt)
     dcoor(:,jvt) = dcoor(:,jvt) + sv(:,ivt,iel)
     iaux(jvt) = iaux(jvt) + 1d0
    END DO
   END DO

   END DO

  DO I=1,nvt
   IF (iAux(i).GT.0) THEN
    dcoor(1,i) = dcoor(1,i)/DBLE(iAux(i))
    dcoor(2,i) = dcoor(2,i)/DBLE(iAux(i))
    dcoor(3,i) = dcoor(3,i)/DBLE(iAux(i))
   END IF
  END DO

!   OPEN (UNIT=555,FILE='mesh.gmv')
!    DO i=1,nvt
!     WRITE(555,'(I8,A,99I8)') i,' | ',coarse%pVERTLINK(:,I)
!    END DO
!   CLOSE(555)

!    OPEN (UNIT=555,FILE='aaa_0.gmv')
!    DO iel=1,pnE
!     DO ivt = 1,nV
!      write(555,'(3D16.8,I)') sv(:,ivt,iel)
!     END DO
!    END DO
!    CLOSE(555)
! 
!   OPEN (UNIT=555,FILE='meshb.gmv')
! 
!   WRITE(555,'(A)')'gmvinput ascii'
!   WRITE(555,*)'nodes ',nvt
! 
!   DO i=1,nvt
!    WRITE(555,'(E13.5)') REAL(dcoor(1,i))
!   END DO
!   DO i=1,nvt
!    WRITE(555,'(E13.5)') REAL(dcoor(2,i))
!   END DO
!   DO i=1,nvt
!    WRITE(555,'(E13.5)') REAL(dcoor(3,i))
!   END DO
! 
!   WRITE(555,*)'cells ',nel
!   DO i=1,nel
!    WRITE(555,*)'hex 8'
!    WRITE(555,'(8I8)') (kvert(j,i),j=1,8)
!   END DO
! 
!   WRITE(555,*)  'endgmv'
! 
!   CLOSE  (555)

END IF

!  pause
END
!
! --------------------------------------------------------------
!
SUBROUTINE EdgeRunner(f,x,y,z,w,v,dcorvg,kvert,kedge,nel,nvt,net,nProjStep,myTime)
USE Parametrization, ONLY: ParametrizeBndr
USE var_QuadScalar, ONLY : myALE,distamce,distance
IMPLICIT NONE

REAL*8 myTime
REAL*8 f(*),x(*),y(*),z(*),w(*),v(*),dcorvg(3,*)
INTEGER kedge(12,*),kvert(8,*),nel,nvt,net,nProjStep
INTEGER NeighE(2,12)
DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
INTEGER i,j,k,ivt1,ivt2,iProjStep,iaux,iel
REAL*8 WeightE,P1(3),P2(3),daux2,daux1,PX,PY,PZ,dScale1,dScale2
REAL*8 :: dOmega = 0.25d0
REAL*8 DIST,dIII,www,mydist,ipc
REAL*4, ALLOCATABLE :: myVol(:)

ipc=0

DO k=nvt+1,nvt+net
 v(k) = 1d0
END DO

CALL E013SUM(v)

ALLOCATE(myVol(nel))

DO iProjStep=1,nProjStep

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
 CALL E013SUM(f)
 CALL E013SUM(w)

 DO i=1,nvt
  f(i) = f(i)/w(i)
 END DO

 DO i=1,nvt
   PX = dcorvg(1,i)
   PY = dcorvg(2,i)
   PZ = dcorvg(3,i)
   !call GetMonitor()
   call getdistanceid(px,py,pz,myALE%Monitor(i),ipc);
   IF (myALE%Monitor(i).GT.0d0) dIII = 75.1d0 * myALE%Monitor(i)
   IF (myALE%Monitor(i).LT.0d0) dIII = 75.1d0 * -2d0*myALE%Monitor(i)
   f(i) = f(i) * MAX((1.0d0-dIII),0.5d0)**6.5d0
   !IF (distamce(i).GT.0d0) dIII = distamce(i)
   !IF (distamce(i).LT.0d0) dIII = -2d0*distamce(i)   
   !f(i) = MAX((0.0333d0-dIII),0.016666d0)**1.25d0   
   myALE%Monitor(i) = f(i)
   !distance=f(i)

!    f(i) = MAX(2.5d0-MIN(ABS(distamce(i)-0.3d0),ABS(SQRT(PX*PX+PY*PY)-22.7)),1d0)
!    IF (distamce(i).lt.-2.5d0) f(i) =  0.333d0
!    f(i) = MIN(MAX(2.5d0,myALE%Monitor(i)),10d0)
!    CALL DistanceEstimation(PX,PY,PZ,f(i),myTime)
 END DO
! return
!  DO i=1,nvt
!   f(i) = 1d0
!  END DO

 x(1:nvt) = 0d0
 y(1:nvt) = 0d0
 z(1:nvt) = 0d0
 w(1:nvt) = 0d0

 k=1
 DO i=1,nel
  DO j=1,12
   IF (k.eq.kedge(j,i)) THEN
    ivt1 = kvert(NeighE(1,j),i)
    ivt2 = kvert(NeighE(2,j),i)
    P1(:) = dcorvg(:,ivt1)
    P2(:) = dcorvg(:,ivt2)

    daux1 = ABS(f(ivt1))
    daux2 = ABS(f(ivt2))
!     DIST = SQRT((P1(1)-P2(1))**2d0 + (P1(2)-P2(2))**2d0 + (P1(3)-P2(3))**2d0)
!     daux1 = dIST**0.25d0
!     daux2 = dIST**0.25d0
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
   END IF
  END DO
 END DO

 CALL E013SUM(x)
 CALL E013SUM(y)
 CALL E013SUM(z)
 CALL E013SUM(w)

 DO i=1,nvt
   PX = x(i)/w(i)
   PY = y(i)/w(i)
   PZ = z(i)/w(i)
   dcorvg(1,i) = MAX(0d0,(1d0-dOmega))*dcorvg(1,i) + dOmega*PX
   dcorvg(2,i) = MAX(0d0,(1d0-dOmega))*dcorvg(2,i) + dOmega*PY
!   IF (dcorvg(3,i).GT.1d0.OR.dcorvg(3,i).LT.-15d0) THEN
   dcorvg(3,i) = MAX(0d0,(1d0-dOmega))*dcorvg(3,i) + dOmega*PZ
!   END IF
 END DO

! TODO: Adjust to mesh data structure
! CALL ParametrizeBndr()

END DO

 CONTAINS

SUBROUTINE DistanceEstimation(x1,y1,z1,f1,t)
! USE Sigma_User, ONLY: mySigma
INTEGER :: myCase = 0
REAL*8 x1,y1,z1,f1,t
REAL*8 a1,a2,xd1,xd2,xl1,xl2
REAL*8 :: Width = 0.010d0,Base = 1d0, Height = 6d0
REAL*8 :: xC1=1.5d0,yC1=0.3d0,zC1=0.2d0,Radius1 = 0.050d0
REAL*8 :: xC2=2.0d0,yC2=0.1d0,zC2=0.2d0,Radius2 = 0.075d0
REAL*8 :: xC3=0.5d0,yC3=0.2d0,xC,yC

IF (myCASE.EQ.0) THEN
 f1 = 1d0
END IF

IF (myCASE.EQ.1) THEN
xd1 = SQRT((x1-xC1)**2d0+(y1-yC1)**2d0+(z1-zC1)**2d0) - Radius1

xl1 = MAX(0d0,0.5d0*(2d0-ABS(xd1)))

xd1=min(Width,ABS(xd1))

xd1 = (Width-xd1)*(xd1+Width)/(Width**2d0)

xd1 = Base + xl1 + (Height-Base)*xd1

! --------------
xd2 = SQRT((x1-xC2)**2d0+(y1-yC2)**2d0+(z1-zC2)**2d0) - Radius2

xl2 = MAX(0d0,0.5d0*(2d0-ABS(xd2)))

xd2=min(Width,ABS(xd2))

xd2 = (Width-xd2)*(xd2+Width)/(Width**2d0)

xd2 = Base +  xl2 + (Height-Base)*xd2

f1 = max(xd1,xd2)
END IF

IF (myCASE.EQ.2) THEN

xd1 = SQRT((x1-xC3)**2d0+(y1-yC3)**2d0) - Radius1

xl1 = MAX(0d0,1.5d0*(2.5d0-ABS(xd1)))

xd1=min(Width,ABS(xd1))

xd1 = (Width-xd1)*(xd1+Width)/(Width**2d0)

xd1 = Base + xl1 + (Height-Base)*xd1

f1 = xd1

END IF

IF (myCASE.EQ.3) THEN

xC = 1.75 + 0.1d0*SIN(t*6.2832)
yC = 0.2  + 0.1d0*COS(t*6.2832)

xd1 = SQRT((x1-xC)**2d0+(y1-yC)**2d0+(z1-zC1)**2d0) - Radius1

xl1 = MAX(0d0,0.5d0*(2d0-ABS(xd1)))

xd1=min(Width,ABS(xd1))

xd1 = (Width-xd1)*(xd1+Width)/(Width**2d0)

xd1 = Base + xl1 + (Height-Base)*xd1

f1 = xd1

END IF

END SUBROUTINE DistanceEstimation

END
!
! --------------------------------------------------------------
!
SUBROUTINE ProlongateCoordinates(dcorvg,dcorvg2,karea,kvert,kedge,nel,nvt,net,nat)
IMPLICIT NONE
REAL*8  dcorvg(3,*)
REAL*8  dcorvg2(3,*)
INTEGER kvert(8,*),kedge(12,*),karea(6,*),nel,nvt,net,nat
REAL*8 PX,PY,PZ
INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4
INTEGER NeighE(2,12),NeighA(4,6)
DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/

k=1
DO i=1,nel
 DO j=1,12
  IF (k.eq.kedge(j,i)) THEN
   ivt1 = kvert(NeighE(1,j),i)
   ivt2 = kvert(NeighE(2,j),i)
   PX = 0.5d0*(dcorvg(1,ivt1)+dcorvg(1,ivt2))
   PY = 0.5d0*(dcorvg(2,ivt1)+dcorvg(2,ivt2))
   PZ = 0.5d0*(dcorvg(3,ivt1)+dcorvg(3,ivt2))
   dcorvg2(:,nvt+k)=[PX,PY,PZ]
   k = k + 1
  END IF
 END DO
END DO
return

k=1
DO i=1,nel
 DO j=1,6
  IF (k.eq.karea(j,i)) THEN
   ivt1 = kvert(NeighA(1,j),i)
   ivt2 = kvert(NeighA(2,j),i)
   ivt3 = kvert(NeighA(3,j),i)
   ivt4 = kvert(NeighA(4,j),i)
   PX = 0.25d0*(dcorvg(1,ivt1)+dcorvg(1,ivt2)+dcorvg(1,ivt3)+dcorvg(1,ivt4))
   PY = 0.25d0*(dcorvg(2,ivt1)+dcorvg(2,ivt2)+dcorvg(2,ivt3)+dcorvg(2,ivt4))
   PZ = 0.25d0*(dcorvg(3,ivt1)+dcorvg(3,ivt2)+dcorvg(3,ivt3)+dcorvg(3,ivt4))
   dcorvg2(:,nvt+net+k)=[PX,PY,PZ]
   k = k + 1
  END IF
 END DO
END DO

DO i=1,nel
 PX = 0d0
 PY = 0d0
 PZ = 0d0
 DO j=1,8
  PX = PX + 0.125d0*(dcorvg(1,kvert(j,i)))
  PY = PY + 0.125d0*(dcorvg(2,kvert(j,i)))
  PZ = PZ + 0.125d0*(dcorvg(3,kvert(j,i)))
 END DO
 dcorvg2(:,nvt+net+nat+i)=[PX,PY,PZ]
END DO

END SUBROUTINE ProlongateCoordinates
!
! ----------------------------------------------
!
SUBROUTINE ExchangeNodeValuesOnCoarseLevel(dcorvg,kvert,nvt,nel)
USE PP3D_MPI

IMPLICIT NONE
REAL*8 dcorvg(3,*)
INTEGER kvert(8,*),nvt,nel
INTEGER pN,i,j,pID
INTEGER iAUX(nvt)
REAL*8 dAUX(3,nvt)

 IF (myid.NE.0) THEN
  pN = NVT
  CALL SENDI_myMPI(pN,0)
  CALL SENDD_myMPI(dcorvg,3*pN,0)
 ELSE

  DO i=1,nvt
   dcorvg(:,i) = 0d0
  END DO
  iAux = 0

  DO pID=1,subnodes
   CALL RECVI_myMPI(pN ,pID)
   CALL RECVD_myMPI(dAux,3*pN,pID)

   DO I=1,pN
    J = coarse%pVERTLINK(pID,I)
    iAux(J)     = iAux(J)     + 1
    dcorvg(:,j) = dcorvg(:,j) + dAux(:,i)
   END DO

  END DO

  DO I=1,nvt
   dcorvg(1,i) = dcorvg(1,i)/DBLE(iAux(i))
   dcorvg(2,i) = dcorvg(2,i)/DBLE(iAux(i))
   dcorvg(3,i) = dcorvg(3,i)/DBLE(iAux(i))
  END DO

!   OPEN (UNIT=555,FILE='mesh.gmv')
!    DO i=1,nvt
!     WRITE(555,'(I8,A,99I8)') i,' | ',coarse%pVERTLINK(:,I)
!    END DO
!   CLOSE(555)

!   OPEN (UNIT=555,FILE='mesh.gmv')
! 
!   WRITE(555,'(A)')'gmvinput ascii'
!   WRITE(555,*)'nodes ',nvt
! 
!   DO i=1,nvt
!    WRITE(555,'(E13.5)') REAL(dcorvg(1,i))
!   END DO
!   DO i=1,nvt
!    WRITE(555,'(E13.5)') REAL(dcorvg(2,i))
!   END DO
!   DO i=1,nvt
!    WRITE(555,'(E13.5)') REAL(dcorvg(3,i))
!   END DO
! 
!   WRITE(555,*)'cells ',nel
!   DO i=1,nel
!    WRITE(555,*)'hex 8'
!    WRITE(555,'(8I8)') (kvert(j,i),j=1,8)
!   END DO
! 
!   WRITE(555,*)  'endgmv'
! 
!   CLOSE  (555)
 
END IF

 CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)


END SUBROUTINE ExchangeNodeValuesOnCoarseLevel
