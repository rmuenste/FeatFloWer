!-----------------------------------------------------------------------------------------
SUBROUTINE CreateSimpleInterface(myIfc,dTotArea)
USE PP3D_MPI, ONLY:myid
USE var_QuadScalar, ONLY: tTetInterface
IMPLICIT NONE
REAL*8 dTotArea
TYPE(tTetInterface) myIfc
INTEGEr i,j,iP
REAL*8 PM(3),P1(3),P2(3),P3(3),D1(3),D2(3),P(3),DN(3)
REAL*8 dAux,dArea

dTotArea = 0d0

myIfc%L(1) = 1
DO i=2,myIfc%nG
 myIfc%L(i) = myIfc%L(i-1) + myIfc%L(i)
END DO

IF (ALLOCATED(myIfc%Y)) DEALLOCATE(myIfc%Y)
ALLOCATE(myIfc%Y(3,myIfc%nG-1))

DO i=1,myIfc%nG-1
 PM    = 0d0
 dArea = 0d0
 DO j=myIfc%L(i),myIfc%L(i+1)-1
  iP = 3*(j-1)
  P1(:) = myIfc%x(:,iP+1)
  P2(:) = myIfc%x(:,iP+2)
  P3(:) = myIfc%x(:,iP+3)

  D1(:) = P2(:)-P1(:)
  D2(:) = P3(:)-P1(:)
  DN(1) = (D2(2)*D1(3))-(D2(3)*D1(2))
  DN(2) = (D2(3)*D1(1))-(D2(1)*D1(3))
  DN(3) = (D2(1)*D1(2))-(D2(2)*D1(1)) !CrossProduct(b-a, p1-a)
  dAux = 0.5d0*SQRT(DN(1)*DN(1) + DN(2)*DN(2) + DN(3)*DN(3))


  P(:) = (P1(:) + P2(:) + P3(:))/3d0
  dArea = dArea + dAux

  PM(:) = PM(:) + dAux*P(:)
 END DO
 myIfc%Y(:,i) = PM(:)/dArea
 dTotArea = dTotArea + dArea

END DO

END
!-----------------------------------------------------------------------------------------
SUBROUTINE GetClosestGroups(P,iSet,d)
USE PP3D_MPI, ONLY :myid
USE var_QuadScalar, ONLY: gInterface
IMPLICIT NONE
REAL*8 P(3),Q(3),d,dist,dMindist,dMinDAux
INTEGER i,j,k,kk,iMinDist,iP,iG,iAux,iTT
REAL*8 R(4,3),dTriang(3,3),dTestVector(4),PP(3)
REAL*8 AX2,AY2,AZ2,AX3,AY3,AZ3,dAux,ProjPdist,ProjP(3),dSgn,dSaux
!!!!!
INTEGER nXX
REAL*8 dEps
PARAMETER (nXX = 10, dEps=1d-10)
REAL*8  :: dMonitor(nXX+1),ddSign(nXX+1)
INTEGER :: iMonitor(nXX+1),iSet(nXX)
LOGICAL bProblem

dMonitor = 1d30
iMonitor = 0

DO i=1,gInterface%nG-1
 Q = gInterface%Y(:,i)
 dist=sqrt((P(1)-Q(1))**2d0 + (P(2)-Q(2))**2d0 + (P(3)-Q(3))**2d0)
 IF (dist.lt.dMonitor(nXX)) THEN
  dMonitor(nXX) = dist
  iMonitor(nXX) = i
  DO kk=nXX-1,1,-1
   IF (dMonitor(kk+1).lt.dMonitor(kk)) THEN
    iAux = iMonitor(kk)
    dAux = dMonitor(kk)
    iMonitor(kk) = iMonitor(kk+1)
    dMonitor(kk) = dMonitor(kk+1)
    iMonitor(kk+1) = iAux
    dMonitor(kk+1) = dAux
   END IF
  END DO
 END IF
END DO

iSet(:) = iMonitor(1:nXX)
d = dMonitor(1)

END
!-----------------------------------------------------------------------------------------
SUBROUTINE GetDistanceToInterfaceInGivenGroups(P,iSet,d)
USE PP3D_MPI, ONLY :myid
USE var_QuadScalar, ONLY: gInterface
IMPLICIT NONE
REAL*8 P(3),Q(3),d,dist,dMindist,dMinDAux
INTEGER i,j,k,kk,iMinDist,iP,iG,iAux,iTT
REAL*8 R(4,3),dTriang(3,3),dTestVector(4),PP(3)
REAL*8 AX2,AY2,AZ2,AX3,AY3,AZ3,dAux,ProjPdist,ProjP(3),dSgn,dSaux
!!!!!
INTEGER nXX
REAL*8 dEps
PARAMETER (nXX = 10, dEps=1d-10)
REAL*8  :: dMonitor(nXX+1),ddSign(nXX+1)
INTEGER :: iMonitor(nXX+1),iSet(nXX),jSet(nXX)
LOGICAL bProblem

dMonitor = 1d30
iMonitor = 0

DO k=1,nXX
 iG = iSet(k)
 IF (iG.EQ.0) CYCLE
 DO i=gInterface%L(iG),gInterface%L(iG+1)-1
  iP = 3*(i-1)
  R(1,:) = gInterface%x(:,iP+1)
  R(2,:) = gInterface%x(:,iP+2)
  R(3,:) = gInterface%x(:,iP+3)
  R(4,:) = (R(1,:) + R(2,:) + R(3,:))/3d0
  R(1,:) = 0.5d0*(R(1,:) + R(4,:))
  R(2,:) = 0.5d0*(R(2,:) + R(4,:))
  R(3,:) = 0.5d0*(R(3,:) + R(4,:))

  dMinDAux = 1d30
  DO j=1,4
   Q = R(j,:)
   dist=sqrt((P(1)-Q(1))**2d0 + (P(2)-Q(2))**2d0 + (P(3)-Q(3))**2d0)
   IF (dist.lt.dMinDAux) dMinDaux = dist
  END DO

  dist = dMinDaux
  IF (dist.lt.dMonitor(nXX)) THEN
   dMonitor(nXX) = dist
   iMonitor(nXX) = iP
   DO kk=nXX-1,1,-1
    IF (dMonitor(kk+1).lt.dMonitor(kk)) THEN
     iAux = iMonitor(kk)
     dAux = dMonitor(kk)
     iMonitor(kk) = iMonitor(kk+1)
     dMonitor(kk) = dMonitor(kk+1)
     iMonitor(kk+1) = iAux
     dMonitor(kk+1) = dAux
     END IF
    END DO
   END IF

 END DO
END DO

jSet(:) = iMonitor(1:nXX)

dMonitor = 1d30
ddSign   = 0d0
iMonitor = 0

dMinDist = 1d30

DO k=1,nXX

iMinDist = jSet(k)

dTriang(:,1) = gInterface%x(:,iMinDist+1)
dTriang(:,2) = gInterface%x(:,iMinDist+2)
dTriang(:,3) = gInterface%x(:,iMinDist+3)

AX2=dTriang(1,2)-dTriang(1,1)
AY2=dTriang(2,2)-dTriang(2,1)
AZ2=dTriang(3,2)-dTriang(3,1)
AX3=dTriang(1,3)-dTriang(1,1)
AY3=dTriang(2,3)-dTriang(2,1)
AZ3=dTriang(3,3)-dTriang(3,1)
dTestVector(1)=(AY3*AZ2)-(AZ3*AY2)
dTestVector(2)=(AZ3*AX2)-(AX3*AZ2)
dTestVector(3)=(AX3*AY2)-(AY3*AX2)
dTestVector(4) = -(dTestVector(1)*dTriang(1,1) + dTestVector(2)*dTriang(2,1) + dTestVector(3)*dTriang(3,1))

dAux = SQRT(dTestVector(1)**2d0 + dTestVector(2)**2d0+dTestVector(3)**2d0)
ProjPdist = (dTestVector(4) + dTestVector(1)*P(1) + dTestVector(2)*P(2) + dTestVector(3)*P(3))/dAux

ProjP(1) = P(1) - ProjPdist*dTestVector(1)/dAux
ProjP(2) = P(2) - ProjPdist*dTestVector(2)/dAux
ProjP(3) = P(3) - ProjPdist*dTestVector(3)/dAux

IF (SameSide(ProjP,dTriang(:,1), dTriang(:,2),dTriang(:,3)).and.&
    SameSide(ProjP,dTriang(:,2), dTriang(:,1),dTriang(:,3)).and.&
    SameSide(ProjP,dTriang(:,3), dTriang(:,1),dTriang(:,2))) THEN
    dist = ABS(ProjPdist)
    IF (dist.LT.dMonitor(nXX)) THEN
     dMinDist = dist 
     PP = ProjP
     iTT = iMinDist
     dSgn = ReturnSign(dTriang(:,1),dTriang(:,2),dTriang(:,3),P)

   dMonitor(nXX) = dist
   ddSign(nXX)   = dSgn
   iMonitor(nXX) = iTT
  
   DO kk=nXX-1,1,-1
    IF (dMonitor(kk+1).lt.dMonitor(kk)) THEN
     iAux = iMonitor(kk)
     dAux = dMonitor(kk)
     dSaux= ddSign(kk)
     iMonitor(kk) = iMonitor(kk+1)
     dMonitor(kk) = dMonitor(kk+1)
     ddSign(kk)   = ddSign(kk+1)
     iMonitor(kk+1) = iAux
     dMonitor(kk+1) = dAux
     ddSign(kk+1)   = dSaux
     END IF
    END DO

!      pause
    END IF
    GOTO 1
ELSE
 ProjP = ProjPointToLine(P,dTriang(:,1),dTriang(:,2))
 IF (.NOT.ISNAN(ProjP(1))) THEN
  dist= SQRT((ProjP(1)-P(1))**2d0 + (ProjP(2)-P(2))**2d0+(ProjP(3)-P(3))**2d0)
  IF (dist.LT.dMonitor(nXX)) THEN
   dMinDist = dist
   PP = ProjP
   iTT = iMinDist
   dSgn = ReturnSign(dTriang(:,1),dTriang(:,2),dTriang(:,3),P)

   dMonitor(nXX) = dist
   ddSign(nXX)   = dSgn
   iMonitor(nXX) = iTT
  
   DO kk=nXX-1,1,-1
    IF (dMonitor(kk+1).lt.dMonitor(kk)) THEN
     iAux = iMonitor(kk)
     dAux = dMonitor(kk)
     dSaux= ddSign(kk)
     iMonitor(kk) = iMonitor(kk+1)
     dMonitor(kk) = dMonitor(kk+1)
     ddSign(kk)   = ddSign(kk+1)
     iMonitor(kk+1) = iAux
     dMonitor(kk+1) = dAux
     ddSign(kk+1)   = dSaux
     END IF
    END DO

  END IF
  GOTO 1
 END IF
 ProjP = ProjPointToLine(P,dTriang(:,2),dTriang(:,3))
 IF (.NOT.ISNAN(ProjP(1))) THEN
  dist= SQRT((ProjP(1)-P(1))**2d0 + (ProjP(2)-P(2))**2d0+(ProjP(3)-P(3))**2d0)
  IF (dist.LT.dMonitor(nXX)) THEN
   dMinDist = dist
   PP = ProjP
   iTT = iMinDist
   dSgn = ReturnSign(dTriang(:,1),dTriang(:,2),dTriang(:,3),P)

   dMonitor(nXX) = dist
   ddSign(nXX)   = dSgn
   iMonitor(nXX) = iTT
  
   DO kk=nXX-1,1,-1
    IF (dMonitor(kk+1).lt.dMonitor(kk)) THEN
     iAux = iMonitor(kk)
     dAux = dMonitor(kk)
     dSaux= ddSign(kk)
     iMonitor(kk) = iMonitor(kk+1)
     dMonitor(kk) = dMonitor(kk+1)
     ddSign(kk)   = ddSign(kk+1)
     iMonitor(kk+1) = iAux
     dMonitor(kk+1) = dAux
     ddSign(kk+1)   = dSaux
     END IF
    END DO

  END IF
  GOTO 1
 END IF
 ProjP = ProjPointToLine(P,dTriang(:,3),dTriang(:,1))
 IF (.NOT.ISNAN(ProjP(1))) THEN
  dist= SQRT((ProjP(1)-P(1))**2d0 + (ProjP(2)-P(2))**2d0+(ProjP(3)-P(3))**2d0)
  IF (dist.LT.dMonitor(nXX)) THEN
   dMinDist = dist
   PP = ProjP
   iTT = iMinDist
   dSgn = ReturnSign(dTriang(:,1),dTriang(:,2),dTriang(:,3),P)

   dMonitor(nXX) = dist
   ddSign(nXX)   = dSgn
   iMonitor(nXX) = iTT
  
   DO kk=nXX-1,1,-1
    IF (dMonitor(kk+1).lt.dMonitor(kk)) THEN
     iAux = iMonitor(kk)
     dAux = dMonitor(kk)
     dSaux= ddSign(kk)
     iMonitor(kk) = iMonitor(kk+1)
     dMonitor(kk) = dMonitor(kk+1)
     ddSign(kk)   = ddSign(kk+1)
     iMonitor(kk+1) = iAux
     dMonitor(kk+1) = dAux
     ddSign(kk+1)   = dSaux
     END IF
    END DO

  END IF
  GOTO 1
 END IF

 DO i=1,3
  dist = SQRT((dTriang(1,i)-P(1))**2d0 + (dTriang(2,i)-P(2))**2d0+(dTriang(3,i)-P(3))**2d0)
  IF (dist.lt.dMonitor(nXX)) THEN
   dMinDist = dist
   PP = dTriang(:,i)
   iTT = iMinDist
   dSgn = ReturnSign(dTriang(:,1),dTriang(:,2),dTriang(:,3),P)

   dMonitor(nXX) = dist
   ddSign(nXX)   = dSgn
   iMonitor(nXX) = iTT
  
   DO kk=nXX-1,1,-1
    IF (dMonitor(kk+1).lt.dMonitor(kk)) THEN
     iAux = iMonitor(kk)
     dAux = dMonitor(kk)
     dSaux= ddSign(kk)
     iMonitor(kk) = iMonitor(kk+1)
     dMonitor(kk) = dMonitor(kk+1)
     ddSign(kk)   = ddSign(kk+1)
     iMonitor(kk+1) = iAux
     dMonitor(kk+1) = dAux
     ddSign(kk+1)   = dSaux
     END IF
    END DO

  END IF
 END DO
 GOTO 1
END IF

1 CONTINUE

END DO

dSgn = 0d0
bProblem = .FALSE.
DO kk=2,nXX
 IF (ABS(dMonitor(1)-dMonitor(kk)).LT.dEps) THEN
  IF (ddSign(1)*ddSign(kk).LT.0d0) bProblem = .TRUE.
 ELSE
  EXIT
 END IF
END DO

IF (bProblem) THEN

  Q = 0d0
  Q = Q + gInterface%x(:,iMonitor(1)+1)
  Q = Q + gInterface%x(:,iMonitor(1)+2)
  Q = Q + gInterface%x(:,iMonitor(1)+3)

  iAux = 3
  DO kk=2,nXX
   IF (ABS(dMonitor(1)-dMonitor(kk)).LT.dEps) THEN
    iAux = iAux + 3
    Q = Q + gInterface%x(:,iMonitor(kk)+1)
    Q = Q + gInterface%x(:,iMonitor(kk)+2)
    Q = Q + gInterface%x(:,iMonitor(kk)+3)
   ELSE
    EXIT
   END IF
  END DO

  Q = Q/DBLE(iAux)

  dTriang(:,1) = gInterface%x(:,iMonitor(1)+1)
  dTriang(:,2) = gInterface%x(:,iMonitor(1)+2)
  dTriang(:,3) = gInterface%x(:,iMonitor(1)+3)
  dSgn = ReturnSign(dTriang(:,1),dTriang(:,2),dTriang(:,3),Q)

  d = -dSgn*dMonitor(1)

ELSE

 d = ddSign(1)*dMonitor(1)

END IF

!   IF (myid.eq.1.and.p(2).lt.-0.078d0.AND.d.GT.0d0) THEN
!    write(*,'(A,10ES14.4)') 'o-o', P,dSgn,ddSign(1)
!    WRITE(*,'(10ES14.4)') dMonitor(1:10)
!    WRITE(*,'(10I14)') iMonitor(1:10)
!    WRITE(*,'(A)')  "---------------------"
!    write(*,'(10ES14.4)') gInterface%x(:,iMonitor(1)+1)
!    write(*,'(10ES14.4)') gInterface%x(:,iMonitor(1)+2)
!    write(*,'(10ES14.4)') gInterface%x(:,iMonitor(1)+3)
!    DO kk=2,nXX
!     IF (ABS(dMonitor(1)-dMonitor(kk)).LT.dEps) THEN
!      write(*,'(10ES14.4)') gInterface%x(:,iMonitor(kk)+1)
!      write(*,'(10ES14.4)') gInterface%x(:,iMonitor(kk)+2)
!      write(*,'(10ES14.4)') gInterface%x(:,iMonitor(kk)+3)
!     ELSE
!      EXIT
!     END IF
!    END DO
!    write(*,'(10ES14.4)') Q
! 
!    pause
! !    dSgn = ddSign(1)
!   END IF



 CONTAINS

 LOGICAL FUNCTION SameSide(p1,p2, a,b)
 IMPLICIT NONE
 REAL*8 p1(3),p2(3), a(3),b(3)
 REAL*8 cp1(3),cp2(3),ddd
 REAL*8 AX2,AY2,AZ2,AX3,AY3,AZ3

 AX2=b(1)  - a(1)
 AY2=b(2)  - a(2)
 AZ2=b(3)  - a(3)
 AX3=p1(1) - a(1)
 AY3=p1(2) - a(2)
 AZ3=p1(3) - a(3)
 cp1 = [(AY3*AZ2)-(AZ3*AY2),(AZ3*AX2)-(AX3*AZ2),(AX3*AY2)-(AY3*AX2)] !CrossProduct(b-a, p1-a)
 AX2=b(1)  - a(1)
 AY2=b(2)  - a(2)
 AZ2=b(3)  - a(3)
 AX3=p2(1) - a(1)
 AY3=p2(2) - a(2)
 AZ3=p2(3) - a(3)
 cp2 = [(AY3*AZ2)-(AZ3*AY2),(AZ3*AX2)-(AX3*AZ2),(AX3*AY2)-(AY3*AX2)] !CrossProduct(b-a, p1-a)
 ddd = cp1(1)*cp2(1) + cp1(2)*cp2(2) + cp1(3)*cp2(3)
 IF (ddd.GE.0) THEN
  SameSide = .true.
 ELSE
  SameSide = .false.
 END IF

 END FUNCTION SameSide

 FUNCTION ProjPointToLine(p,a,b)
 IMPLICIT NONE
 REAL*8 ProjPointToLine(3)
 REAL*8 p(3), a(3),b(3)
 REAL*8 ab(3), ap(3),ddd,t,dot

  ab = b-a
  ddd = ab(1)**2d0 + ab(2)**2d0 + ab(3)**2d0
  ap = p-a
  dot = ab(1)*ap(1) + ab(2)*ap(2) + ab(3)*ap(3)
  t = dot/ddd;
  ProjPointToLine = 0d0/0d0
  IF     (t.ge.0d0.AND.t.le.1d0) THEN
    ProjPointToLine = A + t * AB
  END IF

 END FUNCTION ProjPointToLine

 FUNCTION ReturnSign(P1,P2,P3,PQ)
  IMPLICIT NONE
  REAL*8 ReturnSign
  REAL*8 P1(3),P2(3),P3(3),PQ(3),D,R(3)
  REAL*8 AX2,AY2,AZ2,AX3,AY3,AZ3,dC(4),dAux

  AX2=P2(1) - P1(1)
  AY2=P2(2) - P1(2)
  AZ2=P2(3) - P1(3)
  AX3=P3(1) - P1(1)
  AY3=P3(2) - P1(2)
  AZ3=P3(3) - P1(3)

  dC(1)=(AY3*AZ2)-(AZ3*AY2)
  dC(2)=(AZ3*AX2)-(AX3*AZ2)
  dC(3)=(AX3*AY2)-(AY3*AX2)
  dC(4) = -(dC(1)*P1(1) + dC(2)*P1(2) + dC(3)*P1(3))

  dAux = dC(4) + dC(1)*PQ(1) + dC(2)*PQ(2) + dC(3)*PQ(3)
  !write(*,*) daux
  IF (daux.EQ.0d0) daux = 1d0
  ReturnSign = dAux/ABS(dAux)

 END FUNCTION ReturnSign

END
!-----------------------------------------------------------------------------------------
SUBROUTINE TetraInterfacePoints(Frac,LS,dcorvg,kvert,kedge,karea,kadj,&
           nel,nvt,nat,net,dIntegralVolume,bCounterRun)
USE PP3D_MPI, ONLY:myid,master,myMPI_Barrier
USE var_QuadScalar, ONLY:lInterface,NLMAX,NLMIN

IMPLICIT NONE

REAL*8 dIntegralVolume(3)
INTEGER nel,nvt,nat,net,nPoints,nElems
REAL*8  dcorvg(3,*),LS(*),Frac(*)
INTEGER kvert(8,*),kedge(12,*),karea(6,*),kadj(6,*)
REAL*8 GlobalElemCoor(3,27),LocalElemCoor(3,27),dLocLS(27),dFaceLocLS(9,6)
INTEGER iC,jC,kC,i8,iSubElements(3,8),indice(3,27)
INTEGER i,j,k,l,kk,jj,minus,plus
INTEGER jvt1,jvt2,jvt3
CHARACTER cMyid*2,cMyTime*3
REAL*8 P(3,3),R(3,3),V(3)
REAL*8 dJacComp(8,3)
INTEGER iRecursivity
INTEGER EdgeToFace(2,12),FaceToEdge(4,6),VertToEdge(2,12),iEdge(3,10),iEdgeToFace(4,6),Q2ToFace(9,6),iAxis(3,6)
INTEGER EdgeToVert(3,8)
REAL*8 tnt1,tnt0
!------------------------------------------------------
REAL*8 fVal_Q2_3D,fVal_Q2_2D
!------------------------------------------------------
DATA LocalElemCoor/-1.0,-1.0,-1.0,1.0,-1.0,-1.0,1.0, 1.0,-1.0,-1.0, 1.0,-1.0,&
                   -1.0,-1.0, 1.0,1.0,-1.0, 1.0,1.0, 1.0, 1.0,-1.0, 1.0, 1.0,& ! last nvt
                    0.0,-1.0,-1.0,1.0, 0.0,-1.0,0.0, 1.0,-1.0,-1.0, 0.0,-1.0,&
                   -1.0,-1.0, 0.0,1.0,-1.0, 0.0,1.0, 1.0, 0.0,-1.0, 1.0, 0.0,&
                    0.0,-1.0, 1.0,1.0, 0.0, 1.0,0.0, 1.0, 1.0,-1.0, 0.0, 1.0,& ! last net
                    0.0, 0.0,-1.0,0.0,-1.0, 0.0,1.0, 0.0, 0.0,0.0, 1.0, 0.0,-1.0, 0.0, 0.0,0.0, 0.0, 1.0,& ! last nat
                    0.0, 0.0, 0.0 / ! last nel
DATA VertToEdge/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
DATA EdgeToFace/1,2,1,3,1,4,1,5, 5,2,2,3,3,4,4,5, 2,6,3,6,4,6,5,6/
DATA FaceToEdge/1,2,3,4, 1,6,9,5, 2,7,10,6, 3,8,11,7, 4,5,12,8, 9,10,11,12/
DATA iEdgeToFace/1,3,2,4,5,6,1,9,6,7,2,10,7,8,3,11,8,5,4,12,9,11,10,12/
DATA Q2ToFace/1,2,3,4,  9,10,11,12, 21, &
              1,2,6,5,  9,14,17,13, 22, &
              2,3,7,6, 10,15,18,14, 23, &
              4,3,7,8, 11,15,19,16, 24, &
              1,4,8,5, 12,16,20,13, 25, &
              5,6,7,8, 17,18,19,20, 26/
DATA iAxis/1,2,3, 1,3,2, 2,3,1, 1,3,2, 2,3,1, 1,2,3/
DATA EdgeToVert/12,9,13, 9,10,14, 10,11,15, 11,12,16, 20,17,13, 17,18,14, 18,19,15, 19,20,16/

DATA iSubElements/-1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1, -1,-1,1, 1,-1,1, 1,1,1, -1,1,1/
TYPE tElement
 REAL*8 FineGlobalElemCoor(3,5,5,5),FineLocalElemCoor(3,5,5,5),FinedLocLS(5,5,5)
 REAL*8 dBox(2,3),dX(5),dY(5),dZ(5)
 REAL*8 GlobalElemCoor(3,27),LocalElemCoor(3,27),dLocLS(27)
 REAL*8 Q1(3,24),Q2(3,24)
 INTEGER iEdgeInfo(12)
END TYPE tElement
TYPE(tElement) myEle

REAL*8 dWV(10),dWC(3,10),DFaceVal,dCoorMinMax(2,2),TrialPoint(2)
logical bOK(20,20)
REAL*8 :: DEPS = 1d-10
INTEGER jSide
!
REAL*8 dMinVertSTR(4)
INTEGER jMinVertSTR
!
REAL*8 :: dMinEdgeLength=0.0005d0
LOGICAL bExit(2)

INTEGER Tetra(4,40),ium,iup
CHARACTER cFile*21
INTEGER iTriangles,iiTriangles,iPoints,iGroups
REAL*8  dVolM,dVolP,dVolCM(3),dFullVol
REAL*8 dAux,dHexVolume,dHexCenter(3),dVolElem(2)

LOGICAL bCounterRun
!!! ----------------------------------------------------

! WRITE(cFile(1:21),('(A,I4.4,A)')) '_gmv/Surface_',myid,'.gmv'
! OPEN(FILE = cFile,UNIT=919)
! WRITE(919,'(A)')'gmvinput ascii'
! WRITE(919,*)  'polygons'

CALL myMPI_Barrier()
CALL ZTIME (tnt0)

 tetra(:, 1) = [ 1, 9,21,22]
 tetra(:, 2) = [ 1,21,25,22]
 tetra(:, 3) = [ 1,21,12,25]
 tetra(:, 4) = [ 1,22,25,13]
 tetra(:, 5) = [21,25,22,27]
 tetra(:, 6) = [ 2,10,21,23]
 tetra(:, 7) = [ 2,21,22,23]
 tetra(:, 8) = [ 2,21, 9,22]
 tetra(:, 9) = [ 2,23,22,14]
 tetra(:,10) = [21,22,23,27]
 tetra(:,11) = [ 3,11,21,24]
 tetra(:,12) = [ 3,21,23,24]
 tetra(:,13) = [ 3,21,10,23]
 tetra(:,14) = [ 3,24,23,15]
 tetra(:,15) = [21,23,24,27]
 tetra(:,16) = [ 4,12,21,25]
 tetra(:,17) = [ 4,21,24,25]
 tetra(:,18) = [ 4,21,11,24]
 tetra(:,19) = [ 4,25,24,16]
 tetra(:,20) = [21,24,25,27] 
 tetra(:,21) = [ 5,17,26,22]
 tetra(:,22) = [ 5,26,25,22]
 tetra(:,23) = [ 5,26,20,25]
 tetra(:,24) = [ 5,22,25,13]
 tetra(:,25) = [26,25,22,27]
 tetra(:,26) = [ 6,18,26,23]
 tetra(:,27) = [ 6,26,22,23]
 tetra(:,28) = [ 6,26,17,22]
 tetra(:,29) = [ 6,23,22,14]
 tetra(:,30) = [26,22,23,27]
 tetra(:,31) = [ 7,19,26,24]
 tetra(:,32) = [ 7,26,23,24]
 tetra(:,33) = [ 7,26,18,23]
 tetra(:,34) = [ 7,24,23,15]
 tetra(:,35) = [26,23,24,27]
 tetra(:,36) = [ 8,20,26,25] 
 tetra(:,37) = [ 8,26,24,25]
 tetra(:,38) = [ 8,26,19,24]
 tetra(:,39) = [ 8,25,24,16]
 tetra(:,40) = [26,24,25,27]

iPoints = 0
iGroups = 0
iiTriangles = 0
dVolM = 0d0
dVolP = 0d0
!!!!
dVolCM = 0d0
!!!
dFullVol = 0d0

DO i=1,nel

 ! Extract the Level Set values and the Global Coordinates
 DO j=1,8
   k = kvert(j,i)
   dLocLS(j) = LS(k)
   GlobalElemCoor(:,j) = dcorvg(:,k)
 END DO
 DO j=1,12
   k = nvt+kedge(j,i)
   dLocLS(8+j) = LS(k)
   GlobalElemCoor(:,8+j) = dcorvg(:,k)
 END DO
 DO j=1,6
   k = nvt + net + karea(j,i)
   dLocLS(8+12+j) = LS(k)
   GlobalElemCoor(:,8+12+j) = dcorvg(:,k)
 END DO
 dLocLS(27) = LS(nvt + net + nat + i)
 GlobalElemCoor(:,27) = dcorvg(:,nvt + net + nat + i)

 DO j=1,6
   DO k=1,9
    dFaceLocLS(k,j) = dLocLS(Q2ToFace(k,j))
  END DO
 END DO

 CALL ElementJacobianSetUp(GlobalElemCoor)

 ! interpolation of LevSet onto the finer mesh
 myEle%dBox(1,:) = -1d0
 myEle%dBox(2,: )=  1d0
 DO iC =1,5
  myEle%dX(iC) = myEle%dBox(1,1)+DBLE(iC-1)*(myEle%dBox(2,1)-myEle%dBox(1,1))/DBLE(4)
  myEle%dY(iC) = myEle%dBox(1,2)+DBLE(iC-1)*(myEle%dBox(2,2)-myEle%dBox(1,2))/DBLE(4)
  myEle%dZ(iC) = myEle%dBox(1,3)+DBLE(iC-1)*(myEle%dBox(2,3)-myEle%dBox(1,3))/DBLE(4)
 END DO

 minus = 0
 plus = 0
 DO iC=1,5
  DO jC=1,5
   DO kC=1,5
    myEle%FineLocalElemCoor(:,iC,jC,kC) = [myEle%dX(iC),myEle%dY(jC),myEle%dZ(kC)]
    myEle%FinedLocLS(iC,jC,kC) = fVal_Q2_3D(myEle%dX(iC),myEle%dY(jC),myEle%dZ(kC),dLocLS)
    IF (myEle%FinedLocLS(iC,jC,kC).LT.0d0) THEN
     minus = minus + 1
    END IF
    IF (myEle%FinedLocLS(iC,jC,kC).GT.0d0) THEN
     plus = plus + 1
    END IF
   END DO
  END DO
 END DO

 dHexVolume = GetHexaVolume(GlobalElemCoor(:,1:8))
 CALL GetHexaCenter(GlobalElemCoor,dHexCenter)

 iup = plus; ium = minus

 IF ((plus.eq.0).or.(minus.eq.0)) THEN

  IF (bCounterRun) THEN
  ELSE
   IF (plus.eq.0) THEN
    dVolM  = dVolM  + dHexVolume
    dVolCM = dVolCM + dHexCenter
    Frac(i) = 0d0
   END IF
   IF (minus.eq.0) THEN
    dVolP = dVolP + dHexVolume
    Frac(i) = 1d0
   END IF
  END IF

  CYCLE

 END IF

 ! recompute real coordinates for the finer level mesh
 CALL SPoint_Q2_3D(myEle%FineLocalElemCoor,myEle%FineGlobalElemCoor,125)
 
 dVolElem = 0d0

 DO i8=1,8

    ! Extract the indices ....
   indice(:,1 ) = [3+iSubElements(1,i8)-1, 3+iSubElements(2,i8)-1, 3+iSubElements(3,i8)-1]
   indice(:,2 ) = [3+iSubElements(1,i8)+1, 3+iSubElements(2,i8)-1, 3+iSubElements(3,i8)-1]
   indice(:,3 ) = [3+iSubElements(1,i8)+1, 3+iSubElements(2,i8)+1, 3+iSubElements(3,i8)-1]
   indice(:,4 ) = [3+iSubElements(1,i8)-1, 3+iSubElements(2,i8)+1, 3+iSubElements(3,i8)-1]
   indice(:,5 ) = [3+iSubElements(1,i8)-1, 3+iSubElements(2,i8)-1, 3+iSubElements(3,i8)+1]
   indice(:,6 ) = [3+iSubElements(1,i8)+1, 3+iSubElements(2,i8)-1, 3+iSubElements(3,i8)+1]
   indice(:,7 ) = [3+iSubElements(1,i8)+1, 3+iSubElements(2,i8)+1, 3+iSubElements(3,i8)+1]
   indice(:,8 ) = [3+iSubElements(1,i8)-1, 3+iSubElements(2,i8)+1, 3+iSubElements(3,i8)+1]
 
   indice(:,9 ) = [3+iSubElements(1,i8)  , 3+iSubElements(2,i8)-1, 3+iSubElements(3,i8)-1]
   indice(:,10) = [3+iSubElements(1,i8)+1, 3+iSubElements(2,i8)  , 3+iSubElements(3,i8)-1]
   indice(:,11) = [3+iSubElements(1,i8)  , 3+iSubElements(2,i8)+1, 3+iSubElements(3,i8)-1]
   indice(:,12) = [3+iSubElements(1,i8)-1, 3+iSubElements(2,i8)  , 3+iSubElements(3,i8)-1]
   indice(:,13) = [3+iSubElements(1,i8)-1, 3+iSubElements(2,i8)-1, 3+iSubElements(3,i8)  ]
   indice(:,14) = [3+iSubElements(1,i8)+1, 3+iSubElements(2,i8)-1, 3+iSubElements(3,i8)  ]
   indice(:,15) = [3+iSubElements(1,i8)+1, 3+iSubElements(2,i8)+1, 3+iSubElements(3,i8)  ]
   indice(:,16) = [3+iSubElements(1,i8)-1, 3+iSubElements(2,i8)+1, 3+iSubElements(3,i8)  ]
   indice(:,17) = [3+iSubElements(1,i8)  , 3+iSubElements(2,i8)-1, 3+iSubElements(3,i8)+1]
   indice(:,18) = [3+iSubElements(1,i8)+1, 3+iSubElements(2,i8)  , 3+iSubElements(3,i8)+1]
   indice(:,19) = [3+iSubElements(1,i8)  , 3+iSubElements(2,i8)+1, 3+iSubElements(3,i8)+1]
   indice(:,20) = [3+iSubElements(1,i8)-1, 3+iSubElements(2,i8)  , 3+iSubElements(3,i8)+1]

   indice(:,21) = [3+iSubElements(1,i8)  , 3+iSubElements(2,i8)  , 3+iSubElements(3,i8)-1]
   indice(:,22) = [3+iSubElements(1,i8)  , 3+iSubElements(2,i8)-1, 3+iSubElements(3,i8)  ]
   indice(:,23) = [3+iSubElements(1,i8)+1, 3+iSubElements(2,i8)  , 3+iSubElements(3,i8)  ]
   indice(:,24) = [3+iSubElements(1,i8)  , 3+iSubElements(2,i8)+1, 3+iSubElements(3,i8)  ]
   indice(:,25) = [3+iSubElements(1,i8)-1, 3+iSubElements(2,i8)  , 3+iSubElements(3,i8)  ]
   indice(:,26) = [3+iSubElements(1,i8)  , 3+iSubElements(2,i8)  , 3+iSubElements(3,i8)+1]

   indice(:,27) = [3+iSubElements(1,i8)  , 3+iSubElements(2,i8)  , 3+iSubElements(3,i8)  ]

   DO j=1,27
    myEle%GlobalElemCoor(:,j) = myEle%FineGlobalElemCoor(:,indice(1,j),indice(2,j),indice(3,j))
    myEle%LocalElemCoor (:,j) = myEle%FineLocalElemCoor (:,indice(1,j),indice(2,j),indice(3,j))
    myEle%dLocLS        (  j) = myEle%FinedLocLS        (  indice(1,j),indice(2,j),indice(3,j))
   END DO

   plus = 0; minus = 0
   DO j=1,27
    IF (myEle%dLocLS(j) .LT.0d0) THEN
     minus = minus + 1
    ELSE
     plus = plus + 1
    END IF
   END DO

   IF ((plus.eq.0).or.(minus.eq.0)) THEN

    IF (bCounterRun) THEN
    ELSE
     dAux = GetHexaVolume(myEle%GlobalElemCoor(:,1:8))
     CALL GetHexaCenter(myEle%GlobalElemCoor,dHexCenter)

     IF (plus.eq.0) THEN
      dVolM = dVolM + dAux
      dVolCM = dVolCM + dHexCenter
      dVolElem(2) = dVolElem(2) + dAux
     END IF
     IF (minus.eq.0) THEN
      dVolP = dVolP + dAux
      dVolElem(1) = dVolElem(1) + dAux
     END IF
    END IF

   ELSE

    CALL ProcessHexa()

   END IF

 END DO

 IF (bCounterRun) THEN
 ELSE
  IF (ABS((dHexVolume-(dVolElem(1)+dVolElem(2)))/dHexVolume).GT.1d-2) THEN
!    WRITE(*,*) 'o o ....',100*(dHexVolume-(dVolElem(1)+dVolElem(2)))/dHexVolume,iup,ium
  END IF
  Frac(i) = dVolElem(1)/(dVolElem(1)+dVolElem(2))
 END IF

END DO

99 CONTINUE

lInterface%nG = iGroups
lInterface%nT = iiTriangles

 IF (bCounterRun) THEN
 ELSE
!   WRITE(*,*) myid,dVolCM
!   pause
 END IF

dIntegralVolume = [dVolP,dVolM,dVolCM(3)]

CALL myMPI_Barrier()
CALL ZTIME (tnt1)
CALL myMPI_Barrier()

 CONTAINS
! --------------------------------------------------------------------------------------------
SUBROUTINE ProcessHexa()
INTEGER iiP,iiM,iTet
REAL*8 coor(3,4),dist(4),dP(3),Point(3,6),xPoint(3,6),ePoint(3,6),dS
REAL*8 xTetcoor(3,4),yTetcoor(3,4)
REAL*8 dVol,dVolTet,dVolTetThird,dCenterTetLarge(3),dCenterTetSmall(3)
LOGICAL bEdge(6),bTurn
INTEGER ii,iOmax,ivt(4)

iTriangles = 0

DO iTet=1,40

 iiP = 0; iiM = 0
 DO j=1,4
  IF (myEle%dLocLS(tetra(j,iTet)).LT.0d0) THEN
   iiM = iiM + 1
  ELSE
   iiP = iiP + 1
  END IF
 END DO

 DO j=1,4 
  dist    (j  ) = myEle%dLocLS        (  tetra(j,iTet))
  coor    (:,j) = myEle%LocalElemCoor (:,tetra(j,iTet))
  xTetcoor(:,j) = myEle%GlobalElemCoor(:,tetra(j,iTet))
 END DO

 dAux = 0d0
 DO j=1,4
  IF (ABS(dist(j)).GT.dAux) THEN
   iOmax = j
   dAux = ABS(dist(j))
  END IF
 END DO

 dVolTet = GetTetraVolume(xTetCoor)
 dCenterTetLarge = 0.25d0*(xTetcoor(:,1)+xTetcoor(:,2)+xTetcoor(:,3)+xTetcoor(:,4))*dVolTet

 IF ((iiP.NE.4).AND.(iiM.NE.4)) THEN

  IF (bCounterRun) THEN

   IF (iiP.EQ.2) THEN 
     iTriangles = iTriangles + 2
   ELSE
     iTriangles = iTriangles + 1
   END IF

  ELSE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 j= 0
 bEdge=.FALSE.
 IF (dist(1)*dist(2).LT.0d0) THEN
   bEdge(1)=.TRUE.
   CALL GetPoint(dist(1),dist(2),coor(:,1),coor(:,2),dP)
   j = j + 1
    Point(:,j) = dP
   ePoint(:,1) = dP
 END IF

 IF (dist(2)*dist(3).LT.0d0) THEN
   bEdge(2)=.TRUE.
   CALL GetPoint(dist(2),dist(3),coor(:,2),coor(:,3),dP)
   j = j + 1
   Point(:,j) = dP
   ePoint(:,2) = dP
 END IF

 IF (dist(1)*dist(3).LT.0d0) THEN
   bEdge(3)=.TRUE.
   CALL GetPoint(dist(3),dist(1),coor(:,3),coor(:,1),dP)
   j = j + 1
   Point(:,j) = dP
   ePoint(:,3) = dP
 END IF

 IF (dist(1)*dist(4).LT.0d0) THEN
   bEdge(4)=.TRUE.
   CALL GetPoint(dist(1),dist(4),coor(:,1),coor(:,4),dP)
   j = j + 1
   Point(:,j) = dP
   ePoint(:,4) = dP
 END IF

 IF (dist(2)*dist(4).LT.0d0) THEN
   bEdge(5)=.TRUE.
   CALL GetPoint(dist(2),dist(4),coor(:,2),coor(:,4),dP)
   j = j + 1
   Point(:,j) = dP
   ePoint(:,5) = dP
 END IF

 IF (dist(3)*dist(4).LT.0d0) THEN
   bEdge(6)=.TRUE.
   CALL GetPoint(dist(3),dist(4),coor(:,3),coor(:,4),dP)
   j = j + 1
   Point(:,j) = dP
   ePoint(:,6) = dP
 END IF

 IF (j.eq.3) THEN
  IF (iiP.EQ.1) THEN
   DO ii=1,4
    IF (dist(ii).GT.0d0) Point(:,4) = coor(:,ii)
   END DO
   dS =  1d0 
  END IF
  IF (iiM.EQ.1) THEN
   DO ii=1,4
    IF (dist(ii).LT.0d0) Point(:,4) = coor(:,ii)
   END DO
   dS = -1d0 
  END IF

  CALL SPoint_Q2_3D(Point,xPoint,4)

  yTetcoor(:,1) = [xPoint(1,1),xPoint(2,1),xPoint(3,1)]
  yTetcoor(:,2) = [xPoint(1,2),xPoint(2,2),xPoint(3,2)]
  yTetcoor(:,3) = [xPoint(1,3),xPoint(2,3),xPoint(3,3)]
  yTetcoor(:,4) = [xPoint(1,4),xPoint(2,4),xPoint(3,4)]
  dVol = GetTetraVolume(yTetCoor)
  dCenterTetSmall = 0.25d0*(yTetcoor(:,1)+yTetcoor(:,2)+yTetcoor(:,3)+yTetcoor(:,4))*dVol

  IF (dS.gt.0d0) THEN
   dVolP  = dVolP  + dVol
   dVolM  = dVolM  + (dVolTet-dVol)
   dVolCM = dVolCM + (dCenterTetLarge-dCenterTetSmall)
   dVolElem(1) = dVolElem(1) + dVol
   dVolElem(2) = dVolElem(2) + (dVolTet-dVol)
  END IF
  IF (dS.lt.0d0) THEN
   dVolP  = dVolP  + (dVolTet-dVol)
   dVolM  = dVolM  + dVol
   dVolCM = dVolCM + dCenterTetSmall
   dVolElem(1) = dVolElem(1) + (dVolTet-dVol)
   dVolElem(2) = dVolElem(2) + dVol
  END IF

  CALL PointOrdering(xPoint(:,1),xPoint(:,2),xPoint(:,3),xTetcoor(:,iOmax),dist(iOmax),bTurn)
  IF (bTurn) THEN
!    WRITE(919,'(2i6,12E16.8)') 2,3,(xPoint(ii,1),xPoint(ii,2),xPoint(ii,3),ii=1,3)
   lInterface%X(:,iPoints+1) = xPoint(:,1)
   lInterface%X(:,iPoints+2) = xPoint(:,3)
   lInterface%X(:,iPoints+3) = xPoint(:,2)
  ELSE
!    WRITE(919,'(2i6,12E16.8)') 2,3,(xPoint(ii,1),xPoint(ii,3),xPoint(ii,2),ii=1,3)
   lInterface%X(:,iPoints+1) = xPoint(:,1)
   lInterface%X(:,iPoints+2) = xPoint(:,2)
   lInterface%X(:,iPoints+3) = xPoint(:,3)
  END IF
  iPoints    = iPoints + 3
  iTriangles = iTriangles + 1

 END IF

 IF (j.eq.4) THEN
   IF (bEdge(1).AND.bEdge(2)) THEN
    Point(:,1) = ePoint(:,1)
    Point(:,2) = ePoint(:,2)
    Point(:,3) = ePoint(:,6)
    Point(:,4) = ePoint(:,4)
    Point(:,5) = coor(:,1)
    Point(:,6) = coor(:,3)
    IF (dist(1).LT.0d0) dS=-1d0
    IF (dist(1).GT.0d0) dS=+1d0
   END IF

   IF (bEdge(1).AND.bEdge(3)) THEN
    Point(:,1) = ePoint(:,3)
    Point(:,2) = ePoint(:,1)
    Point(:,3) = ePoint(:,5)
    Point(:,4) = ePoint(:,6)
    Point(:,5) = coor(:,3)
    Point(:,6) = coor(:,2)
    IF (dist(3).LT.0d0) dS=-1d0
    IF (dist(3).GT.0d0) dS=+1d0
   END IF
   IF (bEdge(2).AND.bEdge(3)) THEN
    Point(:,1) = ePoint(:,2)
    Point(:,2) = ePoint(:,3)
    Point(:,3) = ePoint(:,4)
    Point(:,4) = ePoint(:,5)
    Point(:,5) = coor(:,2)
    Point(:,6) = coor(:,1)
    IF (dist(2).LT.0d0) dS=-1d0
    IF (dist(2).GT.0d0) dS=+1d0
   END IF

   CALL SPoint_Q2_3D(Point,xPoint,6)

   dVol = 0d0
   dCenterTetSmall = 0d0

   yTetcoor(:,1) = [xPoint(1,1),xPoint(2,1),xPoint(3,1)]
   yTetcoor(:,2) = [xPoint(1,4),xPoint(2,4),xPoint(3,4)]
   yTetcoor(:,3) = [xPoint(1,5),xPoint(2,5),xPoint(3,5)]
   yTetcoor(:,4) = [xPoint(1,6),xPoint(2,6),xPoint(3,6)]
   dVolTetThird = GetTetraVolume(yTetCoor) 
   dVol = dVol + dVolTetThird
   dCenterTetSmall = dCenterTetSmall + 0.25d0*(yTetcoor(:,1)+yTetcoor(:,2)+yTetcoor(:,3)+yTetcoor(:,4))*dVolTetThird

   yTetcoor(:,1) = [xPoint(1,4),xPoint(2,4),xPoint(3,4)]
   yTetcoor(:,2) = [xPoint(1,2),xPoint(2,2),xPoint(3,2)]
   yTetcoor(:,3) = [xPoint(1,3),xPoint(2,3),xPoint(3,3)]
   yTetcoor(:,4) = [xPoint(1,6),xPoint(2,6),xPoint(3,6)]
   dVolTetThird = GetTetraVolume(yTetCoor) 
   dVol = dVol + dVolTetThird
   dCenterTetSmall = dCenterTetSmall + 0.25d0*(yTetcoor(:,1)+yTetcoor(:,2)+yTetcoor(:,3)+yTetcoor(:,4))*dVolTetThird

   yTetcoor(:,1) = [xPoint(1,1),xPoint(2,1),xPoint(3,1)]
   yTetcoor(:,2) = [xPoint(1,2),xPoint(2,2),xPoint(3,2)]
   yTetcoor(:,3) = [xPoint(1,4),xPoint(2,4),xPoint(3,4)]
   yTetcoor(:,4) = [xPoint(1,6),xPoint(2,6),xPoint(3,6)]
   dVolTetThird = GetTetraVolume(yTetCoor) 
   dVol = dVol + dVolTetThird
   dCenterTetSmall = dCenterTetSmall + 0.25d0*(yTetcoor(:,1)+yTetcoor(:,2)+yTetcoor(:,3)+yTetcoor(:,4))*dVolTetThird

   IF (dS.gt.0d0) THEN
    dVolP  = dVolP  + dVol
    dVolM  = dVolM  + (dVolTet-dVol)
    dVolCM = dVolCM + (dCenterTetLarge-dCenterTetSmall)
    dVolElem(1) = dVolElem(1) + dVol
    dVolElem(2) = dVolElem(2) + (dVolTet-dVol)
   END IF
   IF (dS.lt.0d0) THEN
    dVolP  = dVolP  + (dVolTet-dVol)
    dVolM  = dVolM  + dVol
    dVolCM = dVolCM + dCenterTetSmall
    dVolElem(1) = dVolElem(1) + (dVolTet-dVol)
    dVolElem(2) = dVolElem(2) + dVol
   END IF

   CALL PointOrdering(xPoint(:,1),xPoint(:,2),xPoint(:,3),xTetcoor(:,iOmax),dist(iOmax),bTurn)
   IF (bTurn) THEN
!     WRITE(919,'(2i6,12E16.8)') 2,3,(xPoint(ii,1),xPoint(ii,3),xPoint(ii,2),ii=1,3)
    lInterface%X(:,iPoints+1) = xPoint(:,1)
    lInterface%X(:,iPoints+2) = xPoint(:,3)
    lInterface%X(:,iPoints+3) = xPoint(:,2)
   ELSE
!    WRITE(919,'(2i6,12E16.8)') 2,3,(xPoint(ii,1),xPoint(ii,2),xPoint(ii,3),ii=1,3)
    lInterface%X(:,iPoints+1) = xPoint(:,1)
    lInterface%X(:,iPoints+2) = xPoint(:,2)
    lInterface%X(:,iPoints+3) = xPoint(:,3)
   END IF
   CALL PointOrdering(xPoint(:,1),xPoint(:,4),xPoint(:,3),xTetcoor(:,iOmax),dist(iOmax),bTurn)
!    CALL PointOrdering(Point(:,1),Point(:,4),Point(:,3),coor(:,ivt(iOmax)),dist(ivt(iOmax)),bTurn)
   IF (bTurn) THEN
!     WRITE(919,'(2i6,12E16.8)') 2,3,(xPoint(ii,1),xPoint(ii,3),xPoint(ii,4),ii=1,3)
    lInterface%X(:,iPoints+4) = xPoint(:,1)
    lInterface%X(:,iPoints+5) = xPoint(:,3)
    lInterface%X(:,iPoints+6) = xPoint(:,4)
   ELSE
!     WRITE(919,'(2i6,12E16.8)') 2,3,(xPoint(ii,1),xPoint(ii,4),xPoint(ii,3),ii=1,3)
    lInterface%X(:,iPoints+4) = xPoint(:,1)
    lInterface%X(:,iPoints+5) = xPoint(:,4)
    lInterface%X(:,iPoints+6) = xPoint(:,3)
   END IF
  iPoints    = iPoints + 6
  iTriangles = iTriangles + 2

 END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  END IF

 ELSE

  IF (bCounterRun) THEN
  ELSE
   IF (dist(iOmax).LT.0d0) THEN
    dVolM  = dVolM  + dVolTet
    dVolCM = dVolCM + dCenterTetLarge
    dVolElem(2) = dVolElem(2) + dVolTet
   END IF
   IF (dist(iOmax).GT.0d0) THEN
    dVolP = dVolP + dVolTet
    dVolElem(1) = dVolElem(1) + dVolTet
   END IF
  END IF

 END IF


END DO

IF (iTriangles.NE.0) then
 iGroups = iGroups + 1
 IF (bCounterRun) THEN
 ELSE
  lInterface%L(iGroups) = iTriangles
 END IF
END IF
iiTriangles = iiTriangles + iTriangles

END SUBROUTINE ProcessHexa
! --------------------------------------------------------------------------------------------
REAL*8 FUNCTION GetHexaVolume(DC)
REAL*8 DC(3,8),DV,dX
REAL*8 X1,X2,X3,X4,X5,X6,X7,X8
REAL*8 Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8
REAL*8 Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8

X1=DC(1,1); Y1=DC(2,1); Z1=DC(3,1)
X2=DC(1,2); Y2=DC(2,2); Z2=DC(3,2)
X3=DC(1,3); Y3=DC(2,3); Z3=DC(3,3)
X4=DC(1,4); Y4=DC(2,4); Z4=DC(3,4)
X5=DC(1,5); Y5=DC(2,5); Z5=DC(3,5)
X6=DC(1,6); Y6=DC(2,6); Z6=DC(3,6)
X7=DC(1,7); Y7=DC(2,7); Z7=DC(3,7)
X8=DC(1,8); Y8=DC(2,8); Z8=DC(3,8)

dX = ((DABS((X4-X1)*(Y4-Y3)*(Z4-Z8)+(Y4-Y1)*(Z4-Z3)*(X4-X8)+(Z4-Z1)*(X4-X3)*(Y4-Y8)   - &
            (X4-X8)*(Y4-Y3)*(Z4-Z1)-(Y4-Y8)*(Z4-Z3)*(X4-X1)-(Z4-Z8)*(X4-X3)*(Y4-Y1))) + &
      (DABS((X2-X3)*(Y2-Y1)*(Z2-Z6)+(Y2-Y3)*(Z2-Z1)*(X2-X6)+(Z2-Z3)*(X2-X1)*(Y2-Y6)   - &
            (X2-X6)*(Y2-Y1)*(Z2-Z3)-(Y2-Y6)*(Z2-Z1)*(X2-X3)-(Z2-Z6)*(X2-X1)*(Y2-Y3))) + &
      (DABS((X5-X8)*(Y5-Y6)*(Z5-Z1)+(Y5-Y8)*(Z5-Z6)*(X5-X1)+(Z5-Z8)*(X5-X6)*(Y5-Y1)   - &
            (X5-X1)*(Y5-Y6)*(Z5-Z8)-(Y5-Y1)*(Z5-Z6)*(X5-X8)-(Z5-Z1)*(X5-X6)*(Y5-Y8))) + &
      (DABS((X7-X6)*(Y7-Y8)*(Z7-Z3)+(Y7-Y6)*(Z7-Z8)*(X7-X3)+(Z7-Z6)*(X7-X8)*(Y7-Y3)   - &
            (X7-X3)*(Y7-Y8)*(Z7-Z6)-(Y7-Y3)*(Z7-Z8)*(X7-X6)-(Z7-Z3)*(X7-X8)*(Y7-Y6))) + &
      (DABS((X1-X3)*(Y1-Y8)*(Z1-Z6)+(Y1-Y3)*(Z1-Z8)*(X1-X6)+(Z1-Z3)*(X1-X8)*(Y1-Y6)   - &
            (X1-X6)*(Y1-Y8)*(Z1-Z3)-(Y1-Y6)*(Z1-Z8)*(X1-X3)-(Z1-Z6)*(X1-X8)*(Y1-Y3))))
GetHexaVolume = dX/6d0

END FUNCTION GetHexaVolume
! --------------------------------------------------------------------------------------------
Subroutine GetHexaCenter(DVV,DCC)
REAL*8 DVV(3,27),DCC(3)
REAL*8 DCCi(3),dSubTetVol
INTEGER iSubTet,iVert

DCC = 0d0
DO iSubTet=1,40
 dSubTetVol = GetTetraVolume(DVV(:,tetra(:,iSubTet)))
 DCCi = 0d0
 DO iVert = 1,4
  DCCi = DCCi + 0.250*DVV(:,tetra(iVert,iSubTet))
 END DO
 DCC = DCC + DCCi*dSubTetVol
END DO

END Subroutine GetHexaCenter
! --------------------------------------------------------------------------------------------
REAL*8 FUNCTION GetTetraVolume(DC)
REAL*8 DC(3,4),DV,dX
REAL*8 X1,X2,X3,X4
REAL*8 Y1,Y2,Y3,Y4
REAL*8 Z1,Z2,Z3,Z4

X1=DC(1,1); Y1=DC(2,1); Z1=DC(3,1)
X2=DC(1,2); Y2=DC(2,2); Z2=DC(3,2)
X3=DC(1,3); Y3=DC(2,3); Z3=DC(3,3)
X4=DC(1,4); Y4=DC(2,4); Z4=DC(3,4)

dX = ((DABS((X2-X3)*(Y2-Y1)*(Z2-Z4)+(Y2-Y3)*(Z2-Z1)*(X2-X4)+(Z2-Z3)*(X2-X1)*(Y2-Y4)   - &
            (X2-X4)*(Y2-Y1)*(Z2-Z3)-(Y2-Y4)*(Z2-Z1)*(X2-X3)-(Z2-Z4)*(X2-X1)*(Y2-Y3))))

GetTetraVolume = dX/6d0

END FUNCTION GetTetraVolume
! --------------------------------------------------------------------------------------------
SUBROUTINE GetPoint(dV1,dV2,dP1,dP2,Point)
IMPLICIT NONE
REAL*8 dV1,dV2,dP1(3),dP2(3),Point(3)
REAL*8 daux,dL

 daux = abs(dV1)/abs(dV1-dV2)
 dL = SQRT((dP1(1)-dP2(1))**2d0 + (dP1(2)-dP2(2))**2d0 + (dP1(3)-dP2(3))**2d0)
 Point(:) = dP1(:) + daux*(dP2(:)-dP1(:))

!  WRITE(*,'(9ES12.4)') dP1,dP2,Point
END SUBROUTINE GetPoint
! --------------------------------------------------------------------------------------------

SUBROUTINE ElementJacobianSetUp(dEC)
REAL*8 dEC(3,8)
REAL*8 :: Q8=0.125D0

dJacComp(1,1)=( dEC(1,1)+dEC(1,2)+dEC(1,3)+dEC(1,4)+dEC(1,5)+dEC(1,6)+dEC(1,7)+dEC(1,8))*Q8
dJacComp(1,2)=( dEC(2,1)+dEC(2,2)+dEC(2,3)+dEC(2,4)+dEC(2,5)+dEC(2,6)+dEC(2,7)+dEC(2,8))*Q8
dJacComp(1,3)=( dEC(3,1)+dEC(3,2)+dEC(3,3)+dEC(3,4)+dEC(3,5)+dEC(3,6)+dEC(3,7)+dEC(3,8))*Q8
dJacComp(2,1)=(-dEC(1,1)+dEC(1,2)+dEC(1,3)-dEC(1,4)-dEC(1,5)+dEC(1,6)+dEC(1,7)-dEC(1,8))*Q8
dJacComp(2,2)=(-dEC(2,1)+dEC(2,2)+dEC(2,3)-dEC(2,4)-dEC(2,5)+dEC(2,6)+dEC(2,7)-dEC(2,8))*Q8
dJacComp(2,3)=(-dEC(3,1)+dEC(3,2)+dEC(3,3)-dEC(3,4)-dEC(3,5)+dEC(3,6)+dEC(3,7)-dEC(3,8))*Q8
dJacComp(3,1)=(-dEC(1,1)-dEC(1,2)+dEC(1,3)+dEC(1,4)-dEC(1,5)-dEC(1,6)+dEC(1,7)+dEC(1,8))*Q8
dJacComp(3,2)=(-dEC(2,1)-dEC(2,2)+dEC(2,3)+dEC(2,4)-dEC(2,5)-dEC(2,6)+dEC(2,7)+dEC(2,8))*Q8
dJacComp(3,3)=(-dEC(3,1)-dEC(3,2)+dEC(3,3)+dEC(3,4)-dEC(3,5)-dEC(3,6)+dEC(3,7)+dEC(3,8))*Q8
dJacComp(4,1)=(-dEC(1,1)-dEC(1,2)-dEC(1,3)-dEC(1,4)+dEC(1,5)+dEC(1,6)+dEC(1,7)+dEC(1,8))*Q8
dJacComp(4,2)=(-dEC(2,1)-dEC(2,2)-dEC(2,3)-dEC(2,4)+dEC(2,5)+dEC(2,6)+dEC(2,7)+dEC(2,8))*Q8
dJacComp(4,3)=(-dEC(3,1)-dEC(3,2)-dEC(3,3)-dEC(3,4)+dEC(3,5)+dEC(3,6)+dEC(3,7)+dEC(3,8))*Q8
dJacComp(5,1)=( dEC(1,1)-dEC(1,2)+dEC(1,3)-dEC(1,4)+dEC(1,5)-dEC(1,6)+dEC(1,7)-dEC(1,8))*Q8
dJacComp(5,2)=( dEC(2,1)-dEC(2,2)+dEC(2,3)-dEC(2,4)+dEC(2,5)-dEC(2,6)+dEC(2,7)-dEC(2,8))*Q8
dJacComp(5,3)=( dEC(3,1)-dEC(3,2)+dEC(3,3)-dEC(3,4)+dEC(3,5)-dEC(3,6)+dEC(3,7)-dEC(3,8))*Q8
dJacComp(6,1)=( dEC(1,1)-dEC(1,2)-dEC(1,3)+dEC(1,4)-dEC(1,5)+dEC(1,6)+dEC(1,7)-dEC(1,8))*Q8
dJacComp(6,2)=( dEC(2,1)-dEC(2,2)-dEC(2,3)+dEC(2,4)-dEC(2,5)+dEC(2,6)+dEC(2,7)-dEC(2,8))*Q8
dJacComp(6,3)=( dEC(3,1)-dEC(3,2)-dEC(3,3)+dEC(3,4)-dEC(3,5)+dEC(3,6)+dEC(3,7)-dEC(3,8))*Q8
dJacComp(7,1)=( dEC(1,1)+dEC(1,2)-dEC(1,3)-dEC(1,4)-dEC(1,5)-dEC(1,6)+dEC(1,7)+dEC(1,8))*Q8
dJacComp(7,2)=( dEC(2,1)+dEC(2,2)-dEC(2,3)-dEC(2,4)-dEC(2,5)-dEC(2,6)+dEC(2,7)+dEC(2,8))*Q8
dJacComp(7,3)=( dEC(3,1)+dEC(3,2)-dEC(3,3)-dEC(3,4)-dEC(3,5)-dEC(3,6)+dEC(3,7)+dEC(3,8))*Q8
dJacComp(8,1)=(-dEC(1,1)+dEC(1,2)-dEC(1,3)+dEC(1,4)+dEC(1,5)-dEC(1,6)+dEC(1,7)-dEC(1,8))*Q8
dJacComp(8,2)=(-dEC(2,1)+dEC(2,2)-dEC(2,3)+dEC(2,4)+dEC(2,5)-dEC(2,6)+dEC(2,7)-dEC(2,8))*Q8
dJacComp(8,3)=(-dEC(3,1)+dEC(3,2)-dEC(3,3)+dEC(3,4)+dEC(3,5)-dEC(3,6)+dEC(3,7)-dEC(3,8))*Q8

END SUBROUTINE ElementJacobianSetUp
! --------------------------------------------------------------------------------------------
SUBROUTINE SPoint_Q2_3D(lP,gP,nP)
IMPLICIT NONE
INTEGER iP,nP
REAL*8 lP(3,nP),gP(3,nP)
REAL*8 DJAC(3)

DO iP = 1,nP
 DJAC(1)=dJacComp(2,1)+dJacComp(5,1)*lP(2,iP)+dJacComp(6,1)*lP(3,iP)+dJacComp(8,1)*lP(2,iP)*lP(3,iP)
 DJAC(2)=dJacComp(3,2)+dJacComp(5,2)*lP(1,iP)+dJacComp(7,2)*lP(3,iP)+dJacComp(8,2)*lP(1,iP)*lP(3,iP)
 DJAC(3)=dJacComp(4,3)+dJacComp(6,3)*lP(1,iP)+dJacComp(7,3)*lP(2,iP)+dJacComp(8,3)*lP(1,iP)*lP(2,iP)
 gP(1,iP)=dJacComp(1,1)+DJAC(1)*lP(1,iP)+dJacComp(3,1)*lP(2,iP)+dJacComp(4,1)*lP(3,iP)+dJacComp(7,1)*lP(2,iP)*lP(3,iP)
 gP(2,iP)=dJacComp(1,2)+dJacComp(2,2)*lP(1,iP)+DJAC(2)*lP(2,iP)+dJacComp(4,2)*lP(3,iP)+dJacComp(6,2)*lP(1,iP)*lP(3,iP)
 gP(3,iP)=dJacComp(1,3)+dJacComp(2,3)*lP(1,iP)+dJacComp(3,3)*lP(2,iP)+DJAC(3)*lP(3,iP)+dJacComp(5,3)*lP(1,iP)*lP(2,iP)
END DO

END SUBROUTINE SPoint_Q2_3D
! --------------------------------------------------------------------------------------------
SUBROUTINE PointOrdering(P1,P2,P3,PQ,D,bTurn)
IMPLICIT NONE
LOGICAL bTurn
REAL*8 P1(3),P2(3),P3(3),PQ(3),D,R(3)
REAL*8 AX2,AY2,AZ2,AX3,AY3,AZ3,dC(4),daux

 AX2=P2(1) - P1(1)
 AY2=P2(2) - P1(2)
 AZ2=P2(3) - P1(3)
 AX3=P3(1) - P1(1)
 AY3=P3(2) - P1(2)
 AZ3=P3(3) - P1(3)

 dC(1)=(AY3*AZ2)-(AZ3*AY2)
 dC(2)=(AZ3*AX2)-(AX3*AZ2)
 dC(3)=(AX3*AY2)-(AY3*AX2)
 dC(4) = -(dC(1)*P1(1) + dC(2)*P1(2) + dC(3)*P1(3))

 daux = dC(4) + dC(1)*PQ(1) + dC(2)*PQ(2) + dC(3)*PQ(3)

 IF (daux*D.GT.0d0) THEN
  bTurn =.FALSE. 
 ELSE
  bTurn =.TRUE. 
 END IF

END SUBROUTINE PointOrdering
! --------------------------------------------------------------------------------------------
END SUBROUTINE TetraInterfacePoints
! --------------------------------------------------------------------------------------------
DOUBLE PRECISION FUNCTION fVal_Q2_2D(x1,x2,dPhi)
REAL*8 DH(9),dPhi(9),x1,x2
REAL*8 Q2,Q4
INTEGER iPP
PARAMETER (Q4=.25D0,Q2=.5D0)

DH(1)= Q4*(1D0-X1)*(1D0-X2)*X1*X2
DH(2)=-Q4*(1D0+X1)*(1D0-X2)*X1*X2
DH(3)= Q4*(1D0+X1)*(1D0+X2)*X1*X2
DH(4)=-Q4*(1D0-X1)*(1D0+X2)*X1*X2
DH(5)=-Q2*(1D0-X1*X1)*(1D0-X2)*X2
DH(6)= Q2*(1D0+X1)*(1D0-X2*X2)*X1
DH(7)= Q2*(1D0-X1*X1)*(1D0+X2)*X2
DH(8)=-Q2*(1D0-X1)*(1D0-X2*X2)*X1
DH(9)= (1D0-X1*X1)*(1D0-X2*X2)

fVal_Q2_2D = 0d0
DO ipp=1,9
 fVal_Q2_2D = fVal_Q2_2D + dPhi(ipp) * DH(ipp)
END DO

RETURN

END FUNCTION fVal_Q2_2D
! --------------------------------------------------------------------------------------------
DOUBLE PRECISION FUNCTION fVal_Q2_3D(x1,x2,x3,dPhi)
IMPLICIT NONE
REAL*8 DH(27),dPhi(27),x1,x2,x3
REAL*8 Q8
INTEGER iPP
PARAMETER (Q8=.125D0)

DH(1) =-Q8*X1*(1D0-X1)*X2*(1D0-X2)*X3*(1D0-X3)
DH(2) = Q8*X1*(1D0+X1)*X2*(1D0-X2)*X3*(1D0-X3)
DH(3) =-Q8*X1*(1D0+X1)*X2*(1D0+X2)*X3*(1D0-X3)
DH(4) = Q8*X1*(1D0-X1)*X2*(1D0+X2)*X3*(1D0-X3)
DH(5) = Q8*X1*(1D0-X1)*X2*(1D0-X2)*X3*(1D0+X3)
DH(6) =-Q8*X1*(1D0+X1)*X2*(1D0-X2)*X3*(1D0+X3)
DH(7) = Q8*X1*(1D0+X1)*X2*(1D0+X2)*X3*(1D0+X3)
DH(8) =-Q8*X1*(1D0-X1)*X2*(1D0+X2)*X3*(1D0+X3)
DH(9) = Q8*(2D0-2D0*X1**2D0)*X2*(1D0-X2)*X3*(1D0-X3)
DH(10)=-Q8*X1*(1D0+X1)*(2D0-2D0*X2**2D0)*X3*(1D0-X3)
DH(11)=-Q8*(2D0-2D0*X1**2D0)*X2*(1D0+X2)*X3*(1D0-X3)
DH(12)= Q8*X1*(1D0-X1)*(2D0-2D0*X2**2D0)*X3*(1D0-X3)
DH(13)= Q8*X1*(1D0-X1)*X2*(1D0-X2)*(2D0-2D0*X3**2D0)
DH(14)=-Q8*X1*(1D0+X1)*X2*(1D0-X2)*(2D0-2D0*X3**2D0)
DH(15)= Q8*X1*(1D0+X1)*X2*(1D0+X2)*(2D0-2D0*X3**2D0)
DH(16)=-Q8*X1*(1D0-X1)*X2*(1D0+X2)*(2D0-2D0*X3**2D0)
DH(17)=-Q8*(2D0-2D0*X1**2D0)*X2*(1D0-X2)*X3*(1D0+X3)
DH(18)= Q8*X1*(1D0+X1)*(2D0-2D0*X2**2D0)*X3*(1D0+X3)
DH(19)= Q8*(2D0-2D0*X1**2D0)*X2*(1D0+X2)*X3*(1D0+X3)
DH(20)=-Q8*X1*(1D0-X1)*(2D0-2D0*X2**2D0)*X3*(1D0+X3)
DH(21)=-Q8*(2D0-2D0*X1**2D0)*(2D0-2D0*X2**2D0)*X3*(1D0-X3)
DH(22)=-Q8*(2D0-2D0*X1**2D0)*X2*(1D0-X2)*(2D0-2D0*X3**2D0)
DH(23)= Q8*X1*(1D0+X1)*(2D0-2D0*X2**2D0)*(2D0-2D0*X3**2D0)
DH(24)= Q8*(2D0-2D0*X1**2D0)*X2*(1D0+X2)*(2D0-2D0*X3**2D0)
DH(25)=-Q8*X1*(1D0-X1)*(2D0-2D0*X2**2D0)*(2D0-2D0*X3**2D0)
DH(26)= Q8*(2D0-2D0*X1**2D0)*(2D0-2D0*X2**2D0)*X3*(1D0+X3)
DH(27)= Q8*(2D0-2D0*X1**2D0)*(2D0-2D0*X2**2D0)*(2D0-2D0*X3**2D0)

fVal_Q2_3D = 0d0
DO ipp=1,27
 fVal_Q2_3D = fVal_Q2_3D + dPhi(ipp) * DH(ipp)
END DO

RETURN

END FUNCTION fVal_Q2_3D
!
! --------------------------------------------------------------------------------------------
!
SUBROUTINE GetQ1Distance(DG,dLS,DMM,NU,&
           KVERT,KAREA,KEDGE,DCORVG,ICUB,ELE,dEps)
USE PP3D_MPI, ONLY:myid
IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
CHARACTER SUB*6,FMT*15,CPARAM*120

PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,&
           NNDIM=3,NNCOF=10)
PARAMETER (Q2=0.5D0,Q8=0.125D0)

DIMENSION DG(*),DMM(*),dLS(*)
DIMENSION DCORVG(NNDIM,*)
DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)

DIMENSION KDFG(NNBAS),KDFL(NNBAS)
REAL*8 pTet(4,3),pMat(4,4),dP(3),dLSBas(4)
REAL*8 dTetra(4,3),dCoeff(4,4),dMid(3),dBasX(4),dBasY(4),dBasZ(4)
INTEGEr iGroups(10)

LOGICAL bInside,bPositive,bNegative

COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
COMMON /ERRCTL/ IER,ICHECK
COMMON /CHAR/   SUB,FMT(3),CPARAM
COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,&
                DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),&
                IEL,NDIM
COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,&
                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
COMMON /COAUX1/ KDFG,KDFL,IDFL

!user COMMON blocks
INTEGER  VIPARM 
DIMENSION VIPARM(100)
EQUIVALENCE (IAUSAV,VIPARM)
COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,&
               IMASS,IMASSL,IUPW,IPRECA,IPRECB,&
               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,&
               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,&
               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,&
               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA

SAVE

DO I= 1,NNDER
 BDER(I)=.FALSE.
end do

BDER(1:4)=.TRUE.

IELTYP=-1
CALL ELE(0D0,0D0,0D0,IELTYP)
IDFL=NDFL(IELTYP)

CALL CB3H(ICUB)
IF (IER.NE.0) RETURN

!******************************************************************
!Calculation of the matrix - storage technique 7 or 8
!******************************************************************
ICUBP=ICUB
CALL ELE(0D0,0D0,0D0,-2)

!Loop over all elements
DO IEL=1,NEL

CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
IF (IER.LT.0) RETURN

!Evaluation of coordinates of the vertices
dP = 0d0
DO IVE=1,NVE
JP=KVERT(IVE,IEL)
KVE(IVE)=JP
DX(IVE)=DCORVG(1,JP)
DY(IVE)=DCORVG(2,JP)
DZ(IVE)=DCORVG(3,JP)
dP = dP + [DX(IVE),DY(IVE),DZ(IVE)]
END DO
dP = 0.125d0*dP

bPositive=.TRUE.
bNegative=.TRUE.
DO JDFL=1,IDFL
 IG=KDFG(JDFL)
 IF (dLS(IG).LE.0d0) bPositive = .FALSE.
 IF (dLS(IG).GT.0d0) bNegative = .FALSE.
END DO

CALL GetClosestGroups(dP,iGroups,dLSEstimate)

DJ11=( DX(1)+DX(2)+DX(3)+DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
DJ12=( DY(1)+DY(2)+DY(3)+DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
DJ13=( DZ(1)+DZ(2)+DZ(3)+DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
DJ21=(-DX(1)+DX(2)+DX(3)-DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
DJ22=(-DY(1)+DY(2)+DY(3)-DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
DJ23=(-DZ(1)+DZ(2)+DZ(3)-DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
DJ31=(-DX(1)-DX(2)+DX(3)+DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
DJ32=(-DY(1)-DY(2)+DY(3)+DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
DJ33=(-DZ(1)-DZ(2)+DZ(3)+DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
DJ41=(-DX(1)-DX(2)-DX(3)-DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
DJ42=(-DY(1)-DY(2)-DY(3)-DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
DJ43=(-DZ(1)-DZ(2)-DZ(3)-DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
DJ51=( DX(1)-DX(2)+DX(3)-DX(4)+DX(5)-DX(6)+DX(7)-DX(8))*Q8
DJ52=( DY(1)-DY(2)+DY(3)-DY(4)+DY(5)-DY(6)+DY(7)-DY(8))*Q8
DJ53=( DZ(1)-DZ(2)+DZ(3)-DZ(4)+DZ(5)-DZ(6)+DZ(7)-DZ(8))*Q8
DJ61=( DX(1)-DX(2)-DX(3)+DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
DJ62=( DY(1)-DY(2)-DY(3)+DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
DJ63=( DZ(1)-DZ(2)-DZ(3)+DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
DJ71=( DX(1)+DX(2)-DX(3)-DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
DJ72=( DY(1)+DY(2)-DY(3)-DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
DJ73=( DZ(1)+DZ(2)-DZ(3)-DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
DJ81=(-DX(1)+DX(2)-DX(3)+DX(4)+DX(5)-DX(6)+DX(7)-DX(8))*Q8
DJ82=(-DY(1)+DY(2)-DY(3)+DY(4)+DY(5)-DY(6)+DY(7)-DY(8))*Q8
DJ83=(-DZ(1)+DZ(2)-DZ(3)+DZ(4)+DZ(5)-DZ(6)+DZ(7)-DZ(8))*Q8

!Loop over all cubature points
DO ICUBP=1,NCUBP

XI1=DXI(ICUBP,1)
XI2=DXI(ICUBP,2)
XI3=DXI(ICUBP,3)

!Jacobian of the bilinear mapping onto the reference element
DJAC(1,1)=DJ21+DJ51*XI2+DJ61*XI3+DJ81*XI2*XI3
DJAC(1,2)=DJ31+DJ51*XI1+DJ71*XI3+DJ81*XI1*XI3
DJAC(1,3)=DJ41+DJ61*XI1+DJ71*XI2+DJ81*XI1*XI2
DJAC(2,1)=DJ22+DJ52*XI2+DJ62*XI3+DJ82*XI2*XI3
DJAC(2,2)=DJ32+DJ52*XI1+DJ72*XI3+DJ82*XI1*XI3
DJAC(2,3)=DJ42+DJ62*XI1+DJ72*XI2+DJ82*XI1*XI2
DJAC(3,1)=DJ23+DJ53*XI2+DJ63*XI3+DJ83*XI2*XI3
DJAC(3,2)=DJ33+DJ53*XI1+DJ73*XI3+DJ83*XI1*XI3
DJAC(3,3)=DJ43+DJ63*XI1+DJ73*XI2+DJ83*XI1*XI2
DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3)) &
     -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3)) &
     +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
OM=DOMEGA(ICUBP)*ABS(DETJ)

XX=DJ11+DJAC(1,1)*XI1+DJ31*XI2+DJ41*XI3+DJ71*XI2*XI3
YY=DJ12+DJ22*XI1+DJAC(2,2)*XI2+DJ42*XI3+DJ62*XI1*XI3
ZZ=DJ13+DJ23*XI1+DJ33*XI2+DJAC(3,3)*XI3+DJ53*XI1*XI2

dP = [XX,YY,ZZ]
CALL GetDistanceToInterfaceInGivenGroups(dP,iGroups,dLSdist)

CALL ELE(XI1,XI2,XI3,-3)
IF (IER.LT.0) RETURN

dLevSet = 0d0
DO JDFL=1,IDFL
  IG=KDFG(JDFL)
  IL=KDFL(JDFL)
  HBAS=DBAS(1,IL,1)
  dLevSet = dLevSet + dLS(IG)*HBAS
END DO

dLSsign = dLevSet/ABS(dLevSet)
dLSdist = abs(dLSdist)*dLSsign

!Summing up over all pairs of multiindices
DO JDFL=1,IDFL
  IG=KDFG(JDFL)
  IL=KDFL(JDFL)
  HBAS=DBAS(1,IL,1)
  DG(IG) = DG(IG) + OM*dLSdist*HBAS
  DMM(IG) = DMM(IG) + OM*HBAS
END DO

END DO

END DO

END
