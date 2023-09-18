!========================================================================================
!                           Sub: Move_Particle
!========================================================================================
subroutine Move_xParticle_Fast(dcorvg,kvert,kedge,karea,kel_LdA,kel_ColA,&
                         Vel0U,Vel0V,Vel0W,Vel1U,Vel1V,Vel1W,SHELL,nvt,net,nat,nel,tDelta,tStart)
USE PP3D_MPI, ONLY : myid
USE types, ONLY : myActiveSet,myLostSet,nActiveSet,nLostSet
USE OctTreeSearch
USE xPart_def

! XXXXXXXXXXXX
IMPLICIT NONE
REAL*8 dcorvg(3,*),point(3),tLevel,tDelta,tStart
REAL*8 Vel0U(*),Vel0V(*),Vel0W(*),Vel1U(*),Vel1V(*),Vel1W(*),SHELL(*)
INTEGER nvt,net,nat,nel
INTEGER kel_LdA(*),kel_ColA(*)
INTEGER kvert(8,*),kedge(12,*),karea(6,*)
INTEGER KDFG(27),KDFL(27)
INTEGER i,iPoint,jPoint,ivt,jel,iel,iFoundElem
REAL*8 dist,pointR(3),LocVel0(3,27),LocVel1(3,27)
LOGICAL :: bFound,bOut=.FALSE.

INTEGER :: nBuffElem
INTEGER nXX,kk,iMon
!!!!!!!!!!!!!!!!!!!!!!!!!!!! USER params !!!!!!!!!!!!!!!!!!!!!!!!
PARAMETER(nBuffElem=100)
PARAMETER (nXX = 12)
!!!!!!!!!!!!!!!!!!!!!!!!!!!! USER params !!!!!!!!!!!!!!!!!!!!!!!!

INTEGER :: iIter, iLoc,iParticel,iLostParticel,iActiveParticel,idynType,indice
REAL*8 cpx,cpy,cpz,cnormal(3)

! LOGICAL :: bExit
INTEGER iBuffElem,jBuffElem,BuffElem(nBuffElem),maxBufferSize,maxBufferElem

INTEGER iMonitor(nXX),iAux
REAL*8  dMonitor(nXX),dAux

iLostParticel = 0
iActiveParticel = 0
maxBufferSize = 0
maxBufferElem = 0

DO iParticel = 1,nActiveSet

iIter = 0

point  = myActiveSet(iParticel)%coor
tLevel = myActiveSet(iParticel)%time
dParticleVelo = myActiveSet(iParticel)%velo

55 CONTINUE

bFound = .FALSE.

CALL FindRankingInOctTree(dcorvg,nvt,point,iMonitor,dMonitor,nxx)
call CreateElemBuffer()

! if (iX.eq.1) then
!  write(*,*) 'fsdfdfs',iParticel
! !  pause
! end if
!  bExit = .false.

 DO jBuffElem = 1,iBuffElem
  jel = BuffElem(jBuffElem)
  CALL TraceParticleInTheRightElement()
  maxBufferElem = Max(maxBufferElem,jBuffElem)
  if (bFound) goto 1
 END DO

IF (.not.bFound) THEN
 iLostParticel = iLostParticel + 1
 indice = myActiveSet(iParticel)%indice
 
 myLostSet(indice)%coor   = point
 myLostSet(indice)%time   = tLevel
 myLostSet(indice)%indice = myActiveSet(iParticel)%indice
 myLostSet(indice)%id     = myid
 myLostSet(indice)%velo   = dParticleVelo
 GOTO 2
END IF

1 CONTINUE

iIter = iIter + 1
IF (tLevel.lt.tDelta.AND.iIter.LT.nIter) GOTO 55

IF (bFound.AND.bOut) WRITE(*,'(A,I0,3E14.4,A)') '------------   point ', myActiveSet(iParticel)%indice, point, '  ------------'
IF (bFound.AND.bOut) WRITE(*,'(A,I0,I0,3E14.4,A)') '-------  point ', myActiveSet(iParticel)%indice,iIter,point, '  has been successfully treated -----------'

iActiveParticel = iActiveParticel + 1
myActiveSet(iActiveParticel)%coor   = point
myActiveSet(iActiveParticel)%velo   = dParticleVelo
myActiveSet(iActiveParticel)%time   = tLevel
myActiveSet(iActiveParticel)%indice = myActiveSet(iParticel)%indice
myActiveSet(iActiveParticel)%id     = myActiveSet(iParticel)%id

2 CONTINUE

! pause
END DO

! IF (myid.eq.1) WRITE(*,'(A,2(I0,("/")),I0)') 'BufferManagement, (maxsize,usedsize, maxelem) : ', nBuffElem,maxBufferSize,maxBufferElem

nActiveSet      = iActiveParticel

IF (bOut) WRITE(*,*) 'nActiveSet,nLostSet: ',nActiveSet,iLostParticel

 contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 subroutine GetNeighborsOfFoundElem()
 real*8 ddd
 
 DO i=1,8
  ivt = kvert(i,iFoundElem)
  ddd = (dcorvg(1,ivt)-point(1))**2d0 + (dcorvg(2,ivt)-point(2))**2d0 + (dcorvg(3,ivt)-point(3))**2d0
  dMonitor(i) = ddd
  iMonitor(i) = ivt
 END DO
 
 DO kk=8-1,1,-1
  IF (dMonitor(kk+1).lt.dMonitor(kk)) THEN
   iAux = iMonitor(kk)
   dAux = dMonitor(kk)
   iMonitor(kk) = iMonitor(kk+1)
   dMonitor(kk) = dMonitor(kk+1)
   iMonitor(kk+1) = iAux
   dMonitor(kk+1) = dAux
  END IF
 END DO
 
 iBuffElem = 0
 DO i=1,8
  iPoint = iMonitor(i)
  DO iel = kel_LdA(iPoint),kel_LdA(iPoint+1)-1
    jel = kel_ColA(iel)
    IF (jel.ne.0) then
     DO jBuffElem=1,iBuffElem
      if (jel.eq.BuffElem(jBuffElem)) GOTO 1
     END DO
    else
      GOTO 1
    end if
    
    iBuffElem = iBuffElem + 1
    IF (iBuffElem.eq.nBuffElem) THEN
     write(*,*) 'buffer is too small ... ', myid
     pause
    end if
    BuffElem(iBuffElem) = jel
    
1   CONTINUE    

  END DO
 END DO
 
 maxBufferSize = Max(maxBufferSize,iBuffElem)
 
 end subroutine GetNeighborsOfFoundElem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 subroutine CreateElemBuffer()
 
 iBuffElem = 0
 DO iMon = 1,nXX
  iPoint = iMonitor(iMon)
  if (iPoint.ge.1.and.iPoint.le.nvt) then
   DO iel = kel_LdA(iPoint),kel_LdA(iPoint+1)-1
    jel = kel_ColA(iel)
    IF (jel.ne.0) then
     DO jBuffElem=1,iBuffElem
      if (jel.eq.BuffElem(jBuffElem)) GOTO 1
     END DO
    else
      GOTO 1
    end if
    
    iBuffElem = iBuffElem + 1
    IF (iBuffElem.eq.nBuffElem) THEN
     write(*,*) 'buffer is too small ... ', myid
     pause
    end if
    BuffElem(iBuffElem) = jel
1   CONTINUE    
    
   end do
  end if
 end do
  
 maxBufferSize = Max(maxBufferSize,iBuffElem)
 
 end subroutine CreateElemBuffer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 subroutine TraceParticleInTheRightElement()

    DO i=1,8
     ivt = kvert(i,jel)
     P8(:,i) = dcorvg(:,ivt)
    END DO
   
    dist = 1d30
    jPoint = 0

    DO i=1,8
     daux = sqrt((P8(1,i)-point(1))**2d0 + (P8(2,i)-point(2))**2d0 + (P8(3,i)-point(3))**2d0)
     IF (daux.lt.dist) THEN
      dist = daux
      jPoint = i
     END IF
    END DO

    CALL GetPointFromElement(point,jPoint,bFound,pointR,iParticel)
    IF (Point(3).gt.myParticleParam%OutflowZPos) bFound=.false.

    IF (bFound) THEN
    
     CALL NDFGL(JEL,1,13,KVERT,KEDGE,KAREA,KDFG,KDFL)
     DO i=1,27
      ivt  = KDFG(i)
      iLoc = KDFL(i)
      LocVel0(:,iLoc) = [Vel0U(ivt),Vel0V(ivt),Vel0W(ivt)]
      LocVel1(:,iLoc) = [Vel1U(ivt),Vel1V(ivt),Vel1W(ivt)]
      shellE(iLoc)    = shell(ivt)
     END DO

     IF (bOut) WRITE(*,'(A,I0,3E14.4,A)') '------------   point ', myActiveSet(iParticel)%indice, point, '  ------------'

     DV0 = LocVel0
     DV1 = LocVel1
     P   = point 
     PR  = pointR
     
     CALL MovePointThroughElement(tLevel,tDelta,tStart,iParticel,jel)
     
     point = P
     
     CALL RETURN_Distance()
     
     ! scaling [mm] to [cm]
     distance = 1d1*distance
    
     if (distance.lt. d_CorrDist*0.5d0) then
     
       ! d_CorrDist is in [mm], distance as well is in [mm] / at the end we have to come back to [cm]
       
      cdx = 1d1*P(1) 
      cdy = 1d1*P(2)
      cdz = 1d1*P(3)
      
      CALL RETURN_Normal()
      daux = SQRT(dnormal(1)**2d0 + dnormal(2)**2d0 + dnormal(3)**2d0)
      dnormal = dnormal/daux
      
      cpx = cdx + dnormal(1)*(-distance + d_CorrDist*0.5d0)
      cpy = cdy + dnormal(2)*(-distance + d_CorrDist*0.5d0)
      cpz = cdz + dnormal(3)*(-distance + d_CorrDist*0.5d0)

      P=0.1d0*[cpx,cpy,cpz]
      
      point = [P(1),P(2),P(3)]

     end if
    
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !     CALL GetDistAndProjPToAllSTLs(cdx,cdy,cdz,cpx,cpy,cpz,dist_CGAL)
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !      if (dist_CGAL.lt. d_CorrDist*0.5d0) then
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       cnormal = [cpx-cdx,cpy-cdy,cpz-cdz]
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       daux = SQRT(cnormal(1)**2d0 + cnormal(2)**2d0 + cnormal(3)**2d0)
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       if (dist_CGAL.gt.0d0) then
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !        cnormal = -cnormal/daux
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       else
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !        cnormal = cnormal/daux
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       end if
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       cdx = cpx + cnormal(1)*d_CorrDist*0.5d0
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       cdy = cpy + cnormal(2)*d_CorrDist*0.5d0
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       cdz = cpz + cnormal(3)*d_CorrDist*0.5d0
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       P=0.1d0*[cdx,cdy,cdz]
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       point = [P(1),P(2),P(3)]
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !      end if
     
     iFoundElem = jel
   
!      bExit = .true.
     RETURN
    END IF

 END SUBROUTINE TraceParticleInTheRightElement

END SUBROUTINE Move_xParticle_Fast
