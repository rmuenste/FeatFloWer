!========================================================================================
!                           Sub: Move_Particle
!========================================================================================
subroutine Move_xParticle_Fast(dcorvg,kvert,kedge,karea,kel_LdA,kel_ColA,&
                         Vel0U,Vel0V,Vel0W,Vel1U,Vel1V,Vel1W,SHELL,nvt,net,nat,nel,&
                         tDelta,tStart,iMinParticleID,iMaxParticleID,iActiveParticel)
USE PP3D_MPI, ONLY : myid
USE types, ONLY : myActiveSet,myLostSet,nActiveSet,nLostSet
USE OctTreeSearch, ONLY : FindRankingInOctTreeOMP
USE xPart_def
use omp_lib
use fbm, only : myFBM

! XXXXXXXXXXXX
IMPLICIT NONE
REAL*8 dcorvg(3,*),tDelta,tStart
REAL*8 Vel0U(*),Vel0V(*),Vel0W(*),Vel1U(*),Vel1V(*),Vel1W(*),SHELL(*)
INTEGER nvt,net,nat,nel
integer iMinParticleID,iMaxParticleID,iActiveParticel
INTEGER kel_LdA(*),kel_ColA(*)
INTEGER kvert(8,*),kedge(12,*),karea(6,*)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

INTEGER KDFG(27),KDFL(27)
INTEGER i,iPoint,jPoint,ivt,jel,iel,iFoundElem
REAL*8 dist,point(3),pointR(3),LocVel0(3,27),LocVel1(3,27),tLevel
LOGICAL :: bFound,bOut=.FALSE.

INTEGER :: nBuffElem
INTEGER nXX
INTEGER kk,iMon
!!!!!!!!!!!!!!!!!!!!!!!!!!!! USER params !!!!!!!!!!!!!!!!!!!!!!!!
PARAMETER(nBuffElem=100)
PARAMETER (nXX = 12)
!!!!!!!!!!!!!!!!!!!!!!!!!!!! USER params !!!!!!!!!!!!!!!!!!!!!!!!

INTEGER :: iIter, iLoc,iParticel,iLostParticel,idynType,indice
REAL*8 cpx,cpy,cpz,cnormal(3)
REAL*8 local_P8(3,8)

! LOGICAL :: bExit
INTEGER iBuffElem,jBuffElem,BuffElem(nBuffElem),maxBufferSize,maxBufferElem

INTEGER iMonitor(nXX),iAux
REAL*8  dMonitor(nXX),dAux
INTEGER :: iLastElem,iExitKind,iExitFaceLocal,iExitVertexLocal
INTEGER :: iExitEdgeVerts(2),iExitFaceVerts(4)
LOGICAL :: bHaveLocalHint,bTriedOctree

iLostParticel = 0
iActiveParticel = 0
maxBufferSize = 0
maxBufferElem = 0

DO iParticel = iMinParticleID,iMaxParticleID

iIter = 0
iLastElem = 0
iExitKind = 0
iExitFaceLocal = 0
iExitVertexLocal = 0
iExitEdgeVerts = 0
iExitFaceVerts = 0
bHaveLocalHint = .FALSE.

if (myActiveSet(iParticel)%id.eq.0) then

point  = myActiveSet(iParticel)%coor
tLevel = myActiveSet(iParticel)%time
dParticleVelo = myActiveSet(iParticel)%velo

55 CONTINUE

bFound = .FALSE.
bTriedOctree = .FALSE.

if (bHaveLocalHint) then
 call CreateLocalElemBuffer()
else
 call FindRankingInOctTreeOMP(dcorvg,nvt,point,iMonitor,dMonitor,nxx)
 call CreateElemBuffer()
 bTriedOctree = .TRUE.
end if

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

IF ((.not.bFound).and.(.not.bTriedOctree)) THEN
 CALL FindRankingInOctTreeOMP(dcorvg,nvt,point,iMonitor,dMonitor,nxx)
 call CreateElemBuffer()
 bTriedOctree = .TRUE.

 DO jBuffElem = 1,iBuffElem
  jel = BuffElem(jBuffElem)
  CALL TraceParticleInTheRightElement()
  maxBufferElem = Max(maxBufferElem,jBuffElem)
  if (bFound) goto 1
 END DO
END IF

IF (.not.bFound) THEN
 iLostParticel = iLostParticel + 1
 indice = myActiveSet(iParticel)%indice
 
 myLostSet(indice)%coor   = point
 myLostSet(indice)%time   = tLevel
 myLostSet(indice)%velo   = dParticleVelo
 myLostSet(indice)%indice = indice
 myLostSet(indice)%id     = myid

 myActiveSet(iParticel)%id   = -1
 bHaveLocalHint = .FALSE.
 
 GOTO 2
END IF

1 CONTINUE

iLastElem = iFoundElem
bHaveLocalHint = .FALSE.
if (tLevel.lt.tDelta) then
 call ClassifyElementExit(iExitKind,iExitFaceLocal,iExitVertexLocal,iExitEdgeVerts,iExitFaceVerts)
 if (iExitKind.gt.0) bHaveLocalHint = .TRUE.
end if

iIter = iIter + 1
IF (tLevel.lt.tDelta.AND.iIter.LT.nIter) GOTO 55

IF (bFound.AND.bOut) WRITE(*,'(A,I0,3E14.4,A)') '------------   point ', myActiveSet(iParticel)%indice, point, '  ------------'
IF (bFound.AND.bOut) WRITE(*,'(A,I0,I0,3E14.4,A)') '-------  point ', myActiveSet(iParticel)%indice,iIter,point, '  has been successfully treated -----------'

iActiveParticel = iActiveParticel + 1
myActiveSet(iParticel)%coor   = point
myActiveSet(iParticel)%velo   = dParticleVelo
myActiveSet(iParticel)%time   = tLevel

2 CONTINUE

end if

END DO

! IF (myid.eq.1) WRITE(*,'(A,2(I0,("/")),I0)') 'BufferManagement, (maxsize,usedsize, maxelem) : ', nBuffElem,maxBufferSize,maxBufferElem

! nActiveSet      = iActiveParticel

! IF (bOut) WRITE(*,*) 'nActiveSet,nLostSet: ',nActiveSet,iLostParticel

 contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 subroutine GetNeighborsOfFoundElem()

 DO i=1,8
  ivt = kvert(i,iFoundElem)
  dMonitor(i) = (dcorvg(1,ivt)-point(1))**2d0 + (dcorvg(2,ivt)-point(2))**2d0 + (dcorvg(3,ivt)-point(3))**2d0
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
 subroutine AddElementToBuffer(iElem)
 integer, intent(in) :: iElem

 if (iElem.le.0.or.iElem.gt.nel) return

 do jBuffElem=1,iBuffElem
  if (iElem.eq.BuffElem(jBuffElem)) return
 end do

 iBuffElem = iBuffElem + 1
 IF (iBuffElem.gt.nBuffElem) THEN
  write(*,*) 'buffer is too small ... ', myid
  pause
 end if
 BuffElem(iBuffElem) = iElem

 end subroutine AddElementToBuffer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 logical function ElementHasVertex(iElem,iVertex)
 integer, intent(in) :: iElem,iVertex

 ElementHasVertex = .FALSE.
 do i=1,8
  if (kvert(i,iElem).eq.iVertex) then
   ElementHasVertex = .TRUE.
   return
  end if
 end do

 end function ElementHasVertex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine GetFaceLocalVertices(iFaceLocal,faceVerts)
 integer, intent(in)  :: iFaceLocal
 integer, intent(out) :: faceVerts(4)

 select case (iFaceLocal)
 case (1)
  faceVerts = (/1,4,8,5/)
 case (2)
  faceVerts = (/2,3,7,6/)
 case (3)
  faceVerts = (/1,2,6,5/)
 case (4)
  faceVerts = (/4,3,7,8/)
 case (5)
  faceVerts = (/1,2,3,4/)
 case (6)
  faceVerts = (/5,6,7,8/)
 case default
  faceVerts = 0
 end select

 end subroutine GetFaceLocalVertices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine AddIncidentElementsOfVertex(localVertex)
 integer, intent(in) :: localVertex

 if (localVertex.le.0.or.localVertex.gt.8) return

 iPoint = kvert(localVertex,iLastElem)
 do iel = kel_LdA(iPoint),kel_LdA(iPoint+1)-1
  jel = kel_ColA(iel)
  if (jel.ne.0) call AddElementToBuffer(jel)
 end do

 end subroutine AddIncidentElementsOfVertex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine AddIncidentElementsOfEdge(localEdgeVerts)
 integer, intent(in) :: localEdgeVerts(2)
 integer :: v1,v2

 if (minval(localEdgeVerts).le.0) return

 v1 = kvert(localEdgeVerts(1),iLastElem)
 v2 = kvert(localEdgeVerts(2),iLastElem)
 do iel = kel_LdA(v1),kel_LdA(v1+1)-1
  jel = kel_ColA(iel)
  if (jel.ne.0) then
   if (ElementHasVertex(jel,v2)) call AddElementToBuffer(jel)
  end if
 end do

 end subroutine AddIncidentElementsOfEdge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine AddIncidentElementsOfFaceVertices(faceVerts)
 integer, intent(in) :: faceVerts(4)

 do kk=1,4
  call AddIncidentElementsOfVertex(faceVerts(kk))
 end do

 end subroutine AddIncidentElementsOfFaceVertices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine AddNeighborAcrossFace(iFaceLocal)
 integer, intent(in) :: iFaceLocal
 integer :: faceVerts(4),sharedArea,globalVerts(4)
 logical :: bAllVerts

 if (iFaceLocal.le.0.or.iFaceLocal.gt.6) return

 sharedArea = karea(iFaceLocal,iLastElem)
 call GetFaceLocalVertices(iFaceLocal,faceVerts)
 do kk=1,4
  globalVerts(kk) = kvert(faceVerts(kk),iLastElem)
 end do

 do iel = kel_LdA(globalVerts(1)),kel_LdA(globalVerts(1)+1)-1
  jel = kel_ColA(iel)
  if (jel.eq.0.or.jel.eq.iLastElem) cycle
  if (.not.ElementHasVertex(jel,globalVerts(2))) cycle
  if (.not.ElementHasVertex(jel,globalVerts(3))) cycle
  if (.not.ElementHasVertex(jel,globalVerts(4))) cycle
  bAllVerts = .FALSE.
  do i=1,6
   if (karea(i,jel).eq.sharedArea) then
    bAllVerts = .TRUE.
    exit
   end if
  end do
  if (bAllVerts) then
   call AddElementToBuffer(jel)
   return
  end if
 end do

 end subroutine AddNeighborAcrossFace
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine CreateLocalElemBuffer()
 
 iBuffElem = 0

 if (iLastElem.le.0) return

 select case (iExitKind)
 case (1)
  call AddNeighborAcrossFace(iExitFaceLocal)
 case (2)
  call AddIncidentElementsOfEdge(iExitEdgeVerts)
 case (3)
  call AddIncidentElementsOfVertex(iExitVertexLocal)
 end select

 if (minval(iExitFaceVerts).gt.0) then
  call AddIncidentElementsOfFaceVertices(iExitFaceVerts)
 end if

 maxBufferSize = Max(maxBufferSize,iBuffElem)
 
 end subroutine CreateLocalElemBuffer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine ClassifyElementExit(iKind,iFaceLocal,iVertexLocal,iEdgeVerts,iFaceVerts)
 integer, intent(out) :: iKind,iFaceLocal,iVertexLocal
 integer, intent(out) :: iEdgeVerts(2),iFaceVerts(4)
 integer :: iViolationCount,violationDim(3),signDim(3),dominantDim
 real*8  :: overshoot(3),localXi(3)

 iKind = 0
 iFaceLocal = 0
 iVertexLocal = 0
 iEdgeVerts = 0
 iFaceVerts = 0

 localXi = (/XI1,XI2,XI3/)
 overshoot = 0d0
 violationDim = 0
 signDim = 0
 iViolationCount = 0

 do i=1,3
  if (localXi(i).gt.1d0) then
   iViolationCount = iViolationCount + 1
   violationDim(iViolationCount) = i
   signDim(iViolationCount) = 1
   overshoot(i) = localXi(i) - 1d0
  else if (localXi(i).lt.-1d0) then
   iViolationCount = iViolationCount + 1
   violationDim(iViolationCount) = i
   signDim(iViolationCount) = -1
   overshoot(i) = -1d0 - localXi(i)
  end if
 end do

 if (iViolationCount.le.0) return

 dominantDim = violationDim(1)
 do i=2,iViolationCount
  if (overshoot(violationDim(i)).gt.overshoot(dominantDim)) dominantDim = violationDim(i)
 end do

 iFaceLocal = GetFaceIndex(dominantDim,localXi(dominantDim))
 call GetFaceLocalVertices(iFaceLocal,iFaceVerts)

 select case (iViolationCount)
 case (1)
  iKind = 1
 case (2)
  iKind = 2
  call GetEdgeLocalVertices(violationDim(1),signDim(1),violationDim(2),signDim(2),iEdgeVerts)
 case default
  iKind = 3
  iVertexLocal = GetVertexLocalIndex(signDim(1),signDim(2),signDim(3))
 end select

 end subroutine ClassifyElementExit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 integer function GetFaceIndex(iDim,xiValue)
 integer, intent(in) :: iDim
 real*8,  intent(in) :: xiValue

 select case (iDim)
 case (1)
  if (xiValue.lt.0d0) then
   GetFaceIndex = 1
  else
   GetFaceIndex = 2
  end if
 case (2)
  if (xiValue.lt.0d0) then
   GetFaceIndex = 3
  else
   GetFaceIndex = 4
  end if
 case (3)
  if (xiValue.lt.0d0) then
   GetFaceIndex = 5
  else
   GetFaceIndex = 6
  end if
 case default
  GetFaceIndex = 0
 end select

 end function GetFaceIndex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 integer function GetVertexLocalIndex(sign1,sign2,sign3)
 integer, intent(in) :: sign1,sign2,sign3

 GetVertexLocalIndex = 0
 do i=1,8
  if (GetVertexXiSign(i,1).ne.sign1) cycle
  if (GetVertexXiSign(i,2).ne.sign2) cycle
  if (GetVertexXiSign(i,3).ne.sign3) cycle
  GetVertexLocalIndex = i
  return
 end do

 end function GetVertexLocalIndex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 integer function GetVertexXiSign(iVertex,iDim)
 integer, intent(in) :: iVertex,iDim

 select case (iVertex)
 case (1)
  if (iDim.eq.1) GetVertexXiSign = -1
  if (iDim.eq.2) GetVertexXiSign = -1
  if (iDim.eq.3) GetVertexXiSign = -1
 case (2)
  if (iDim.eq.1) GetVertexXiSign =  1
  if (iDim.eq.2) GetVertexXiSign = -1
  if (iDim.eq.3) GetVertexXiSign = -1
 case (3)
  if (iDim.eq.1) GetVertexXiSign =  1
  if (iDim.eq.2) GetVertexXiSign =  1
  if (iDim.eq.3) GetVertexXiSign = -1
 case (4)
  if (iDim.eq.1) GetVertexXiSign = -1
  if (iDim.eq.2) GetVertexXiSign =  1
  if (iDim.eq.3) GetVertexXiSign = -1
 case (5)
  if (iDim.eq.1) GetVertexXiSign = -1
  if (iDim.eq.2) GetVertexXiSign = -1
  if (iDim.eq.3) GetVertexXiSign =  1
 case (6)
  if (iDim.eq.1) GetVertexXiSign =  1
  if (iDim.eq.2) GetVertexXiSign = -1
  if (iDim.eq.3) GetVertexXiSign =  1
 case (7)
  if (iDim.eq.1) GetVertexXiSign =  1
  if (iDim.eq.2) GetVertexXiSign =  1
  if (iDim.eq.3) GetVertexXiSign =  1
 case (8)
  if (iDim.eq.1) GetVertexXiSign = -1
  if (iDim.eq.2) GetVertexXiSign =  1
  if (iDim.eq.3) GetVertexXiSign =  1
 case default
  GetVertexXiSign = 0
 end select

 end function GetVertexXiSign
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine GetEdgeLocalVertices(dimA,signA,dimB,signB,edgeVerts)
 integer, intent(in)  :: dimA,signA,dimB,signB
 integer, intent(out) :: edgeVerts(2)
 integer :: nMatch

 edgeVerts = 0
 nMatch = 0

 do i=1,8
  if (GetVertexXiSign(i,dimA).ne.signA) cycle
  if (GetVertexXiSign(i,dimB).ne.signB) cycle
  nMatch = nMatch + 1
  if (nMatch.le.2) edgeVerts(nMatch) = i
 end do

 end subroutine GetEdgeLocalVertices
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
     local_P8(:,i) = dcorvg(:,ivt)
    END DO
   
    CALL SetLocal_P8(local_P8)

    dist = 1d30
    jPoint = 0

    DO i=1,8
     daux = sqrt((local_P8(1,i)-point(1))**2d0 + (local_P8(2,i)-point(2))**2d0 + (local_P8(3,i)-point(3))**2d0)
     IF (daux.lt.dist) THEN
      dist = daux
      jPoint = i
     END IF
    END DO

    CALL GetPointFromElement(point,jPoint,bFound,pointR,iParticel)

    IF (bFound) THEN
    
     !CALL NDFGL(JEL,1,13,KVERT,KEDGE,KAREA,KDFG,KDFL)
     do i=1,27
      kdfl(i) = i
     end do
     do i=1,8
      kdfg(i) = kvert(i,jel)
     end do
     do i=1,12
      kdfg(8+i) = nvt+kedge(i,jel)
     end do
     do i=1,6
      kdfg(20+i) = nvt+net+karea(i,jel)
     end do
     kdfg(27) = nvt+net+nat+jel
     
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
     
     cdx = 1d1*P(1)
     cdy = 1d1*P(2)
     cdz = 1d1*P(3)

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

!      end if

    end if
    
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !      if (distance.lt. d_CorrDist*0.5d0) then
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !      myFBM%nParticles = 1
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !      CALL GetDistAndProjPToAllSTLs(cdx,cdy,cdz,cpx,cpy,cpz,dist_CGAL)
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !      dist_CGAL = - dist_CGAL
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !      if (dist_CGAL.lt. d_CorrDist*0.5d0) then
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       write (*,*) 'yes',dist_CGAL,d_CorrDist*0.5d0
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       cnormal = [cpx-cdx,cpy-cdy,cpz-cdz]
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       daux = SQRT(cnormal(1)**2d0 + cnormal(2)**2d0 + cnormal(3)**2d0)
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       if (dist_CGAL.gt.0d0) then
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !        cnormal = -cnormal/daux
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       else
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !        cnormal = cnormal/daux
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       end if
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       cdx = cpx + cnormal(1)*(-dist_CGAL + d_CorrDist*0.5d0)
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       cdy = cpy + cnormal(2)*(-dist_CGAL + d_CorrDist*0.5d0)
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       cdz = cpz + cnormal(3)*(-dist_CGAL + d_CorrDist*0.5d0)
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       P=0.1d0*[cdx,cdy,cdz]
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       point = [P(1),P(2),P(3)]
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !      end if

! !       CALL GetDistAndProjPToAllSTLs(cdx,cdy,cdz,cpx,cpy,cpz,dist_CGAL)
! !
! !      if (dist_CGAL.lt. d_CorrDist*0.5d0) then
! !
! !       cnormal = [cpx-cdx,cpy-cdy,cpz-cdz]
! !       daux = SQRT(cnormal(1)**2d0 + cnormal(2)**2d0 + cnormal(3)**2d0)
! !       if (dist_CGAL.gt.0d0) then
! !        cnormal = -cnormal/daux
! !       else
! !        cnormal = cnormal/daux
! !       end if
! !
! !       cdx = cpx + cnormal(1)*d_CorrDist*0.5d0
! !       cdy = cpy + cnormal(2)*d_CorrDist*0.5d0
! !       cdz = cpz + cnormal(3)*d_CorrDist*0.5d0
! !
! !       P=0.1d0*[cdx,cdy,cdz]
! !       point = [P(1),P(2),P(3)]
! !      end if

     iFoundElem = jel

!      bExit = .true.
     RETURN
    END IF

 END SUBROUTINE TraceParticleInTheRightElement

END SUBROUTINE Move_xParticle_Fast


!========================================================================================
!                           Sub: ExtractMatrialProperties
!========================================================================================
subroutine ExtractMatrialProperties(Mat,MatO,S,dcorvg,kvert,kedge,karea,kel_LdA,kel_ColA,nvt,net,nat,nel,iTooFar)
USE PP3D_MPI, ONLY : myid
USE xPart_def, ONLY : minDist

IMPLICIT NONE
INTEGER Mat(*),MatO(*)
REAL*8 dcorvg(3,*),S(*)
INTEGER nvt,net,nat,nel,iTooFar
INTEGER kel_LdA(*),kel_ColA(*)
INTEGER kvert(8,*),kedge(12,*),karea(6,*)
REAL*8 Xmindist,P(3)
INTEGER nList
PARAMETER (nList=256)
integer list(2,nList)
logical bDoX,bIN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER :: iRec,MaxRec = 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer i,ivt,iel,jel,iMat,nnel,mmel,iFoundEl,jj, iCheck,nCheck

iTooFar = 0

do iel=1,nel
 iMat = Mat(iel)
  
 if (iMat.eq.-myid) then
 
!   bIN = .true.
  bIN = .false.
  do i = 1,8
   ivt = kvert(i,iel)
   if (S(ivt).gt.minDist) bIN = .true.
  end do
  
  if (.not.bIn) then
   MatO(iel) = 0
  else

   P = dcorvg(:,nvt+net+nat+iel)
   Xmindist = 1d8
   nnel = 0
   mmel = 0
   iFoundEl = -1
   list = 0
   
   do i=1,8
    ivt = kvert(i,iel)
    DO jj=kel_LdA(ivt),kel_LdA(ivt+1)-1
     jel = kel_ColA(jj)
     iRec = 0
     CALL check_kel(jel,bDoX,iRec)
     if (bDoX) then
      CALL CheckElem(jel,iRec)
     end if
    END DO
   end do
    
!    if (iFoundEl.eq.-1) write(*,*) 'problem!!!',myid
   if (iFoundEl.gt.0) THEN
    MatO(iel) = -Mat(iFoundEl)
   else
    iTooFar = iTooFar + 1
   end if
   
!    if (myid.eq.1) pause
  end if
 else
  if (iMat.gt.0) MatO(iel) = iMat
 end if
 
end do

 contains
 
 RECURSIVE SUBROUTINE CheckElem(kel,jRec)
 integer kel,jRec
 integer j,jvt,lel,ll
 real*8 dist,Q(3)
 logical bDo

 IF (Mat(kel).le.0) THEN
!   write(*,*) "kel:",kel
  do j=1,8
   jvt = kvert(j,kel)
   DO ll=kel_LdA(jvt),kel_LdA(jvt+1)-1
    lel = kel_ColA(ll)
   
    CALL check_kel(lel,bDo,jRec)
    if (bDo) then
     Q = dcorvg(:,nvt+net+nat+lel)
     dist = sqrt((P(1)-Q(1))**2d0 + (P(2)-Q(2))**2d0 + (P(3)-Q(3))**2d0)
     if (dist.lt.Xmindist.and.jRec.le.MaxRec) then
!       if (myid.eq.1) write(*,*) "xx: ",lel,dist,jRec
      jRec = jRec + 1
      CALL CheckElem(lel,jRec)
      jRec = jRec - 1
     end if
    end if
   END DO
  end do
 ELSE
!   iFoundEl = kel
  
  Q = dcorvg(:,nvt+net+nat+kel)
  if (nnel.eq.0) then
   Xmindist=sqrt((P(1)-Q(1))**2d0 + (P(2)-Q(2))**2d0 + (P(3)-Q(3))**2d0)
   iFoundEl = kel
   nnel = nnel + 1
   
!    write(*,*) 'nnel',nnel,kel
  else
   dist=sqrt((P(1)-Q(1))**2d0 + (P(2)-Q(2))**2d0 + (P(3)-Q(3))**2d0)
   if (dist.lt.Xmindist) then
    Xmindist = dist
    iFoundEl = kel
    nnel = nnel + 1
   end if
  end if
 end if
  
 END 
 
 subroutine check_kel(mel,bAdded,jRec)
 integer jRec
 integer mel
 logical bAdded
 integer j
 
 bAdded = .false.
 
 if (mmel.ge.nList) return
 
 do j=1,mmel
  if (list(1,j).eq.mel.and.list(2,j).le.jRec) return
  if (list(1,j).eq.mel.and.list(2,j).gt.jRec) THEN
   list(2,j) = jRec
   bAdded = .true.
   return
  end if
 end do
 
 
 if (S(nvt+net+nat+mel).gt.minDist) then
  bAdded = .true.
 
  mmel = mmel + 1 
  list(1,mmel) = mel
  list(2,mmel) = jRec

!  if (myid.eq.1)  write(*,*) 'added:', mmel,mel,nnel,Mat(mel),iFoundEl
  
 end if
 
 end
end subroutine ExtractMatrialProperties
