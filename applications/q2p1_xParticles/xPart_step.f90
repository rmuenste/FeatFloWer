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

! LOGICAL :: bExit
INTEGER iBuffElem,jBuffElem,BuffElem(nBuffElem),maxBufferSize,maxBufferElem

INTEGER iMonitor(nXX),iAux
REAL*8  dMonitor(nXX),dAux

save KDFG,KDFL,i,iPoint,jPoint,ivt,jel,iel,iFoundElem,kk,iMon
save dist,point,pointR,LocVel0,LocVel1,tLevel,bFound
save iIter, iLoc,iParticel,iLostParticel,idynType,indice,cpx,cpy,cpz,cnormal
save iBuffElem,jBuffElem,BuffElem,maxBufferSize,maxBufferElem
save iMonitor,dMonitor,iAux,dAux

!$OMP   THREADPRIVATE(KDFG,KDFL,i,iPoint,jPoint,ivt,jel,iel,iFoundElem,kk,iMon)
!$OMP   THREADPRIVATE(dist,point,pointR,LocVel0,LocVel1,tLevel,bFound)
!$OMP   THREADPRIVATE(iIter, iLoc,iParticel,iLostParticel,idynType,indice,cpx,cpy,cpz,cnormal)
!$OMP   THREADPRIVATE(iBuffElem,jBuffElem,BuffElem,maxBufferSize,maxBufferElem)
!$OMP   THREADPRIVATE(iMonitor,dMonitor,iAux,dAux)

iLostParticel = 0
iActiveParticel = 0
maxBufferSize = 0
maxBufferElem = 0

DO iParticel = iMinParticleID,iMaxParticleID

iIter = 0

if (myActiveSet(iParticel)%id.eq.0) then

point  = myActiveSet(iParticel)%coor
tLevel = myActiveSet(iParticel)%time
dParticleVelo = myActiveSet(iParticel)%velo

55 CONTINUE

bFound = .FALSE.

CALL FindRankingInOctTreeOMP(dcorvg,nvt,point,iMonitor,dMonitor,nxx)
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
 myLostSet(indice)%velo   = dParticleVelo
 myLostSet(indice)%indice = indice
 myLostSet(indice)%id     = myid
 
 myActiveSet(iParticel)%id   = -1
 
 GOTO 2
END IF

1 CONTINUE

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