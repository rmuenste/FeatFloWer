module geometry_processing
  use Sigma_User, only: mySigma, myProcess,KNET_elem,SKNET_elem,EKNET_elem,&
  FOERD_elem, SME_elem, ZME_elem , TrueKNET_elem, DistTolerance
 
  use var_QuadScalar, only: Shell,Screw,MixerKNPR,ScrewDist,myHeatObjects
  !------------------------------------------------------------------------------------------------
  ! A module for functions that perform operations on the 
  ! geometry immersed into the simulation domain or boundary geometry
  !------------------------------------------------------------------------------------------------
  implicit none

  ! a variable for counting the outputs
  integer :: ifile = 0

  real*8 :: dEpsDist

contains
!
! ----------------------------------------------
!
SUBROUTINE QuadScalar_MixerKnpr(dcorvg,kvert,kedge,karea,nel,nvt,nat,net,dist1,dist2,dist3)
REAL*8  dcorvg(3,*),dist1(*),dist2(*),dist3(*)
integer nel,nvt,nat,net
INTEGER kvert(8,*),kedge(12,*),karea(6,*)
REAL*8 PX,PY,PZ,DIST,d12,d4,d5,dS0,dS1,dS2,timeLevel
INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4
INTEGER NeighE(2,12),NeighA(4,6)
DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/

timeLevel = (myProcess%Angle/360d0)/(myProcess%Umdr/60d0)

DO i=1,nvt
 PX = dcorvg(1,I)
 PY = dcorvg(2,I)
 PZ = dcorvg(3,I)
 CALL GetMixerKnpr(PX,PY,PZ,MixerKNPR(i),dS0,dS1,dS2,timeLevel)
 Dist1(i) = dS1
 Dist2(i) = dS2
 Dist3(i) = dS0
 ScrewDist(:,i) = [dS1,dS2]
END DO

k=1
DO i=1,nel
 DO j=1,12
  IF (k.eq.kedge(j,i)) THEN
   ivt1 = kvert(NeighE(1,j),i)
   ivt2 = kvert(NeighE(2,j),i)
   PX = 0.5d0*(dcorvg(1,ivt1)+dcorvg(1,ivt2))
   PY = 0.5d0*(dcorvg(2,ivt1)+dcorvg(2,ivt2))
   PZ = 0.5d0*(dcorvg(3,ivt1)+dcorvg(3,ivt2))
   CALL GetMixerKnpr(PX,PY,PZ,MixerKNPR(nvt+k),dS0,dS1,dS2,timeLevel)
   Dist1(nvt+k) = dS1
   Dist2(nvt+k) = dS2
   Dist3(nvt+k) = dS0
   ScrewDist(:,nvt+k) = [dS1,dS2]
   k = k + 1
  END IF
 END DO
END DO

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
   CALL GetMixerKnpr(PX,PY,PZ,MixerKNPR(nvt+net+k),dS0,dS1,dS2,timeLevel)
   Dist1(nvt+net+k) = dS1
   Dist2(nvt+net+k) = dS2
   Dist3(nvt+net+k) = dS0
   ScrewDist(:,nvt+net+k) = [dS1,dS2]
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
 CALL GetMixerKnpr(PX,PY,PZ,MixerKNPR(nvt+net+nat+i),dS0,dS1,dS2,timeLevel)
 Dist1(nvt+net+nat+i) = dS1
 Dist2(nvt+net+nat+i) = dS2
 Dist3(nvt+net+nat+i) = dS0
 ScrewDist(:,nvt+net+nat+i) = [dS1,dS2]
END DO

do i=1,nvt+net+nat+nel
 PX = dcorvg(1,I)
 PY = dcorvg(2,I)
 PZ = dcorvg(3,I)
 Screw(i) = min(Dist1(i),Dist2(i))
 Shell(i) = Dist3(i)
! Screw(i) = sqrt(PX*PX + PY*PY) - 75d0
end do

! if (myid.eq.2) WRITE(*,*) dcorvg(:,17146),Distamce(17146),nvt+net+nat

END SUBROUTINE QuadScalar_MixerKnpr
!
!********************GEOMETRY****************************
!
SUBROUTINE GetMixerKnpr(X,Y,Z,inpr,D0,D1,D2,t)
IMPLICIT NONE
REAL*8 X,Y,Z,t,d0,d1,d2
REAL*8 :: dKnetMin,dKnetMax
INTEGER :: inpr, l, k, nmbr,lKnet
REAL*8 :: dAlpha,XT,YT,ZT,XB,YB,ZB
REAL*8 :: dBeta,XTT,YTT,ZTT,dist
REAL*8 dScale, dist1, dist2,dEps,dCut
REAL*8 dSeg0,dSeg1,dSeg2,tt
REAL*8 myPI
LOGICAL :: bProjection=.FALSE.
REAL*8 :: ProjP(3)

dEps = mySigma%a/1d5

myPI = dATAN(1d0)*4d0

D0 = DistTolerance
D1 = DistTolerance
D2 = DistTolerance
inpr = 0

tt = t 
!----------------------------------------------------------
DO k=1, mySigma%NumberOfSeg
 dSeg0=DistTolerance
 dSeg1=DistTolerance
 dSeg2=DistTolerance
 
 IF (mySigma%mySegment(k)%ART.EQ.'STL_R')  CALL STLR_elem(X,Y,Z,tt,k,dSeg1,dSeg2,inpr,bProjection,ProjP)
 IF (mySigma%mySegment(k)%ART.EQ.'STL_L')  CALL STLL_elem(X,Y,Z,tt,k,dSeg1,dSeg2,inpr,bProjection,ProjP)
 IF (mySigma%mySegment(k)%ART.EQ.'STL'  )  CALL STL_elem(X,Y,Z,tt,k,dSeg0,dSeg1,dSeg2,inpr,bProjection,ProjP)
 IF (mySigma%mySegment(k)%ART.EQ.'KNET' ) CALL KNET_elem(X,Y,Z,tt,k,dSeg1,dSeg2,inpr)
 IF (mySigma%mySegment(k)%ART.EQ.'TKNET') CALL TrueKNET_elem(X,Y,Z,tt,k,dSeg1,dSeg2,inpr)
 IF (mySigma%mySegment(k)%ART.EQ.'EKNET') CALL EKNET_elem(X,Y,Z,tt,k,dSeg1,dSeg2,inpr)
 IF (mySigma%mySegment(k)%ART.EQ.'SKNET') CALL SKNET_elem(X,Y,Z,tt,k,dSeg1,dSeg2,inpr)
 IF (mySigma%mySegment(k)%ART.EQ.'FOERD') CALL FOERD_elem(X,Y,Z,tt,k,dSeg1,dSeg2,inpr)
 IF (mySigma%mySegment(k)%ART.EQ.'SME'  ) CALL SME_elem(X,Y,Z,tt,k,dSeg1,dSeg2,inpr)
 IF (mySigma%mySegment(k)%ART.EQ.'ZME'  ) CALL ZME_elem(X,Y,Z,tt,k,dSeg1,dSeg2,inpr)
 D0 = max(-DistTolerance,min(dSeg0,D0))
 D1 = max(-DistTolerance,min(dSeg1,D1))
 D2 = max(-DistTolerance,min(dSeg2,D2))
END DO
!-------------------------------------------------------

return

END SUBROUTINE GetMixerKnpr
!
!********************GEOMETRY****************************
!
SUBROUTINE STL_elem(X,Y,Z,t,iSeg,d0,d1,d2,inpr,bProjection,ProjP)
IMPLICIT NONE
INTEGER inpr
REAL*8 X,Y,Z,d1,d2,d0,daux
REAL*8 dAlpha, XT,YT,ZT, XP,YP,ZP, XB,YB,ZB,XT_STL,YT_STL,ZT_STL
REAL*8 t,myPI,dUnitScale
INTEGER iSTL,iSeg,iFile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
LOGICAL :: bProjection
REAL*8 :: ProjP(3),DDDD

! d1 = 5d0
! d2 = 5d0
IF (mySigma%mySegment(iSeg)%Unit.eq.'MM') dUnitScale = 1d+1
IF (mySigma%mySegment(iSeg)%Unit.eq.'CM') dUnitScale = 1d0
IF (mySigma%mySegment(iSeg)%Unit.eq.'DM') dUnitScale = 1d-1

IF (mySigma%mySegment(iSeg)%ObjectType.eq.'SCREW') THEN

  t = (myProcess%Angle/360d0)/(myProcess%Umdr/60d0)

  myPI = dATAN(1d0)*4d0

  IF (ADJUSTL(TRIM(mySigma%RotationAxis)).EQ."PARALLEL") THEN
  XB = X
  YB = Y-mySigma%a/2d0
  ZB = Z
  ELSE
  CALL TransformPointToNonparallelRotAxis(x,y,z,XB,YB,ZB,-1d0)
  END IF


  ! First the point needs to be transformed back to time = 0
  dAlpha = mySigma%mySegment(iSeg)%StartAlpha + (-t*myPI*(myProcess%Umdr/3d1) + 0d0*myPI/2d0)*DBLE(myProcess%iInd*myProcess%ind)
  XT = XB*cos(dAlpha) - YB*sin(dAlpha)
  YT = XB*sin(dAlpha) + YB*cos(dAlpha)
  ZT = ZB

  XT_STL = dUnitScale*XT
  YT_STL = dUnitScale*YT
  ZT_STL = dUnitScale*ZT - dUnitScale*mySigma%mySegment(iSeg)%Min

  DO iFile=1,mySigma%mySegment(iSeg)%nOFFfiles
  iSTL = mySigma%mySegment(iSeg)%idxCgal(iFile)
  CALL GetDistToSTL(XT_STL,YT_STL,ZT_STL,iSTL,d1,.TRUE.)
  if (bProjection.and.d1.lt.DistTolerance) THEN
   call getclosestpointid(XT_STL,YT_STL,ZT_STL,ProjP(1),ProjP(2),ProjP(3),DDDD,iSTL-1)
   ProjP(1) =  ProjP(1)/dUnitScale
   ProjP(2) =  ProjP(2)/dUnitScale
   ProjP(3) = (ProjP(3) + mySigma%mySegment(iSeg)%Min*dUnitScale)/dUnitScale
   XT = ProjP(1)*cos(-dAlpha) - ProjP(2)*sin(-dAlpha)
   YT = ProjP(1)*sin(-dAlpha) + ProjP(2)*cos(-dAlpha)
   ZT = ProjP(3)
   ProjP(1) =  XT 
   ProjP(2) =  YT+mySigma%a/2d0
   ProjP(3) =  ZT
  end if
  d1 = max(-DistTolerance,min(DistTolerance,d1/dUnitScale))
  END DO

  if (d1.lt.0d0) inpr = 101

  IF (ADJUSTL(TRIM(mySigma%RotationAxis)).EQ."PARALLEL") THEN
  XB = X
  YB = Y+mySigma%a/2d0
  ZB = Z
  ELSE
  CALL TransformPointToNonparallelRotAxis(x,y,z,XB,YB,ZB,+1d0)
  END IF

  ! First the point needs to be transformed back to time = 0
  IF (mySigma%GANGZAHL .EQ. 1) dAlpha = mySigma%mySegment(iSeg)%StartAlpha -t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
  IF (mySigma%GANGZAHL .EQ. 2) dAlpha = mySigma%mySegment(iSeg)%StartAlpha + (-t*myPI*(myProcess%Umdr/3d1)+myPI/2d0)*myProcess%ind
  IF (mySigma%GANGZAHL .EQ. 3) dAlpha = mySigma%mySegment(iSeg)%StartAlpha -t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
  IF (mySigma%GANGZAHL .EQ. 4) dAlpha = mySigma%mySegment(iSeg)%StartAlpha + (-t*myPI*(myProcess%Umdr/3d1)+myPI/4d0)*myProcess%ind
  !dAlpha = 0.d0 + (-t*myPI*(myProcess%Umdr/3d1) + 0d0*myPI/2d0)*myProcess%ind

  XT = XB*cos(dAlpha) - YB*sin(dAlpha)
  YT = XB*sin(dAlpha) + YB*cos(dAlpha)
  ZT = ZB

  XT_STL = dUnitScale*XT
  YT_STL = dUnitScale*YT
  ZT_STL = dUnitScale*ZT  - dUnitScale*mySigma%mySegment(iSeg)%Min

  DO iFile=1,mySigma%mySegment(iSeg)%nOFFfiles
  iSTL = mySigma%mySegment(iSeg)%idxCgal(iFile)
  CALL GetDistToSTL(XT_STL,YT_STL,ZT_STL,iSTL,d2,.TRUE.)
  if (bProjection.and.d2.lt.DistTolerance) THEN
   call getclosestpointid(XT_STL,YT_STL,ZT_STL,ProjP(1),ProjP(2),ProjP(3),DDDD,iSTL-1)
   ProjP(1) =  ProjP(1)/dUnitScale
   ProjP(2) =  ProjP(2)/dUnitScale
   ProjP(3) = (ProjP(3) + mySigma%mySegment(iSeg)%Min*dUnitScale)/dUnitScale
   XT = ProjP(1)*cos(-dAlpha) - ProjP(2)*sin(-dAlpha)
   YT = ProjP(1)*sin(-dAlpha) + ProjP(2)*cos(-dAlpha)
   ZT = ProjP(3)
   ProjP(1) =  XT 
   ProjP(2) =  YT-mySigma%a/2d0
   ProjP(3) =  ZT
  end if
  d2 = max(-DistTolerance,min(DistTolerance,d2/dUnitScale))
  END DO

  if (d2.lt.0d0) inpr = 102
END IF


IF (mySigma%mySegment(iSeg)%ObjectType.eq.'OBSTACLE') THEN

  XT_STL = dUnitScale*X
  YT_STL = dUnitScale*Y
  ZT_STL = dUnitScale*Z - dUnitScale*mySigma%mySegment(iSeg)%Min
  
  DO iFile=1,mySigma%mySegment(iSeg)%nOFFfiles
   iSTL = mySigma%mySegment(iSeg)%idxCgal(iFile)
   CALL GetDistToSTL(XT_STL,YT_STL,ZT_STL,iSTL,daux,.TRUE.)
   if (bProjection.and.daux.lt.DistTolerance) THEN
    call getclosestpointid(XT_STL,YT_STL,ZT_STL,ProjP(1),ProjP(2),ProjP(3),DDDD,iSTL-1)
    ProjP(1) =  ProjP(1)/dUnitScale
    ProjP(2) =  ProjP(2)/dUnitScale
    ProjP(3) = (ProjP(3) + mySigma%mySegment(iSeg)%Min*dUnitScale)/dUnitScale
   end if
   daux = max(-DistTolerance,min(DistTolerance,daux/dUnitScale))
   d0   = min(daux,d0)
  END DO
  
  if (d0.lt.0d0) inpr = 100
  
END IF

IF (mySigma%mySegment(iSeg)%ObjectType.eq.'DIE') THEN

  XT_STL = dUnitScale*X
  YT_STL = dUnitScale*Y
  ZT_STL = dUnitScale*Z - dUnitScale*mySigma%mySegment(iSeg)%Min
  
  DO iFile=1,mySigma%mySegment(iSeg)%nOFFfiles
   iSTL = mySigma%mySegment(iSeg)%idxCgal(iFile)
   CALL GetDistToSTL(XT_STL,YT_STL,ZT_STL,iSTL,daux,.TRUE.)
   daux = -daux
   if (bProjection.and.daux.lt.DistTolerance) THEN
    call getclosestpointid(XT_STL,YT_STL,ZT_STL,ProjP(1),ProjP(2),ProjP(3),DDDD,iSTL-1)
    ProjP(1) =  ProjP(1)/dUnitScale
    ProjP(2) =  ProjP(2)/dUnitScale
    ProjP(3) = (ProjP(3) + mySigma%mySegment(iSeg)%Min*dUnitScale)/dUnitScale
!    write(*,*) 'P',XT_STL,YT_STL,ZT_STL,ProjP(1),ProjP(2),ProjP(3),daux,DDDD
   end if
   daux = max(-DistTolerance,min(DistTolerance,daux/dUnitScale))
   d0   = min(daux,d0)
  END DO
  
  if (d0.lt.0d0) inpr = 100
  
END IF

END SUBROUTINE STL_elem
!
!********************GEOMETRY****************************
!
SUBROUTINE STLL_elem(X,Y,Z,t,iSeg,d1,d2,inpr,bProjection,ProjP)
IMPLICIT NONE
INTEGER inpr
REAL*8 X,Y,Z,d1,d2
REAL*8 dAlpha, XT,YT,ZT, XP,YP,ZP, XB,YB,ZB,XT_STL,YT_STL,ZT_STL
REAL*8 t,myPI,dUnitScale
INTEGER iSTL,iSeg,iFile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
LOGICAL :: bProjection
REAL*8 :: ProjP(3),daux

! d1 = 5d0
! d2 = 5d0
IF (mySigma%mySegment(iSeg)%Unit.eq.'MM') dUnitScale = 1d+1
IF (mySigma%mySegment(iSeg)%Unit.eq.'CM') dUnitScale = 1d0
IF (mySigma%mySegment(iSeg)%Unit.eq.'DM') dUnitScale = 1d-1

t = (myProcess%Angle/360d0)/(myProcess%Umdr/60d0)

myPI = dATAN(1d0)*4d0

IF (ADJUSTL(TRIM(mySigma%RotationAxis)).EQ."PARALLEL") THEN
 XB = X 
 YB = Y - mySigma%a/2d0
 ZB = Z
ELSE
 CALL TransformPointToNonparallelRotAxis(x,y,z,XB,YB,ZB,-1d0)
END IF


! First the point needs to be transformed back to time = 0
dAlpha = mySigma%mySegment(iSeg)%StartAlpha + (-t*myPI*(myProcess%Umdr/3d1) + 0d0*myPI/2d0)*(DBLE(myProcess%iInd*myProcess%ind))
! write(*,*) 'L',dAlpha
XT = XB*cos(dAlpha) - YB*sin(dAlpha)
YT = XB*sin(dAlpha) + YB*cos(dAlpha)
ZT = ZB

XT_STL = dUnitScale*XT
YT_STL = dUnitScale*YT
ZT_STL = dUnitScale*ZT - dUnitScale*mySigma%mySegment(iSeg)%Min

DO iFile=1,mySigma%mySegment(iSeg)%nOFFfiles
 iSTL = mySigma%mySegment(iSeg)%idxCgal(iFile)
 CALL GetDistToSTL(XT_STL,YT_STL,ZT_STL,iSTL,d1,.TRUE.)
 if (bProjection.and.d1.lt.DistTolerance) THEN
  call getclosestpointid(XT_STL,YT_STL,ZT_STL,ProjP(1),ProjP(2),ProjP(3),daux,iSTL-1)
  ProjP(1) =  ProjP(1)/dUnitScale
  ProjP(2) =  ProjP(2)/dUnitScale
  ProjP(3) = (ProjP(3) + mySigma%mySegment(iSeg)%Min*dUnitScale)/dUnitScale
  XT = ProjP(1)*cos(-dAlpha) - ProjP(2)*sin(-dAlpha)
  YT = ProjP(1)*sin(-dAlpha) + ProjP(2)*cos(-dAlpha)
  ZT = ProjP(3)
  ProjP(1) =  XT 
  ProjP(2) =  YT+mySigma%a/2d0
  ProjP(3) =  ZT
 end if
 d1 = max(-DistTolerance,min(DistTolerance,d1/dUnitScale))
END DO
d1 = d1 + 0.0d0

if (d1.lt.0d0) inpr = 101

d2 = DistTolerance
if (d2.lt.0d0) inpr = 102

END SUBROUTINE STLL_elem
!
!********************GEOMETRY****************************
!
SUBROUTINE STLR_elem(X,Y,Z,t,iSeg,d1,d2,inpr,bProjection,ProjP)
IMPLICIT NONE
INTEGER inpr
REAL*8 X,Y,Z,d1,d2
REAL*8 dAlpha, XT,YT,ZT, XP,YP,ZP, XB,YB,ZB,XT_STL,YT_STL,ZT_STL
REAL*8 t,myPI,dUnitScale
INTEGER iSTL,iSeg,iFile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
LOGICAL :: bProjection
REAL*8 :: ProjP(3),daux

! d1 = 5d0
! d2 = 5d0
IF (mySigma%mySegment(iSeg)%Unit.eq.'MM') dUnitScale = 1d+1
IF (mySigma%mySegment(iSeg)%Unit.eq.'CM') dUnitScale = 1d0
IF (mySigma%mySegment(iSeg)%Unit.eq.'DM') dUnitScale = 1d-1

t = (myProcess%Angle/360d0)/(myProcess%Umdr/60d0)

myPI = dATAN(1d0)*4d0

IF (ADJUSTL(TRIM(mySigma%RotationAxis)).EQ."PARALLEL") THEN
 XB = X 
 YB = Y +  mySigma%a/2d0
 ZB = Z
ELSE
 CALL TransformPointToNonparallelRotAxis(x,y,z,XB,YB,ZB,+1d0)
END IF

! dAlpha = mySigma%mySegment(iSeg)%StartAlpha + (-t*myPI*(myProcess%Umdr/3d1) + 0d0*myPI/2d0)*myProcess%ind
! write(*,*) 'R',dAlpha
! ! First the point needs to be transformed back to time = 0
IF (mySigma%GANGZAHL .EQ. 1) dAlpha = mySigma%mySegment(iSeg)%StartAlpha -t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
IF (mySigma%GANGZAHL .EQ. 2) dAlpha = mySigma%mySegment(iSeg)%StartAlpha + (-t*myPI*(myProcess%Umdr/3d1)+myPI/2d0)*myProcess%ind
IF (mySigma%GANGZAHL .EQ. 3) dAlpha = mySigma%mySegment(iSeg)%StartAlpha -t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
IF (mySigma%GANGZAHL .EQ. 4) dAlpha = mySigma%mySegment(iSeg)%StartAlpha + (-t*myPI*(myProcess%Umdr/3d1)+myPI/4d0)*myProcess%ind
!dAlpha = 0.d0 + (-t*myPI*(myProcess%Umdr/3d1) + 0d0*myPI/2d0)*myProcess%ind

XT = XB*cos(dAlpha) - YB*sin(dAlpha)
YT = XB*sin(dAlpha) + YB*cos(dAlpha)
ZT = ZB

XT_STL = dUnitScale*XT
YT_STL = dUnitScale*YT
ZT_STL = dUnitScale*ZT  - dUnitScale*mySigma%mySegment(iSeg)%Min

DO iFile=1,mySigma%mySegment(iSeg)%nOFFfiles
 iSTL = mySigma%mySegment(iSeg)%idxCgal(iFile)
 CALL GetDistToSTL(XT_STL,YT_STL,ZT_STL,iSTL,d2,.TRUE.)
 if (bProjection.and.d2.lt.DistTolerance) THEN
  call getclosestpointid(XT_STL,YT_STL,ZT_STL,ProjP(1),ProjP(2),ProjP(3),daux,iSTL-1)
  ProjP(1) =  ProjP(1)/dUnitScale
  ProjP(2) =  ProjP(2)/dUnitScale
  ProjP(3) = (ProjP(3) + mySigma%mySegment(iSeg)%Min*dUnitScale)/dUnitScale
  XT = ProjP(1)*cos(-dAlpha) - ProjP(2)*sin(-dAlpha)
  YT = ProjP(1)*sin(-dAlpha) + ProjP(2)*cos(-dAlpha)
  ZT = ProjP(3)
  ProjP(1) =  XT 
  ProjP(2) =  YT-mySigma%a/2d0
  ProjP(3) =  ZT
 end if
 d2 = max(-DistTolerance,min(DistTolerance,d2/dUnitScale))
END DO
d2 = d2 + 0.0d0

if (d2.lt.0d0) inpr = 102

d1 = DistTolerance
if (d1.lt.0d0) inpr = 101

END SUBROUTINE STLR_elem
!------------------------------------------------------------------------------------------------
! The subroutine calculates a distance function on the mesh given 
! by dcorvg, etc. The distance is computed to the list of screw geometries
! maintained in the global structure mySigma. This function for distance computations
! is intented for screw setups only.
!------------------------------------------------------------------------------------------------
subroutine calcDistanceFunction_sse(dcorvg,kvert,kedge,karea,nel,nvt,nat,net,dst1,dst2,dVaux)

  real*8, dimension(:,:), intent(inout) :: dcorvg

  ! output array for the distance to a die or obstacles
  real*8, dimension(:), intent(inout) :: dst1

  ! output array for the distance to the screw
  real*8, dimension(:), intent(inout) :: dst2

  real*8, dimension(:), intent(inout) :: dVaux

  integer :: nel,nvt,nat,net

  integer, dimension(:,:), intent(inout) :: karea, kedge, kvert

  ! local variables
  integer ElemInd(27)

  real*8 ElemCoor(3,27),ElemDist(27),ElemRad,ElemSign,PointSign,dist

  real*8 d,PX,PY,PZ,dS,dUnitScale,dLocEpsDist,dRotAngle

  real*8, dimension(3) :: point

  integer :: i,j,ndof,iel,ipc

  integer iSTL,iSeg,iFile

  ndof = nel+nvt+nat+net

  ! Loop over all segments of the screw
  DO iSeg=1,mySigma%NumberOfSeg

  ! If we find a discrete geometry segment
  IF (adjustl(trim(mySigma%mySegment(iSeg)%ART)).eq."STL") THEN

    IF (mySigma%mySegment(iSeg)%ObjectType.eq.'SCREW') THEN
     dRotAngle = myProcess%Angle
    END IF
    IF (mySigma%mySegment(iSeg)%ObjectType.eq.'DIE'.or.mySigma%mySegment(iSeg)%ObjectType.eq.'OBSTACLE') THEN
     dRotAngle  = 0d0
    END IF

    IF (mySigma%mySegment(iSeg)%Unit.eq.'MM') dUnitScale = 1d+1
    IF (mySigma%mySegment(iSeg)%Unit.eq.'CM') dUnitScale = 1d0
    IF (mySigma%mySegment(iSeg)%Unit.eq.'DM') dUnitScale = 1d-1

    dLocEpsDist = dUnitScale*dEpsDist

    DO iFile=1,mySigma%mySegment(iSeg)%nOFFfiles

    dVaux(1:ndof) = 0d0  
    iSTL = mySigma%mySegment(iSeg)%idxCgal(iFile)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! iSTL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO iel=1,nel

    ! Indices of the vertices
    ElemInd( 1:8 )  = kvert(1:8,iel)

    ! Indices of the edge mid points
    ElemInd( 9:20)  = nvt+kedge(1:12,iel)

    ! Indices of the face mid points
    ElemInd(21:26)  = nvt+net+karea(1:6,iel)

    ! Index of the element mid point
    ElemInd(   27)  = nvt+net+nat+iel

    ! Coordinates of the element dofs
    ElemCoor(1,:) = dcorvg(1,ElemInd(:))
    ElemCoor(2,:) = dcorvg(2,ElemInd(:))
    ElemCoor(3,:) = dcorvg(3,ElemInd(:))

    CALL TransformPoints()

    ! Coordinates of the element mid point
    PX = ElemCoor(1,27)
    PY = ElemCoor(2,27)
    PZ = ElemCoor(3,27)

    ! Compute the element radius
    CALL GetElemRad()

    CALL GetDistToSTL(PX,PY,PZ,iSTL,d,.true.)

    ! Set the distance in the element mid point
    dVaux(ElemInd(27))= d

    ! Do not divide by zero
    IF (d .ne. 0d0)THEN
      ElemSign = d/abs(d)
    ELSE
      ElemSign = -1d0
    END IF

    IF (abs(d).lt.dLocEpsDist) THEN 

      ! If the distance to the surface is greater than
      ! the element radius, we do not need to compute
      ! the FBM-function for the other dofs
      IF (abs(d).gt.ElemRad) THEN

        DO i=1,26
        PX = ElemCoor(1,i)
        PY = ElemCoor(2,i)
        PZ = ElemCoor(3,i)

        ! check if dist is already computed
        IF (dVaux(ElemInd(i)).eq.0d0) THEN
          CALL GetDistToSTL(PX,PY,PZ,iSTL,d,.false.)
          ! Set distance to minium of eps and d
          dVaux(ElemInd(i))= ElemSign*MIN(ABS(d),dLocEpsDist)
        END IF
        END DO

      ELSE ! abs(d).gt.ElemRad

        DO i=1,26
        PX = ElemCoor(1,i)
        PY = ElemCoor(2,i)
        PZ = ElemCoor(3,i)

        IF (dVaux(ElemInd(i)).eq.0d0) THEN

          CALL GetDistToSTL(PX,PY,PZ,iSTL,d,.true.)

          IF (d.ne.0d0)THEN
            PointSign = d/abs(d)
          ELSE
            PointSign  = -1d0
          END IF

          ! Set distance to minium of eps and d
          dVaux(ElemInd(i))= PointSign*MIN(ABS(d),dLocEpsDist)

        END IF
        END DO
      END IF

    ELSE ! abs(d).lt.dLocEpsDist

      ! The distance to element the mid point
      ! is greater than eps: we can set
      ! all other dofs to eps too

      DO i=1,27

      IF (dVaux(ElemInd(i)).eq.0d0) THEN
        dVaux(ElemInd(i))= ElemSign*dLocEpsDist
      end if

      end do

    end if ! abs(d).lt.dLocEpsDist

    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! iSTL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (adjustl(trim(mySigma%mySegment(iSeg)%ObjectType)).eq."DIE") THEN
      do i=1,ndof
      dst1(i) = min(dst1(i),-dVaux(i)/dUnitScale)
      end do
    end if  

!  ipc = 0
!  do i=1,ndof
!
!    Px = dcorvg(1,i)
!    Py = dcorvg(2,i)
!    Pz = dcorvg(3,i)
!
!    call TransformPoint(px,py,pz,point(1),point(2),point(3))
!
!    call getdistanceid(point(1), point(2), point(3), dist, ipc);        
!    distance(i) = dist
!  end do

    if (adjustl(trim(mySigma%mySegment(iSeg)%ObjectType)).eq."SCREW") THEN
      do i=1,ndof
      dst2(i) = min(dst2(i),+dVaux(i)/dUnitScale)
      end do
    end if  

    if (adjustl(trim(mySigma%mySegment(iSeg)%ObjectType)).eq."OBSTACLE") THEN
      do i=1,ndof
      dst1(i) = min(dst1(i),+dVaux(i)/dUnitScale)
      end do
    end if  

    end do

  end if
  end do

  DO i=1,ndof
   Shell(i) = dst1(i)
   IF (Shell(i).le.0d0) THEN
    MixerKNPR(i) = 100
   END IF

   Screw(i) = dst2(i)
   IF (Screw(i).le.0d0) THEN
    MixerKNPR(i) = 103
   END IF
  END DO

contains
subroutine GetElemRad
  REAL*8 daux

  ElemRad = 0d0

  do i=1,8
  dist = (ElemCoor(1,i)-PX)**2d0 + (ElemCoor(2,i)-PY)**2d0 + (ElemCoor(3,i)-PZ)**2d0
  ElemRad = max(dist,ElemRad)
  end do
  ElemRad = 2.1d0*SQRT(ElemRad)

End SUBROUTINE GetElemRad

SUBROUTINE TransformPoints()
  REAL*8 t,XB,YB,ZB,XT,YT,ZT,dAlpha
  REAL*8 :: myPI = dATAN(1d0)*4d0
  INTEGER iP

  t = (dRotAngle/360d0)/(myProcess%Umdr/60d0)

  DO iP=1,27

  ! Point needs to be offseted to its segment position
  XB = dUnitScale*ElemCoor(1,iP)
  YB = dUnitScale*ElemCoor(2,iP)
  ZB = dUnitScale*(ElemCoor(3,iP) - mySigma%mySegment(iSeg)%Min)

  ! Point needs to be transformed back to a time=0 due to the rotation
  dAlpha = 0d0- t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
  XT = XB*cos(dAlpha) - YB*sin(dAlpha)
  YT = XB*sin(dAlpha) + YB*cos(dAlpha)
  ZT = ZB

  ElemCoor(1,iP) = XT
  ElemCoor(2,iP) = YT
  ElemCoor(3,iP) = ZT

  END DO

  END 

SUBROUTINE TransformPoint(x, y, z, xb, yb, zb)
  REAL*8 t,XB,YB,ZB,X,Y,Z,dAlpha
  REAL*8 :: myPI = dATAN(1d0)*4d0
  INTEGER iP

  t = (myProcess%Angle/360d0)/(myProcess%Umdr/60d0)

  ! Point needs to be offseted to its segment position
  XB = dUnitScale * x
  YB = dUnitScale * y
  ZB = dUnitScale * z - dUnitScale*mySigma%mySegment(1)%Min

  ! Point needs to be transformed back to a time=0 due to the rotation
  dAlpha = 0d0- t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
  XB = XB*cos(dAlpha) - YB*sin(dAlpha)
  YB = XB*sin(dAlpha) + YB*cos(dAlpha)
  ZB = ZB

  END 

end subroutine calcDistanceFunction_sse
!------------------------------------------------------------------------------------------------
! The subroutine calculates a distance function on the mesh given 
! by dcorvg, etc. The distance is computed to the list of screw geometries
! maintained in the global structure mySigma. This function for distance computations
! is intented for screw setups only.
!------------------------------------------------------------------------------------------------
subroutine calcDistanceFunction_netzsch(dcorvg,kvert,kedge,karea,nel,nvt,nat,net,dst1,dst2,dVaux)

  real*8, dimension(:,:), intent(inout) :: dcorvg

  ! output array for the distance to a die or obstacles
  real*8, dimension(:), intent(inout) :: dst1

  ! output array for the distance to the screw
  real*8, dimension(:), intent(inout) :: dst2

  real*8, dimension(:), intent(inout) :: dVaux

  integer :: nel,nvt,nat,net

  integer, dimension(:,:), intent(inout) :: karea, kedge, kvert

  ! local variables
  integer ElemInd(27)

  real*8 ElemCoor(3,27),ElemDist(27),ElemRad,ElemSign,PointSign,dist

  real*8 d,PX,PY,PZ,dS,dUnitScale,dLocEpsDist,dRotAngle,dYcShift,dYtShift

  real*8, dimension(3) :: point

  integer :: i,j,ndof,iel,ipc

  integer iSTL,iSeg,iFile

  ndof = nel+nvt+nat+net

  ! Loop over all segments of the screw
  DO iSeg=1,mySigma%NumberOfSeg

  ! If we find a discrete geometry segment
  IF (adjustl(trim(mySigma%mySegment(iSeg)%ART)).eq."STL") THEN

    IF (mySigma%mySegment(iSeg)%ObjectType.eq.'SCREW') THEN
     dRotAngle = myProcess%Angle
     dYtShift   = 2d0*8.8d0*dcos(myProcess%Angle*(dATAN(1d0)*8d0)/360d0)
     dYcShift   = 2d0*8.8d0
    END IF

    IF (mySigma%mySegment(iSeg)%Unit.eq.'MM') dUnitScale = 1d+1
    IF (mySigma%mySegment(iSeg)%Unit.eq.'CM') dUnitScale = 1d0
    IF (mySigma%mySegment(iSeg)%Unit.eq.'DM') dUnitScale = 1d-1

    dLocEpsDist = dUnitScale*dEpsDist

    DO iFile=1,mySigma%mySegment(iSeg)%nOFFfiles

    dVaux(1:ndof) = 0d0  
    iSTL = mySigma%mySegment(iSeg)%idxCgal(iFile)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! iSTL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO iel=1,nel

    ! Indices of the vertices
    ElemInd( 1:8 )  = kvert(1:8,iel)

    ! Indices of the edge mid points
    ElemInd( 9:20)  = nvt+kedge(1:12,iel)

    ! Indices of the face mid points
    ElemInd(21:26)  = nvt+net+karea(1:6,iel)

    ! Index of the element mid point
    ElemInd(   27)  = nvt+net+nat+iel

    ! Coordinates of the element dofs
    ElemCoor(1,:) = dcorvg(1,ElemInd(:))
    ElemCoor(2,:) = dcorvg(2,ElemInd(:))
    ElemCoor(3,:) = dcorvg(3,ElemInd(:))

    CALL TransformPoints()

    ! Coordinates of the element mid point
    PX = ElemCoor(1,27)
    PY = ElemCoor(2,27)
    PZ = ElemCoor(3,27)

    ! Compute the element radius
    CALL GetElemRad()

    CALL GetDistToSTL(PX,PY,PZ,iSTL,d,.true.)

    ! Set the distance in the element mid point
    dVaux(ElemInd(27))= d

    ! Do not divide by zero
    IF (d .ne. 0d0)THEN
      ElemSign = d/abs(d)
    ELSE
      ElemSign = -1d0
    END IF

    IF (abs(d).lt.dLocEpsDist) THEN 

      ! If the distance to the surface is greater than
      ! the element radius, we do not need to compute
      ! the FBM-function for the other dofs
      IF (abs(d).gt.ElemRad) THEN

        DO i=1,26
        PX = ElemCoor(1,i)
        PY = ElemCoor(2,i)
        PZ = ElemCoor(3,i)

        ! check if dist is already computed
        IF (dVaux(ElemInd(i)).eq.0d0) THEN
          CALL GetDistToSTL(PX,PY,PZ,iSTL,d,.false.)
          ! Set distance to minium of eps and d
          dVaux(ElemInd(i))= ElemSign*MIN(ABS(d),dLocEpsDist)
        END IF
        END DO

      ELSE ! abs(d).gt.ElemRad

        DO i=1,26
        PX = ElemCoor(1,i)
        PY = ElemCoor(2,i)
        PZ = ElemCoor(3,i)

        IF (dVaux(ElemInd(i)).eq.0d0) THEN

          CALL GetDistToSTL(PX,PY,PZ,iSTL,d,.true.)

          IF (d.ne.0d0)THEN
            PointSign = d/abs(d)
          ELSE
            PointSign  = -1d0
          END IF

          ! Set distance to minium of eps and d
          dVaux(ElemInd(i))= PointSign*MIN(ABS(d),dLocEpsDist)

        END IF
        END DO
      END IF

    ELSE ! abs(d).lt.dLocEpsDist

      ! The distance to element the mid point
      ! is greater than eps: we can set
      ! all other dofs to eps too

      DO i=1,27

      IF (dVaux(ElemInd(i)).eq.0d0) THEN
        dVaux(ElemInd(i))= ElemSign*dLocEpsDist
      end if

      end do

    end if ! abs(d).lt.dLocEpsDist

    end do

    if (adjustl(trim(mySigma%mySegment(iSeg)%ObjectType)).eq."DIE") THEN
      do i=1,ndof
      dst1(i) = min(dst1(i),-dVaux(i)/dUnitScale)
      end do
    end if  

    if (adjustl(trim(mySigma%mySegment(iSeg)%ObjectType)).eq."SCREW") THEN
      do i=1,ndof
      dst2(i) = min(dst2(i),+dVaux(i)/dUnitScale)
      end do
    end if  

    if (adjustl(trim(mySigma%mySegment(iSeg)%ObjectType)).eq."OBSTACLE") THEN
      do i=1,ndof
      dst1(i) = min(dst1(i),+dVaux(i)/dUnitScale)
      end do
    end if  

    end do

  end if
  end do

  DO i=1,ndof
   Shell(i) = dst1(i)
   IF (Shell(i).le.0d0) THEN
    MixerKNPR(i) = 100
   END IF

   Screw(i) = dst2(i)
   IF (Screw(i).le.0d0) THEN
    MixerKNPR(i) = 104
   END IF
  END DO

contains
subroutine GetElemRad
  REAL*8 daux

  ElemRad = 0d0

  do i=1,8
  dist = (ElemCoor(1,i)-PX)**2d0 + (ElemCoor(2,i)-PY)**2d0 + (ElemCoor(3,i)-PZ)**2d0
  ElemRad = max(dist,ElemRad)
  end do
  ElemRad = 2.1d0*SQRT(ElemRad)

End SUBROUTINE GetElemRad

SUBROUTINE TransformPoints()
  REAL*8 t,XB,YB,ZB,XT,YT,ZT,dAlpha
  REAL*8 :: myPI = dATAN(1d0)*4d0
  INTEGER iP

  t = (dRotAngle/360d0)/(myProcess%Umdr/60d0)

  DO iP=1,27

  ! Point needs to be offseted to its segment position
  XB = dUnitScale*ElemCoor(1,iP)
  YB = dUnitScale*ElemCoor(2,iP) + dYtShift
  ZB = dUnitScale*ElemCoor(3,iP) - mySigma%mySegment(iSeg)%Min

  ! Point needs to be transformed back to a time=0 due to the rotation
  dAlpha = 0d0- t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
  XT = XB*cos(dAlpha) - YB*sin(dAlpha)
  YT = XB*sin(dAlpha) + YB*cos(dAlpha)
  ZT = ZB

  YT = YT - dYcShift

  ElemCoor(1,iP) = XT
  ElemCoor(2,iP) = YT
  ElemCoor(3,iP) = ZT

  END DO

  END 

SUBROUTINE TransformPoint(x, y, z, xb, yb, zb)
  REAL*8 t,XB,YB,ZB,X,Y,Z,dAlpha
  REAL*8 :: myPI = dATAN(1d0)*4d0
  INTEGER iP

  t = (myProcess%Angle/360d0)/(myProcess%Umdr/60d0)

  ! Point needs to be offseted to its segment position
  XB = dUnitScale * x
  YB = dUnitScale * y
  ZB = dUnitScale * z - dUnitScale*mySigma%mySegment(1)%Min

  ! Point needs to be transformed back to a time=0 due to the rotation
  dAlpha = 0d0- t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
  XB = XB*cos(dAlpha) - YB*sin(dAlpha)
  YB = XB*sin(dAlpha) + YB*cos(dAlpha)
  ZB = ZB

  END 

end subroutine calcDistanceFunction_netzsch
!------------------------------------------------------------------------------------------------
! The subroutine calculates a distance function on the mesh given 
! by dcorvg, etc. The distance is computed to the list of screw geometries
! maintained in the global structure mySigma. This function for distance computations
! is intented for screw setups only.
!------------------------------------------------------------------------------------------------
subroutine calcDistanceFunction_heat(dcorvg,kvert,kedge,karea,nel,nvt,nat,net,dst1,dst2,dst3,dVaux)

  real*8, dimension(:,:), intent(inout) :: dcorvg

  ! output array for the distance to a die or obstacles
  real*8, dimension(:), intent(inout) :: dst1

  ! output array for the distance to the screw
  real*8, dimension(:), intent(inout) :: dst2

  ! output array for the distance to flowchannel
  real*8, dimension(:), intent(inout) :: dst3

  real*8, dimension(:), intent(inout) :: dVaux

  integer :: nel,nvt,nat,net

  integer, dimension(:,:), intent(inout) :: karea, kedge, kvert

  ! local variables
  integer ElemInd(27)

  real*8 ElemCoor(3,27),ElemDist(27),ElemRad,ElemSign,PointSign,dist

  real*8 d,PX,PY,PZ,dS,dUnitScale,dLocEpsDist,dRotAngle,dAuxDist

  real*8, dimension(3) :: point

  integer :: i,j,ndof,iel,ipc

  integer iSTL,iSeg,iFile

  ndof = nel+nvt+nat+net

  ! Loop over all segments of the screw
  DO iSeg=1,mySigma%NumberOfSeg

  ! If we find a discrete geometry segment
  IF (adjustl(trim(mySigma%mySegment(iSeg)%ART)).eq."STL") THEN

    IF (mySigma%mySegment(iSeg)%Unit.eq.'MM') dUnitScale = 1d+1
    IF (mySigma%mySegment(iSeg)%Unit.eq.'CM') dUnitScale = 1d0
    IF (mySigma%mySegment(iSeg)%Unit.eq.'DM') dUnitScale = 1d-1

    dLocEpsDist = dUnitScale*dEpsDist

    DO iFile=1,mySigma%mySegment(iSeg)%nOFFfiles

    dVaux(1:ndof) = 0d0  
    iSTL = mySigma%mySegment(iSeg)%idxCgal(iFile)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! iSTL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO iel=1,nel

    ! Indices of the vertices
    ElemInd( 1:8 )  = kvert(1:8,iel)

    ! Indices of the edge mid points
    ElemInd( 9:20)  = nvt+kedge(1:12,iel)

    ! Indices of the face mid points
    ElemInd(21:26)  = nvt+net+karea(1:6,iel)

    ! Index of the element mid point
    ElemInd(   27)  = nvt+net+nat+iel

    ! Coordinates of the element dofs
    ElemCoor(1,:) = dcorvg(1,ElemInd(:))
    ElemCoor(2,:) = dcorvg(2,ElemInd(:))
    ElemCoor(3,:) = dcorvg(3,ElemInd(:))

    CALL TransformPoints()

    ! Coordinates of the element mid point
    PX = ElemCoor(1,27)
    PY = ElemCoor(2,27)
    PZ = ElemCoor(3,27)

    ! Compute the element radius
    CALL GetElemRad()

    CALL GetDistToSTL(PX,PY,PZ,iSTL,d,.true.)

    ! Set the distance in the element mid point
    dVaux(ElemInd(27))= d

    ! Do not divide by zero
    IF (d .ne. 0d0)THEN
      ElemSign = d/abs(d)
    ELSE
      ElemSign = -1d0
    END IF

    IF (abs(d).lt.dLocEpsDist) THEN 

      ! If the distance to the surface is greater than
      ! the element radius, we do not need to compute
      ! the FBM-function for the other dofs
      IF (abs(d).gt.ElemRad) THEN

        DO i=1,26
        PX = ElemCoor(1,i)
        PY = ElemCoor(2,i)
        PZ = ElemCoor(3,i)

        ! check if dist is already computed
        IF (dVaux(ElemInd(i)).eq.0d0) THEN
          CALL GetDistToSTL(PX,PY,PZ,iSTL,d,.false.)
          ! Set distance to minium of eps and d
          dVaux(ElemInd(i))= ElemSign*MIN(ABS(d),dLocEpsDist)
        END IF
        END DO

      ELSE ! abs(d).gt.ElemRad

        DO i=1,26
        PX = ElemCoor(1,i)
        PY = ElemCoor(2,i)
        PZ = ElemCoor(3,i)

        IF (dVaux(ElemInd(i)).eq.0d0) THEN

          CALL GetDistToSTL(PX,PY,PZ,iSTL,d,.true.)

          IF (d.ne.0d0)THEN
            PointSign = d/abs(d)
          ELSE
            PointSign  = -1d0
          END IF

          ! Set distance to minium of eps and d
          dVaux(ElemInd(i))= PointSign*MIN(ABS(d),dLocEpsDist)

        END IF
        END DO
      END IF

    ELSE ! abs(d).lt.dLocEpsDist

      ! The distance to element the mid point
      ! is greater than eps: we can set
      ! all other dofs to eps too

      DO i=1,27

      IF (dVaux(ElemInd(i)).eq.0d0) THEN
        dVaux(ElemInd(i))= ElemSign*dLocEpsDist
      end if

      end do

    end if ! abs(d).lt.dLocEpsDist

    end do

    do i=1,ndof
      IF (dVaux(i).lt.0d0) myHeatObjects%Segment(i) = iSeg
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! iSTL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (adjustl(trim(mySigma%mySegment(iSeg)%ObjectType)).eq."BLOCK") THEN
      do i=1,ndof
       dst1(i) = min(dst1(i),dVaux(i)/dUnitScale)
      end do
    end if  

    if (adjustl(trim(mySigma%mySegment(iSeg)%ObjectType)).eq."WIRE") THEN
      do i=1,ndof
      dst2(i) = min(dst2(i),dVaux(i)/dUnitScale)
      end do
    end if  

    if (adjustl(trim(mySigma%mySegment(iSeg)%ObjectType)).eq."MELT") THEN
      do i=1,ndof
       dst3(i) = min(dst3(i),dVaux(i)/dUnitScale)
      end do
    end if  
    
    end do

  end if
  end do

  DO i=1,ndof
    myHeatObjects%Block(i) = -dst1(i)
!    IF (Shell(i).le.0d0) THEN
! !     MixerKNPR(i) = 100
!    END IF

    myHeatObjects%Wire(i) = -dst2(i)
!    IF (Screw(i).le.0d0) THEN
! !     MixerKNPR(i) = 103
!    END IF

    myHeatObjects%Channel(i) = -dst3(i)
!     if (myHeatObjects%Channel(i).gt.0d0) then
!      myHeatObjects%Channel(i) = max(myHeatObjects%Channel(i),myHeatObjects%Block(i))
!     end if

 END DO

contains
subroutine GetElemRad
  REAL*8 daux

  ElemRad = 0d0

  do i=1,8
  dist = (ElemCoor(1,i)-PX)**2d0 + (ElemCoor(2,i)-PY)**2d0 + (ElemCoor(3,i)-PZ)**2d0
  ElemRad = max(dist,ElemRad)
  end do
  ElemRad = 2.1d0*SQRT(ElemRad)

End SUBROUTINE GetElemRad

SUBROUTINE TransformPoints()
  REAL*8 :: myPI = dATAN(1d0)*4d0
  INTEGER iP

  DO iP=1,27

   ElemCoor(1,iP) = dUnitScale*ElemCoor(1,iP)
   ElemCoor(2,iP) = dUnitScale*ElemCoor(2,iP)
   ElemCoor(3,iP) = dUnitScale*ElemCoor(3,iP)

  END DO

  END 

end subroutine calcDistanceFunction_heat
!------------------------------------------------------------------------------------------------
! Subroutine that calculates the distance of point px,py,pz to a single geometry
! If neccessary also the FBM function is computed
!------------------------------------------------------------------------------------------------
SUBROUTINE GetDistToSTL(Px,Py,Pz,iSTL,dist,bRay)
implicit none
LOGICAL bRay
INTEGER iSTL
REAL*8 Px,Py,Pz,dist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
REAL*8,dimension(3) :: P

REAL*8 :: dS

integer :: isin, ipc 


! P=[1.5461,-1.8252,17.412]
P=[Px,Py,Pz]

! If bRay is set to true then
! we also compute the FBM-function for
! this point

if (bRay) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! Determine the sign of the resulting distance !!!!!!!!!!!!!!!!!

  isin = 0
  ipc = iSTL-1
  ! This routine computes the FBM-function 
  ! and sets the parameter isin to 1 if
  ! the point is inside of the geometry
  call isinelementid(P(1),P(2),P(3),ipc,isin)

  ! If the point is outside of the
  ! geometry we set dS to -1
  if(isin .gt. 0)then
    dS = -1d0
  else
    dS = +1d0
  end if

else

  dS = -1d0

end if

!!!!!!!!!!!!!!!!! Determine the sign of the resulting distance !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ipc = iSTL-1
call getdistanceid(P(1),P(2),P(3),dist,ipc);        

! The signed distance is the distance to the geometry * dS * geometry_modifier
! dist = mySTL(iSTL)%dS * dS * dist
dist = dS * dist

END SUBROUTINE GetDistToSTL

end module geometry_processing
