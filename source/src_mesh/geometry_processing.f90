module geometry_processing
  use Sigma_User, only: mySigma, myProcess
  use var_QuadScalar, only: distance
  !------------------------------------------------------------------------------------------------
  ! A module for functions that perform operations on the 
  ! geometry immersed into the simulation domain or boundary geometry
  !------------------------------------------------------------------------------------------------
  implicit none

  ! a variable for counting the outputs
  integer :: ifile = 0

  real*8 :: dEpsDist

contains
!------------------------------------------------------------------------------------------------
! The subroutine calculates a distance function on the mesh given 
! by dcorvg, etc. The distance is computed to the list of screw geometries
! maintained in the global structure mySigma. This function for distance computations
! is intented for screw setups only.
!------------------------------------------------------------------------------------------------
subroutine calcDistanceFunction(dcorvg,kvert,kedge,karea,nel,nvt,nat,net,dst1,dst2,dVaux)

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
    ELSE
     dRotAngle = 0d0
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
  ZB = dUnitScale*ElemCoor(3,iP) - dUnitScale*mySigma%mySegment(iSeg)%Min

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

end subroutine calcDistanceFunction
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
