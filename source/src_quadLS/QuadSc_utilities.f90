!=========================================================================
! QuadSc_utilities.f90
!
! Miscellaneous utility functions for QuadScalar module
! Extracted from QuadSc_main.f90 for better code organization
!=========================================================================
!
!=========================================================================
SUBROUTINE TESTER(DX,DXP)
  REAL*8 DX(*),DXP(*)

  CALL GetParPressure(DX,DXP)

END SUBROUTINE TESTER
!=========================================================================
!
!=========================================================================
SUBROUTINE Analyzer
  INTEGER I,J

  J=0
  DO I=1,QuadSc%ndof
  IF (ABS(QuadSC%auxU(I)).GT.1d-10) THEN
    J = J + 1
    !   WRITE(*,*) I,J,QuadSC%auxU(I)
  END IF
  END DO

END SUBROUTINE Analyzer
!=========================================================================
!
!=========================================================================
SUBROUTINE  GetMonitor()
  implicit none
  INTEGER i
  REAL*8 daux,px,py,pz
  return
  ILEV = NLMAX
  CALL SETLEV(2)

  !
  !
  !  DO i=1,nvt
  !
  !    PX = dcorvg(1,i)
  !    PY = dcorvg(2,i)
  !    PZ = dcorvg(3,i)
  !    getdistanceid(px,py,pz,daux,i);
  !
  !   myALE%Monitor(i) = sqrt(daux) !HogenPowerlaw(daux)
  !
  !  END DO

END SUBROUTINE  GetMonitor
!=========================================================================
!
!=========================================================================
SUBROUTINE DetermineIfGoalsWereReached(bGoalsReached)
use, intrinsic :: ieee_arithmetic
REAL*8 myinf
LOGICAL bGoalsReached

if(ieee_support_inf(myInf))then
  myInf = ieee_value(myInf, ieee_negative_inf)
endif

bGoalsReached = .true.

! IF (ADJUSTL(TRIM(mySigma%cType)).EQ."TSE".or.ADJUSTL(TRIM(mySigma%cType)).EQ."SSE".or.ADJUSTL(TRIM(mySigma%cType)).EQ."XSE") THEN
!  if (myProcess%FillingDegree.eq.myInf .or. myProcess%FillingDegree .eq. 1d0) then
!     if (itns.ge.nitns) bGoalsReached=.false.
!  end if
! END IF

IF (ADJUSTL(TRIM(mySigma%cType)).EQ."DIE") THEN
 if (istart.eq.1) then
   if (itns.ge.nitns) bGoalsReached=.false.
 end if
END IF


if (.not.bGoalsReached) THEN
 if (myid.eq.1) write(*,*) 'max time steps have been reached // the simulation has - most probably - not converged! '
end if

END SUBROUTINE DetermineIfGoalsWereReached
