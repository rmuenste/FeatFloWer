!=========================================================================
! QuadSc_utilities.f90
!
! Miscellaneous utility functions for QuadScalar module
! Extracted from QuadSc_main.f90 for better code organization
!=========================================================================

!=========================================================================
! TESTER - Test subroutine for parallel pressure calculation
!=========================================================================
subroutine TESTER(DX, DXP)
  implicit none
  real(8), intent(inout) :: DX(*)
  real(8), intent(inout) :: DXP(*)

  call GetParPressure(DX, DXP)

end subroutine TESTER

!=========================================================================
! Analyzer - Analyze auxiliary velocity field
! Counts non-zero entries in QuadSc%auxU array
!=========================================================================
subroutine Analyzer
  implicit none
  integer :: i, j

  j = 0
  do i = 1, QuadSc%ndof
    if (abs(QuadSC%auxU(i)) > 1.0d-10) then
      j = j + 1
    end if
  end do

end subroutine Analyzer

!=========================================================================
! GetMonitor - Get monitor function for mesh adaptation (currently disabled)
!=========================================================================
subroutine GetMonitor()
  implicit none

  ! This subroutine is currently disabled (early return)
  ! Originally used for computing distance-based monitor function
  return

end subroutine GetMonitor

!=========================================================================
! DetermineIfGoalsWereReached - Check if simulation goals are reached
! Determines if the simulation should continue or terminate
!=========================================================================
subroutine DetermineIfGoalsWereReached(bGoalsReached)
  use, intrinsic :: ieee_arithmetic
  implicit none
  logical, intent(out) :: bGoalsReached
  real(8) :: myInf

  ! Initialize infinity value if supported
  if (ieee_support_inf(myInf)) then
    myInf = ieee_value(myInf, ieee_negative_inf)
  end if

  bGoalsReached = .true.

  ! Check for DIE simulation type
  if (adjustl(trim(mySigma%cType)) == "DIE") then
    if (istart == 1) then
      if (itns >= nitns) bGoalsReached = .false.
    end if
  end if

  ! Report if goals were not reached
  if (.not. bGoalsReached) then
    if (myid == 1) then
      write (*, *) 'Maximum time steps reached - simulation may not have converged!'
    end if
  end if

end subroutine DetermineIfGoalsWereReached
