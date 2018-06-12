! An optional callback function that allows
! the user to customize the handling of fictitious boundary velocity 
! boundary condition update.
interface 
  subroutine fbm_velBC_handler(x,y,z,valu,valv,valw,ip,t)

  use var_QuadScalar, only : myFBM,bRefFrame
  implicit none

  ! The FBM value of the boundary vertex
  integer, intent(in) :: ip

  ! The coordinates of the boundary vertex
  ! and the simulation time
  real*8 , intent(in) :: x,y,z,t

  ! The velocitiy values of the boundary vertex
  real*8 , intent(inout) :: valu,valv,valw

  end subroutine
end interface
