! An optional callback function that allows
! the user to customize the handling of the FBM geom update.
interface 
  subroutine fbm_geom_handler(px, py, pz, bndryId, fictId, dist)

    implicit none
   
    ! Coordinates of the query point 
    real*8, intent(in) :: px, py, pz 

    ! Id of the boundary component
    integer, intent(inout) :: bndryId

    ! fictId
    integer(kind=16), intent(inout) :: fictId

    ! Distance solution in the query point 
    real*8, intent(inout) :: dist 

  end subroutine
end interface
