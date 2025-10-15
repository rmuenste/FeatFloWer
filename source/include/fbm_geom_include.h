! An optional callback function that allows
! the user to customize the handling of the FBM geom update.
interface 
  subroutine fbm_geom_handler(px, py, pz, bndryId, fictId, dist, vidx, longFictId)
    use types
    use iso_c_binding, only: c_short

    implicit none
   
    ! Coordinates of the query point 
    real*8, intent(in) :: px, py, pz 

    ! Id of the boundary component
    integer, intent(inout) :: bndryId

    ! fictId
    integer, intent(inout) :: fictId

    ! Distance solution in the query point 
    real*8, intent(inout) :: dist 

    ! vidx
    integer, intent(in) :: vidx

    ! long FictId
    type(tUint64), intent(inout) :: longFictId

  end subroutine
end interface
