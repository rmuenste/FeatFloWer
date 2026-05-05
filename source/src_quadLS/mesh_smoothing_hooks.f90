module mesh_smoothing_hooks
!-------------------------------------------------------------------------------------------------
! Thin app-initialization hooks for optional mesh smoothing features.
!-------------------------------------------------------------------------------------------------
  use attraction_smoother, only : tAttractionConfig, AttractionSmoother_GetConfig, &
                                  AttractionSmoother_Apply
  use PP3D_MPI, only : myid
  use var_QuadScalar, only : tMultiMesh

  implicit none

  private
  public :: ApplyAttractionSmootherHook

contains

  subroutine ApplyAttractionSmootherHook(mgMesh, ilev)
    implicit none
    type(tMultiMesh), intent(inout) :: mgMesh
    integer, intent(in) :: ilev

    interface
      subroutine UmbrellaSmoother_STRCT(myTime, nSteps)
        implicit none
        real*8 :: myTime
        integer :: nSteps
      end subroutine UmbrellaSmoother_STRCT
    end interface

    type(tAttractionConfig) :: cfg
    integer :: iUmbrella, iApplyLev

    call AttractionSmoother_GetConfig(cfg)
    if (.not. cfg%enabled) return
    if (trim(adjustl(cfg%mode)) .eq. "distance_weight_only") return

    iApplyLev = ilev
    if (cfg%applyOnLevel .gt. 0) iApplyLev = cfg%applyOnLevel

    if (myid .ne. 0) then
      call AttractionSmoother_Apply(mgMesh, iApplyLev, cfg)
    end if

    do iUmbrella = 1, cfg%nPostUmbrellaSteps
      call UmbrellaSmoother_STRCT(0d0, 1)
    end do

  end subroutine ApplyAttractionSmootherHook

end module mesh_smoothing_hooks
