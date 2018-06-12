module types
!-------------------------------------------------------------------------------------------------
! A module that contains several variants of the Laplacian 
! smoother 'Umbrella'. Besides the standard version of this 
! smoother type, versions with user-defined weighting functions
! are available
!-------------------------------------------------------------------------------------------------
! No implicit variables in this module
implicit none

type tBoundingBox
  double precision, dimension(3,2)  :: vertices
end type tBoundingBox

type tmgBoundingBox
  type(tBoundingBox), dimension(:), allocatable :: bb
end type tmgBoundingBox

type tPostprocessingParams
  real*8 :: U_mean = 0.2d0
  ! This is the setting for 2D FAC, for
  ! the full 3D FAC the value is H = 0.205d0
  real*8 :: H = 0.05d0
  real*8 :: D = 0.1d0
  real*8 :: Sc_U = 1d0
  real*8 :: Sc_Mu = 1d0
  real*8 :: Sc_a = 1d0
end type tPostprocessingParams

end module types
