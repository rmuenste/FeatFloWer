module dem_query
!================================================================================================
! Module USE
!================================================================================================
use iso_c_binding

implicit none
!================================================================================================
! Variables and Types
!================================================================================================
type, bind(C) :: tParticleData
  real(c_double), dimension(3)  :: position
  real(c_double), dimension(3)  :: velocity
  real(c_double), dimension(3)  :: force = (/0d0, 0d0, 0d0/)
  real(c_double), dimension(3)  :: torque = (/0d0, 0d0, 0d0/)
  real(c_double) :: time
  integer(c_int) :: localIdx
  integer(c_int) :: uniqueIdx
  integer(c_int) :: systemIdx
  integer(c_short), dimension(8) :: bytes
end type tParticleData

!================================================================================================
! Interfaces
!================================================================================================
interface
integer(c_int) function getNumParticles() bind(C, name="getNumParticles")
  use iso_c_binding, only: c_int
  end function
end interface

interface
integer(c_int) function getNumRemParticles() bind(C, name="getNumRemParticles")
  use iso_c_binding, only: c_int
  end function
end interface

interface
subroutine getParticle(idx, lidx, uidx, time, pos, vel) bind(C, name="getParticle")
  use iso_c_binding, only: c_int, c_double
  integer(c_int), value :: idx
  integer(c_int)        :: lidx
  integer(c_int)        :: uidx
  real(c_double)        :: time
  real(c_double)        :: pos(*)
  real(c_double)        :: vel(*)
  end subroutine
end interface

interface
subroutine getParticle2(idx, particle) bind(C, name="getParticle2")
  use iso_c_binding, only: c_int, c_double
  import tParticleData
  integer(c_int), value :: idx
  type(tParticleData) :: particle
  end subroutine
end interface

interface
subroutine getRemoteParticle(idx, lidx, uidx, time, pos, vel) bind(C, name="getRemoteParticle")
  use iso_c_binding, only: c_int, c_double
  integer(c_int) :: idx
  integer(c_int) :: lidx
  integer(c_int) :: uidx
  real(c_double) :: time
  real(c_double) :: pos(*)
  real(c_double) :: vel(*)
  end subroutine
end interface

interface
subroutine getRemoteParticle2(idx, particle) bind(C, name="getRemoteParticle2")
  use iso_c_binding, only: c_int, c_double
  import tParticleData
  integer(c_int), value :: idx
  type(tParticleData) :: particle
  end subroutine
end interface

interface
subroutine setParticle(idx, lidx, uidx, time, pos, vel) bind(C, name="setParticle")
  use iso_c_binding, only: c_int, c_double
  integer(c_int) :: idx
  integer(c_int) :: lidx
  integer(c_int) :: uidx
  real(c_double) :: time
  real(c_double) :: pos(*)
  real(c_double) :: vel(*)
  end subroutine
end interface

interface
subroutine setForces(idx, lidx, uidx, force, torque) bind(C, name="setForces")
  use iso_c_binding, only: c_int, c_double
  integer(c_int) :: idx
  integer(c_int) :: lidx
  integer(c_int) :: uidx
  real(c_double) :: force(*)
  real(c_double) :: torque(*)
  end subroutine
end interface

interface
subroutine setRemoteForces(idx, lidx, uidx, force, torque) bind(C, name="setRemoteForces")
  use iso_c_binding, only: c_int, c_double
  integer(c_int) :: idx
  integer(c_int) :: lidx
  integer(c_int) :: uidx
  real(c_double) :: force(*)
  real(c_double) :: torque(*)
  end subroutine
end interface

interface
logical(c_bool) function isSphere(idx) bind(C, name="isTypeSphere")
  use iso_c_binding, only: c_int, c_bool
  integer(c_int) :: idx
  end function
end interface

interface
real(c_double) function getParticleRadius(idx) bind(C, name="getParticleRadius")
  use iso_c_binding, only: c_int, c_double
  integer(c_int) :: idx
  end function
end interface

interface
logical(c_bool) function pointInsideObject(idx, pos) bind(C, name="pointInsideObject")
  use iso_c_binding, only: c_int, c_bool, c_double
  integer(c_int), value :: idx
  real(c_double) :: pos(*)
  end function
end interface

interface
logical(c_bool) function pointInsideRemObject(idx, pos) bind(C, name="pointInsideRemObject")
  use iso_c_binding, only: c_int, c_bool, c_double
  integer(c_int), value :: idx
  real(c_double) :: pos(*)
  end function
end interface

interface
logical(c_bool) function check_rem_id(fbmid, id) bind(C, name="check_rem_id")
  use iso_c_binding, only: c_int, c_bool, c_double
  integer(c_int), value :: fbmid
  integer(c_int), value :: id
  end function
end interface

interface
subroutine rem_particles_index_map(idxMap) bind(C, name="rem_particles_index_map")
  use iso_c_binding, only: c_int, c_double
  integer(c_int) :: idxMap(*)
  end subroutine
end interface

interface
subroutine particles_index_map(idxMap) bind(C, name="particles_index_map")
  use iso_c_binding, only: c_int, c_double
  integer(c_int) :: idxMap(*)
  end subroutine
end interface

interface
subroutine get_bytes(bytes) bind(C, name="get_bytes")
  use iso_c_binding, only: c_short
  integer(c_short), dimension(8) :: bytes 
end subroutine
end interface

interface
logical(c_bool) function map_local_to_system(lidx, vidx) bind(C, name="map_local_to_system")
  use iso_c_binding, only: c_int, c_bool 
  integer(c_int), value :: lidx
  integer(c_int), value :: vidx
  end function
end interface

interface
logical(c_bool) function map_local_to_system2(lidx, vidx) bind(C, name="map_local_to_system2")
  use iso_c_binding, only: c_int, c_bool 
  integer(c_int), value :: lidx
  integer(c_int), value :: vidx
  end function
end interface

!integer :: numLocalParticles

contains

!================================================================================================
!                              Function numLocalParticles
!================================================================================================

integer function numLocalParticles()
  implicit none

  numLocalParticles = getNumParticles()

end function numLocalParticles

!================================================================================================
!                              Function numLocalParticles
!================================================================================================

integer function numRemParticles()
  implicit none

  numRemParticles = getNumRemParticles()

end function numRemParticles
!================================================================================================
!                              Function numLocalParticles
!================================================================================================

integer function numTotalParticles()
  implicit none

  numTotalParticles = numRemParticles() + numLocalParticles()

end function numTotalParticles
!================================================================================================
!                              Function getParticlesIndexMap
!================================================================================================

subroutine getParticlesIndexMap(indexMap)
  implicit none
  integer, dimension(:), intent(inout) :: indexMap
  
  call particles_index_map(indexMap)

end subroutine getParticlesIndexMap
!================================================================================================
!                              Function getRemoteParticlesIndexMap
!================================================================================================

subroutine getRemoteParticlesIndexMap(indexMap)
  implicit none
  integer, dimension(:), intent(inout) :: indexMap
  
  call rem_particles_index_map(indexMap)

end subroutine getRemoteParticlesIndexMap
!================================================================================================
!                              Function getAllParticles
!================================================================================================

subroutine getAllParticles(theParticles)
  implicit none
  type(tParticleData), dimension(:), intent(inout) :: theParticles

  type(tParticleData) :: temp
  integer, allocatable, dimension(:) :: indexMap
  integer :: a,i,idx
  
  a = size(theParticles) 

  allocate(indexMap(a))

  call getParticlesIndexMap(indexMap)
  
  do i=1,a
    idx = indexMap(i)
    call getParticle2(idx,&
                      temp)

    theParticles(i) = temp
  end do

end subroutine getAllParticles
!================================================================================================
!                              Function getAllRemoteParticles
!================================================================================================

subroutine getAllRemoteParticles(theParticles)
  implicit none
  type(tParticleData), dimension(:), intent(inout) :: theParticles

  type(tParticleData) :: temp
  integer, allocatable, dimension(:) :: indexMap
  integer :: a,i,idx
  
  a = size(theParticles) 

  allocate(indexMap(a))

  call getRemoteParticlesIndexMap(indexMap)
!  write(*,*)'IndexMap>', indexMap(:)
  
  do i=1,a
    idx = indexMap(i)
    call getRemoteParticle2(idx, temp)  
!    call getRemoteParticle(idx,&
!                           temp%localIdx,& 
!                           temp%uniqueIdx,& 
!                           temp%time,& 
!                           temp%position,& 
!                           temp%velocity)

    theParticles(i) = temp
  end do

end subroutine getAllRemoteParticles
!================================================================================================
!                              Function setForcesMapped
!================================================================================================

subroutine setForcesMapped(particle)
  implicit none

  type(tParticleData), intent(inout) :: particle
  
  call setForces(particle%uniqueIdx,&
                 particle%localIdx,& 
                 particle%uniqueIdx,& 
                 particle%force,& 
                 particle%torque)

end subroutine setForcesMapped
!================================================================================================
!                              Function setRemoteForcesMapped
!================================================================================================

subroutine setRemoteForcesMapped(particle)
  implicit none

  type(tParticleData), intent(inout) :: particle
  
  call setRemoteForces(particle%uniqueIdx,&
                       particle%localIdx,& 
                       particle%uniqueIdx,& 
                       particle%force,& 
                       particle%torque)

end subroutine setRemoteForcesMapped
!================================================================================================
!                              Function getLocalParticle
!================================================================================================

type(tParticleData) function getLocalParticle(idx)
  implicit none
  integer, intent(in) :: idx

  type(tParticleData) :: temp
  integer :: a
  
  a = numLocalParticles()
  
  call getParticle(idx,&
                   temp%localIdx,& 
                   temp%uniqueIdx,& 
                   temp%time,& 
                   temp%position,& 
                   temp%velocity)

  getLocalParticle = temp

end function getLocalParticle
!================================================================================================
!                              Function getLocalParticle2
!================================================================================================

type(tParticleData) function getLocalParticle2(idx)
  implicit none
  integer, intent(in) :: idx

  type(tParticleData) :: temp
  integer :: a
  
  a = numLocalParticles()
  
  call getParticle2(idx, temp)

  getLocalParticle2 = temp

end function getLocalParticle2

!================================================================================================
!                              Function getRadius
!================================================================================================

double precision function getRadius(idx)
  implicit none
  integer, intent(in) :: idx

  getRadius = getParticleRadius(idx)

end function getRadius

!================================================================================================
!                              Function objectContainsPoint
!================================================================================================

logical function objectContainsPoint(idx, point)
  implicit none
  integer, intent(in) :: idx
  double precision, dimension(3), intent(inout) :: point

  objectContainsPoint = pointInsideObject(idx, point) 

end function objectContainsPoint

!================================================================================================
!                              Function remObjectContainsPoint
!================================================================================================

logical function remObjectContainsPoint(idx, point)
  implicit none
  integer, intent(in) :: idx
  double precision, dimension(3), intent(inout) :: point

  remObjectContainsPoint = pointInsideRemObject(idx, point) 

end function remObjectContainsPoint

!================================================================================================
!                              Function isSphereType
!================================================================================================

logical function isSphereType(idx)
  implicit none
  integer, intent(in) :: idx

  isSphereType = isSphere(idx)

end function isSphereType

!================================================================================================
!                              Subroutine testParticleGet 
!================================================================================================

subroutine testParticleGet(idx)
  implicit none

  type(tParticleData) :: temp

  integer,intent(in), optional :: idx

  integer :: psize, idx2
  
  if(.not. present(idx))then
    idx2 = 0
  else
    idx2 = idx
  end if

  psize = numLocalParticles()


  if (isSphereType(idx2))then
    temp = getLocalParticle2(idx2)
    write(*,*)"Particle(",idx2,") pos:",temp%position, " vel:", temp%velocity
  else
    write(*,*)"Particle(",idx2,") is not a sphere. "
  end if


end subroutine testParticleGet

!================================================================================================
!                              Subroutine testParticleRadius 
!================================================================================================

subroutine testParticleRadius()
  implicit none

  double precision :: temp

  integer :: idx, psize
  
  idx = 5

  temp = getRadius(idx)

end subroutine testParticleRadius
  
!================================================================================================
!                              Subroutine testMapParticles 
!================================================================================================

subroutine testMapParticles()
  implicit none

  
  call map_particles()

end subroutine testMapParticles

!================================================================================================
!                              checkLongId
!================================================================================================

logical function longIdMatch(idx, longFictId)
  use iso_c_binding, only: c_short
  use var_QuadScalar, ONLY : FictKNPR_uint64
  implicit none
  integer, intent(in) :: idx
  integer(c_short), dimension(8) :: longFictId

  longIdMatch = .false.

  if( (FictKNPR_uint64(idx)%bytes(1) .eq. longFictId(1)) .and. &
      (FictKNPR_uint64(idx)%bytes(2) .eq. longFictId(2)) .and. &
      (FictKNPR_uint64(idx)%bytes(3) .eq. longFictId(3)) .and. &
      (FictKNPR_uint64(idx)%bytes(4) .eq. longFictId(4)) .and. &
      (FictKNPR_uint64(idx)%bytes(5) .eq. longFictId(5)) .and. &
      (FictKNPR_uint64(idx)%bytes(6) .eq. longFictId(6)) .and. &
      (FictKNPR_uint64(idx)%bytes(7) .eq. longFictId(7)) .and. &
      (FictKNPR_uint64(idx)%bytes(8) .eq. longFictId(8)) )then
      longIdMatch = .true.
  end if

end function longIdMatch

!================================================================================================
!                              Function convertSystemId
!================================================================================================

!integer(kind=16) function convertSystemId(bytes)
!  use iso_c_binding, only: c_short
!  implicit none
!  integer(c_short), dimension(8), intent(inout) :: bytes 
!  ! locals
!  integer(kind = 16) :: result, temp 
!  integer :: i
!
!
!  result = huge(temp)
!!!  call get_bytes(bytes)
!!  do i = 1, 8
!!    result = ishft(result, 8)
!!    temp = int(bytes(i), 16)
!!    result = ior(result, temp)
!!  end do
!!
!  convertSystemId = result
!
!end function convertSystemId

end module dem_query

