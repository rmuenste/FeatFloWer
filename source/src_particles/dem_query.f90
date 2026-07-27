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
  real(c_double), dimension(3)  :: angvel
  real(c_double), dimension(3)  :: force = (/0d0, 0d0, 0d0/)
  real(c_double), dimension(3)  :: torque = (/0d0, 0d0, 0d0/)
  real(c_double) :: time
  real(c_double) :: density
  real(c_double), dimension(3)  :: aabb
  real(c_double) :: radius
  integer(c_int) :: typeId
  integer(c_int) :: localIdx
  integer(c_int) :: uniqueIdx
  integer(c_int) :: systemIdx
  integer(c_short), dimension(8) :: bytes
end type tParticleData

#ifdef HAVE_PE 
!================================================================================================
! Interfaces
!================================================================================================


!================================================================================================
!                              Function getTotalParticles
!================================================================================================
! C++ implementation: getTotalParts() in libs/pe/src/interface/object_queries.cpp
interface
integer(c_int) function getTotalParticles() bind(C, name="getTotalParticles")
  use iso_c_binding, only: c_int
  end function
end interface


!================================================================================================
!                              Function getNumParticles
!================================================================================================
! C++ implementation: getNumParts() in libs/pe/src/interface/object_queries.cpp
interface
integer(c_int) function getNumParticles() bind(C, name="getNumParticles")
  use iso_c_binding, only: c_int
  end function
end interface

!================================================================================================
!                              Function getNumRemParticles
!================================================================================================
! C++ implementation: getNumRemParts() in libs/pe/src/interface/object_queries.cpp
interface
integer(c_int) function getNumRemParticles() bind(C, name="getNumRemParticles")
  use iso_c_binding, only: c_int
  end function
end interface

!================================================================================================
!                              Function getElSubsteps
!================================================================================================
! C++ implementation: getElSubsteps() in libs/pe/pe/interface/c_interface_queries.h
! Returns the PE sub-step count n_sub = max(1, substeps_), so the Euler-Lagrange
! feedback can reproduce PE's sub-cycled semi-implicit drag exactly.
interface
integer(c_int) function getElSubsteps() bind(C, name="getElSubsteps")
  use iso_c_binding, only: c_int
  end function
end interface

!================================================================================================
!                              Function getElStepsize
!================================================================================================
! C++ implementation: getElStepsize() in libs/pe/pe/interface/c_interface_queries.h
interface
real(c_double) function getElStepsize() bind(C, name="getElStepsize")
  use iso_c_binding, only: c_double
  end function
end interface

!================================================================================================
!                              Function getElContactCount
!================================================================================================
! C++ implementation: getElContactCount() in libs/pe/pe/interface/c_interface_queries.h
interface
integer(c_int) function getElContactCount() bind(C, name="getElContactCount")
  use iso_c_binding, only: c_int
  end function
end interface

!================================================================================================
!                              Function getElMaxPenetration
!================================================================================================
! C++ implementation: getElMaxPenetration() in libs/pe/pe/interface/c_interface_queries.h
interface
real(c_double) function getElMaxPenetration() bind(C, name="getElMaxPenetration")
  use iso_c_binding, only: c_double
  end function
end interface

!================================================================================================
!                              Subroutine getElLubricationVirial
!================================================================================================
! C++ implementation: getElLubricationVirial() in libs/pe/pe/interface/c_interface_queries.h
! Rank-local impulse virial of the current macro step (9 components, row-major
! force x separation); MPI-SUM counts each pair exactly once (0.5 weight per
! locally-owned member). Particle-phase stress: sigma = -virial/(dt*V).
interface
subroutine getElLubricationVirial(sigma) bind(C, name="getElLubricationVirial")
  use iso_c_binding, only: c_double
  real(c_double) :: sigma(*)
  end subroutine
end interface

!================================================================================================
!                              Function getElLubricationPairs
!================================================================================================
! C++ implementation: getElLubricationPairs() in libs/pe/pe/interface/c_interface_queries.h
interface
integer(c_int) function getElLubricationPairs() bind(C, name="getElLubricationPairs")
  use iso_c_binding, only: c_int
  end function
end interface

!================================================================================================
!                              Subroutine getElLubricationImpulse
!================================================================================================
! C++ implementation: getElLubricationImpulse() in libs/pe/pe/interface/c_interface_queries.h
! Rank-local net momentum folded by the lubrication sweep this macro step; the
! MPI sum across ranks measures pair one-sidedness (must be ~0).
interface
subroutine getElLubricationImpulse(dp) bind(C, name="getElLubricationImpulse")
  use iso_c_binding, only: c_double
  real(c_double) :: dp(*)
  end subroutine
end interface

!================================================================================================
!                              Subroutine getElContactVirial
!================================================================================================
! C++ implementation: getElContactVirial() in libs/pe/pe/interface/c_interface_queries.h
! Rank-local impulse virial of the converged sphere-sphere PGS contact impulses
! this macro step (9 components, row-major, Sum p x r12); the PE contact filter
! treats each contact on exactly one rank, so MPI-SUM counts each contact once.
! Same normalization as the lubrication virial: sigma = -virial/(dt*V).
interface
subroutine getElContactVirial(sigma) bind(C, name="getElContactVirial")
  use iso_c_binding, only: c_double
  real(c_double) :: sigma(*)
  end subroutine
end interface

!================================================================================================
!                              Function getElContactVirialPairs
!================================================================================================
! C++ implementation: getElContactVirialPairs() in libs/pe/pe/interface/c_interface_queries.h
interface
integer(c_int) function getElContactVirialPairs() bind(C, name="getElContactVirialPairs")
  use iso_c_binding, only: c_int
  end function
end interface

!================================================================================================
!                              Subroutine getParticle
! C++ implementation: getObjByIdx() in libs/pe/src/interface/object_queries.cpp
!================================================================================================
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

!================================================================================================
!                              Subroutine getParticle2
! C++ implementation: getPartStructByIdx() in libs/pe/src/interface/object_queries.cpp
!================================================================================================
interface
subroutine getParticle2(idx, particle) bind(C, name="getParticle2")
  use iso_c_binding, only: c_int, c_double
  import tParticleData
  integer(c_int), value :: idx
  type(tParticleData) :: particle
  end subroutine
end interface

!================================================================================================
!                              Subroutine setParticle2
! C++ implementation: setPartStruct() in libs/pe/src/interface/object_queries.cpp
!================================================================================================
interface
subroutine setParticle2(particle) bind(C, name="setParticle2")
  use iso_c_binding, only: c_int, c_double
  import tParticleData
  type(tParticleData) :: particle
  end subroutine
end interface

!================================================================================================
!                              Subroutine getRemoteParticle
! C++ implementation: getRemoteObjByIdx() in libs/pe/src/interface/object_queries.cpp
!================================================================================================
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

!================================================================================================
!                              Subroutine getRemoteParticle2
! C++ implementation: getRemPartStructByIdx() in libs/pe/src/interface/object_queries.cpp
!================================================================================================
interface
subroutine getRemoteParticle2(idx, particle) bind(C, name="getRemoteParticle2")
  use iso_c_binding, only: c_int, c_double
  import tParticleData
  integer(c_int), value :: idx
  type(tParticleData) :: particle
  end subroutine
end interface

!================================================================================================
!                              Subroutine setRemoteParticle2
! C++ implementation: setRemPartStruct() in libs/pe/src/interface/object_queries.cpp
!================================================================================================
interface
subroutine setRemoteParticle2(particle) bind(C, name="setRemoteParticle2")
  use iso_c_binding, only: c_int, c_double
  import tParticleData
  type(tParticleData) :: particle
  end subroutine
end interface

!================================================================================================
!                              Subroutine setParticle
! C++ implementation: setObjByIdx() in libs/pe/src/interface/object_queries.cpp
!================================================================================================
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

!================================================================================================
!                              Subroutine setForces
! C++ implementation: setForcesByIdx() in libs/pe/src/interface/object_queries.cpp
!================================================================================================
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

!================================================================================================
!                              Subroutine setRemoteForces
! C++ implementation: setRemoteForcesByIdx() in libs/pe/src/interface/object_queries.cpp
!================================================================================================
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

!================================================================================================
!                              Function isSphere
! C++ implementation: isSphereType() in libs/pe/src/interface/object_queries.cpp
!================================================================================================
interface
logical(c_bool) function isSphere(idx) bind(C, name="isTypeSphere")
  use iso_c_binding, only: c_int, c_bool
  integer(c_int) :: idx
  end function
end interface

!================================================================================================
!                              Subroutine getParticleRadius
! C++ implementation: getObjRadius() in libs/pe/src/interface/object_queries.cpp
!================================================================================================
interface
real(c_double) function getParticleRadius(idx) bind(C, name="getParticleRadius")
  use iso_c_binding, only: c_int, c_double
  integer(c_int) :: idx
  end function
end interface

!================================================================================================
!                              Function pointInsideObject
! C++ implementation: isInsideObject() in libs/pe/src/interface/object_queries.cpp
!================================================================================================
interface
logical(c_bool) function pointInsideObject(idx, pos) bind(C, name="pointInsideObject")
  use iso_c_binding, only: c_int, c_bool, c_double
  integer(c_int), value :: idx
  real(c_double) :: pos(*)
  end function
end interface

!================================================================================================
!                              Function pointInsideRemObject
! C++ implementation: isInsideRemObject() in libs/pe/src/interface/object_queries.cpp
!================================================================================================
interface
logical(c_bool) function pointInsideRemObject(idx, pos) bind(C, name="pointInsideRemObject")
  use iso_c_binding, only: c_int, c_bool, c_double
  integer(c_int), value :: idx
  real(c_double) :: pos(*)
  end function
end interface

!================================================================================================
!                              Function check_rem_id
! C++ implementation: checkRemoteFBM() in libs/pe/src/interface/object_queries.cpp
!================================================================================================
interface
logical(c_bool) function check_rem_id(fbmid, id) bind(C, name="check_rem_id")
  use iso_c_binding, only: c_int, c_bool, c_double
  integer(c_int), value :: fbmid
  integer(c_int), value :: id
  end function
end interface

! C++ implementation: getRemoteParticlesIndexMap() in libs/pe/src/interface/object_queries.cpp
interface
subroutine rem_particles_index_map(idxMap) bind(C, name="rem_particles_index_map")
  use iso_c_binding, only: c_int, c_double
  integer(c_int) :: idxMap(*)
  end subroutine
end interface

! C++ implementation: getParticlesIndexMap() in libs/pe/src/interface/object_queries.cpp
interface
subroutine particles_index_map(idxMap) bind(C, name="particles_index_map")
  use iso_c_binding, only: c_int, c_double
  integer(c_int) :: idxMap(*)
  end subroutine
end interface

! C++ implementation: uint64toByteArray() in libs/pe/src/interface/object_queries.h
interface
subroutine get_bytes(bytes) bind(C, name="get_bytes")
  use iso_c_binding, only: c_short
  integer(c_short), dimension(8) :: bytes 
end subroutine
end interface

! C++ implementation: not found
interface
logical(c_bool) function map_local_to_system(lidx, vidx) bind(C, name="map_local_to_system")
  use iso_c_binding, only: c_int, c_bool 
  integer(c_int), value :: lidx
  integer(c_int), value :: vidx
  end function
end interface

! C++ implementation: not found
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
!                              Function numTotalParticles
!================================================================================================
! Returns total number of particles
! - In PE_SERIAL_MODE: All domains have access to all particles
! - In parallel mode: Sum of local and remote particles for this domain
!================================================================================================

integer function numTotalParticles()
  implicit none

#ifdef PE_SERIAL_MODE
  ! Serial PE mode: all ranks have access to all particles
  numTotalParticles = getNumParticles()
#else
  ! Parallel PE mode: sum local and remote particles
  numTotalParticles = numRemParticles() + numLocalParticles()
#endif

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

  type(tParticleData) :: temp, zeroParticle
  integer, allocatable, dimension(:) :: indexMap
  integer :: a,i,idx,nLocal

  a = size(theParticles)

  ! Number of particles actually owned by this rank. In parallel PE a particle
  ! may live on a single rank only, so the caller-provided array (sized by the
  ! GLOBAL particle count, e.g. in ComputeCFL) can be larger than this rank's
  ! local index map. Reading the surplus slots from an uninitialised index map
  ! yields a garbage body index and aborts getParticle2 ("Body index ... out of
  ! range"). Read only the local slots and zero-fill the rest. This is a no-op
  ! when local == global (the usual multi-particle, every-rank-populated case).
  nLocal = getNumParticles()
  if (nLocal .gt. a) nLocal = a

  allocate(indexMap(a))

  if (nLocal .gt. 0) call getParticlesIndexMap(indexMap)

  zeroParticle%position = 0d0
  zeroParticle%velocity = 0d0
  zeroParticle%angvel   = 0d0
  zeroParticle%force    = 0d0
  zeroParticle%torque   = 0d0
  zeroParticle%time     = 0d0
  zeroParticle%density  = 0d0
  zeroParticle%aabb     = 0d0
  zeroParticle%radius   = 0d0
  zeroParticle%typeId   = 0
  zeroParticle%localIdx = 0
  zeroParticle%uniqueIdx= 0
  zeroParticle%systemIdx= 0
  zeroParticle%bytes    = 0

  do i=1,a
    if (i .le. nLocal) then
      idx = indexMap(i)
      call getParticle2(idx,&
                        temp)
      theParticles(i) = temp
    else
      theParticles(i) = zeroParticle
    end if
  end do

  deallocate(indexMap)

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
  
  do i=1,a
    idx = indexMap(i)
    call getRemoteParticle2(idx, temp)  
    theParticles(i) = temp
  end do

end subroutine getAllRemoteParticles
!================================================================================================
!                              Function setForcesMapped
!================================================================================================

subroutine setForcesMapped(particle)
  implicit none

  type(tParticleData), intent(inout) :: particle

  ! Canonical CFD->PE force path: write the full particle struct so
  ! force, torque, and identity stay coupled in one interface.
  call setLocalParticle2(particle)

end subroutine setForcesMapped
!================================================================================================
!                              Function setRemoteForcesMapped
!================================================================================================

subroutine setRemoteForcesMapped(particle)
  implicit none

  type(tParticleData), intent(inout) :: particle

  ! Remote/shadow-particle equivalent of the canonical struct-based
  ! CFD->PE force path.
  call setRemoteParticle2(particle)

end subroutine setRemoteForcesMapped
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
!                              Subroutine setLocalParticle2
!================================================================================================

subroutine setLocalParticle2(particle)
  implicit none

  type(tParticleData), intent(inout) :: particle
  
  call setParticle2(particle)

end subroutine setLocalParticle2
  

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

#endif 
end module dem_query
