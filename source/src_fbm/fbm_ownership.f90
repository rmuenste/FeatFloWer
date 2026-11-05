!##############################################################################
!# fbm_ownership - multiple ownership of FBM DOFs, as a CSR exception table
!#
!# The canonical owner of a solid DOF stays where it always was: FictKNPR (the
!# flag) and FictKNPR_uint64 (the pe system id), dense, one per DOF, and since
!# pe 8ea90af the MINIMUM id among the bodies containing the point. Overlapping
!# bodies are real (contact), and a DOF inside two of them has two owners that
!# the velocity BC and the force attribution cannot both honour; today they
!# silently take one. This module makes the ambiguity VISIBLE and COUNTABLE
!# without touching the hot path: only DOFs inside two or more bodies get a
!# row, so a run with no overlap allocates no payload at all.
!#
!#   FictOwn_row(i)      0 = fluid or single owner, else the row r of DOF i
!#   FictOwn_dof(r)      the DOF of row r (ascending, the sweep order)
!#   FictOwn_ptr(r)      first index of row r in FictOwn_id; ptr(nexc+1) = nnz+1
!#   FictOwn_id(k)       owner system ids, ASCENDING within a row
!#
!# Invariant the gate asserts: FictOwn_id(FictOwn_ptr(r)) equals the canonical
!# owner FictKNPR_uint64(FictOwn_dof(r)) - the first id of a row IS the minimum.
!# A consumer that ignores the table therefore sees exactly today's behaviour.
!#
!# Filled only when SimPar@FbmOwnerTable = Yes (default No), by the
!# classification handler, from pe's multi-owner query. The first use is
!# diagnostic: how many DOFs have several owners in a real run is the number
!# that decides whether a consumer-side rule is worth writing. Nothing here is
!# collective except fbm_own_report, which must run on a rank-symmetric path.
!##############################################################################
module fbm_ownership

  use iso_c_binding, only: c_int, c_double, c_long_long
  use var_QuadScalar, only: FictKNPR_uint64, myErrorCode

  implicit none
  private
  public :: fbm_own_reset, fbm_own_record, fbm_own_finalize, fbm_own_report
  public :: pointOwners, pointOwnersAccelerated, FBM_OWN_MAX
  public :: FictOwn_row, FictOwn_dof, FictOwn_ptr, FictOwn_id
  public :: fbm_own_nexc, fbm_own_nnz

  ! Owners per DOF the query buffer holds. More than this is counted as overflow and
  ! reported; it is not an error, it is a number worth knowing.
  integer, parameter :: FBM_OWN_MAX = 16

  integer :: own_ndof = 0
  integer :: own_nexc = 0        ! rows = DOFs with two or more owners
  integer :: own_nnz  = 0        ! stored ids
  integer :: own_maxper = 0      ! most owners seen on one DOF (before truncation)
  integer :: own_overflow = 0    ! DOFs with more than FBM_OWN_MAX owners
  integer :: own_inv_viol = 0    ! rows whose first id is not the canonical owner

  integer,         allocatable :: FictOwn_row(:)
  integer,         allocatable :: FictOwn_dof(:)
  integer,         allocatable :: FictOwn_ptr(:)
  integer(kind=8), allocatable :: FictOwn_id(:)

  ! pe: every body containing the point, system ids ascending; nOwners is the TOTAL count,
  ! ids(1:min(nOwners,maxOwners)) are filled. libs/pe/src/interface/object_queries.cpp.
  interface
    subroutine pointOwners(pos, maxOwners, ids, nOwners) bind(C, name="pointOwners")
      import :: c_double, c_int, c_long_long
      real(c_double),       intent(in)  :: pos(3)
      integer(c_int),       intent(in)  :: maxOwners
      integer(c_long_long), intent(out) :: ids(*)
      integer(c_int),       intent(out) :: nOwners
    end subroutine
    subroutine pointOwnersAccelerated(pos, maxOwners, ids, nOwners) bind(C, name="pointOwnersAccelerated")
      import :: c_double, c_int, c_long_long
      real(c_double),       intent(in)  :: pos(3)
      integer(c_int),       intent(in)  :: maxOwners
      integer(c_long_long), intent(out) :: ids(*)
      integer(c_int),       intent(out) :: nOwners
    end subroutine
  end interface

contains

  integer function fbm_own_nexc()
    fbm_own_nexc = own_nexc
  end function
  integer function fbm_own_nnz()
    fbm_own_nnz = own_nnz
  end function

  ! Called before a classification sweep. The dense row index is the one array sized by
  ! the mesh (one default integer per DOF, the same order as FictKNPR itself); the
  ! payload starts small and grows only if overlaps exist.
  subroutine fbm_own_reset(ndof)
    integer, intent(in) :: ndof
    own_ndof = ndof
    if (allocated(FictOwn_row)) then
      if (size(FictOwn_row) /= ndof) deallocate(FictOwn_row)
    end if
    if (.not. allocated(FictOwn_row)) allocate(FictOwn_row(ndof))
    FictOwn_row = 0
    if (.not. allocated(FictOwn_dof)) allocate(FictOwn_dof(64))
    if (.not. allocated(FictOwn_ptr)) allocate(FictOwn_ptr(65))
    if (.not. allocated(FictOwn_id))  allocate(FictOwn_id(256))
    own_nexc = 0; own_nnz = 0; own_maxper = 0; own_overflow = 0; own_inv_viol = 0
  end subroutine fbm_own_reset

  ! One DOF's owner list from pe. Nothing is stored for zero or one owner.
  subroutine fbm_own_record(vidx, nOwners, ids)
    integer,              intent(in) :: vidx
    integer(c_int),       intent(in) :: nOwners
    integer(c_long_long), intent(in) :: ids(*)
    integer :: m
    if (nOwners < 2) return
    if (vidx < 1 .or. vidx > own_ndof) return
    m = min(int(nOwners), FBM_OWN_MAX)
    if (nOwners > FBM_OWN_MAX) own_overflow = own_overflow + 1
    own_maxper = max(own_maxper, int(nOwners))
    if (own_nexc + 1 > size(FictOwn_dof)) call grow_i4(FictOwn_dof, 2*size(FictOwn_dof))
    if (own_nexc + 2 > size(FictOwn_ptr)) call grow_i4(FictOwn_ptr, 2*size(FictOwn_ptr))
    if (own_nnz + m > size(FictOwn_id))   call grow_i8(FictOwn_id, max(2*size(FictOwn_id), own_nnz + m))
    own_nexc = own_nexc + 1
    FictOwn_dof(own_nexc) = vidx
    FictOwn_ptr(own_nexc) = own_nnz + 1
    FictOwn_id(own_nnz+1:own_nnz+m) = ids(1:m)
    own_nnz = own_nnz + m
    FictOwn_row(vidx) = own_nexc
  end subroutine fbm_own_record

  ! Called after the sweep: close the CSR and check the invariant against the canonical
  ! owner the sweep just stored. Rank-local, no collectives.
  subroutine fbm_own_finalize()
    integer :: r, dof
    integer(kind=8) :: canonical
    if (.not. allocated(FictOwn_ptr)) return
    FictOwn_ptr(own_nexc+1) = own_nnz + 1
    own_inv_viol = 0
    if (.not. allocated(FictKNPR_uint64)) return
    do r = 1, own_nexc
      dof = FictOwn_dof(r)
      if (dof > size(FictKNPR_uint64)) cycle
      canonical = bytes_to_i8(FictKNPR_uint64(dof)%bytes)
      if (FictOwn_id(FictOwn_ptr(r)) /= canonical) own_inv_viol = own_inv_viol + 1
    end do
  end subroutine fbm_own_finalize

  ! The collective report and verdict for one stage. RANK-SYMMETRIC path only.
  ! expectMulti: the stage was built to produce overlapping bodies on real DOFs (stage E),
  ! so a count of zero means its construction reached no DOF - a failure of the harness,
  ! reported as such.
  subroutine fbm_own_report(simTime, stageName, expectMulti)
    use PP3D_MPI, only: myid, COMM_SUMMN, MPI_COMM_WORLD
    real*8,           intent(in) :: simTime
    character(len=*), intent(in) :: stageName
    logical,          intent(in) :: expectMulti
    real*8  :: buf(5)
    integer :: ierr
    logical :: failed
    character(len=4) :: verdictWord

    buf = 0d0
    if (myid /= 0) then
      buf(1) = dble(own_nexc)
      buf(2) = dble(own_nnz)
      buf(3) = dble(own_overflow)
      buf(4) = dble(own_inv_viol)
      buf(5) = dble(own_maxper)   ! summed over ranks: serial PE sees the same world on each
    end if
    call COMM_SUMMN(buf, 5)

    failed = (nint(buf(4)) > 0) .or. (expectMulti .and. nint(buf(1)) == 0)
    verdictWord = merge('FAIL', 'pass', failed)

    if (myid == 1) then
      write(*,'(A,ES16.8,A,A,A,I0,A,I0,A,I0,A,I0,A,I0,A,A,A,A)') &
        'DNS_QUERY_OWNERS time= ', simTime, ' stage= ', trim(stageName), &
        ' multi_owned_dofs= ', nint(buf(1)), ' nnz= ', nint(buf(2)), &
        ' max_owners_this_rank= ', own_maxper, ' overflow= ', nint(buf(3)), &
        ' canonical_violations= ', nint(buf(4)), &
        ' expected_multi= ', merge('>0 ', 'any', expectMulti), ' verdict= ', verdictWord
      if (failed) then
        if (nint(buf(4)) > 0) then
          write(*,'(A)') 'DNS_QUERY_OWNERS FAILED: the first id of an owner row is not the'
          write(*,'(A)') '  canonical owner in FictKNPR_uint64 - the minimum-id rule and the'
          write(*,'(A)') '  sorted owner list disagree, which should be impossible.'
        else
          write(*,'(A)') 'DNS_QUERY_OWNERS FAILED (control): this stage was built to overlap'
          write(*,'(A)') '  bodies on real DOFs and no DOF has two owners. The overlap'
          write(*,'(A)') '  construction reached no DOF, so it tests nothing.'
        end if
      end if
      flush(6)
    end if
    if (failed) call MPI_Abort(MPI_COMM_WORLD, myErrorCode%FBM_QUERY_GATE_FAILED, ierr)
  end subroutine fbm_own_report

  ! pe writes byte_array(i) = id >> (7-i)*8: bytes(1) is the MOST significant byte.
  pure function bytes_to_i8(bytes) result(id)
    integer(kind=2), intent(in) :: bytes(8)
    integer(kind=8) :: id
    integer :: i
    id = 0
    do i = 1, 8
      id = ior(ishft(id, 8), int(iand(int(bytes(i)), 255), kind=8))
    end do
  end function bytes_to_i8

  subroutine grow_i4(a, n)
    integer, allocatable, intent(inout) :: a(:)
    integer, intent(in) :: n
    integer, allocatable :: t(:)
    allocate(t(n)); t = 0
    t(1:size(a)) = a
    call move_alloc(t, a)
  end subroutine grow_i4

  subroutine grow_i8(a, n)
    integer(kind=8), allocatable, intent(inout) :: a(:)
    integer, intent(in) :: n
    integer(kind=8), allocatable :: t(:)
    allocate(t(n)); t = 0
    t(1:size(a)) = a
    call move_alloc(t, a)
  end subroutine grow_i8

end module fbm_ownership
