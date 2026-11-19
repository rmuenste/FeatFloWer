!=========================================================================
! QuadSc_pco_gather.f90 - collection of the coarse pressure operator of the workers on the master
!
! Part of the periodic coarse pressure operator (overview: QuadSc_pco_operator.f90).
!
! What the workers apply on the coarse level lev (E012_DAX) is, row by row for the rows of
! their own elements,  (C x)_r + (CP xp)_r : C couples own elements, CP couples ghost
! elements of other ranks (periodic neighbours included), xp being the ghost values of the
! pressure exchange GetParPressure. Both matrices are already assembled on the workers with
! the actual coefficients (weighted EL mass, slip projection, constraints). Collecting them
! into the master numbering therefore gives the master exactly the operator the workers
! apply - the reason this mode is the default and the oracle of the wrapped one.
!
! Global numbering: master element of a worker element on level lev from the coarse element
! link (coarse%myELEMLINK) and the refinement numbering (child 1 keeps the parent number,
! children 2..8 of parent p of a level with N parents follow at N + 7(p-1) + 1..7, on the
! worker and on the master alike): the mapping of E012DISTR_L1/L2/L3. Global P1 unknown
! 4(E-1)+m. The ghost columns get their global numbers by exchanging the global numbers of
! the own unknowns with GetParPressure, as the MUMPS exporter (MUMPS_SetUpP1_SLAVE) does.
!
! Plan (built once per pattern, PCO_GATHER_BUILD, collective): every worker sends the global
! (row, column) of each entry of C and CP in their storage order; the master finds the slot of
! each in its layout (pco_coarse_solver) and keeps the map. Validation on the master:
!   missing    an entry whose (row, column) is not in the master layout (a worker coupling
!              the periodic classes of the master do not know: missing partner, wrong period)
!   ownership  a master row filled by two workers, or by none (unowned)
!   uncovered  layout slots no worker fills (structural zeros; reported, not refused)
! Values (PCO_GATHER_VALUES, every renewal, collective): the workers send their C and CP value
! arrays as they are; the master sums them into its slots (duplicates add, which is what the
! worker product does).
! Immutability (PCO_GATHER_UNCHANGED, every renewal, collective): the plan is valid only for the
! pattern it was built for. Every worker recomputes the global (row, column) of its C and CP
! entries (one pressure exchange for the ghost numbers) and compares them with the ones of the
! plan; a changed pattern, ghost numbering or ownership is reported on every rank (the
! coordinator stops: a changed coarse pattern needs a new layout, classes and gauge).
!=========================================================================
MODULE pco_gather

  USE PP3D_MPI, ONLY: myid, master, numnodes, subnodes, coarse, MPI_COMM_WORLD, MPI_INTEGER, &
                      MPI_DOUBLE_PRECISION, MPI_LAND, MPI_LOGICAL
  USE var_QuadScalar, ONLY: mg_mesh, mg_lMat, mg_lPMat, mg_CMat, mg_CPMat
  USE def_FEAT, ONLY: ILEV
  USE pco_periodic_topology, ONLY: PCO_GHOST_BUFFER_SIZE

  IMPLICIT NONE
  PRIVATE

  TYPE, PUBLIC :: tPcoGatherPlan
    LOGICAL :: built = .FALSE.
    INTEGER :: lev = 0
    INTEGER :: nC = 0, nCP = 0                 ! value counts of this worker at plan time
    INTEGER, ALLOCATABLE :: nelm(:)            ! master elements per level 1..lev
    INTEGER, ALLOCATABLE :: irow(:), icol(:)   ! worker: global row/column of every entry
    INTEGER, ALLOCATABLE :: counts(:), displs(:)
    INTEGER, ALLOCATABLE :: slot(:)            ! master: received value -> raw slot
    REAL*8,  ALLOCATABLE :: buf(:)
    INTEGER :: nmissing = 0, nownership = 0, nunowned = 0, nuncovered = 0
  END TYPE tPcoGatherPlan

  PUBLIC :: PCO_GATHER_BUILD, PCO_GATHER_VALUES, PCO_GATHER_UNCHANGED, PCO_GATHER_FREE
  PUBLIC :: PCO_MASTER_ELEMENT_IDS

CONTAINS

  !-----------------------------------------------------------------------
  ! Master element numbers (level lev) of the elements of this worker of level lev.
  ! nelm(l): number of master elements on level l = 1..lev.
  !-----------------------------------------------------------------------
  SUBROUTINE PCO_MASTER_ELEMENT_IDS(lev, nelm, gid)
    INTEGER, INTENT(IN) :: lev, nelm(:)
    INTEGER, INTENT(OUT) :: gid(:)
    INTEGER :: iel
    DO iel = 1, mg_mesh%level(lev)%nel
      gid(iel) = master_elem(iel, lev)
    END DO
  CONTAINS
    RECURSIVE INTEGER FUNCTION master_elem(e, l) RESULT(g)
      INTEGER, INTENT(IN) :: e, l
      INTEGER :: nwp, k
      IF (l.EQ.1) THEN
        g = coarse%myELEMLINK(e)
        RETURN
      END IF
      nwp = mg_mesh%level(l-1)%nel
      IF (e.LE.nwp) THEN
        g = master_elem(e, l-1)
      ELSE
        k = e - nwp
        g = nelm(l-1) + (master_elem((k-1)/7+1, l-1) - 1)*7 + MOD(k-1,7) + 1
      END IF
    END FUNCTION master_elem
  END SUBROUTINE PCO_MASTER_ELEMENT_IDS

  !-----------------------------------------------------------------------
  ! Build the plan for level lev into the master layout ldA/colA (master arguments; dummies
  ! elsewhere). Collective. The counters of the validation are broadcast to every rank.
  !-----------------------------------------------------------------------
  SUBROUTINE PCO_GATHER_BUILD(plan, lev, ldA, colA)
    TYPE(tPcoGatherPlan), INTENT(INOUT) :: plan
    INTEGER, INTENT(IN) :: lev
    INTEGER, INTENT(IN) :: ldA(:), colA(:)
    INTEGER, ALLOCATABLE :: irow(:), icol(:), grow(:), gcol(:), owner(:)
    LOGICAL, ALLOCATABLE :: covered(:)
    INTEGER :: l, k, n, ierr, pid, ntot, nrow, slot, cnt(4)

    CALL PCO_GATHER_FREE(plan)
    plan%lev = lev
    ALLOCATE(plan%nelm(lev))
    IF (myid.EQ.master) THEN
      DO l = 1, lev
        plan%nelm(l) = mg_mesh%level(l)%nel
      END DO
    END IF
    CALL MPI_Bcast(plan%nelm, lev, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL worker_entry_ids(lev, plan%nelm, irow, icol, plan%nC, plan%nCP)
    n = plan%nC + plan%nCP
    IF (myid.NE.master) THEN
      ALLOCATE(plan%irow(MAX(n,1)), plan%icol(MAX(n,1)))
      plan%irow = irow
      plan%icol = icol
    END IF

    ALLOCATE(plan%counts(0:numnodes-1), plan%displs(0:numnodes-1))
    CALL MPI_Gather(n, 1, MPI_INTEGER, plan%counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    IF (myid.EQ.master) THEN
      plan%displs(0) = 0
      DO pid = 1, numnodes-1
        plan%displs(pid) = plan%displs(pid-1) + plan%counts(pid-1)
      END DO
      ntot = plan%displs(numnodes-1) + plan%counts(numnodes-1)
      ALLOCATE(grow(MAX(ntot,1)), gcol(MAX(ntot,1)))
    ELSE
      ALLOCATE(grow(1), gcol(1))
      ntot = 0
    END IF
    CALL MPI_Gatherv(irow, n, MPI_INTEGER, grow, plan%counts, plan%displs, MPI_INTEGER, 0, &
                     MPI_COMM_WORLD, ierr)
    CALL MPI_Gatherv(icol, n, MPI_INTEGER, gcol, plan%counts, plan%displs, MPI_INTEGER, 0, &
                     MPI_COMM_WORLD, ierr)

    cnt = 0
    IF (myid.EQ.master) THEN
      nrow = SIZE(ldA) - 1
      ALLOCATE(plan%slot(MAX(ntot,1)), owner(nrow), covered(MAX(SIZE(colA),1)))
      plan%slot = 0
      owner = 0
      covered = .FALSE.
      DO pid = 1, numnodes-1
        DO k = plan%displs(pid)+1, plan%displs(pid)+plan%counts(pid)
          slot = 0
          IF (grow(k).GE.1 .AND. grow(k).LE.nrow .AND. gcol(k).GE.1 .AND. gcol(k).LE.nrow) THEN
            slot = find_slot(ldA, colA, grow(k), gcol(k))
            IF (owner(grow(k)).EQ.0) THEN
              owner(grow(k)) = pid
            ELSE IF (owner(grow(k)).NE.pid) THEN
              cnt(2) = cnt(2) + 1
            END IF
          END IF
          IF (slot.EQ.0) THEN
            cnt(1) = cnt(1) + 1
          ELSE
            covered(slot) = .TRUE.
          END IF
          plan%slot(k) = slot
        END DO
      END DO
      cnt(3) = COUNT(owner.EQ.0)
      cnt(4) = COUNT(.NOT.covered(1:SIZE(colA)))
      ALLOCATE(plan%buf(MAX(ntot,1)))
      DEALLOCATE(owner, covered)
    ELSE
      ALLOCATE(plan%buf(MAX(n,1)))
    END IF
    CALL MPI_Bcast(cnt, 4, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    plan%nmissing = cnt(1)
    plan%nownership = cnt(2)
    plan%nunowned = cnt(3)
    plan%nuncovered = cnt(4)
    plan%built = .TRUE.
    DEALLOCATE(irow, icol, grow, gcol)

  END SUBROUTINE PCO_GATHER_BUILD

  INTEGER FUNCTION find_slot(ldA, colA, row, col) RESULT(slot)
    INTEGER, INTENT(IN) :: ldA(:), colA(:), row, col
    INTEGER :: k
    slot = 0
    DO k = ldA(row), ldA(row+1)-1
      IF (colA(k).EQ.col) THEN
        slot = k
        RETURN
      END IF
    END DO
  END FUNCTION find_slot

  !-----------------------------------------------------------------------
  ! Global (row, column) of every C and CP entry of this worker on level lev, in storage
  ! order (C entries first); collective (the pressure exchange of the ghost numbers). On the
  ! master: empty, nC = nCP = 0.
  !-----------------------------------------------------------------------
  SUBROUTINE worker_entry_ids(lev, nelm, irow, icol, nC, nCP)
    INTEGER, INTENT(IN) :: lev, nelm(:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: irow(:), icol(:)
    INTEGER, INTENT(OUT) :: nC, nCP
    INTEGER, ALLOCATABLE :: gid(:)
    REAL*8, ALLOCATABLE :: x(:), xp(:)
    INTEGER :: nel, r, k, ilev_save

    ilev_save = ILEV
    IF (myid.NE.master) THEN
      nel = mg_mesh%level(lev)%nel
      ALLOCATE(gid(nel), x(4*nel), xp(PCO_GHOST_BUFFER_SIZE(lev)))
      CALL PCO_MASTER_ELEMENT_IDS(lev, nelm, gid)
      DO r = 1, 4*nel
        x(r) = DBLE(4*(gid((r-1)/4+1)-1) + MOD(r-1,4) + 1)
      END DO
      xp = 0d0
      nC = mg_lMat(lev)%LdA(4*nel+1) - 1
      nCP = mg_lPMat(lev)%LdA(4*nel+1) - 1
      ALLOCATE(irow(MAX(nC+nCP,1)), icol(MAX(nC+nCP,1)))
    ELSE
      ALLOCATE(x(4), xp(4), irow(1), icol(1))
      nC = 0
      nCP = 0
    END IF
    ILEV = lev
    CALL GetParPressure(x, xp)
    ILEV = ilev_save
    IF (myid.NE.master) THEN
      DO r = 1, 4*nel
        DO k = mg_lMat(lev)%LdA(r), mg_lMat(lev)%LdA(r+1)-1
          irow(k) = NINT(x(r))
          icol(k) = NINT(x(mg_lMat(lev)%ColA(k)))
        END DO
        DO k = mg_lPMat(lev)%LdA(r), mg_lPMat(lev)%LdA(r+1)-1
          irow(nC+k) = NINT(x(r))
          icol(nC+k) = NINT(xp(mg_lPMat(lev)%ColA(k)))
        END DO
      END DO
      DEALLOCATE(gid)
    END IF
    DEALLOCATE(x, xp)
  END SUBROUTINE worker_entry_ids

  !-----------------------------------------------------------------------
  ! .TRUE. on every rank when every worker still has the pattern, the ghost numbering and the
  ! ownership the plan was built for (see the module header). Collective.
  !-----------------------------------------------------------------------
  LOGICAL FUNCTION PCO_GATHER_UNCHANGED(plan) RESULT(same)
    TYPE(tPcoGatherPlan), INTENT(IN) :: plan
    INTEGER, ALLOCATABLE :: irow(:), icol(:)
    INTEGER :: nC, nCP, ierr
    LOGICAL :: mine
    CALL worker_entry_ids(plan%lev, plan%nelm, irow, icol, nC, nCP)
    mine = plan%built
    IF (myid.NE.master .AND. mine) THEN
      mine = (nC.EQ.plan%nC .AND. nCP.EQ.plan%nCP)
      IF (mine .AND. nC+nCP.GT.0) mine = ALL(irow(1:nC+nCP).EQ.plan%irow(1:nC+nCP)) .AND. &
                                         ALL(icol(1:nC+nCP).EQ.plan%icol(1:nC+nCP))
    END IF
    DEALLOCATE(irow, icol)
    CALL MPI_Allreduce(mine, same, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
  END FUNCTION PCO_GATHER_UNCHANGED

  !-----------------------------------------------------------------------
  ! Collect the current values of C and CP (level plan%lev) into raw (master: zeroed and
  ! summed; elsewhere unused). Collective.
  !-----------------------------------------------------------------------
  SUBROUTINE PCO_GATHER_VALUES(plan, raw)
    TYPE(tPcoGatherPlan), INTENT(INOUT) :: plan
    REAL*8, INTENT(INOUT) :: raw(:)
    INTEGER :: n, k, ierr, ntot
    REAL*8 :: none(1)
    IF (myid.NE.master) THEN
      n = plan%nC + plan%nCP
      IF (plan%nC.GT.0) plan%buf(1:plan%nC) = mg_CMat(plan%lev)%a(1:plan%nC)
      IF (plan%nCP.GT.0) plan%buf(plan%nC+1:n) = mg_CPMat(plan%lev)%a(1:plan%nCP)
      CALL MPI_Gatherv(plan%buf, n, MPI_DOUBLE_PRECISION, none, plan%counts, plan%displs, &
                       MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    ELSE
      none = 0d0
      CALL MPI_Gatherv(none, 0, MPI_DOUBLE_PRECISION, plan%buf, plan%counts, plan%displs, &
                       MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      ntot = plan%displs(numnodes-1) + plan%counts(numnodes-1)
      raw = 0d0
      DO k = 1, ntot
        IF (plan%slot(k).GT.0) raw(plan%slot(k)) = raw(plan%slot(k)) + plan%buf(k)
      END DO
    END IF
  END SUBROUTINE PCO_GATHER_VALUES

  SUBROUTINE PCO_GATHER_FREE(plan)
    TYPE(tPcoGatherPlan), INTENT(INOUT) :: plan
    IF (ALLOCATED(plan%counts)) DEALLOCATE(plan%counts)
    IF (ALLOCATED(plan%displs)) DEALLOCATE(plan%displs)
    IF (ALLOCATED(plan%slot)) DEALLOCATE(plan%slot)
    IF (ALLOCATED(plan%buf)) DEALLOCATE(plan%buf)
    IF (ALLOCATED(plan%nelm)) DEALLOCATE(plan%nelm)
    IF (ALLOCATED(plan%irow)) DEALLOCATE(plan%irow)
    IF (ALLOCATED(plan%icol)) DEALLOCATE(plan%icol)
    plan%built = .FALSE.
  END SUBROUTINE PCO_GATHER_FREE

END MODULE pco_gather
