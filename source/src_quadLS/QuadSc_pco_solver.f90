!=========================================================================
! QuadSc_pco_solver.f90 - storage, gauge and solves of the periodic coarse pressure operator
! on the master
!
! Part of the periodic coarse pressure operator (overview: QuadSc_pco_operator.f90). Owns
! the master-side system of a non-legacy mode; the coordinator writes the raw values, the
! multigrid (QuadSc_mg.f90) solves with it. Never entered in legacy mode (PCO_ACTIVE false).
!
! Numbering and layout (the conventions of the legacy master matrix, Get_CMatStruct iP=1):
!   unknown 4(e-1)+m is mode m of the P1 pressure of master element e of the coarse level
!   (m = 1 the element mean, 2..4 the coefficients of the physical linear basis x-x0, y-y0,
!   z-z0); the row block of e lists the own block first, then the blocks of the neighbour
!   elements in ascending order (element graph through vertex classes,
!   pco_periodic_topology); row 4(e-1)+m starts with its diagonal, then the other modes of e
!   ascending, then 4 entries per neighbour. Hence the mean-mode system is read with stride 4,
!   as the legacy type-3/4 extraction does.
!
! raw  : the operator as the workers apply it (gathered) or its master-side equivalent
!        (wrapped); ungauged and never modified by a solve. It is the operator of the type-3/4
!        smoothing (PCO_RAW_SOR: the SOR sweep of crs_E012_SOR_Smoother; a row with a zero
!        diagonal is left unchanged instead of divided by zero) and of the checks.
! gauge: decided ONCE per layout, at the first refresh, from the structure of the operator,
!        with the matrix as a consistency check:
!        - components: elements e, f are coupled when any entry of their 4x4 block is nonzero
!          (a weak but real coupling joins them; two regions that touch only through
!          constrained dofs have an all-zero block and are two components);
!        - structure (authoritative): s_c = the largest structural flux of a constant pressure
!          over the free velocity dofs of the elements of c (PCO_FREE_FLUX on the workers:
!          |sum| / sum|.| of the mean-mode entries of each B row with a free velocity
!          component; rounding level where the dof support closes, order one at a free
!          boundary, whatever the cell aspect ratio), handed over by PCO_SET_STRUCTURAL_FLUX.
!          s_c <= PCO_FLUX_CLOSED (256 eps = 2^-44 = 5.7e-14): closed; s_c >= PCO_FLUX_OPEN
!          (1e-6): open; in between the decision is undecidable. The closed bound sits just
!          above the measured rounding population of closed components (1.1e-15 to 1.4e-14,
!          i.e. 5 to 63 eps over every gate; factor 4 to the largest): a physical flux below
!          it cannot be told from the rounding of the assembled entries;
!        - matrix (consistency): right(c) = max over EVERY row i of |(C e_c)_i| / sum_j |C_ij|
!          (all columns: the normwise scale). A closed component must have right(c) <=
!          PCO_CLOSURE_TOL (512 eps = 2^-43 = 1.1e-13; measured closed 1.3e-15 to 2.0e-14,
!          factor 5.7); an open one right(c) > PCO_CLOSURE_TOL, otherwise the constant is
!          numerically in the kernel of an operator that is structurally regular.
!          A weak outlet (cells stretched across it) gives an open component with a small
!          right(c) (2.5e-10 in the unit test): it stays open, no reference;
!        - undecidable or inconsistent: status PCO_GAUGE_UNDECIDED, nothing is pinned (the
!          coordinator refuses an explicit key or falls back to legacy, logged);
!        - one reference per closed component, the mean mode of its lowest element, set to
!          zero in the SOLVER COPY only: row and column cleared except the diagonal, the
!          right-hand-side entry zeroed (symmetric elimination of a zero value); a zero
!          diagonal of a reference row (a closed component whose mean row vanishes) is
!          replaced on the copy by the largest diagonal modulus of the copy, so the reference
!          equation reads d x = 0 and not 0 x = 0.
!        Numerical limit: a component whose physical boundary flux is below the closed bound
!        (relative to the B entries around the element) is treated as closed and pinned. At
!        that level the ungauged system is singular to working precision (unit test 20: with
!        the right-hand side perturbed at the backward-error level, the direct solve of the
!        ungauged ring at outlet strength 2^-50 has a mean error of 0.74; 5e-15 with the exact one);
!        between the bound and 1e-6 nothing is pinned (fallback / refusal).
!        Later refreshes keep the decision (no flipping between pinned and unpinned), check the
!        matrix measure (a closed component above PCO_CLOSURE_TOL or an open one at or below it
!        is reported, pco%drift) and keep the CURRENT smallest measure of the open components
!        (pco%closure_open) for the admission of the P0 types. The left measure e_c^T C (column sums of
!        moduli as scale) is computed at the decision and reported; the right-hand-side
!        compatibility measure below assumes the left null vector e (symmetric operator).
!        A right-hand side whose sum over a closed component does not vanish is reported
!        (component, absolute net value, relative imbalance), not changed.
! solver copy: full (types 1, 2) or mean modes only (types 3, 4), pattern fixed (explicit zeros
!        kept), values refreshed from raw. Types 2/4: UMFPACK in explicit column storage of this
!        matrix (tUmfCSC: it solves A x = b, not A^T x = b); symbolic analysis once per pattern,
!        numeric factorisation at every refresh. Types 1/3: BiCGStab of this module with an
!        SSOR preconditioner (every diagonal of the copy must be nonzero, checked at the
!        refresh) that scans every off-diagonal entry with a lower/upper predicate
!        (the legacy PARID117_MASTER stops at the first upper column of a row, which in this
!        diagonal-first layout skips lower neighbours); a zero right-hand side returns the zero
!        correction, convergence is tested before every division, a breakdown is reported.
! Status: factorisation failures are returned to the coordinator (collective stop there); a
! failing solve inside the multigrid runs on the master only and aborts with a message.
! Diagnostics (PCO_SOLVE_CHECK, PCO_PROBE_IMBALANCE) never touch the production statistics.
!=========================================================================
MODULE pco_coarse_solver

  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  USE PP3D_MPI, ONLY: myid, master, MPI_COMM_WORLD
  USE def_FEAT, ONLY: MTERM
  USE UMFPackSolver, ONLY: tUmfCSC, UmfCSC_SetStructure, UmfCSC_SetValues, UmfCSC_Analyse, &
                           UmfCSC_Factorize, UmfCSC_Solve, UmfCSC_Release, UMF_STATUS_OK

  IMPLICIT NONE
  PRIVATE

  TYPE, PUBLIC :: tPcoSystem
    LOGICAL :: active = .FALSE.          ! a non-legacy operator serves the master coarse solve
    INTEGER :: mode = 0                  ! pressure_coarse_key value
    INTEGER :: lev = 0, nel = 0, n = 0
    INTEGER :: solver_type = 0
    ! layout and raw operator
    INTEGER, ALLOCATABLE :: eptr(:), enbr(:)
    INTEGER, ALLOCATABLE :: ldA(:), colA(:)
    REAL*8,  ALLOCATABLE :: a(:)
    ! gauge (decided once per layout)
    LOGICAL :: gauge_decided = .FALSE.
    INTEGER :: ncomp = 0, nclosed = 0
    INTEGER, ALLOCATABLE :: comp_of(:), comp_ref(:)
    LOGICAL, ALLOCATABLE :: comp_closed(:)
    REAL*8  :: closure_open = 0d0        ! smallest closure measure of the open components
    REAL*8  :: closure_closed = 0d0      ! largest closure measure of the closed components
    REAL*8  :: left_closed = 0d0         ! largest |e^T A| / e^T |A| of the closed components
    LOGICAL :: drift = .FALSE.           ! a later refresh contradicts the gauge decision
    REAL*8, ALLOCATABLE :: elem_flux(:)  ! structural flux per element (PCO_SET_STRUCTURAL_FLUX)
    REAL*8  :: flux_closed = 0d0         ! largest structural flux of the closed components
    REAL*8  :: flux_open = 0d0           ! smallest structural flux of the open components
    ! solver copy
    LOGICAL :: reduced = .FALSE.
    INTEGER :: ns = 0
    INTEGER, ALLOCATABLE :: sldA(:), scolA(:), sfrom(:)
    REAL*8,  ALLOCATABLE :: sa(:)
    INTEGER, ALLOCATABLE :: scomp(:)     ! component of a solver unknown
    LOGICAL, ALLOCATABLE :: sref(:)      ! reference unknown (x = 0, rhs zeroed)
    TYPE(tUmfCSC) :: lu
    LOGICAL :: analysed = .FALSE.
    REAL*8,  ALLOCATABLE :: rhs(:), sol(:), w(:,:)
    ! production statistics (diagnostic solves do not count)
    INTEGER :: nsolve = 0, nincompat = 0, nbreakdown = 0
    REAL*8  :: max_incompat = 0d0        ! largest relative imbalance |sum b| / sum |b|
    REAL*8  :: max_incompat_net = 0d0    ! the net sum at that solve
    INTEGER :: max_incompat_comp = 0     ! and its component
  END TYPE tPcoSystem

  TYPE(tPcoSystem), PUBLIC, SAVE, TARGET :: pco

  ! closed bounds: multiples of the unit roundoff just above the measured rounding population
  REAL*8, PARAMETER, PUBLIC :: PCO_FLUX_CLOSED = 256d0*EPSILON(1d0)   ! 2^-44
  REAL*8, PARAMETER, PUBLIC :: PCO_CLOSURE_TOL = 512d0*EPSILON(1d0)   ! 2^-43
  REAL*8, PARAMETER, PUBLIC :: PCO_FLUX_OPEN = 1d-6
  INTEGER, PARAMETER, PUBLIC :: PCO_GAUGE_UNDECIDED = 5
  ! types 3/4 (P0 coarse system, correction over-relaxed by 3 in crs_P0_Correction) diverge
  ! with a weakly open component (an outlet across cells stretched in z): measured on a 12 x 12
  ! x 120 box with a z outlet - matrix measure of the open component 0.90 (unstretched): 2
  ! cycles per solve; 1.2e-2 (stretch 10) and 1.2e-3 (stretch 30): divergence with types 3
  ! and 4, while type 2 converges and legacy stagnates. Below PCO_WEAK_OPEN the P0 types are
  ! refused (status PCO_P0_WEAK_OPEN), checked with the CURRENT measure at every refresh and
  ! before a type reload (PCO_P0_ADMISSIBLE): at the first renewal fallback to legacy (key
  ! absent) or refusal (explicit key); later (a weight change crossing it) a collective stop
  ! with the same message (both keys); at a reload the previous type is kept, its factors
  ! untouched. The cutoff is empirical: working at 0.90, failing at 1.2e-2; the interval
  ! between is not measured (the ring at 0.11 in the unit test checks the refresh only).
  REAL*8, PARAMETER, PUBLIC :: PCO_WEAK_OPEN = 1d-1
  INTEGER, PARAMETER, PUBLIC :: PCO_P0_WEAK_OPEN = 6
  ! reported right-hand-side imbalance |sum b| / sum |b| over a closed component. The coarse
  ! right-hand side of a converging cycle is a small defect whose sum is rounding of the
  ! divergence cancellation: relative imbalances up to ~2e-2 occur in compatible single-phase
  ! runs at initial defects of 1e-14; a real imbalance (a source that does not integrate to
  ! zero) is O(1). Reported above INCOMPAT_TOL with the absolute net value; the largest value
  ! is kept for the PCO_TIMING lines.
  REAL*8, PARAMETER :: INCOMPAT_TOL = 1d-1

  PUBLIC :: PCO_ACTIVE, PCO_INSTALL_LAYOUT, PCO_CHANGE_SOLVER_TYPE, PCO_REFRESH_SOLVER
  PUBLIC :: PCO_SOLVE_FULL, PCO_SOLVE_P0, PCO_RAW_SOR, PCO_RAW_MATVEC, PCO_BLOCK_SLOT
  PUBLIC :: PCO_NEIGHBOUR_POS, PCO_RELEASE, PCO_SOLVE_CHECK, PCO_PROBE_IMBALANCE
  PUBLIC :: PCO_SET_STRUCTURAL_FLUX, PCO_P0_ADMISSIBLE

CONTAINS

  LOGICAL FUNCTION PCO_ACTIVE()
    PCO_ACTIVE = pco%active
  END FUNCTION PCO_ACTIVE

  !-----------------------------------------------------------------------
  ! Position (1-based) of element f among the neighbours of e, 0 if f is not one.
  !-----------------------------------------------------------------------
  INTEGER FUNCTION PCO_NEIGHBOUR_POS(e, f) RESULT(pos)
    INTEGER, INTENT(IN) :: e, f
    INTEGER :: l, r, mid
    pos = 0
    l = pco%eptr(e)
    r = pco%eptr(e+1) - 1
    DO WHILE (l.LE.r)
      mid = (l + r)/2
      IF (pco%enbr(mid).EQ.f) THEN
        pos = mid - pco%eptr(e) + 1
        RETURN
      ELSE IF (pco%enbr(mid).LT.f) THEN
        l = mid + 1
      ELSE
        r = mid - 1
      END IF
    END DO
  END FUNCTION PCO_NEIGHBOUR_POS

  !-----------------------------------------------------------------------
  ! Slot of entry (mode mp of element e, mode mq of element f) in the raw layout; bpos is
  ! PCO_NEIGHBOUR_POS(e,f) for f /= e and 0 for f = e.
  !-----------------------------------------------------------------------
  INTEGER FUNCTION PCO_BLOCK_SLOT(e, mp, bpos, mq) RESULT(slot)
    INTEGER, INTENT(IN) :: e, mp, bpos, mq
    INTEGER :: row
    row = 4*(e-1) + mp
    IF (bpos.EQ.0) THEN
      IF (mq.EQ.mp) THEN
        slot = pco%ldA(row)
      ELSE IF (mq.LT.mp) THEN
        slot = pco%ldA(row) + mq
      ELSE
        slot = pco%ldA(row) + mq - 1
      END IF
    ELSE
      slot = pco%ldA(row) + 4*bpos + mq - 1
    END IF
  END FUNCTION PCO_BLOCK_SLOT

  !-----------------------------------------------------------------------
  ! Install the layout of level lev from the element graph (eptr/enbr, moved in) and the
  ! solver copy pattern for coarse solver type crs (1-4). The gauge is decided at the first
  ! refresh. Master only.
  !-----------------------------------------------------------------------
  SUBROUTINE PCO_INSTALL_LAYOUT(lev, nel, eptr, enbr, crs, mode)
    INTEGER, INTENT(IN) :: lev, nel, crs, mode
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: eptr(:), enbr(:)
    INTEGER :: e, mp, mq, j, row, k, nnb, f

    CALL PCO_RELEASE()
    pco%lev = lev
    pco%nel = nel
    pco%n = 4*nel
    pco%mode = mode
    CALL MOVE_ALLOC(eptr, pco%eptr)
    CALL MOVE_ALLOC(enbr, pco%enbr)

    ALLOCATE(pco%ldA(pco%n+1))
    pco%ldA(1) = 1
    DO e = 1, nel
      nnb = pco%eptr(e+1) - pco%eptr(e)
      DO mp = 1, 4
        row = 4*(e-1) + mp
        pco%ldA(row+1) = pco%ldA(row) + 4*(nnb+1)
      END DO
    END DO
    ALLOCATE(pco%colA(pco%ldA(pco%n+1)-1), pco%a(pco%ldA(pco%n+1)-1))
    pco%a = 0d0
    DO e = 1, nel
      DO mp = 1, 4
        row = 4*(e-1) + mp
        k = pco%ldA(row)
        pco%colA(k) = row
        DO mq = 1, 4
          IF (mq.EQ.mp) CYCLE
          k = k + 1
          pco%colA(k) = 4*(e-1) + mq
        END DO
        DO j = pco%eptr(e), pco%eptr(e+1)-1
          f = pco%enbr(j)
          DO mq = 1, 4
            k = k + 1
            pco%colA(k) = 4*(f-1) + mq
          END DO
        END DO
      END DO
    END DO
    ALLOCATE(pco%comp_of(nel))
    pco%gauge_decided = .FALSE.
    pco%drift = .FALSE.
    CALL install_copy(crs)

  END SUBROUTINE PCO_INSTALL_LAYOUT

  ! solver copy pattern of type crs (full or mean modes)
  SUBROUTINE install_copy(crs)
    INTEGER, INTENT(IN) :: crs
    INTEGER :: e, j, k, nel

    CALL UmfCSC_Release(pco%lu)
    IF (ALLOCATED(pco%sldA)) DEALLOCATE(pco%sldA, pco%scolA, pco%sfrom, pco%sa, pco%scomp, pco%sref, &
                                        pco%rhs, pco%sol, pco%w)
    nel = pco%nel
    pco%solver_type = crs
    pco%reduced = (crs.EQ.3 .OR. crs.EQ.4)
    IF (pco%reduced) THEN
      pco%ns = nel
      ALLOCATE(pco%sldA(nel+1))
      pco%sldA(1) = 1
      DO e = 1, nel
        pco%sldA(e+1) = pco%sldA(e) + 1 + pco%eptr(e+1) - pco%eptr(e)
      END DO
      ALLOCATE(pco%scolA(pco%sldA(nel+1)-1), pco%sfrom(pco%sldA(nel+1)-1))
      DO e = 1, nel
        k = pco%sldA(e)
        pco%scolA(k) = e
        pco%sfrom(k) = PCO_BLOCK_SLOT(e, 1, 0, 1)
        DO j = pco%eptr(e), pco%eptr(e+1)-1
          k = k + 1
          pco%scolA(k) = pco%enbr(j)
          pco%sfrom(k) = PCO_BLOCK_SLOT(e, 1, j - pco%eptr(e) + 1, 1)
        END DO
      END DO
    ELSE
      pco%ns = pco%n
      ALLOCATE(pco%sldA(pco%n+1), pco%scolA(SIZE(pco%colA)), pco%sfrom(SIZE(pco%colA)))
      pco%sldA = pco%ldA
      pco%scolA = pco%colA
      DO k = 1, SIZE(pco%colA)
        pco%sfrom(k) = k
      END DO
    END IF
    ALLOCATE(pco%sa(SIZE(pco%scolA)), pco%scomp(pco%ns), pco%sref(pco%ns), pco%rhs(pco%ns), &
             pco%sol(pco%ns), pco%w(pco%ns,6))
    pco%sref = .FALSE.
    pco%scomp = 0
    IF (pco%gauge_decided) CALL map_gauge_to_copy()
    pco%analysed = .FALSE.
    IF (crs.EQ.2 .OR. crs.EQ.4) CALL UmfCSC_SetStructure(pco%lu, pco%ns, pco%sldA, pco%scolA)
  END SUBROUTINE install_copy

  ! components and references of the gauge decision onto the unknowns of the copy
  SUBROUTINE map_gauge_to_copy()
    INTEGER :: e, c
    pco%sref = .FALSE.
    IF (pco%reduced) THEN
      pco%scomp = pco%comp_of
    ELSE
      DO e = 1, pco%nel
        pco%scomp(4*(e-1)+1:4*e) = pco%comp_of(e)
      END DO
    END IF
    DO c = 1, pco%ncomp
      IF (.NOT.pco%comp_closed(c)) CYCLE
      e = pco%comp_ref(c)
      pco%sref(MERGE(e, 4*(e-1)+1, pco%reduced)) = .TRUE.
    END DO
  END SUBROUTINE map_gauge_to_copy

  !-----------------------------------------------------------------------
  ! A new coarse solver type 1-4 on the installed layout (runtime parameter reload): new copy
  ! pattern, same gauge decision, refreshed values and factors. Master only.
  !-----------------------------------------------------------------------
  SUBROUTINE PCO_CHANGE_SOLVER_TYPE(crs, status, msg)
    INTEGER, INTENT(IN) :: crs
    INTEGER, INTENT(OUT) :: status
    CHARACTER(len=*), INTENT(OUT) :: msg
    CALL install_copy(crs)
    CALL PCO_REFRESH_SOLVER(status, msg)
  END SUBROUTINE PCO_CHANGE_SOLVER_TYPE

  !-----------------------------------------------------------------------
  ! Structural flux per element (master; the coordinator collects PCO_FREE_FLUX of the
  ! workers at the first renewal, unit tests set it directly). Required by the gauge decision.
  !-----------------------------------------------------------------------
  SUBROUTINE PCO_SET_STRUCTURAL_FLUX(flux)
    REAL*8, INTENT(IN) :: flux(:)
    IF (ALLOCATED(pco%elem_flux)) DEALLOCATE(pco%elem_flux)
    ALLOCATE(pco%elem_flux(SIZE(flux)))
    pco%elem_flux = flux
  END SUBROUTINE PCO_SET_STRUCTURAL_FLUX

  !-----------------------------------------------------------------------
  ! Matrix measures per component c (indicator e_c = 1 on the mean modes of its elements), one
  ! sparse traversal, O(nnz):
  !   right(c) = max over the rows i of c of |(A e_c)_i| / sum_j |A_ij| (all columns),
  !   left(c)  = max over the columns j of c of |(e_c^T A)_j| / sum_i |A_ij| (all rows).
  ! Rows and columns outside c: with components from nonzero connectivity, every nonzero
  ! entry joins two elements of one component, so (A e_c) and (e_c^T A) vanish outside c
  ! exactly; a nonzero entry between two components would contradict the decomposition and is
  ! counted in nleak (the caller treats it as an inconsistency). left is optional (the drift
  ! check needs only right).
  !-----------------------------------------------------------------------
  SUBROUTINE closure_measures(right, nleak, left)
    REAL*8, INTENT(OUT) :: right(:)
    INTEGER, INTENT(OUT) :: nleak
    REAL*8, INTENT(OUT), OPTIONAL :: left(:)
    REAL*8, ALLOCATABLE :: cs(:), ca(:)
    INTEGER :: row, k, c, cr
    REAL*8 :: s, rs
    right = 0d0
    nleak = 0
    IF (PRESENT(left)) THEN
      left = 0d0
      ALLOCATE(cs(pco%n), ca(pco%n))
      cs = 0d0
      ca = 0d0
    END IF
    DO row = 1, pco%n
      cr = pco%comp_of((row-1)/4+1)
      s = 0d0
      rs = 0d0
      DO k = pco%ldA(row), pco%ldA(row+1)-1
        IF (pco%a(k).EQ.0d0) CYCLE
        c = pco%comp_of((pco%colA(k)-1)/4+1)
        IF (c.NE.cr) nleak = nleak + 1
        s = s + ABS(pco%a(k))
        IF (MOD(pco%colA(k)-1,4).EQ.0) rs = rs + pco%a(k)
        IF (PRESENT(left)) THEN
          ca(pco%colA(k)) = ca(pco%colA(k)) + ABS(pco%a(k))
          IF (MOD(row-1,4).EQ.0) cs(pco%colA(k)) = cs(pco%colA(k)) + pco%a(k)
        END IF
      END DO
      IF (s.GT.0d0) right(cr) = MAX(right(cr), ABS(rs)/s)
    END DO
    IF (PRESENT(left)) THEN
      DO k = 1, pco%n
        c = pco%comp_of((k-1)/4+1)
        IF (ca(k).GT.0d0) left(c) = MAX(left(c), ABS(cs(k))/ca(k))
      END DO
      DEALLOCATE(cs, ca)
    END IF
  END SUBROUTINE closure_measures

  !-----------------------------------------------------------------------
  ! Gauge decision (see the module header): components by union-find over the nonzero
  ! element blocks, structure (authoritative) and matrix (consistency) per component,
  ! references. status 0, or PCO_GAUGE_UNDECIDED with msg (nothing pinned).
  !-----------------------------------------------------------------------
  SUBROUTINE decide_gauge(status, msg)
    INTEGER, INTENT(OUT) :: status
    CHARACTER(len=*), INTENT(OUT) :: msg
    INTEGER, ALLOCATABLE :: parent(:), root2c(:)
    REAL*8, ALLOCATABLE :: right(:), left(:), sflux(:)
    INTEGER :: e, f, j, mp, mq, c, nel, pos, nleak, cbad
    REAL*8 :: bmax
    CHARACTER(len=16) :: why

    nel = pco%nel
    ALLOCATE(parent(nel), root2c(nel))
    DO e = 1, nel
      parent(e) = e
    END DO
    DO e = 1, nel
      DO j = pco%eptr(e), pco%eptr(e+1)-1
        f = pco%enbr(j)
        pos = j - pco%eptr(e) + 1
        bmax = 0d0
        DO mp = 1, 4
          DO mq = 1, 4
            bmax = MAX(bmax, ABS(pco%a(PCO_BLOCK_SLOT(e, mp, pos, mq))))
          END DO
        END DO
        IF (bmax.GT.0d0) CALL union(e, f)
      END DO
    END DO
    root2c = 0
    pco%ncomp = 0
    DO e = 1, nel
      f = find(e)
      IF (root2c(f).EQ.0) THEN
        pco%ncomp = pco%ncomp + 1
        root2c(f) = pco%ncomp
      END IF
      pco%comp_of(e) = root2c(f)
    END DO
    IF (ALLOCATED(pco%comp_ref)) DEALLOCATE(pco%comp_ref, pco%comp_closed)
    ALLOCATE(pco%comp_ref(pco%ncomp), pco%comp_closed(pco%ncomp), right(pco%ncomp), left(pco%ncomp))
    pco%comp_ref = 0
    DO e = nel, 1, -1
      pco%comp_ref(pco%comp_of(e)) = e
    END DO
    status = 0
    msg = ' '
    IF (.NOT.ALLOCATED(pco%elem_flux)) THEN
      status = PCO_GAUGE_UNDECIDED
      msg = 'no structural flux of the constant pressure was provided (internal)'
      DEALLOCATE(parent, root2c, right, left)
      RETURN
    END IF
    IF (SIZE(pco%elem_flux).NE.nel) THEN
      status = PCO_GAUGE_UNDECIDED
      msg = 'structural flux of the wrong length (internal)'
      DEALLOCATE(parent, root2c, right, left)
      RETURN
    END IF
    ALLOCATE(sflux(pco%ncomp))
    sflux = 0d0
    DO e = 1, nel
      sflux(pco%comp_of(e)) = MAX(sflux(pco%comp_of(e)), pco%elem_flux(e))
    END DO
    CALL closure_measures(right, nleak, left)
    cbad = 0
    why = ' '
    DO c = 1, pco%ncomp
      IF (sflux(c).LE.PCO_FLUX_CLOSED) THEN
        pco%comp_closed(c) = .TRUE.
        IF (right(c).GT.PCO_CLOSURE_TOL .AND. cbad.EQ.0) THEN
          cbad = c
          why = 'closed-matrix'
        END IF
      ELSE IF (sflux(c).GE.PCO_FLUX_OPEN) THEN
        pco%comp_closed(c) = .FALSE.
        IF (right(c).LE.PCO_CLOSURE_TOL .AND. cbad.EQ.0) THEN
          cbad = c
          why = 'open-matrix'
        END IF
      ELSE
        pco%comp_closed(c) = .FALSE.
        IF (cbad.EQ.0) THEN
          cbad = c
          why = 'flux'
        END IF
      END IF
    END DO
    IF (nleak.GT.0) THEN
      status = PCO_GAUGE_UNDECIDED
      WRITE(msg,'(I0,A)') nleak, ' nonzero entries between two components (internal inconsistency)'
    ELSE IF (cbad.GT.0) THEN
      status = PCO_GAUGE_UNDECIDED
      SELECT CASE (TRIM(why))
      CASE ('flux')
        WRITE(msg,'(A,I0,A,ES9.2,A,ES8.1,A,ES8.1,A)') 'component ', cbad, ': structural flux of a constant pressure ', &
          sflux(cbad), ' between closed (<=', PCO_FLUX_CLOSED, ') and open (>=', PCO_FLUX_OPEN, '): undecidable'
      CASE ('closed-matrix')
        WRITE(msg,'(A,I0,A,ES9.2,A,ES9.2,A)') 'component ', cbad, ': structurally closed (flux ', sflux(cbad), &
          ') but the operator does not annihilate the constant (', right(cbad), ')'
      CASE DEFAULT
        WRITE(msg,'(A,I0,A,ES9.2,A,ES9.2,A)') 'component ', cbad, ': structurally open (flux ', sflux(cbad), &
          ') but the constant is numerically in the kernel of the operator (', right(cbad), ')'
      END SELECT
      msg = TRIM(msg)//'; no reference is set'
    END IF
    IF (status.NE.0) THEN
      pco%comp_closed = .FALSE.
      pco%nclosed = 0
      DEALLOCATE(parent, root2c, right, left, sflux)
      RETURN
    END IF
    pco%nclosed = COUNT(pco%comp_closed)
    pco%closure_closed = 0d0
    pco%closure_open = HUGE(1d0)
    pco%left_closed = 0d0
    pco%flux_closed = 0d0
    pco%flux_open = HUGE(1d0)
    DO c = 1, pco%ncomp
      IF (pco%comp_closed(c)) THEN
        pco%closure_closed = MAX(pco%closure_closed, right(c))
        pco%left_closed = MAX(pco%left_closed, left(c))
        pco%flux_closed = MAX(pco%flux_closed, sflux(c))
      ELSE
        pco%closure_open = MIN(pco%closure_open, right(c))
        pco%flux_open = MIN(pco%flux_open, sflux(c))
      END IF
    END DO
    DEALLOCATE(sflux)
    pco%gauge_decided = .TRUE.
    CALL map_gauge_to_copy()
    DEALLOCATE(parent, root2c, right, left)

  CONTAINS

    INTEGER FUNCTION find(i) RESULT(r)
      INTEGER, INTENT(IN) :: i
      INTEGER :: k, t
      r = i
      DO WHILE (parent(r).NE.r)
        r = parent(r)
      END DO
      k = i
      DO WHILE (parent(k).NE.r)
        t = parent(k)
        parent(k) = r
        k = t
      END DO
    END FUNCTION find

    SUBROUTINE union(i, k)
      INTEGER, INTENT(IN) :: i, k
      INTEGER :: ri, rk
      ri = find(i)
      rk = find(k)
      IF (ri.LT.rk) THEN
        parent(rk) = ri
      ELSE IF (rk.LT.ri) THEN
        parent(ri) = rk
      END IF
    END SUBROUTINE union

  END SUBROUTINE decide_gauge

  ! the decision of decide_gauge against the current values (reported once per drift)
  SUBROUTINE check_gauge()
    REAL*8, ALLOCATABLE :: right(:)
    INTEGER :: c, nleak
    LOGICAL :: bad
    ALLOCATE(right(pco%ncomp))
    CALL closure_measures(right, nleak)
    bad = (nleak.GT.0)
    DO c = 1, pco%ncomp
      IF (pco%comp_closed(c) .AND. right(c).GT.PCO_CLOSURE_TOL) bad = .TRUE.
      IF (.NOT.pco%comp_closed(c) .AND. right(c).LE.PCO_CLOSURE_TOL) bad = .TRUE.
    END DO
    ! the current smallest measure of the open components (admission of the P0 types)
    IF (pco%nclosed.LT.pco%ncomp) &
      pco%closure_open = MINVAL(right, MASK=.NOT.pco%comp_closed(1:pco%ncomp))
    IF (bad .AND. .NOT.pco%drift) WRITE(MTERM,'(A)') 'PCO gauge drift: the closure of a '// &
      'component no longer matches the decision of the first refresh (kept; see pco%drift)'
    pco%drift = pco%drift .OR. bad
    DEALLOCATE(right)
  END SUBROUTINE check_gauge

  !-----------------------------------------------------------------------
  ! May solver type crs run on the current operator? Types 3/4 (P0 system) not with an open
  ! component whose current matrix measure is below PCO_WEAK_OPEN (see there). Master.
  !-----------------------------------------------------------------------
  LOGICAL FUNCTION PCO_P0_ADMISSIBLE(crs) RESULT(ok)
    INTEGER, INTENT(IN) :: crs
    ok = .TRUE.
    IF (crs.NE.3 .AND. crs.NE.4) RETURN
    IF (.NOT.pco%gauge_decided) RETURN
    IF (pco%nclosed.LT.pco%ncomp .AND. pco%closure_open.LT.PCO_WEAK_OPEN) ok = .FALSE.
  END FUNCTION PCO_P0_ADMISSIBLE

  !-----------------------------------------------------------------------
  ! After the raw values changed: gauge decision (first refresh) or its check, solver copy
  ! values, numeric factorisation (types 2, 4). status 0 = ok; msg says what failed. Master.
  !-----------------------------------------------------------------------
  SUBROUTINE PCO_REFRESH_SOLVER(status, msg)
    INTEGER, INTENT(OUT) :: status
    CHARACTER(len=*), INTENT(OUT) :: msg
    INTEGER :: k, row, ist
    REAL*8 :: dref

    status = 0
    msg = ' '
    IF (.NOT.ALL(IEEE_IS_FINITE(pco%a))) THEN
      status = 3
      msg = 'non-finite entry in the coarse pressure operator'
      RETURN
    END IF
    IF (.NOT.pco%gauge_decided) THEN
      CALL decide_gauge(status, msg)
      IF (status.NE.0) RETURN
    ELSE
      CALL check_gauge()
    END IF
    IF (.NOT.PCO_P0_ADMISSIBLE(pco%solver_type)) THEN
      status = PCO_P0_WEAK_OPEN
      WRITE(msg,'(A,I0,A,ES9.2,A,ES8.1,A)') 'coarse solver type ', pco%solver_type, &
        ' (P0 system) with a weakly open component (matrix measure ', pco%closure_open, ' below ', &
        PCO_WEAK_OPEN, '): its over-relaxed P0 correction diverges there; use type 2 (or 1)'
      RETURN
    END IF

    DO k = 1, SIZE(pco%sa)
      pco%sa(k) = pco%a(pco%sfrom(k))
    END DO
    DO row = 1, pco%ns
      DO k = pco%sldA(row)+1, pco%sldA(row+1)-1
        IF (pco%sref(row) .OR. pco%sref(pco%scolA(k))) pco%sa(k) = 0d0
      END DO
    END DO
    ! a reference equation must not read 0 x = 0 (closed component with a vanishing mean row)
    IF (ANY(pco%sref)) THEN
      dref = 0d0
      DO row = 1, pco%ns
        dref = MAX(dref, ABS(pco%sa(pco%sldA(row))))
      END DO
      IF (dref.LE.0d0) dref = 1d0
      DO row = 1, pco%ns
        IF (pco%sref(row) .AND. pco%sa(pco%sldA(row)).EQ.0d0) pco%sa(pco%sldA(row)) = dref
      END DO
    END IF
    ! the SSOR preconditioner of types 1/3 divides by every diagonal
    IF (pco%solver_type.EQ.1 .OR. pco%solver_type.EQ.3) THEN
      DO row = 1, pco%ns
        IF (pco%sa(pco%sldA(row)).EQ.0d0) THEN
          status = 4
          WRITE(msg,'(A,I0,A,I0,A)') 'zero diagonal in row ', row, ' of the coarse pressure copy (solver type ', &
            pco%solver_type, '): the SSOR preconditioner is undefined; use a direct coarse solver (type 2 or 4)'
          RETURN
        END IF
      END DO
    END IF

    IF (pco%solver_type.EQ.2 .OR. pco%solver_type.EQ.4) THEN
      CALL UmfCSC_SetValues(pco%lu, pco%sa)
      IF (.NOT.pco%analysed) THEN
        CALL UmfCSC_Analyse(pco%lu, ist)
        IF (ist.NE.UMF_STATUS_OK) THEN
          status = 1
          WRITE(msg,'(A,I0)') 'UMFPACK symbolic analysis of the coarse pressure matrix failed, status ', ist
          RETURN
        END IF
        pco%analysed = .TRUE.
      END IF
      CALL UmfCSC_Factorize(pco%lu, ist)
      IF (ist.NE.UMF_STATUS_OK) THEN
        status = 2
        WRITE(msg,'(A,I0,A,I0,A,I0,A)') 'UMFPACK factorisation of the gauged coarse pressure matrix failed, status ', &
          ist, ' (1 = singular); components ', pco%ncomp, ', closed ', pco%nclosed, &
          ' (a singular gauged matrix means a closed component was not recognised)'
        RETURN
      END IF
    END IF

  END SUBROUTINE PCO_REFRESH_SOLVER

  !-----------------------------------------------------------------------
  ! Right-hand side of the solver copy: copy (stride 1 or 4), imbalance per closed component,
  ! zero references. diagnostic: no counters, no output; imb returns the largest relative
  ! imbalance |sum b| / sum |b|.
  !-----------------------------------------------------------------------
  SUBROUTINE prepare_rhs(b, stride, diagnostic, imb)
    REAL*8, INTENT(IN) :: b(*)
    INTEGER, INTENT(IN) :: stride
    LOGICAL, INTENT(IN) :: diagnostic
    REAL*8, INTENT(OUT) :: imb
    REAL*8, ALLOCATABLE :: sm(:), sa(:)
    INTEGER :: i, c, cw
    REAL*8 :: rel

    DO i = 1, pco%ns
      pco%rhs(i) = b(stride*(i-1)+1)
    END DO
    imb = 0d0
    IF (pco%nclosed.EQ.0) RETURN
    ALLOCATE(sm(pco%ncomp), sa(pco%ncomp))
    sm = 0d0
    sa = 0d0
    DO i = 1, pco%ns
      IF (.NOT.pco%reduced .AND. MOD(i-1,4).NE.0) CYCLE
      c = pco%scomp(i)
      sm(c) = sm(c) + pco%rhs(i)
      sa(c) = sa(c) + ABS(pco%rhs(i))
    END DO
    cw = 0
    DO c = 1, pco%ncomp
      IF (.NOT.pco%comp_closed(c) .OR. sa(c).LE.0d0) CYCLE
      rel = ABS(sm(c))/sa(c)
      IF (rel.GT.imb) THEN
        imb = rel
        cw = c
      END IF
    END DO
    IF (.NOT.diagnostic .AND. cw.GT.0) THEN
      IF (imb.GT.pco%max_incompat) THEN
        pco%max_incompat = imb
        pco%max_incompat_net = sm(cw)
        pco%max_incompat_comp = cw
      END IF
      IF (imb.GT.INCOMPAT_TOL) THEN
        pco%nincompat = pco%nincompat + 1
        IF (pco%nincompat.LE.5 .OR. MOD(pco%nincompat,1000).EQ.0) &
          WRITE(MTERM,'(A,I0,A,ES12.4,A,ES12.4,A,ES12.4,A,I0,A,I0)') 'PCO_COARSE_RHS component ', cw, &
          ' net ', sm(cw), ' sum of moduli ', sa(cw), ' relative ', imb, ' solve ', pco%nsolve, &
          ' occurrence ', pco%nincompat
      END IF
    END IF
    DEALLOCATE(sm, sa)
    DO i = 1, pco%ns
      IF (pco%sref(i)) pco%rhs(i) = 0d0
    END DO
  END SUBROUTINE prepare_rhs

  !-----------------------------------------------------------------------
  ! Types 1, 2: x (4 nel) = coarse solution of C x = b. nx: length of x at the caller (checked).
  ! iter: in = maximum iterations (type 1), out = used. Master only.
  !-----------------------------------------------------------------------
  SUBROUTINE PCO_SOLVE_FULL(x, b, nx, iter, dcrit, rlx)
    REAL*8, INTENT(INOUT) :: x(*)
    REAL*8, INTENT(IN) :: b(*)
    INTEGER, INTENT(IN) :: nx
    INTEGER, INTENT(INOUT) :: iter
    REAL*8, INTENT(IN) :: dcrit, rlx
    REAL*8 :: imb
    CALL check_dims(nx, .FALSE.)
    pco%nsolve = pco%nsolve + 1
    CALL prepare_rhs(b, 1, .FALSE., imb)
    CALL solve_copy(iter, dcrit, rlx, .FALSE.)
    x(1:pco%ns) = pco%sol
  END SUBROUTINE PCO_SOLVE_FULL

  !-----------------------------------------------------------------------
  ! Types 3, 4: p0 (nel) = solution of the mean-mode system for the mean modes of b.
  !-----------------------------------------------------------------------
  SUBROUTINE PCO_SOLVE_P0(b, p0, iter, dcrit, rlx)
    REAL*8, INTENT(IN) :: b(*)
    REAL*8, INTENT(OUT) :: p0(:)
    INTEGER, INTENT(INOUT) :: iter
    REAL*8, INTENT(IN) :: dcrit, rlx
    REAL*8 :: imb
    CALL check_dims(SIZE(p0), .TRUE.)
    pco%nsolve = pco%nsolve + 1
    CALL prepare_rhs(b, 4, .FALSE., imb)
    CALL solve_copy(iter, dcrit, rlx, .FALSE.)
    p0(1:pco%ns) = pco%sol
  END SUBROUTINE PCO_SOLVE_P0

  ! the solver kind and length of the caller against the installed copy (a structural reload that
  ! bypassed the coordinator would index past the arrays)
  SUBROUTINE check_dims(n, reduced)
    INTEGER, INTENT(IN) :: n
    LOGICAL, INTENT(IN) :: reduced
    INTEGER :: ierr
    IF (n.EQ.pco%ns .AND. (reduced.EQV.pco%reduced)) RETURN
    WRITE(*,'(A,I0,A,L1,A,I0,A,I0)') 'PCO: coarse solve called with length ', n, ' reduced ', reduced, &
      ' but the installed copy is of solver type ', pco%solver_type, ' and length ', pco%ns
    FLUSH(6)
    CALL MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  END SUBROUTINE check_dims

  SUBROUTINE solve_copy(iter, dcrit, rlx, diagnostic)
    INTEGER, INTENT(INOUT) :: iter
    REAL*8, INTENT(IN) :: dcrit, rlx
    LOGICAL, INTENT(IN) :: diagnostic
    INTEGER :: ist, ierr
    pco%sol = 0d0
    IF (pco%solver_type.EQ.2 .OR. pco%solver_type.EQ.4) THEN
      CALL UmfCSC_Solve(pco%lu, pco%sol, pco%rhs, ist)
      IF (ist.LT.0 .OR. ist.EQ.1) THEN
        WRITE(*,'(A,I0,A,I0)') 'PCO: UMFPACK solve of the coarse pressure system failed on the master, status ', &
          ist, ', solve ', pco%nsolve
        FLUSH(6)
        CALL MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      END IF
      iter = 1
    ELSE
      CALL bicgstab(pco%sol, pco%rhs, pco%ns, iter, dcrit, rlx, ist)
      IF (ist.NE.0) THEN
        IF (.NOT.diagnostic) pco%nbreakdown = pco%nbreakdown + 1
        WRITE(*,'(A,I0,A,I0,A,I0)') 'PCO: BiCGStab breakdown in the coarse pressure solve on the master, type ', &
          pco%solver_type, ', status ', ist, ', solve ', pco%nsolve
        FLUSH(6)
        CALL MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      END IF
    END IF
  END SUBROUTINE solve_copy

  !-----------------------------------------------------------------------
  ! Left-preconditioned BiCGStab on the gauged copy (the iteration of the legacy master
  ! BiCGStab, E012_BiCGStabSolverMaster, with the guards of the module header).
  ! x: in/out (zero start), b: rhs, iter: in = maximum, out = used; crit: relative reduction
  ! of the preconditioned residual; status 0 converged or at the maximum, 1 rho breakdown,
  ! 2 (r0, v) breakdown, 3 non-finite.
  !-----------------------------------------------------------------------
  SUBROUTINE bicgstab(x, b, n, iter, crit, rlx, status)
    INTEGER, INTENT(IN) :: n
    REAL*8, INTENT(INOUT) :: x(n)
    REAL*8, INTENT(IN) :: b(n), crit, rlx
    INTEGER, INTENT(INOUT) :: iter
    INTEGER, INTENT(OUT) :: status
    REAL*8 :: rho0, rho1, alpha, omega, beta, def0, def, denom, tt
    INTEGER :: it, maxit

    ASSOCIATE(r => pco%w(:,1), r0 => pco%w(:,2), p => pco%w(:,3), v => pco%w(:,4), &
              s => pco%w(:,5), t => pco%w(:,6))
    status = 0
    maxit = iter
    CALL copy_matvec(x, r)
    r = b - r
    CALL ssor(r, rlx)
    def0 = NORM2(r)
    iter = 0
    IF (def0.EQ.0d0) RETURN                     ! zero right-hand side: zero correction
    r0 = r
    p = 0d0
    v = 0d0
    rho0 = 1d0
    alpha = 1d0
    omega = 1d0
    DO it = 1, maxit
      iter = it
      rho1 = DOT_PRODUCT(r0, r)
      IF (rho1.EQ.0d0) THEN
        status = 1
        EXIT
      END IF
      beta = (rho1/rho0)*(alpha/omega)
      p = r + beta*(p - omega*v)
      CALL copy_matvec(p, v)
      CALL ssor(v, rlx)
      denom = DOT_PRODUCT(r0, v)
      IF (denom.EQ.0d0) THEN
        status = 2
        EXIT
      END IF
      alpha = rho1/denom
      s = r - alpha*v
      IF (NORM2(s).LE.crit*def0) THEN            ! converged at the half step
        x = x + alpha*p
        EXIT
      END IF
      CALL copy_matvec(s, t)
      CALL ssor(t, rlx)
      tt = DOT_PRODUCT(t, t)
      IF (tt.EQ.0d0) THEN
        x = x + alpha*p
        status = 2
        EXIT
      END IF
      omega = DOT_PRODUCT(t, s)/tt
      x = x + alpha*p + omega*s
      r = s - omega*t
      def = NORM2(r)
      IF (def.NE.def) THEN
        status = 3
        EXIT
      END IF
      IF (def.LE.crit*def0) EXIT
      IF (omega.EQ.0d0) THEN
        status = 1
        EXIT
      END IF
      rho0 = rho1
    END DO
    END ASSOCIATE
  END SUBROUTINE bicgstab

  ! y = gauged copy * x
  SUBROUTINE copy_matvec(x, y)
    REAL*8, INTENT(IN) :: x(:)
    REAL*8, INTENT(OUT) :: y(:)
    INTEGER :: i, k
    DO i = 1, pco%ns
      y(i) = 0d0
      DO k = pco%sldA(i), pco%sldA(i+1)-1
        y(i) = y(i) + pco%sa(k)*x(pco%scolA(k))
      END DO
    END DO
  END SUBROUTINE copy_matvec

  ! SSOR sweep pair of PARID117 on the copy, every off-diagonal entry scanned (lower: column
  ! below the row, upper: above), so the diagonal-first layout loses no lower neighbour
  SUBROUTINE ssor(x, omega)
    REAL*8, INTENT(INOUT) :: x(:)
    REAL*8, INTENT(IN) :: omega
    INTEGER :: i, k, j
    REAL*8 :: aux
    DO i = 1, pco%ns
      aux = 0d0
      DO k = pco%sldA(i)+1, pco%sldA(i+1)-1
        j = pco%scolA(k)
        IF (j.LT.i) aux = aux + pco%sa(k)*x(j)
      END DO
      x(i) = (x(i) - aux*omega)/pco%sa(pco%sldA(i))
    END DO
    DO i = pco%ns-1, 1, -1
      aux = 0d0
      DO k = pco%sldA(i)+1, pco%sldA(i+1)-1
        j = pco%scolA(k)
        IF (j.GT.i) aux = aux + pco%sa(k)*x(j)
      END DO
      x(i) = x(i) - aux*omega/pco%sa(pco%sldA(i))
    END DO
  END SUBROUTINE ssor

  !-----------------------------------------------------------------------
  ! SOR sweeps with the raw (ungauged) operator: the smoothing of the type-3/4 coarse
  ! correction, the sweep of crs_E012_SOR_Smoother on the legacy matrix (2 nsweep forward
  ! sweeps, diagonal first in every row). A row with a zero diagonal (an identically zero row
  ! of a closed component, see the gauge) keeps its value.
  !-----------------------------------------------------------------------
  SUBROUTINE PCO_RAW_SOR(x, b, nsweep, omega)
    REAL*8, INTENT(INOUT) :: x(*)
    REAL*8, INTENT(IN) :: b(*)
    INTEGER, INTENT(IN) :: nsweep
    REAL*8, INTENT(IN) :: omega
    INTEGER :: i, ieq, k
    REAL*8 :: aux, d
    DO i = 1, 2*nsweep
      DO ieq = 1, pco%n
        d = pco%a(pco%ldA(ieq))
        IF (d.EQ.0d0) CYCLE
        aux = b(ieq)
        DO k = pco%ldA(ieq)+1, pco%ldA(ieq+1)-1
          aux = aux - pco%a(k)*x(pco%colA(k))
        END DO
        x(ieq) = omega*(aux/d - x(ieq)) + x(ieq)
      END DO
    END DO
  END SUBROUTINE PCO_RAW_SOR

  ! y = raw * x
  SUBROUTINE PCO_RAW_MATVEC(x, y)
    REAL*8, INTENT(IN) :: x(:)
    REAL*8, INTENT(OUT) :: y(:)
    INTEGER :: i, k
    DO i = 1, pco%n
      y(i) = 0d0
      DO k = pco%ldA(i), pco%ldA(i+1)-1
        y(i) = y(i) + pco%a(k)*x(pco%colA(k))
      END DO
    END DO
  END SUBROUTINE PCO_RAW_MATVEC

  !-----------------------------------------------------------------------
  ! Diagnostic check of the solver copy (no production statistics): for a compatible right-
  ! hand side b = raw * v the UNGAUGED residual |b - raw x| / |b| (types 2, 4: rounding level;
  ! types 1, 3: the iteration tolerance 1e-11 amplified by the conditioning), the deviation of
  ! x - v from a constant per closed component, relative to |v|, the value of the reference
  ! unknowns relative to |x|.
  !-----------------------------------------------------------------------
  SUBROUTINE PCO_SOLVE_CHECK(v, resid, devi, xref)
    REAL*8, INTENT(IN) :: v(:)
    REAL*8, INTENT(OUT) :: resid, devi, xref
    REAL*8, ALLOCATABLE :: b(:), x(:), y(:), cm(:), cn(:), vs(:)
    INTEGER :: i, c, it
    REAL*8 :: d, imb

    ALLOCATE(b(pco%ns), x(pco%ns), y(pco%ns), vs(pco%ns))
    IF (pco%reduced) THEN
      DO i = 1, pco%ns
        vs(i) = v(4*(i-1)+1)
      END DO
    ELSE
      vs = v(1:pco%ns)
    END IF
    CALL copy_raw_matvec(vs, b)
    it = 2000
    CALL prepare_rhs(b, 1, .TRUE., imb)
    CALL solve_copy(it, 1d-11, 1d0, .TRUE.)
    x = pco%sol
    CALL copy_raw_matvec(x, y)
    resid = NORM2(b - y)/MAX(NORM2(b), 1d-300)
    ALLOCATE(cm(pco%ncomp), cn(pco%ncomp))
    cm = 0d0
    cn = 0d0
    DO i = 1, pco%ns
      IF (.NOT.pco%reduced .AND. MOD(i-1,4).NE.0) CYCLE
      c = pco%scomp(i)
      cm(c) = cm(c) + (x(i) - vs(i))
      cn(c) = cn(c) + 1d0
    END DO
    cm = cm/MAX(cn, 1d0)
    devi = 0d0
    DO i = 1, pco%ns
      d = x(i) - vs(i)
      IF ((pco%reduced .OR. MOD(i-1,4).EQ.0) .AND. pco%comp_closed(pco%scomp(i))) d = d - cm(pco%scomp(i))
      devi = MAX(devi, ABS(d))
    END DO
    devi = devi/MAX(MAXVAL(ABS(vs)), 1d-300)
    xref = 0d0
    IF (ANY(pco%sref)) xref = MAXVAL(ABS(x), MASK=pco%sref)/MAX(MAXVAL(ABS(x)), 1d-300)
    DEALLOCATE(b, x, y, cm, cn, vs)
  END SUBROUTINE PCO_SOLVE_CHECK

  !-----------------------------------------------------------------------
  ! Diagnostic: largest relative imbalance of b (full length 4 nel) over the closed
  ! components, as the production solve would measure it (no solve, no statistics).
  !-----------------------------------------------------------------------
  REAL*8 FUNCTION PCO_PROBE_IMBALANCE(b) RESULT(imb)
    REAL*8, INTENT(IN) :: b(:)
    CALL prepare_rhs(b, MERGE(4, 1, pco%reduced), .TRUE., imb)
  END FUNCTION PCO_PROBE_IMBALANCE

  ! raw operator restricted to the solver unknowns (mean modes for the reduced system)
  SUBROUTINE copy_raw_matvec(x, y)
    REAL*8, INTENT(IN) :: x(:)
    REAL*8, INTENT(OUT) :: y(:)
    INTEGER :: i, k
    DO i = 1, pco%ns
      y(i) = 0d0
      DO k = pco%sldA(i), pco%sldA(i+1)-1
        y(i) = y(i) + pco%a(pco%sfrom(k))*x(pco%scolA(k))
      END DO
    END DO
  END SUBROUTINE copy_raw_matvec

  SUBROUTINE PCO_RELEASE()
    CALL UmfCSC_Release(pco%lu)
    IF (ALLOCATED(pco%eptr)) DEALLOCATE(pco%eptr)
    IF (ALLOCATED(pco%enbr)) DEALLOCATE(pco%enbr)
    IF (ALLOCATED(pco%ldA)) DEALLOCATE(pco%ldA)
    IF (ALLOCATED(pco%colA)) DEALLOCATE(pco%colA)
    IF (ALLOCATED(pco%a)) DEALLOCATE(pco%a)
    IF (ALLOCATED(pco%comp_of)) DEALLOCATE(pco%comp_of)
    IF (ALLOCATED(pco%comp_ref)) DEALLOCATE(pco%comp_ref)
    IF (ALLOCATED(pco%comp_closed)) DEALLOCATE(pco%comp_closed)
    IF (ALLOCATED(pco%elem_flux)) DEALLOCATE(pco%elem_flux)
    IF (ALLOCATED(pco%sldA)) DEALLOCATE(pco%sldA)
    IF (ALLOCATED(pco%scolA)) DEALLOCATE(pco%scolA)
    IF (ALLOCATED(pco%sfrom)) DEALLOCATE(pco%sfrom)
    IF (ALLOCATED(pco%sa)) DEALLOCATE(pco%sa)
    IF (ALLOCATED(pco%scomp)) DEALLOCATE(pco%scomp)
    IF (ALLOCATED(pco%sref)) DEALLOCATE(pco%sref)
    IF (ALLOCATED(pco%rhs)) DEALLOCATE(pco%rhs)
    IF (ALLOCATED(pco%sol)) DEALLOCATE(pco%sol)
    IF (ALLOCATED(pco%w)) DEALLOCATE(pco%w)
    pco%active = .FALSE.
    pco%analysed = .FALSE.
    pco%gauge_decided = .FALSE.
  END SUBROUTINE PCO_RELEASE

END MODULE pco_coarse_solver
