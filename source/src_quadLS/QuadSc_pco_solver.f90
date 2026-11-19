!=========================================================================
! QuadSc_pco_solver.f90 - storage, gauge and factorisation of the periodic coarse pressure
! operator on the master
!
! Part of the periodic coarse pressure operator (overview: QuadSc_pco_operator.f90). Owns
! the master-side system of a non-legacy mode; the coordinator writes the raw values, the
! multigrid (QuadSc_mg.f90) solves with it. Never entered in legacy mode (PCO_ACTIVE false).
!
! Numbering and layout (the conventions of the legacy master matrix, Get_CMatStruct iP=1):
!   unknown 4(e-1)+m is mode m of the P1 pressure of master element e of the coarse level
!   (m = 1 the element mean, 2..4 the slopes); the row block of e lists the own block first,
!   then the blocks of the neighbour elements in ascending order (element graph through
!   vertex classes, pco_periodic_topology); row 4(e-1)+m starts with its diagonal, then the
!   other modes of e ascending, then 4 entries per neighbour. Hence the P0 (mode 1) system
!   is read with stride 4, as the legacy type-3/4 extraction does.
!
! raw  : the operator as the workers apply it (gathered) or its master-side equivalent
!        (wrapped); ungauged and never modified by a solve. It is the operator of the type-3/4
!        smoothing (PCO_RAW_SOR) and of the checks.
! gauge: on every connected component of the element graph whose operator annihilates the
!        constant pressure (C e = 0, e = 1 on the mean modes: a closed fluid domain) one
!        reference unknown, the mean mode of the lowest element of the component, is set to
!        zero in the SOLVER COPY only: its row and column are cleared except the diagonal and
!        its right-hand side entry is zeroed (symmetric elimination of a zero value). The
!        decision is taken from the operator at every refresh (closure measure below 1e-8),
!        not from SimPar@NoOutflow; a disagreement with NoOutflow is reported by the
!        coordinator. A right-hand side whose component sums do not vanish (incompatible with
!        the closed component) is reported, not changed: the gauge cannot repair a net
!        continuity imbalance.
! solver copy: full (types 1, 2) or mean modes only (types 3, 4), same diagonal-first layout,
!        pattern fixed (explicit zeros kept), values refreshed from raw. Types 2/4: UMFPACK in
!        explicit column storage of this matrix (tUmfCSC: it solves A x = b, not A^T x = b);
!        the symbolic analysis is done once per pattern, the numeric factorisation at every
!        refresh. Types 1/3: BiCGStab of the master with the SSOR preconditioner of the legacy
!        path, on the gauged copy.
! Status: factorisation failures are returned to the coordinator (collective stop there); a
! failing solve inside the multigrid runs on the master only and aborts with a message.
!=========================================================================
MODULE pco_coarse_solver

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
    ! components of the element graph and the gauge
    INTEGER :: ncomp = 0, nclosed = 0
    INTEGER, ALLOCATABLE :: comp_of(:), comp_ref(:)
    LOGICAL, ALLOCATABLE :: comp_closed(:)
    REAL*8  :: closure_open = 0d0        ! smallest closure measure of the open components
    REAL*8  :: closure_closed = 0d0      ! largest closure measure of the closed components
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
    ! statistics
    INTEGER :: nsolve = 0, nincompat = 0
    REAL*8  :: max_incompat = 0d0
  END TYPE tPcoSystem

  TYPE(tPcoSystem), PUBLIC, SAVE, TARGET :: pco

  REAL*8, PARAMETER, PUBLIC :: PCO_CLOSURE_TOL = 1d-8   ! closed component: |C e| / |C| e below
  ! reported right-hand-side imbalance |sum b| / sum |b| over a closed component. The coarse
  ! right-hand side of a converging cycle is a small defect whose sum is rounding of the
  ! divergence cancellation: relative imbalances up to ~1e-3 occur in compatible runs
  ! (single-phase probe); a real imbalance (a source that does not integrate to zero) is O(1).
  ! The largest value is kept (pco%max_incompat, printed on the PCO_TIMING lines).
  REAL*8, PARAMETER :: INCOMPAT_TOL = 1d-2
  REAL*8, SAVE :: pco_rlx = 1d0                         ! relaxation of the type-1/3 preconditioner

  PUBLIC :: PCO_ACTIVE, PCO_REDUCED, PCO_INSTALL_LAYOUT, PCO_REFRESH_SOLVER, PCO_SOLVE_FULL
  PUBLIC :: PCO_SOLVE_P0, PCO_RAW_SOR, PCO_RAW_MATVEC, PCO_BLOCK_SLOT, PCO_NEIGHBOUR_POS
  PUBLIC :: PCO_RELEASE, PCO_SOLVE_CHECK

CONTAINS

  LOGICAL FUNCTION PCO_ACTIVE()
    PCO_ACTIVE = pco%active
  END FUNCTION PCO_ACTIVE

  LOGICAL FUNCTION PCO_REDUCED()
    PCO_REDUCED = pco%reduced
  END FUNCTION PCO_REDUCED

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
  ! Install the layout of level lev from the element graph (eptr/enbr, moved in), the
  ! components and the solver copy pattern for coarse solver type crs (1-4). Master only.
  !-----------------------------------------------------------------------
  SUBROUTINE PCO_INSTALL_LAYOUT(lev, nel, eptr, enbr, crs, mode)
    INTEGER, INTENT(IN) :: lev, nel, crs, mode
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: eptr(:), enbr(:)
    INTEGER :: e, mp, mq, j, row, k, nnb, ic, f, head, tail
    INTEGER, ALLOCATABLE :: queue(:)

    CALL PCO_RELEASE()
    pco%lev = lev
    pco%nel = nel
    pco%n = 4*nel
    pco%solver_type = crs
    pco%mode = mode
    CALL MOVE_ALLOC(eptr, pco%eptr)
    CALL MOVE_ALLOC(enbr, pco%enbr)

    ! raw layout
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

    ! connected components of the element graph; reference = lowest element
    ALLOCATE(pco%comp_of(nel), queue(nel))
    pco%comp_of = 0
    pco%ncomp = 0
    DO e = 1, nel
      IF (pco%comp_of(e).NE.0) CYCLE
      pco%ncomp = pco%ncomp + 1
      pco%comp_of(e) = pco%ncomp
      head = 1
      tail = 1
      queue(1) = e
      DO WHILE (head.LE.tail)
        ic = queue(head)
        head = head + 1
        DO j = pco%eptr(ic), pco%eptr(ic+1)-1
          f = pco%enbr(j)
          IF (pco%comp_of(f).EQ.0) THEN
            pco%comp_of(f) = pco%ncomp
            tail = tail + 1
            queue(tail) = f
          END IF
        END DO
      END DO
    END DO
    DEALLOCATE(queue)
    ALLOCATE(pco%comp_ref(pco%ncomp), pco%comp_closed(pco%ncomp))
    pco%comp_ref = 0
    DO e = nel, 1, -1
      pco%comp_ref(pco%comp_of(e)) = e
    END DO
    pco%comp_closed = .FALSE.

    ! solver copy pattern
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
      ALLOCATE(pco%scomp(nel))
      pco%scomp = pco%comp_of
    ELSE
      pco%ns = pco%n
      ALLOCATE(pco%sldA(pco%n+1), pco%scolA(SIZE(pco%colA)), pco%sfrom(SIZE(pco%colA)))
      pco%sldA = pco%ldA
      pco%scolA = pco%colA
      DO k = 1, SIZE(pco%colA)
        pco%sfrom(k) = k
      END DO
      ALLOCATE(pco%scomp(pco%n))
      DO e = 1, nel
        pco%scomp(4*(e-1)+1:4*e) = pco%comp_of(e)
      END DO
    END IF
    ALLOCATE(pco%sa(SIZE(pco%scolA)), pco%sref(pco%ns), pco%rhs(pco%ns), pco%sol(pco%ns))
    ALLOCATE(pco%w(pco%ns,5))
    pco%sref = .FALSE.
    pco%analysed = .FALSE.
    IF (crs.EQ.2 .OR. crs.EQ.4) CALL UmfCSC_SetStructure(pco%lu, pco%ns, pco%sldA, pco%scolA)

  END SUBROUTINE PCO_INSTALL_LAYOUT

  !-----------------------------------------------------------------------
  ! After the raw values changed: closure of each component, gauge, solver copy values,
  ! numeric factorisation (types 2, 4). status 0 = ok; msg says what failed. Master only.
  !-----------------------------------------------------------------------
  SUBROUTINE PCO_REFRESH_SOLVER(status, msg)
    INTEGER, INTENT(OUT) :: status
    CHARACTER(len=*), INTENT(OUT) :: msg
    INTEGER :: e, c, k, row, g, ist
    REAL*8, ALLOCATABLE :: worst(:)
    REAL*8 :: r, s

    status = 0
    msg = ' '
    ! closure measure per component: max over the rows of |sum of the mean-mode columns| /
    ! sum of their moduli (C e with e = 1 on the mean modes)
    ALLOCATE(worst(pco%ncomp))
    worst = 0d0
    DO row = 1, pco%n
      r = 0d0
      s = 0d0
      DO k = pco%ldA(row), pco%ldA(row+1)-1
        IF (MOD(pco%colA(k)-1,4).NE.0) CYCLE
        r = r + pco%a(k)
        s = s + ABS(pco%a(k))
      END DO
      c = pco%comp_of((row-1)/4+1)
      IF (s.GT.0d0) worst(c) = MAX(worst(c), ABS(r)/s)
    END DO
    pco%comp_closed = (worst.LE.PCO_CLOSURE_TOL)
    pco%nclosed = COUNT(pco%comp_closed)
    pco%closure_closed = 0d0
    pco%closure_open = HUGE(1d0)
    DO c = 1, pco%ncomp
      IF (pco%comp_closed(c)) THEN
        pco%closure_closed = MAX(pco%closure_closed, worst(c))
      ELSE
        pco%closure_open = MIN(pco%closure_open, worst(c))
      END IF
    END DO
    DEALLOCATE(worst)

    ! solver copy with the gauge
    DO k = 1, SIZE(pco%sa)
      pco%sa(k) = pco%a(pco%sfrom(k))
    END DO
    pco%sref = .FALSE.
    DO c = 1, pco%ncomp
      IF (.NOT.pco%comp_closed(c)) CYCLE
      e = pco%comp_ref(c)
      g = MERGE(e, 4*(e-1)+1, pco%reduced)
      pco%sref(g) = .TRUE.
    END DO
    DO row = 1, pco%ns
      DO k = pco%sldA(row)+1, pco%sldA(row+1)-1     ! off-diagonal entries
        IF (pco%sref(row) .OR. pco%sref(pco%scolA(k))) pco%sa(k) = 0d0
      END DO
    END DO

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
  ! Right-hand side of the solver copy: copy, compatibility measure per closed component
  ! (|sum| / sum of moduli over its mean modes; reported above INCOMPAT_TOL), zero references.
  !-----------------------------------------------------------------------
  SUBROUTINE prepare_rhs(b, stride)
    REAL*8, INTENT(IN) :: b(*)
    INTEGER, INTENT(IN) :: stride
    REAL*8, ALLOCATABLE :: sm(:), sa(:)
    INTEGER :: i, c
    REAL*8 :: worst

    DO i = 1, pco%ns
      pco%rhs(i) = b(stride*(i-1)+1)
    END DO
    IF (pco%nclosed.GT.0) THEN
      ALLOCATE(sm(pco%ncomp), sa(pco%ncomp))
      sm = 0d0
      sa = 0d0
      DO i = 1, pco%ns
        IF (.NOT.pco%reduced .AND. MOD(i-1,4).NE.0) CYCLE
        c = pco%scomp(i)
        sm(c) = sm(c) + pco%rhs(i)
        sa(c) = sa(c) + ABS(pco%rhs(i))
      END DO
      worst = 0d0
      DO c = 1, pco%ncomp
        IF (pco%comp_closed(c) .AND. sa(c).GT.0d0) worst = MAX(worst, ABS(sm(c))/sa(c))
      END DO
      DEALLOCATE(sm, sa)
      pco%max_incompat = MAX(pco%max_incompat, worst)
      IF (worst.GT.INCOMPAT_TOL) THEN
        pco%nincompat = pco%nincompat + 1
        IF (pco%nincompat.LE.5 .OR. MOD(pco%nincompat,1000).EQ.0) &
          WRITE(MTERM,'(A,ES12.4,A,I0,A,I0)') 'PCO_COARSE_RHS incompatible with a closed component: '// &
          'imbalance ', worst, ' solve ', pco%nsolve, ' occurrence ', pco%nincompat
      END IF
      DO i = 1, pco%ns
        IF (pco%sref(i)) pco%rhs(i) = 0d0
      END DO
    END IF
  END SUBROUTINE prepare_rhs

  !-----------------------------------------------------------------------
  ! Types 1, 2: x (4 nel) = coarse solution of C x = b. iter: iterations used (type 1: in =
  ! maximum, out = used). Master only.
  !-----------------------------------------------------------------------
  SUBROUTINE PCO_SOLVE_FULL(x, b, iter, dcrit, rlx)
    REAL*8, INTENT(INOUT) :: x(*)
    REAL*8, INTENT(IN) :: b(*)
    INTEGER, INTENT(INOUT) :: iter
    REAL*8, INTENT(IN) :: dcrit, rlx
    pco%nsolve = pco%nsolve + 1
    CALL prepare_rhs(b, 1)
    CALL solve_copy(iter, dcrit, rlx)
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
    pco%nsolve = pco%nsolve + 1
    CALL prepare_rhs(b, 4)
    CALL solve_copy(iter, dcrit, rlx)
    p0(1:pco%ns) = pco%sol
  END SUBROUTINE PCO_SOLVE_P0

  SUBROUTINE solve_copy(iter, dcrit, rlx)
    INTEGER, INTENT(INOUT) :: iter
    REAL*8, INTENT(IN) :: dcrit, rlx
    INTEGER :: ist, ierr
    EXTERNAL E012_BiCGStabSolverMaster
    pco%sol = 0d0
    pco_rlx = rlx
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
      ! the first iteration multiplies the direction vectors by zero: they must be finite
      pco%w = 0d0
      CALL E012_BiCGStabSolverMaster(pco%sol, pco%rhs, pco%ns, iter, copy_dax, copy_prec, .TRUE., &
        pco%w(:,1), pco%w(:,2), pco%w(:,3), pco%w(:,4), pco%w(:,5), dcrit)
      IF (ANY(pco%sol.NE.pco%sol)) THEN
        WRITE(*,'(A,I0,A,I0)') 'PCO: BiCGStab breakdown (NaN) in the coarse pressure solve on the master, type ', &
          pco%solver_type, ', solve ', pco%nsolve
        FLUSH(6)
        CALL MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      END IF
    END IF
  END SUBROUTINE solve_copy

  ! matrix-vector product and preconditioner of the gauged copy (BiCGStab callbacks)
  SUBROUTINE copy_dax(dx, da, neq, a1, a2)
    REAL*8 :: dx(*), da(*), a1, a2
    INTEGER :: neq, i, k
    DO i = 1, neq
      da(i) = a2*da(i)
      DO k = pco%sldA(i), pco%sldA(i+1)-1
        da(i) = da(i) + a1*pco%sa(k)*dx(pco%scolA(k))
      END DO
    END DO
  END SUBROUTINE copy_dax

  SUBROUTINE copy_prec(dx, neq)
    REAL*8 :: dx(*)
    INTEGER :: neq
    EXTERNAL PARID117_MASTER
    CALL PARID117_MASTER(pco%sa, pco%scolA, pco%sldA, dx, neq, pco_rlx)
  END SUBROUTINE copy_prec

  !-----------------------------------------------------------------------
  ! SOR sweeps with the raw (ungauged) operator: the smoothing of the type-3/4 coarse
  ! correction, the arithmetic of crs_E012_SOR on the legacy matrix.
  !-----------------------------------------------------------------------
  SUBROUTINE PCO_RAW_SOR(x, b, nsweep, omega)
    REAL*8, INTENT(INOUT) :: x(*)
    REAL*8, INTENT(IN) :: b(*)
    INTEGER, INTENT(IN) :: nsweep
    REAL*8, INTENT(IN) :: omega
    INTEGER :: i
    REAL*8 :: om
    EXTERNAL crs_E012_SOR_Smoother
    om = omega
    DO i = 1, 2*nsweep
      CALL crs_E012_SOR_Smoother(b, pco%a, pco%colA, pco%ldA, x, pco%n, om)
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
  ! Check of the solver copy (start-up check of the coordinator): for a compatible right-hand
  ! side b = raw * v the solve must give a solution whose UNGAUGED residual |b - raw x| / |b|
  ! is at rounding level (types 2, 4, direct) or at the iteration tolerance 1e-11 (types 1, 3,
  ! BiCGStab; a tighter tolerance runs the legacy BiCGStab into breakdown) (types 1, 2: full
  ! system; types 3, 4: mean-mode system), and the
  ! solution differs from v by a constant per closed component only (returned as the largest
  ! deviation of x - v from its component mean, relative to |v|).
  !-----------------------------------------------------------------------
  SUBROUTINE PCO_SOLVE_CHECK(v, resid, devi, incompat_seen)
    REAL*8, INTENT(IN) :: v(:)
    REAL*8, INTENT(OUT) :: resid, devi
    LOGICAL, INTENT(OUT) :: incompat_seen
    REAL*8, ALLOCATABLE :: b(:), x(:), y(:), cm(:), cn(:), vs(:)
    INTEGER :: i, c, it, nin
    REAL*8 :: dummy

    ALLOCATE(b(pco%ns), x(pco%ns), y(pco%ns), vs(pco%ns))
    IF (pco%reduced) THEN
      DO i = 1, pco%ns
        vs(i) = v(4*(i-1)+1)
      END DO
    ELSE
      vs = v(1:pco%ns)
    END IF
    CALL copy_raw_matvec(vs, b)
    nin = pco%nincompat
    it = 2000
    CALL prepare_rhs(b, 1)
    CALL solve_copy(it, 1d-11, 1d0)
    x = pco%sol
    incompat_seen = (pco%nincompat.GT.nin)
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
      dummy = x(i) - vs(i)
      IF (pco%reduced .OR. MOD(i-1,4).EQ.0) dummy = dummy - cm(pco%scomp(i))
      devi = MAX(devi, ABS(dummy))
    END DO
    devi = devi/MAX(MAXVAL(ABS(vs)), 1d-300)
    DEALLOCATE(b, x, y, cm, cn, vs)
  END SUBROUTINE PCO_SOLVE_CHECK

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
  END SUBROUTINE PCO_RELEASE

END MODULE pco_coarse_solver
