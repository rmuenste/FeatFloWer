PROGRAM TEST_PCO_COARSE

  ! Unit tests of the periodic coarse pressure operator (QuadSc_pco_*.f90) that need no MPI
  ! topology; one rank, MPI_Init only for the harness. The gathered operator, the link census
  ! and the operator identity against the worker product need the partitioned solver and are
  ! checked in the application (Pres@PeriodicCoarseCheck, Slurm probes; see
  ! docs/md_docs/pressure_periodic_coarse_operator_2026-11-19.md).
  !
  !   1 UMFPACK orientation: a non-symmetric matrix in CSR; the explicit adapter (tUmfCSC)
  !     solves A x = b, the legacy wrapper (myUmfPack_*) solves A^T x = b (documented
  !     convention of the certified path), statuses are returned
  !   2 periodic classes on the structured unit cube (Q2 dofs of level 1 at the vertices of
  !     level 2): sizes per kind for x, xy, xyz periodicity (face 2, edge 4, corner 8 images),
  !     no refusal on a valid mesh; a corner image moved off its place -> incomplete class; a
  !     duplicated dof -> ambiguous; unequal Dirichlet flags over a class -> conflict
  !   3 element graph through the classes: neighbour counts on the xy-periodic cube
  !   4 P0 -> vertex averaging (types 3/4) on a non-constant periodic field: every image of a
  !     periodic vertex gets the class value sum(num)/sum(den); vertices without images are
  !     bit for bit those of INTPVB; INTPVB itself differs between images (the defect)
  !   5 gauge and solver on a ring of elements (mean mode singular, slopes regular), solver
  !     types 1-4: closed component detected (C e = 0), compatible right-hand side solved
  !     with an ungauged residual at rounding level and a reference value of zero, the
  !     solution equal to the exact one up to the kernel; an incompatible right-hand side is
  !     reported; the same ring with an outflow (open) gets no reference and is solved exactly
  !   6 constraint masks: all 64 pairs of three-component Dirichlet masks
  !   7 gauge components from effective couplings: two closed rings that share only all-zero
  !     structural blocks are two components with a reference each and an imbalance of
  !     opposite sign in the two is reported; an open and a closed ring sharing a zero block:
  !     one reference (in the closed one), factorisation and solve succeed
  !   8 BiCGStab (types 1, 3): a zero right-hand side returns zero; exact convergence at the
  !     half step (diagonal operator, exact preconditioner) returns the solution, no 0/0
  !   9 vertex averaging on a 1 x 1 x n strip (more vertices, 4(n+1), than pressure unknowns,
  !     4 n): own nvt-sized workspace, exact values
  !  10 runtime change of the coarse solver type 2 -> 4 -> 1 on the installed layout
  !  11 census consistency classification of a coupled pair: adjacent cells of width 0.6
  !     period share a vertex unshifted, a pair across the plane needs the shift, a pair
  !     without a shared vertex is inconsistent
  !  12 operator comparison (wrapped validation): NaN or Inf in either operand rejected, a
  !     small nonzero entry in a row that is zero in the oracle rejected, rounding passes
  !  13 weak bridge: two closed rings joined by one conservative edge of weight 1e-12 are ONE
  !     component with one reference (each indicator checked on its own); a compatible
  !     right-hand side that needs a pressure difference of 1 across the bridge is solved
  !     with that difference (types 2, 4); also the two separate rings keep two references
  !  14 closed singleton whose mean row vanishes: the reference equation gets a nonzero
  !     diagonal on the copy, full (types 1, 2) and reduced (types 3, 4) solves succeed; the
  !     raw SOR leaves the zero row unchanged; a zero diagonal outside a reference is
  !     refused for the SSOR types
  !  15 link provenance of the census: two cells along a periodic axis, one per rank, share
  !     the ordinary face x = 1/2 and the periodic face x = 0 ~ 1. The E013 acceptance test
  !     (coincident or one period apart, per axis) links the dofs of both faces, and
  !     COUNT_PERIODIC_LINK counts the periodic ones on x (the pair classification alone
  !     says unshifted); cells of width 0.6 period get no periodic link
  !  16 closure scale: a closed ring whose slope rows carry rounding-noise couplings to the
  !     mean modes (1e-17, as on the symmetric 2 x 2 x 2 periodic box) stays closed (the
  !     measure is normwise, not relative to the mean columns of the row alone)
  !  17 gauge from the structure: a weak outlet (open diagonal 2 + 1e-9, matrix closure
  !     2.5e-10, below the matrix tolerance) with a structural flux of order one stays open,
  !     no reference, and the exact solution (one on every mean mode) is returned (type 2);
  !     the P0 types 3/4 are refused for it (their over-relaxed P0 correction diverges at weak
  !     outlets, measured), a well open ring is accepted;
  !     next to it the closed ring with noise stays closed (16). Inconsistent or undecidable
  !     cases set no reference and return PCO_GAUGE_UNDECIDED: structurally closed but the
  !     operator does not annihilate the constant; structurally open but the constant is in
  !     the numerical kernel; a structural flux between the closed and the open thresholds
  !  18 worker-product comparison: a zero row of the gathered operator with a nonzero worker
  !     result fails (it was skipped)
  !  19 coordinate tolerance per coordinate (as the E013 links): an image offset by 0.8
  !     DEpsPrec in two coordinates is coincident for the pair check and the master classes;
  !     1.2 DEpsPrec in one coordinate is not
  !  20 the numerical limit of the gauge decision: a ring with an outlet of relative strength
  !     d (structural flux d, matrix measure about d) for d = 2^-10, 2^-20, 2^-30, 2^-40,
  !     2^-50 and 0: open / undecided / undecided / undecided / closed / closed; for the closed
  !     ones the direct solve of the UNGAUGED system fails (d = 0) or loses the mean under a
  !     perturbation of the right-hand side at the backward-error level eps |A| |x| of a stable
  !     solve (d = 2^-50); with the exact right-hand side of the manufactured case the mean is
  !     still recovered (printed: the data, not the operator, carry the information there)
  !  21 the P0 types follow the CURRENT matrix measure: a ring renewed on the same layout from
  !     2.5 to 2 + 1e-3 at the outlet is no longer admissible for type 4 (admission check and
  !     refresh)
  !  22 the worker-product comparison fails on a non-finite product or scale

  USE MPI
  USE TYPES, ONLY: tMultiMesh, TMatrix
  USE var_QuadScalar, ONLY: mg_mesh
  USE UMFPackSolver, ONLY: tUmfCSC, UmfCSC_SetStructure, UmfCSC_SetValues, UmfCSC_Analyse, &
                           UmfCSC_Factorize, UmfCSC_Solve, UmfCSC_Release, &
                           myUmfPack_Factorize, myUmfPack_Solve, myUmfPack_Free
  USE pco_periodic_topology, ONLY: tPcoClasses, PCO_BUILD_CLASSES, PCO_CLASS_MASK_CONFLICTS, &
                                   PCO_ELEMENT_GRAPH, PCO_MASKS_DIFFER, PCO_PAIR_SHIFT
  USE pco_p0_reconstruct, ONLY: PCO_SET_VERTEX_CLASSES, PCO_P0_TO_VERTEX, pco_vertex_values
  USE pco_coarse_solver, ONLY: pco, PCO_INSTALL_LAYOUT, PCO_REFRESH_SOLVER, PCO_SOLVE_CHECK, &
                               PCO_BLOCK_SLOT, PCO_NEIGHBOUR_POS, PCO_SOLVE_FULL, PCO_SOLVE_P0, &
                               PCO_RELEASE, PCO_CHANGE_SOLVER_TYPE, PCO_RAW_SOR, PCO_RAW_MATVEC, &
                               PCO_SET_STRUCTURAL_FLUX, PCO_GAUGE_UNDECIDED, PCO_P0_WEAK_OPEN, &
                               PCO_P0_ADMISSIBLE, PCO_FLUX_CLOSED, PCO_CLOSURE_TOL
  USE pco_diagnostics, ONLY: PCO_COMPARE_OPERATORS, PCO_COMPARE_PRODUCTS
  USE PP3D_MPI, ONLY: COUNT_PERIODIC_LINK, DEpsPrec
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_VALUE, IEEE_QUIET_NAN, IEEE_POSITIVE_INF
  USE mg_QuadScalar, ONLY: INTPVB
  USE EL_TEST_HELPERS, ONLY: EL_TEST_INIT, EL_TEST_CHECK, EL_TEST_REPORT, BUILD_STRUCTURED_MESH, &
                             BUILD_DOF_COORDINATES

  IMPLICIT NONE

  INTEGER :: rank, nproc, ierr, nf
  TYPE(tMultiMesh) :: tmp

  CALL MPI_Init(ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  CALL EL_TEST_INIT()
  nf = 0

  CALL test_umfpack_orientation()
  CALL test_classes()
  CALL test_reconstruction()
  CALL test_gauge_solver()
  CALL test_masks()
  CALL test_components()
  CALL test_bicgstab_edges()
  CALL test_strip_reconstruction()
  CALL test_solver_type_change()
  CALL test_pair_shift()
  CALL test_compare_operators()
  CALL test_weak_bridge()
  CALL test_zero_reference_row()
  CALL test_link_provenance()
  CALL test_closure_scale()
  CALL test_structural_gauge()
  CALL test_compare_products()
  CALL test_componentwise_tolerance()
  CALL test_numerical_limit()
  CALL test_current_measure()

  CALL EL_TEST_REPORT(nf, rank, nproc, 'pco coarse')
  CALL MPI_Finalize(ierr)
  IF (nf.NE.0) STOP 1

CONTAINS

  !---------------------------------------------------------------------
  SUBROUTINE test_umfpack_orientation()
    ! A = [4 1 0 0; 2 5 1 0; 0 3 6 1; 1 0 2 7] (non-symmetric), CSR
    INTEGER, PARAMETER :: n = 4
    INTEGER :: ld(n+1) = (/1, 3, 6, 9, 12/)
    INTEGER :: col(11) = (/1, 2, 1, 2, 3, 2, 3, 4, 1, 3, 4/)
    REAL*8 :: a(11) = (/4d0, 1d0, 2d0, 5d0, 1d0, 3d0, 6d0, 1d0, 1d0, 2d0, 7d0/)
    REAL*8 :: xt(n) = (/1d0, -2d0, 3d0, 0.5d0/), b(n), x(n), bt(n), xl(n), ad(n,n)
    TYPE(tUmfCSC) :: F
    TYPE(TMatrix) :: M
    REAL*8 :: acopy(11)
    INTEGER :: st, i, k

    WRITE(*,'(A)') '1 UMFPACK orientation'
    ad = 0d0
    DO i = 1, n
      DO k = ld(i), ld(i+1)-1
        ad(i, col(k)) = a(k)
      END DO
    END DO
    b = MATMUL(ad, xt)
    CALL UmfCSC_SetStructure(F, n, ld, col)
    CALL UmfCSC_SetValues(F, a)
    CALL UmfCSC_Analyse(F, st)
    CALL EL_TEST_CHECK(st.EQ.0, 'adapter: symbolic analysis status 0', nf, rank)
    CALL UmfCSC_Factorize(F, st)
    CALL EL_TEST_CHECK(st.EQ.0, 'adapter: numeric factorisation status 0', nf, rank)
    x = 0d0
    CALL UmfCSC_Solve(F, x, b, st)
    CALL EL_TEST_CHECK(st.EQ.0 .AND. MAXVAL(ABS(x-xt)).LT.1d-13, 'adapter solves A x = b', nf, rank)
    CALL UmfCSC_Release(F)
    CALL UmfCSC_Solve(F, x, b, st)
    CALL EL_TEST_CHECK(st.LT.0, 'adapter: solve without factors returns an error status', nf, rank)

    ! legacy wrapper with the same CSR arrays: solves the transpose
    M%nu = n
    M%na = 11
    ALLOCATE(M%LdA(n+1), M%ColA(11))
    M%LdA = ld
    M%ColA = col
    acopy = a
    CALL myUmfPack_Factorize(acopy, M)
    xl = 0d0
    CALL myUmfPack_Solve(xl, b, acopy, M, 1)
    bt = MATMUL(TRANSPOSE(ad), xl)
    CALL EL_TEST_CHECK(MAXVAL(ABS(bt-b)).LT.1d-12 .AND. MAXVAL(ABS(xl-xt)).GT.1d-3, &
      'legacy wrapper solves A^T x = b (row storage handed over as column storage)', nf, rank)
    CALL myUmfPack_Free()
    CALL myUmfPack_Free()     ! second release is a no-op
    CALL EL_TEST_CHECK(.TRUE., 'legacy release twice is safe', nf, rank)

    ! singular matrix: the adapter reports it
    a = 0d0
    a(1) = 1d0
    CALL UmfCSC_SetStructure(F, n, ld, col)
    CALL UmfCSC_SetValues(F, a)
    CALL UmfCSC_Analyse(F, st)
    CALL UmfCSC_Factorize(F, st)
    CALL EL_TEST_CHECK(st.EQ.1, 'adapter reports a singular matrix (status 1)', nf, rank)
    CALL UmfCSC_Release(F)
  END SUBROUTINE test_umfpack_orientation

  !---------------------------------------------------------------------
  ! level 1 = unit cube with n cells, level 2 dcorvg = Q2 dof coordinates of level 1
  SUBROUTINE make_mesh(n, perturb, dup)
    INTEGER, INTENT(IN) :: n
    LOGICAL, INTENT(IN) :: perturb, dup
    INTEGER :: ndof
    CALL BUILD_STRUCTURED_MESH(tmp, n)
    IF (ALLOCATED(mg_mesh%level)) DEALLOCATE(mg_mesh%level)
    ALLOCATE(mg_mesh%level(2))
    mg_mesh%level(1)%nel = tmp%level(1)%nel
    mg_mesh%level(1)%nvt = tmp%level(1)%nvt
    mg_mesh%level(1)%net = tmp%level(1)%net
    mg_mesh%level(1)%nat = tmp%level(1)%nat
    ALLOCATE(mg_mesh%level(1)%dcorvg(3, tmp%level(1)%nvt))
    mg_mesh%level(1)%dcorvg = tmp%level(1)%dcorvg
    ALLOCATE(mg_mesh%level(1)%kvert(8, tmp%level(1)%nel))
    mg_mesh%level(1)%kvert = tmp%level(1)%kvert
    ALLOCATE(mg_mesh%level(1)%dvol(tmp%level(1)%nel))
    mg_mesh%level(1)%dvol = tmp%level(1)%dvol
    ndof = tmp%level(1)%nvt + tmp%level(1)%net + tmp%level(1)%nat + tmp%level(1)%nel
    ALLOCATE(mg_mesh%level(2)%dcorvg(3, ndof))
    CALL BUILD_DOF_COORDINATES(tmp, mg_mesh%level(2)%dcorvg)
    ! vertex (1,1,0) of the cube is the last vertex of the first z layer
    IF (perturb) mg_mesh%level(2)%dcorvg(1, (n+1)*(n+1)) = 1d0 - 1d-3
    ! duplicate: an edge midpoint placed onto vertex 1
    IF (dup) mg_mesh%level(2)%dcorvg(:, tmp%level(1)%nvt + 1) = mg_mesh%level(2)%dcorvg(:, 2)
  END SUBROUTINE make_mesh

  SUBROUTINE count_sizes(cls, nvt, hist)
    TYPE(tPcoClasses), INTENT(IN) :: cls
    INTEGER, INTENT(IN) :: nvt
    INTEGER, INTENT(OUT) :: hist(0:8)
    INTEGER :: k, s
    hist = 0
    DO k = 1, cls%ncls
      s = cls%ptr(k+1) - cls%ptr(k)
      IF (cls%mem(cls%ptr(k)).GT.nvt) CYCLE
      hist(MIN(s,8)) = hist(MIN(s,8)) + 1
    END DO
  END SUBROUTINE count_sizes

  SUBROUTINE test_classes()
    TYPE(tPcoClasses) :: cls
    INTEGER, PARAMETER :: n = 3
    INTEGER :: hist(0:8), nvt, ndof, nconf, e, ncnt(0:30)
    INTEGER, ALLOCATABLE :: ku(:), kv(:), kw(:), eptr(:), enbr(:)
    REAL*8 :: per(3) = (/1d0, 1d0, 1d0/)

    WRITE(*,'(A)') '2 periodic classes'
    CALL make_mesh(n, .FALSE., .FALSE.)
    nvt = (n+1)**3
    ndof = SIZE(mg_mesh%level(2)%dcorvg, 2)

    CALL PCO_BUILD_CLASSES(1, (/.TRUE., .FALSE., .FALSE./), per, cls)
    CALL count_sizes(cls, nvt, hist)
    WRITE(*,'(A,9I5)') '    vertex class sizes 0..8 (x):', hist
    CALL EL_TEST_CHECK(hist(2).EQ.(n+1)**2 .AND. hist(1).EQ.(n-1)*(n+1)**2 .AND. &
      cls%nambiguous+cls%nincomplete+cls%nmixed.EQ.0, 'x periodic: vertex pairs on the x planes', nf, rank)
    CALL EL_TEST_CHECK(cls%nmulti.EQ.(n+1)**2 + 2*n*(n+1) + n*n, &
      'x periodic: vertex, edge and face classes of the x planes', nf, rank)

    CALL PCO_BUILD_CLASSES(1, (/.TRUE., .TRUE., .FALSE./), per, cls)
    CALL count_sizes(cls, nvt, hist)
    WRITE(*,'(A,9I5)') '    vertex class sizes 0..8 (xy):', hist
    CALL EL_TEST_CHECK(hist(4).EQ.n+1 .AND. hist(2).EQ.2*(n-1)*(n+1) .AND. &
      cls%nambiguous+cls%nincomplete+cls%nmixed.EQ.0, 'xy periodic: edge vertices with 4 images', nf, rank)

    CALL PCO_BUILD_CLASSES(1, (/.TRUE., .TRUE., .TRUE./), per, cls)
    CALL count_sizes(cls, nvt, hist)
    WRITE(*,'(A,9I5)') '    vertex class sizes 0..8 (xyz):', hist
    CALL EL_TEST_CHECK(hist(8).EQ.1 .AND. hist(4).EQ.3*(n-1) .AND. hist(2).EQ.3*(n-1)**2 .AND. &
      cls%nambiguous+cls%nincomplete+cls%nmixed.EQ.0, 'xyz periodic: one corner class with 8 images', nf, rank)
    CALL EL_TEST_CHECK(cls%maxsize.EQ.8, 'xyz periodic: largest class 8', nf, rank)

    ! constraint flags: equal over the classes, then one image flagged alone
    ALLOCATE(ku(ndof), kv(ndof), kw(ndof))
    ku = 0
    kv = 0
    kw = 0
    nconf = PCO_CLASS_MASK_CONFLICTS(cls, ku, kv, kw)
    CALL EL_TEST_CHECK(nconf.EQ.0, 'no constraint conflict with equal flags', nf, rank)
    kw(1) = 1
    nconf = PCO_CLASS_MASK_CONFLICTS(cls, ku, kv, kw)
    CALL EL_TEST_CHECK(nconf.EQ.1, 'one image with a Dirichlet flag of its own -> conflict', nf, rank)
    DEALLOCATE(ku, kv, kw)

    ! element graph: xy periodic cube, n = 3 -> 17 neighbours in the end layers, 26 inside
    CALL PCO_BUILD_CLASSES(1, (/.TRUE., .TRUE., .FALSE./), per, cls)
    CALL PCO_ELEMENT_GRAPH(1, cls, eptr, enbr)
    ncnt = 0
    DO e = 1, mg_mesh%level(1)%nel
      ncnt(MIN(eptr(e+1)-eptr(e),30)) = ncnt(MIN(eptr(e+1)-eptr(e),30)) + 1
    END DO
    CALL EL_TEST_CHECK(ncnt(17).EQ.2*n*n .AND. ncnt(26).EQ.n*n*(n-2), &
      'element graph through xy classes: 17 / 26 neighbours', nf, rank)
    DEALLOCATE(eptr, enbr)

    ! a corner image moved off its place: the class of that corner is incomplete
    CALL make_mesh(n, .TRUE., .FALSE.)
    CALL PCO_BUILD_CLASSES(1, (/.TRUE., .TRUE., .FALSE./), per, cls)
    CALL EL_TEST_CHECK(cls%nincomplete.GE.1, 'moved corner image -> incomplete class refused', nf, rank)

    ! two dofs at one place: ambiguous partner
    CALL make_mesh(n, .FALSE., .TRUE.)
    CALL PCO_BUILD_CLASSES(1, (/.TRUE., .TRUE., .FALSE./), per, cls)
    CALL EL_TEST_CHECK(cls%nambiguous.GE.1, 'duplicated dof -> ambiguous partner refused', nf, rank)
  END SUBROUTINE test_classes

  !---------------------------------------------------------------------
  SUBROUTINE test_reconstruction()
    TYPE(tPcoClasses) :: cls
    INTEGER, PARAMETER :: n = 3
    INTEGER :: nel, nvt, e, j, k, m, v
    REAL*8, ALLOCATABLE :: dp(:), dpl(:), dl0(:), da0(:), num(:), den(:)
    REAL*8 :: c(3), pi, worst_img, worst_val, worst_int, legacy_gap, nn, dd
    REAL*8 :: per(3) = (/1d0, 1d0, 1d0/)
    LOGICAL :: on_image

    WRITE(*,'(A)') '4 P0 -> vertex averaging at periodic planes'
    CALL make_mesh(n, .FALSE., .FALSE.)
    nel = mg_mesh%level(1)%nel
    nvt = mg_mesh%level(1)%nvt
    pi = 4d0*ATAN(1d0)
    ALLOCATE(dp(nel), dpl(nvt), dl0(nvt), da0(nvt), num(nvt), den(nvt))
    DO e = 1, nel
      c = 0d0
      DO j = 1, 8
        c = c + mg_mesh%level(1)%dcorvg(:, mg_mesh%level(1)%kvert(j,e))
      END DO
      c = 0.125d0*c
      dp(e) = SIN(2d0*pi*c(1))*COS(2d0*pi*c(2)) + c(3) + 0.1d0*DBLE(MOD(e,5))
      mg_mesh%level(1)%dvol(e) = mg_mesh%level(1)%dvol(e)*(1d0 + 0.05d0*DBLE(MOD(3*e,7)))
    END DO
    CALL PCO_BUILD_CLASSES(1, (/.TRUE., .TRUE., .FALSE./), per, cls)
    CALL PCO_SET_VERTEX_CLASSES(cls%nvt, cls%ncls, cls%ptr, cls%mem)
    CALL PCO_P0_TO_VERTEX(dp, mg_mesh%level(1)%dvol, mg_mesh%level(1)%kvert, nel, nvt)
    dpl = pco_vertex_values(1:nvt)
    CALL INTPVB(dp, dl0, da0, mg_mesh%level(1)%dvol, mg_mesh%level(1)%kvert, nel, nvt)
    ! independent numerators / denominators
    num = 0d0
    den = 0d0
    DO e = 1, nel
      DO j = 1, 8
        v = mg_mesh%level(1)%kvert(j,e)
        num(v) = num(v) + 0.125d0*mg_mesh%level(1)%dvol(e)*dp(e)
        den(v) = den(v) + 0.125d0*mg_mesh%level(1)%dvol(e)
      END DO
    END DO
    worst_img = 0d0
    worst_val = 0d0
    worst_int = 0d0
    legacy_gap = 0d0
    DO k = 1, cls%ncls
      IF (cls%mem(cls%ptr(k)).GT.nvt) CYCLE
      on_image = (cls%ptr(k+1)-cls%ptr(k).GT.1)
      nn = 0d0
      dd = 0d0
      DO m = cls%ptr(k), cls%ptr(k+1)-1
        nn = nn + num(cls%mem(m))
        dd = dd + den(cls%mem(m))
      END DO
      DO m = cls%ptr(k), cls%ptr(k+1)-1
        v = cls%mem(m)
        IF (on_image) THEN
          worst_img = MAX(worst_img, ABS(dpl(v) - dpl(cls%mem(cls%ptr(k)))))
          worst_val = MAX(worst_val, ABS(dpl(v) - nn/dd))
          legacy_gap = MAX(legacy_gap, ABS(dl0(v) - dl0(cls%mem(cls%ptr(k)))))
        ELSE
          worst_int = MAX(worst_int, ABS(dpl(v) - dl0(v))/MAX(ABS(dl0(v)), 1d-300))
        END IF
      END DO
    END DO
    CALL EL_TEST_CHECK(worst_img.EQ.0d0, 'every image of a periodic vertex gets the same value', nf, rank)
    CALL EL_TEST_CHECK(worst_val.LT.1d-14, 'the value is sum(num)/sum(den) over the class', nf, rank)
    CALL EL_TEST_CHECK(worst_int.LT.1d-14, 'vertices without images: INTPVB to rounding', nf, rank)
    CALL EL_TEST_CHECK(legacy_gap.GT.1d-3, 'INTPVB gives different values to the images (defect shown)', nf, rank)
    DEALLOCATE(dp, dpl, dl0, da0, num, den)
  END SUBROUTINE test_reconstruction

  !---------------------------------------------------------------------
  ! ring of nel elements, neighbours e-1 and e+1; mean modes: periodic Laplacian (+ outflow
  ! term on element 1 when open), slopes: diagonally dominant
  SUBROUTINE ring_operator(nel, crs, open)
    INTEGER, INTENT(IN) :: nel, crs
    LOGICAL, INTENT(IN) :: open
    INTEGER, ALLOCATABLE :: eptr(:), enbr(:)
    INTEGER :: e, f(2), j, mp, pos
    ALLOCATE(eptr(nel+1), enbr(2*nel))
    DO e = 1, nel
      eptr(e) = 2*(e-1) + 1
      f(1) = MOD(e-2+nel, nel) + 1
      f(2) = MOD(e, nel) + 1
      enbr(eptr(e)) = MINVAL(f)
      enbr(eptr(e)+1) = MAXVAL(f)
    END DO
    eptr(nel+1) = 2*nel + 1
    CALL PCO_INSTALL_LAYOUT(1, nel, eptr, enbr, crs, 2)
    pco%a = 0d0
    DO e = 1, nel
      DO mp = 1, 4
        pco%a(PCO_BLOCK_SLOT(e, mp, 0, mp)) = MERGE(2d0, 3d0 + 0.1d0*DBLE(e), mp.EQ.1)
      END DO
      ! a non-symmetric coupling between two slopes of the same element (keeps C e = 0 and
      ! e^T C = 0, e = 1 on the means)
      pco%a(PCO_BLOCK_SLOT(e, 3, 0, 2)) = 0.3d0
      f(1) = MOD(e-2+nel, nel) + 1
      f(2) = MOD(e, nel) + 1
      DO j = 1, 2
        pos = PCO_NEIGHBOUR_POS(e, f(j))
        DO mp = 1, 4
          pco%a(PCO_BLOCK_SLOT(e, mp, pos, mp)) = pco%a(PCO_BLOCK_SLOT(e, mp, pos, mp)) - 1d0
        END DO
      END DO
    END DO
    IF (open) pco%a(PCO_BLOCK_SLOT(1, 1, 0, 1)) = 2.5d0
    CALL set_flux(nel, MERGE(1, 0, open))
  END SUBROUTINE ring_operator

  ! structural flux of a constant pressure per element (PCO_FREE_FLUX of the workers in the
  ! application): order one on element 1 when iopen = 1, zero elsewhere
  SUBROUTINE set_flux(nel, iopen)
    INTEGER, INTENT(IN) :: nel, iopen
    REAL*8 :: fl(nel)
    fl = 0d0
    IF (iopen.EQ.1) fl(1) = 1d0
    IF (iopen.EQ.2) fl = 1d0
    CALL PCO_SET_STRUCTURAL_FLUX(fl)
  END SUBROUTINE set_flux

  SUBROUTINE test_gauge_solver()
    INTEGER, PARAMETER :: nel = 6
    INTEGER :: crs, st, i, it
    CHARACTER(len=256) :: msg
    CHARACTER(len=64) :: lab
    REAL*8 :: v(4*nel), resid, devi, b(4*nel), x(4*nel), p0(nel), tol, xref

    WRITE(*,'(A)') '5 gauge and solver (ring of elements)'
    DO i = 1, 4*nel
      v(i) = SIN(1.7d0*DBLE(i)) + 0.2d0*DBLE(MOD(i,3))
    END DO
    DO crs = 1, 4
      WRITE(lab,'(A,I0)') 'type ', crs
      CALL ring_operator(nel, crs, .FALSE.)
      CALL PCO_REFRESH_SOLVER(st, msg)
      CALL EL_TEST_CHECK(st.EQ.0 .AND. pco%ncomp.EQ.1 .AND. pco%nclosed.EQ.1 .AND. &
        pco%closure_closed.LT.1d-14, TRIM(lab)//': closed ring detected (C e = 0), factorised', nf, rank)
      i = pco%nsolve
      CALL PCO_SOLVE_CHECK(v, resid, devi, xref)
      WRITE(*,'(A,3ES11.3)') '    residual, deviation from v + kernel, reference:', resid, devi, xref
      tol = MERGE(1d-12, 1d-9, crs.EQ.2 .OR. crs.EQ.4)
      CALL EL_TEST_CHECK(resid.LT.tol .AND. devi.LT.100d0*tol, &
        TRIM(lab)//': compatible rhs, ungauged residual and solution up to the kernel', nf, rank)
      CALL EL_TEST_CHECK(xref.LE.tol, TRIM(lab)//': reference mean value zero', nf, rank)
      CALL EL_TEST_CHECK(pco%nsolve.EQ.i, TRIM(lab)//': the check does not count as a production solve', nf, rank)
      ! incompatible: a constant on the means is not in the range
      b = 0d0
      DO i = 1, nel
        b(4*(i-1)+1) = 1d0
      END DO
      it = 500
      i = pco%nincompat
      IF (crs.LE.2) THEN
        CALL PCO_SOLVE_FULL(x, b, 4*nel, it, 1d-12, 1d0)
      ELSE
        CALL PCO_SOLVE_P0(b, p0, it, 1d-12, 1d0)
      END IF
      CALL EL_TEST_CHECK(pco%nincompat.EQ.i+1, TRIM(lab)//': incompatible rhs reported', nf, rank)
      ! open ring: no reference, exact solve
      CALL ring_operator(nel, crs, .TRUE.)
      CALL PCO_REFRESH_SOLVER(st, msg)
      CALL EL_TEST_CHECK(st.EQ.0 .AND. pco%nclosed.EQ.0 .AND. .NOT.ANY(pco%sref), &
        TRIM(lab)//': open ring, no reference', nf, rank)
      CALL PCO_SOLVE_CHECK(v, resid, devi, xref)
      WRITE(*,'(A,2ES11.3)') '    open: residual, deviation:', resid, devi
      CALL EL_TEST_CHECK(resid.LT.tol .AND. devi.LT.100d0*tol, TRIM(lab)//': open ring solved exactly', nf, rank)
    END DO
    CALL PCO_RELEASE()
  END SUBROUTINE test_gauge_solver

  !---------------------------------------------------------------------
  SUBROUTINE test_masks()
    INTEGER :: i, j, a(3), b(3), nbad
    TYPE(tPcoClasses) :: cls
    INTEGER, ALLOCATABLE :: ku(:), kv(:), kw(:)
    REAL*8 :: per(3) = (/1d0, 1d0, 1d0/)
    WRITE(*,'(A)') '6 constraint masks'
    nbad = 0
    DO i = 0, 7
      DO j = 0, 7
        a = (/MOD(i,2), MOD(i/2,2), MOD(i/4,2)/)
        b = (/MOD(j,2), MOD(j/2,2), MOD(j/4,2)/)
        IF (PCO_MASKS_DIFFER(a(1), a(2), a(3), b(1), b(2), b(3)) .NEQV. ANY(a.NE.b)) nbad = nbad + 1
      END DO
    END DO
    CALL EL_TEST_CHECK(nbad.EQ.0, 'all 64 mask pairs: differ exactly when a component differs', nf, rank)
    CALL make_mesh(3, .FALSE., .FALSE.)
    CALL PCO_BUILD_CLASSES(1, (/.TRUE., .FALSE., .FALSE./), per, cls)
    ALLOCATE(ku(cls%ndof), kv(cls%ndof), kw(cls%ndof))
    ku = 1
    kv = 1
    kw = 0
    CALL EL_TEST_CHECK(PCO_CLASS_MASK_CONFLICTS(cls, ku, kv, kw).EQ.0, &
      'masks (1,1,0) on every image: no conflict', nf, rank)
    ku(1) = 0
    kv(1) = 0
    CALL EL_TEST_CHECK(PCO_CLASS_MASK_CONFLICTS(cls, ku, kv, kw).EQ.1, &
      'masks (0,0,0) against (1,1,0) on one class: one conflict', nf, rank)
    DEALLOCATE(ku, kv, kw)
  END SUBROUTINE test_masks

  !---------------------------------------------------------------------
  ! two rings A = 1..3, B = 4..6 (ring couplings -1 on every mode), every element structurally
  ! neighbour of every other (blocks between A and B stay zero); open_a: outflow on element 1
  SUBROUTINE two_rings(crs, open_a)
    INTEGER, INTENT(IN) :: crs
    LOGICAL, INTENT(IN) :: open_a
    INTEGER, ALLOCATABLE :: eptr(:), enbr(:)
    INTEGER :: e, f, k, mp, g0
    ALLOCATE(eptr(7), enbr(30))
    k = 0
    DO e = 1, 6
      eptr(e) = k + 1
      DO f = 1, 6
        IF (f.EQ.e) CYCLE
        k = k + 1
        enbr(k) = f
      END DO
    END DO
    eptr(7) = k + 1
    CALL PCO_INSTALL_LAYOUT(1, 6, eptr, enbr, crs, 2)
    pco%a = 0d0
    DO e = 1, 6
      g0 = 3*((e-1)/3)
      DO mp = 1, 4
        pco%a(PCO_BLOCK_SLOT(e, mp, 0, mp)) = MERGE(2d0, 4d0, mp.EQ.1)
      END DO
      DO f = g0+1, g0+3
        IF (f.EQ.e) CYCLE
        DO mp = 1, 4
          pco%a(PCO_BLOCK_SLOT(e, mp, PCO_NEIGHBOUR_POS(e, f), mp)) = -1d0
        END DO
      END DO
    END DO
    IF (open_a) pco%a(PCO_BLOCK_SLOT(1, 1, 0, 1)) = 3d0
    CALL set_flux(6, MERGE(1, 0, open_a))
  END SUBROUTINE two_rings

  SUBROUTINE test_components()
    INTEGER :: crs, st, it, n0
    CHARACTER(len=256) :: msg
    CHARACTER(len=32) :: lab
    REAL*8 :: v(24), b(24), x(24), p0(6), resid, devi, xref
    INTEGER :: i
    WRITE(*,'(A)') '7 gauge components from effective couplings'
    DO i = 1, 24
      v(i) = COS(0.7d0*DBLE(i))
    END DO
    DO crs = 2, 4, 2
      WRITE(lab,'(A,I0)') 'type ', crs
      CALL two_rings(crs, .FALSE.)
      CALL PCO_REFRESH_SOLVER(st, msg)
      CALL EL_TEST_CHECK(st.EQ.0 .AND. pco%ncomp.EQ.2 .AND. pco%nclosed.EQ.2 .AND. &
        COUNT(pco%sref).EQ.2, TRIM(lab)//': two closed rings with a zero structural block: two references', nf, rank)
      CALL PCO_SOLVE_CHECK(v, resid, devi, xref)
      CALL EL_TEST_CHECK(resid.LT.1d-12 .AND. devi.LT.1d-10, TRIM(lab)//': compatible rhs solved', nf, rank)
      b = 0d0
      DO i = 1, 3
        b(4*(i-1)+1) = 1d0
        b(4*(i+2)+1) = -1d0
      END DO
      n0 = pco%nincompat
      it = 100
      IF (crs.EQ.2) THEN
        CALL PCO_SOLVE_FULL(x, b, 24, it, 1d-12, 1d0)
      ELSE
        CALL PCO_SOLVE_P0(b, p0, it, 1d-12, 1d0)
      END IF
      CALL EL_TEST_CHECK(pco%nincompat.EQ.n0+1, TRIM(lab)//': opposite imbalances of the two rings reported', nf, rank)
      CALL two_rings(crs, .TRUE.)
      CALL PCO_REFRESH_SOLVER(st, msg)
      CALL EL_TEST_CHECK(st.EQ.0 .AND. pco%ncomp.EQ.2 .AND. pco%nclosed.EQ.1 .AND. COUNT(pco%sref).EQ.1 &
        .AND. pco%comp_closed(pco%comp_of(4)), TRIM(lab)//': open and closed ring: one reference, in the closed one', &
        nf, rank)
      CALL PCO_SOLVE_CHECK(v, resid, devi, xref)
      CALL EL_TEST_CHECK(resid.LT.1d-12 .AND. devi.LT.1d-10, TRIM(lab)//': open/closed pair solved', nf, rank)
    END DO
    CALL PCO_RELEASE()
  END SUBROUTINE test_components

  !---------------------------------------------------------------------
  ! isolated elements with diagonal blocks (open, no gauge): the SSOR preconditioner is exact
  SUBROUTINE diagonal_operator(crs)
    INTEGER, INTENT(IN) :: crs
    INTEGER, ALLOCATABLE :: eptr(:), enbr(:)
    INTEGER :: e, mp
    ALLOCATE(eptr(5), enbr(1))
    eptr = 1
    CALL PCO_INSTALL_LAYOUT(1, 4, eptr, enbr, crs, 2)
    pco%a = 0d0
    DO e = 1, 4
      DO mp = 1, 4
        pco%a(PCO_BLOCK_SLOT(e, mp, 0, mp)) = DBLE(e + mp)
      END DO
    END DO
    CALL set_flux(4, 2)
  END SUBROUTINE diagonal_operator

  SUBROUTINE test_bicgstab_edges()
    INTEGER :: crs, st, it, i
    CHARACTER(len=256) :: msg
    CHARACTER(len=32) :: lab
    REAL*8 :: b(16), x(16), p0(4), ex(16)
    WRITE(*,'(A)') '8 BiCGStab edge cases'
    DO crs = 1, 3, 2
      WRITE(lab,'(A,I0)') 'type ', crs
      CALL diagonal_operator(crs)
      CALL PCO_REFRESH_SOLVER(st, msg)
      b = 0d0
      it = 50
      IF (crs.EQ.1) THEN
        x = 1d0
        CALL PCO_SOLVE_FULL(x, b, 16, it, 1d-12, 1d0)
        CALL EL_TEST_CHECK(ALL(x.EQ.0d0) .AND. it.EQ.0, TRIM(lab)//': zero rhs -> zero, no iteration', nf, rank)
      ELSE
        p0 = 1d0
        CALL PCO_SOLVE_P0(b, p0, it, 1d-12, 1d0)
        CALL EL_TEST_CHECK(ALL(p0.EQ.0d0) .AND. it.EQ.0, TRIM(lab)//': zero rhs -> zero, no iteration', nf, rank)
      END IF
      DO i = 1, 16
        b(i) = DBLE(i)
        ex(i) = b(i)/DBLE((i-1)/4 + 1 + MOD(i-1,4) + 1)
      END DO
      it = 50
      IF (crs.EQ.1) THEN
        CALL PCO_SOLVE_FULL(x, b, 16, it, 1d-12, 1d0)
        CALL EL_TEST_CHECK(MAXVAL(ABS(x-ex)).LT.1d-14 .AND. it.EQ.1, &
          TRIM(lab)//': exact convergence at the half step, no 0/0', nf, rank)
      ELSE
        CALL PCO_SOLVE_P0(b, p0, it, 1d-12, 1d0)
        CALL EL_TEST_CHECK(MAXVAL(ABS(p0-ex(1:16:4))).LT.1d-14 .AND. it.EQ.1, &
          TRIM(lab)//': exact convergence at the half step, no 0/0', nf, rank)
      END IF
    END DO
    CALL PCO_RELEASE()
  END SUBROUTINE test_bicgstab_edges

  !---------------------------------------------------------------------
  SUBROUTINE test_strip_reconstruction()
    INTEGER, PARAMETER :: n = 5
    INTEGER :: kvert(8,n), e, k, nvt
    REAL*8 :: dvol(n), dp(n)
    INTEGER :: ptr(2), mem(1)
    WRITE(*,'(A)') '9 vertex averaging on a 1 x 1 x n strip'
    nvt = 4*(n+1)
    DO e = 1, n
      DO k = 1, 4
        kvert(k, e) = 4*(e-1) + k
        kvert(k+4, e) = 4*e + k
      END DO
      dvol(e) = 1d0
      dp(e) = 2.5d0
    END DO
    ptr = (/1, 2/)
    mem = 1
    CALL PCO_SET_VERTEX_CLASSES(nvt, 1, ptr, mem)
    CALL PCO_P0_TO_VERTEX(dp, dvol, kvert, n, nvt)
    CALL EL_TEST_CHECK(SIZE(pco_vertex_values).GE.nvt .AND. nvt.GT.4*n .AND. &
      MAXVAL(ABS(pco_vertex_values(1:nvt)-2.5d0)).LT.1d-15, &
      'strip: nvt = 4(n+1) > 4 n values in the own workspace, constant field reproduced', nf, rank)
  END SUBROUTINE test_strip_reconstruction

  !---------------------------------------------------------------------
  SUBROUTINE test_solver_type_change()
    INTEGER :: st
    CHARACTER(len=256) :: msg
    REAL*8 :: v(24), resid, devi, xref
    INTEGER :: i
    WRITE(*,'(A)') '10 runtime change of the coarse solver type'
    DO i = 1, 24
      v(i) = SIN(0.3d0*DBLE(i))
    END DO
    CALL ring_operator(6, 2, .FALSE.)
    CALL PCO_REFRESH_SOLVER(st, msg)
    CALL PCO_CHANGE_SOLVER_TYPE(4, st, msg)
    CALL EL_TEST_CHECK(st.EQ.0 .AND. pco%ns.EQ.6 .AND. pco%solver_type.EQ.4 .AND. COUNT(pco%sref).EQ.1, &
      '2 -> 4: mean-mode copy installed with the gauge', nf, rank)
    CALL PCO_SOLVE_CHECK(v, resid, devi, xref)
    CALL EL_TEST_CHECK(resid.LT.1d-12, '2 -> 4: solve', nf, rank)
    CALL PCO_CHANGE_SOLVER_TYPE(1, st, msg)
    CALL EL_TEST_CHECK(st.EQ.0 .AND. pco%ns.EQ.24 .AND. pco%solver_type.EQ.1, '4 -> 1: full copy installed', nf, rank)
    CALL PCO_SOLVE_CHECK(v, resid, devi, xref)
    CALL EL_TEST_CHECK(resid.LT.1d-9, '4 -> 1: solve', nf, rank)
    CALL PCO_RELEASE()
  END SUBROUTINE test_solver_type_change

  !---------------------------------------------------------------------
  SUBROUTINE test_pair_shift()
    REAL*8 :: va(24), vb(24), vc(24), vd(24), per(3)
    INTEGER :: shift(3), kind
    WRITE(*,'(A)') '11 census classification of coupled pairs'
    per = (/1d0, 1d9, 1d9/)
    CALL cell(0d0, 0.6d0, va)
    CALL cell(0.6d0, 1.2d0, vb)
    CALL PCO_PAIR_SHIFT(va, vb, per, 1d-5, shift, kind)
    CALL EL_TEST_CHECK(kind.EQ.0, 'adjacent cells of width 0.6 period: unshifted (no periodic link)', nf, rank)
    CALL cell(0.4d0, 1.0d0, vc)
    CALL PCO_PAIR_SHIFT(va, vc, per, 1d-5, shift, kind)
    CALL EL_TEST_CHECK(kind.EQ.1 .AND. ALL(shift.EQ.(/1, 0, 0/)), 'pair across the plane x = 0: shift (1,0,0)', &
      nf, rank)
    CALL cell(2.0d0, 2.5d0, vd)
    CALL PCO_PAIR_SHIFT(va, vd, per, 1d-5, shift, kind)
    CALL EL_TEST_CHECK(kind.EQ.-1, 'pair without a shared vertex in any image: inconsistent', nf, rank)
    per = (/1d9, 1d9, 1d9/)
    CALL PCO_PAIR_SHIFT(va, vc, per, 1d-5, shift, kind)
    CALL EL_TEST_CHECK(kind.EQ.-1, 'no finite period: no shift is tried', nf, rank)
  END SUBROUTINE test_pair_shift

  !---------------------------------------------------------------------
  SUBROUTINE test_compare_operators()
    ! two rows: row 1 a normal row, row 2 zero in the oracle
    INTEGER :: ld(3) = (/1, 3, 5/), kmax
    REAL*8 :: ref(4), a(4), nan, inf, d
    WRITE(*,'(A)') '12 operator comparison'
    nan = IEEE_VALUE(nan, IEEE_QUIET_NAN)
    inf = IEEE_VALUE(inf, IEEE_POSITIVE_INF)
    ref = (/2d0, -1d0, 0d0, 0d0/)
    a = ref
    a(1) = a(1)*(1d0 + 1d-15)
    d = PCO_COMPARE_OPERATORS(a, ref, kmax, ld)
    CALL EL_TEST_CHECK(d.LT.1d-14, 'rounding difference passes', nf, rank)
    a = ref
    a(2) = nan
    d = PCO_COMPARE_OPERATORS(a, ref, kmax, ld)
    CALL EL_TEST_CHECK(d.GE.HUGE(1d0) .AND. kmax.EQ.2, 'NaN in the tested operator rejected', nf, rank)
    a = ref
    ref(2) = nan
    d = PCO_COMPARE_OPERATORS(a, ref, kmax, ld)
    CALL EL_TEST_CHECK(d.GE.HUGE(1d0) .AND. kmax.EQ.2, 'NaN in the oracle (gathered) rejected', nf, rank)
    ref = (/2d0, -1d0, 0d0, 0d0/)
    a = ref
    ref(1) = inf
    a(1) = inf
    d = PCO_COMPARE_OPERATORS(a, ref, kmax, ld)
    CALL EL_TEST_CHECK(d.GE.HUGE(1d0) .AND. kmax.EQ.1, 'Inf in both operands (Inf/Inf) rejected', nf, rank)
    ref = (/2d0, -1d0, 0d0, 0d0/)
    a = ref
    a(1) = -inf
    d = PCO_COMPARE_OPERATORS(a, ref, kmax, ld)
    CALL EL_TEST_CHECK(d.GE.HUGE(1d0) .AND. kmax.EQ.1, 'Inf in the tested operator rejected', nf, rank)
    a = ref
    a(4) = 1d-20
    d = PCO_COMPARE_OPERATORS(a, ref, kmax, ld)
    CALL EL_TEST_CHECK(d.GE.1d0 .AND. kmax.EQ.4, 'small nonzero entry in a zero oracle row rejected', nf, rank)
  END SUBROUTINE test_compare_operators

  !---------------------------------------------------------------------
  ! two_rings plus one conservative mean-mode edge of weight w between elements 1 and 4
  SUBROUTINE bridged_rings(crs, w)
    INTEGER, INTENT(IN) :: crs
    REAL*8, INTENT(IN) :: w
    CALL two_rings(crs, .FALSE.)
    pco%a(PCO_BLOCK_SLOT(1, 1, 0, 1)) = pco%a(PCO_BLOCK_SLOT(1, 1, 0, 1)) + w
    pco%a(PCO_BLOCK_SLOT(4, 1, 0, 1)) = pco%a(PCO_BLOCK_SLOT(4, 1, 0, 1)) + w
    pco%a(PCO_BLOCK_SLOT(1, 1, PCO_NEIGHBOUR_POS(1, 4), 1)) = -w
    pco%a(PCO_BLOCK_SLOT(4, 1, PCO_NEIGHBOUR_POS(4, 1), 1)) = -w
  END SUBROUTINE bridged_rings

  SUBROUTINE test_weak_bridge()
    INTEGER :: crs, st, it, i
    CHARACTER(len=256) :: msg
    CHARACTER(len=32) :: lab
    REAL*8 :: xt(24), b(24), x(24), p0(6), r(24), jump, berr
    WRITE(*,'(A)') '13 weak bridge between two closed rings'
    xt = 0d0
    DO i = 4, 6
      xt(4*(i-1)+1) = 1d0                       ! ring B one unit above ring A
    END DO
    DO crs = 2, 4, 2
      WRITE(lab,'(A,I0)') 'type ', crs
      CALL bridged_rings(crs, 1d-12)
      CALL PCO_REFRESH_SOLVER(st, msg)
      CALL EL_TEST_CHECK(st.EQ.0 .AND. pco%ncomp.EQ.1 .AND. pco%nclosed.EQ.1 .AND. COUNT(pco%sref).EQ.1 &
        .AND. pco%closure_closed.LT.1d-14, TRIM(lab)//': bridge 1e-12: one component, one reference', nf, rank)
      CALL PCO_RAW_MATVEC(xt, b)
      it = 10
      IF (crs.EQ.2) THEN
        CALL PCO_SOLVE_FULL(x, b, 24, it, 1d-12, 1d0)
      ELSE
        CALL PCO_SOLVE_P0(b, p0, it, 1d-12, 1d0)
        x = 0d0
        DO i = 1, 6
          x(4*(i-1)+1) = p0(i)
        END DO
        b(2:24:4) = 0d0
        b(3:24:4) = 0d0
        b(4:24:4) = 0d0
      END IF
      CALL PCO_RAW_MATVEC(x, r)
      berr = MAXVAL(ABS(r - b))/(4d0*MAXVAL(ABS(x)))
      jump = (SUM(x(13:24:4)) - SUM(x(1:12:4)))/3d0
      WRITE(*,'(A,2ES11.3)') '    pressure jump across the bridge, backward error:', jump, berr
      CALL EL_TEST_CHECK(ABS(jump - 1d0).LT.1d-2 .AND. berr.LT.1d-12, &
        TRIM(lab)//': compatible rhs solved with the jump 1 across the weak bridge', nf, rank)
      CALL bridged_rings(crs, 0d0)
      CALL PCO_REFRESH_SOLVER(st, msg)
      CALL EL_TEST_CHECK(st.EQ.0 .AND. pco%ncomp.EQ.2 .AND. COUNT(pco%sref).EQ.2, &
        TRIM(lab)//': without the bridge: two components, two references', nf, rank)
    END DO
    CALL PCO_RELEASE()
  END SUBROUTINE test_weak_bridge

  !---------------------------------------------------------------------
  ! one element, no neighbours: block diag(0, d2, d3, d4) (closed: the mean row vanishes)
  SUBROUTINE singleton(crs, d, iopen)
    INTEGER, INTENT(IN) :: crs, iopen
    REAL*8, INTENT(IN) :: d(4)
    INTEGER, ALLOCATABLE :: eptr(:), enbr(:)
    INTEGER :: mp
    ALLOCATE(eptr(2), enbr(1))
    eptr = 1
    enbr = 0
    CALL PCO_INSTALL_LAYOUT(1, 1, eptr, enbr, crs, 2)
    pco%a = 0d0
    DO mp = 1, 4
      pco%a(PCO_BLOCK_SLOT(1, mp, 0, mp)) = d(mp)
    END DO
    CALL set_flux(1, iopen)
  END SUBROUTINE singleton

  SUBROUTINE test_zero_reference_row()
    INTEGER :: crs, st, it
    CHARACTER(len=256) :: msg
    CHARACTER(len=32) :: lab
    REAL*8 :: b(4), x(4), p0(1), xs(4)
    WRITE(*,'(A)') '14 closed singleton with a vanishing mean row'
    b = (/0d0, 3d0, 8d0, 15d0/)
    DO crs = 1, 4
      WRITE(lab,'(A,I0)') 'type ', crs
      CALL singleton(crs, (/0d0, 3d0, 4d0, 5d0/), 0)
      CALL PCO_REFRESH_SOLVER(st, msg)
      CALL EL_TEST_CHECK(st.EQ.0 .AND. pco%nclosed.EQ.1 .AND. pco%sref(1) .AND. &
        pco%sa(pco%sldA(1)).NE.0d0, TRIM(lab)//': reference diagonal nonzero on the copy, refresh ok', nf, rank)
      it = 50
      IF (crs.LE.2) THEN
        x = -1d0
        CALL PCO_SOLVE_FULL(x, b, 4, it, 1d-13, 1d0)
        CALL EL_TEST_CHECK(MAXVAL(ABS(x - (/0d0, 1d0, 2d0, 3d0/))).LT.1d-13, &
          TRIM(lab)//': full solve, reference 0 and the slopes', nf, rank)
      ELSE
        p0 = -1d0
        CALL PCO_SOLVE_P0(b, p0, it, 1d-13, 1d0)
        CALL EL_TEST_CHECK(p0(1).EQ.0d0, TRIM(lab)//': reduced solve (1 x 1 zero system), reference 0', nf, rank)
      END IF
    END DO
    xs = (/7d0, 0d0, 0d0, 0d0/)
    CALL PCO_RAW_SOR(xs, b, 1, 1d0)
    CALL EL_TEST_CHECK(xs(1).EQ.7d0 .AND. MAXVAL(ABS(xs(2:4) - (/1d0, 2d0, 3d0/))).LT.1d-15, &
      'raw SOR: the zero row keeps its value, no division by zero', nf, rank)
    CALL singleton(1, (/1d0, 0d0, 4d0, 5d0/), 1)
    CALL PCO_REFRESH_SOLVER(st, msg)
    CALL EL_TEST_CHECK(st.EQ.4, 'type 1: zero diagonal outside a reference refused at the refresh', nf, rank)
    CALL PCO_RELEASE()
  END SUBROUTINE test_zero_reference_row

  !---------------------------------------------------------------------
  ! the dofs (3 x 3 Q2 positions) of the faces x = xa and x = xb of a unit cross-section cell
  SUBROUTINE face_dofs(xa, xb, p)
    REAL*8, INTENT(IN) :: xa, xb
    REAL*8, INTENT(OUT) :: p(3,18)
    INTEGER :: i, j, k
    k = 0
    DO i = 0, 2
      DO j = 0, 2
        k = k + 1
        p(:,k) = (/xa, 0.5d0*DBLE(i), 0.5d0*DBLE(j)/)
        p(:,k+9) = (/xb, 0.5d0*DBLE(i), 0.5d0*DBLE(j)/)
      END DO
    END DO
  END SUBROUTINE face_dofs

  ! the acceptance test of E013_CreateComm*: every remote dof against the local ones, first hit
  SUBROUTINE e013_links(prem, ploc, per, nlink, nper)
    REAL*8, INTENT(IN) :: prem(:,:), ploc(:,:), per(3)
    INTEGER, INTENT(OUT) :: nlink, nper(3)
    INTEGER :: i, j
    REAL*8, PARAMETER :: eps = 1d-5
    LOGICAL :: hit
    nlink = 0
    nper = 0
    DO i = 1, SIZE(prem,2)
      DO j = 1, SIZE(ploc,2)
        hit = ((ABS(prem(1,i)-ploc(1,j)).LT.eps).OR.(ABS(ABS(prem(1,i)-ploc(1,j))-per(1)).LT.eps)).AND. &
              ((ABS(prem(2,i)-ploc(2,j)).LT.eps).OR.(ABS(ABS(prem(2,i)-ploc(2,j))-per(2)).LT.eps)).AND. &
              ((ABS(prem(3,i)-ploc(3,j)).LT.eps).OR.(ABS(ABS(prem(3,i)-ploc(3,j))-per(3)).LT.eps))
        IF (hit) THEN
          nlink = nlink + 1
          CALL COUNT_PERIODIC_LINK(prem(:,i), ploc(:,j), nper)
          EXIT
        END IF
      END DO
    END DO
  END SUBROUTINE e013_links

  SUBROUTINE test_link_provenance()
    REAL*8 :: p1(3,18), p2(3,18), va(24), vb(24), per(3)
    INTEGER :: nlink, nper(3), shift(3), kind
    WRITE(*,'(A)') '15 census link provenance'
    per = (/1d0, 1d9, 1d9/)
    CALL cell(0d0, 0.5d0, va)
    CALL cell(0.5d0, 1d0, vb)
    CALL PCO_PAIR_SHIFT(va, vb, per, 1d-5, shift, kind)
    CALL EL_TEST_CHECK(kind.EQ.0, 'two cells along a periodic axis: the pair classification says unshifted', &
      nf, rank)
    CALL face_dofs(0d0, 0.5d0, p1)
    CALL face_dofs(0.5d0, 1d0, p2)
    CALL e013_links(p2, p1, per, nlink, nper)
    CALL EL_TEST_CHECK(nlink.EQ.18 .AND. ALL(nper.EQ.(/9, 0, 0/)), &
      'two cells, one per rank: 9 ordinary + 9 periodic links, axis x active from the provenance', nf, rank)
    CALL face_dofs(0d0, 0.6d0, p1)
    CALL face_dofs(0.6d0, 1.2d0, p2)
    CALL e013_links(p2, p1, per, nlink, nper)
    CALL EL_TEST_CHECK(nlink.EQ.9 .AND. ALL(nper.EQ.0), 'cells of width 0.6 period: ordinary links only', nf, rank)
  END SUBROUTINE test_link_provenance

  !---------------------------------------------------------------------
  SUBROUTINE test_closure_scale()
    INTEGER :: st, crs
    CHARACTER(len=256) :: msg
    WRITE(*,'(A)') '16 closure scale with noise in the slope-mean couplings'
    DO crs = 2, 4, 2
      CALL ring_operator(6, crs, .FALSE.)
      pco%a(PCO_BLOCK_SLOT(2, 2, 0, 1)) = 1d-17     ! slope row 2 of element 2 -> its own mean
      pco%a(PCO_BLOCK_SLOT(5, 3, 0, 1)) = -2d-17
      CALL PCO_REFRESH_SOLVER(st, msg)
      CALL EL_TEST_CHECK(st.EQ.0 .AND. pco%ncomp.EQ.1 .AND. pco%nclosed.EQ.1 .AND. pco%closure_closed.LT.1d-16, &
        'noise 1e-17 in slope-mean couplings: the ring stays closed, one reference', nf, rank)
    END DO
    CALL PCO_RELEASE()
  END SUBROUTINE test_closure_scale

  !---------------------------------------------------------------------
  SUBROUTINE test_structural_gauge()
    INTEGER :: st, crs, it
    CHARACTER(len=256) :: msg
    CHARACTER(len=32) :: lab
    REAL*8 :: xt(24), b(24), x(24)
    WRITE(*,'(A)') '17 gauge from the structure, the matrix as a consistency check'
    xt = 0d0
    xt(1:24:4) = 1d0
    crs = 2
    WRITE(lab,'(A,I0)') 'type ', crs
    CALL ring_operator(6, crs, .TRUE.)
    pco%a(PCO_BLOCK_SLOT(1, 1, 0, 1)) = 2d0 + 1d-9          ! weak outlet
    CALL PCO_REFRESH_SOLVER(st, msg)
    WRITE(*,'(A,ES11.3)') '    weak outlet: matrix closure measure', pco%closure_open
    CALL EL_TEST_CHECK(st.EQ.0 .AND. pco%nclosed.EQ.0 .AND. .NOT.ANY(pco%sref) .AND. pco%closure_open.LT.1d-8, &
      TRIM(lab)//': weak outlet (matrix closure below 1e-8) structurally open: no reference', nf, rank)
    CALL PCO_RAW_MATVEC(xt, b)
    it = 10
    CALL PCO_SOLVE_FULL(x, b, 24, it, 1d-12, 1d0)
    WRITE(*,'(A,ES11.3)') '    weak outlet: max |x - 1| on the means', MAXVAL(ABS(x(1:24:4) - 1d0))
    CALL EL_TEST_CHECK(MAXVAL(ABS(x(1:24:4) - 1d0)).LT.1d-5, &
      TRIM(lab)//': weak outlet solved: the exact solution 1 (a reference would return 0)', nf, rank)
    ! the P0 types are refused for a weakly open component
    DO crs = 3, 4
      WRITE(lab,'(A,I0)') 'type ', crs
      CALL ring_operator(6, crs, .TRUE.)
      pco%a(PCO_BLOCK_SLOT(1, 1, 0, 1)) = 2d0 + 1d-9
      CALL PCO_REFRESH_SOLVER(st, msg)
      CALL EL_TEST_CHECK(st.EQ.PCO_P0_WEAK_OPEN .AND. pco%nclosed.EQ.0, &
        TRIM(lab)//': weakly open component with the P0 system: refused (use type 2)', nf, rank)
    END DO
    CALL ring_operator(6, 4, .TRUE.)                          ! well open (0.11): accepted
    CALL PCO_REFRESH_SOLVER(st, msg)
    CALL EL_TEST_CHECK(st.EQ.0, 'type 4: a well open component (matrix measure 0.11) is accepted', nf, rank)
    ! structurally closed, but the operator does not annihilate the constant
    CALL ring_operator(6, 2, .TRUE.)
    CALL set_flux(6, 0)
    CALL PCO_REFRESH_SOLVER(st, msg)
    CALL EL_TEST_CHECK(st.EQ.PCO_GAUGE_UNDECIDED .AND. .NOT.ANY(pco%sref), &
      'structure closed, matrix open: undecided, no reference ('//TRIM(msg(1:60))//')', nf, rank)
    ! structurally open, but the constant is in the numerical kernel
    CALL ring_operator(6, 2, .FALSE.)
    CALL set_flux(6, 1)
    CALL PCO_REFRESH_SOLVER(st, msg)
    CALL EL_TEST_CHECK(st.EQ.PCO_GAUGE_UNDECIDED .AND. .NOT.ANY(pco%sref), &
      'structure open, matrix closed to rounding: undecided, no reference', nf, rank)
    ! structural flux between the thresholds
    CALL ring_operator(6, 2, .FALSE.)
    CALL PCO_SET_STRUCTURAL_FLUX((/1d-8, 0d0, 0d0, 0d0, 0d0, 0d0/))
    CALL PCO_REFRESH_SOLVER(st, msg)
    CALL EL_TEST_CHECK(st.EQ.PCO_GAUGE_UNDECIDED .AND. .NOT.ANY(pco%sref), &
      'structural flux 1e-8 (between 1e-10 and 1e-6): undecided, no reference', nf, rank)
    CALL PCO_RELEASE()
  END SUBROUTINE test_structural_gauge

  !---------------------------------------------------------------------
  SUBROUTINE test_compare_products()
    REAL*8 :: yd(3), y(3), sc(3)
    WRITE(*,'(A)') '18 worker-product comparison on zero rows'
    sc = (/2d0, 0d0, 1d0/)
    y = (/1d0, 0d0, -1d0/)
    yd = y
    yd(1) = yd(1) + 1d-16
    CALL EL_TEST_CHECK(PCO_COMPARE_PRODUCTS(yd, y, sc).LT.1d-15, 'rounding on scaled rows passes', nf, rank)
    yd(2) = 1d-30
    CALL EL_TEST_CHECK(PCO_COMPARE_PRODUCTS(yd, y, sc).GE.1d0, &
      'zero gathered row, nonzero worker result: fails (was skipped)', nf, rank)
    yd = y
    yd(3) = IEEE_VALUE(1d0, IEEE_QUIET_NAN)
    CALL EL_TEST_CHECK(PCO_COMPARE_PRODUCTS(yd, y, sc).GE.HUGE(1d0), '22 NaN worker product: fails', nf, rank)
    yd = y
    sc(1) = IEEE_VALUE(1d0, IEEE_POSITIVE_INF)
    CALL EL_TEST_CHECK(PCO_COMPARE_PRODUCTS(yd, y, sc).GE.HUGE(1d0), '22 infinite scale: fails', nf, rank)
  END SUBROUTINE test_compare_products

  !---------------------------------------------------------------------
  SUBROUTINE test_componentwise_tolerance()
    REAL*8 :: va(24), vb(24), per(3), tol
    INTEGER :: shift(3), kind, n, iv
    TYPE(tPcoClasses) :: cls
    REAL*8 :: pe(3) = (/1d0, 1d0, 1d0/)
    WRITE(*,'(A)') '19 coordinate tolerance per coordinate'
    tol = DEpsPrec
    per = (/1d0, 1d9, 1d9/)
    CALL cell(0d0, 0.5d0, va)
    CALL cell(1d0, 1.5d0, vb)                 ! image of va across x = 0 ~ 1 shares x = 1
    vb(2:24:3) = vb(2:24:3) + 0.8d0*tol      ! offset 0.8 tol in y and z: Euclidean 1.13 tol
    vb(3:24:3) = vb(3:24:3) + 0.8d0*tol
    CALL PCO_PAIR_SHIFT(va, vb, per, tol, shift, kind)
    CALL EL_TEST_CHECK(kind.EQ.1, 'pair check: offset 0.8 tol in two coordinates is coincident (as E013)', nf, rank)
    n = 3
    CALL make_mesh(n, .FALSE., .FALSE.)
    iv = (n+1)*(n+1)                          ! vertex (1,1,0)
    mg_mesh%level(2)%dcorvg(2, iv) = mg_mesh%level(2)%dcorvg(2, iv) - 0.8d0*tol
    mg_mesh%level(2)%dcorvg(3, iv) = mg_mesh%level(2)%dcorvg(3, iv) + 0.8d0*tol
    CALL PCO_BUILD_CLASSES(1, (/.TRUE., .TRUE., .FALSE./), pe, cls)
    CALL EL_TEST_CHECK(cls%nincomplete.EQ.0 .AND. cls%nambiguous.EQ.0, &
      'master classes: an image offset by 0.8 tol in two coordinates stays in its class', nf, rank)
    CALL make_mesh(n, .FALSE., .FALSE.)
    mg_mesh%level(2)%dcorvg(2, iv) = mg_mesh%level(2)%dcorvg(2, iv) - 1.2d0*tol
    CALL PCO_BUILD_CLASSES(1, (/.TRUE., .TRUE., .FALSE./), pe, cls)
    CALL EL_TEST_CHECK(cls%nincomplete.GE.1, 'master classes: 1.2 tol in one coordinate is not coincident', &
      nf, rank)
  END SUBROUTINE test_componentwise_tolerance

  !---------------------------------------------------------------------
  SUBROUTINE test_numerical_limit()
    INTEGER :: k, st, it, ist
    CHARACTER(len=256) :: msg
    CHARACTER(len=12) :: dec
    REAL*8 :: d, xt(24), b(24), x(24), merr, perr, bp(24), anorm
    REAL*8 :: fl(6)
    INTEGER :: i
    INTEGER :: ex(6) = (/10, 20, 30, 40, 50, 0/)
    TYPE(tUmfCSC) :: lu
    WRITE(*,'(A,2ES10.2)') '20 numerical limit of the gauge decision; closed bounds structure, matrix:', &
      PCO_FLUX_CLOSED, PCO_CLOSURE_TOL
    xt = 0d0
    xt(1:24:4) = 1d0
    DO k = 1, 6
      d = 0d0
      IF (ex(k).GT.0) d = 2d0**(-ex(k))
      CALL ring_operator(6, 2, .FALSE.)
      pco%a(PCO_BLOCK_SLOT(1, 1, 0, 1)) = 2d0 + 4d0*d          ! outlet of relative strength d
      fl = 0d0
      fl(1) = d
      CALL PCO_SET_STRUCTURAL_FLUX(fl)
      CALL PCO_REFRESH_SOLVER(st, msg)
      IF (st.EQ.PCO_GAUGE_UNDECIDED) THEN
        dec = 'undecided'
      ELSE IF (pco%nclosed.EQ.1) THEN
        dec = 'closed'
      ELSE
        dec = 'open'
      END IF
      ! the ungauged system, solved directly: b = A xt (mean 1 on every element)
      CALL PCO_RAW_MATVEC(xt, b)
      CALL UmfCSC_SetStructure(lu, 24, pco%ldA, pco%colA)
      CALL UmfCSC_SetValues(lu, pco%a)
      CALL UmfCSC_Analyse(lu, ist)
      CALL UmfCSC_Factorize(lu, ist)
      merr = -1d0
      perr = -1d0
      IF (ist.EQ.0) THEN
        CALL UmfCSC_Solve(lu, x, b, ist)
        merr = ABS(SUM(x(1:24:4))/6d0 - 1d0)
        ! the same with the right-hand side perturbed at the backward-error level of a stable
        ! solve, eps |A| |x| (normwise): the accuracy any computed coarse solution can claim
        anorm = 0d0
        DO i = 1, 24
          anorm = MAX(anorm, SUM(ABS(pco%a(pco%ldA(i):pco%ldA(i+1)-1))))
        END DO
        DO i = 1, 24
          bp(i) = b(i) + EPSILON(1d0)*anorm*MAXVAL(ABS(xt))*MERGE(1d0, -1d0, MOD(i*7,3).EQ.0)
        END DO
        CALL UmfCSC_Solve(lu, x, bp, ist)
        perr = ABS(SUM(x(1:24:4))/6d0 - 1d0)
      END IF
      CALL UmfCSC_Release(lu)
      WRITE(*,'(A,I3,A,ES9.2,A,A10,A,I2,A,ES10.2,A,ES10.2)') '    d = 2^-', ex(k), ' (', d, ') decision ', dec, &
        '   ungauged direct solve: status', ist, ', mean error exact rhs', merr, ', rhs perturbed by eps|A||x|', perr
      SELECT CASE (ex(k))
      CASE (10)
        CALL EL_TEST_CHECK(TRIM(dec).EQ.'open' .AND. merr.GE.0d0 .AND. merr.LT.1d-10, &
          'd = 2^-10: open, the ungauged system determines the mean', nf, rank)
      CASE (20, 30, 40)
        CALL EL_TEST_CHECK(TRIM(dec).EQ.'undecided' .AND. .NOT.ANY(pco%sref), &
          'd between the closed bound and 1e-6: undecided, nothing pinned', nf, rank)
      CASE DEFAULT
        CALL EL_TEST_CHECK(TRIM(dec).EQ.'closed' .AND. (ist.NE.0 .OR. perr.GT.1d-2), &
          'd <= closed bound: closed; the ungauged mean is lost at the backward-error level', nf, rank)
      END SELECT
    END DO
    CALL PCO_RELEASE()
  END SUBROUTINE test_numerical_limit

  !---------------------------------------------------------------------
  SUBROUTINE test_current_measure()
    INTEGER :: st
    CHARACTER(len=256) :: msg
    WRITE(*,'(A)') '21 the P0 admission follows the current matrix measure'
    CALL ring_operator(6, 2, .TRUE.)                        ! outlet 2.5: measure 0.11
    CALL PCO_REFRESH_SOLVER(st, msg)
    CALL EL_TEST_CHECK(st.EQ.0 .AND. PCO_P0_ADMISSIBLE(4), 'outlet 2.5: type 4 admissible', nf, rank)
    pco%a(PCO_BLOCK_SLOT(1, 1, 0, 1)) = 2d0 + 1d-3           ! same layout, weaker outlet
    CALL PCO_REFRESH_SOLVER(st, msg)
    WRITE(*,'(A,ES10.2)') '    current open measure', pco%closure_open
    CALL EL_TEST_CHECK(st.EQ.0 .AND. .NOT.PCO_P0_ADMISSIBLE(4) .AND. pco%closure_open.LT.1d-3, &
      'outlet 2 + 1e-3 on the same layout: type 4 no longer admissible (current measure)', nf, rank)
    CALL PCO_CHANGE_SOLVER_TYPE(4, st, msg)
    CALL EL_TEST_CHECK(st.EQ.PCO_P0_WEAK_OPEN, 'a type-4 copy on it is refused at the refresh', nf, rank)
    CALL PCO_RELEASE()
  END SUBROUTINE test_current_measure

  ! unit-height hexahedron [x0,x1] x [0,1] x [0,1], vertex coordinates x1,y1,z1, x2, ...
  SUBROUTINE cell(x0, x1, v)
    REAL*8, INTENT(IN) :: x0, x1
    REAL*8, INTENT(OUT) :: v(24)
    INTEGER :: i
    REAL*8 :: c(3,8)
    c(:,1) = (/x0, 0d0, 0d0/)
    c(:,2) = (/x1, 0d0, 0d0/)
    c(:,3) = (/x1, 1d0, 0d0/)
    c(:,4) = (/x0, 1d0, 0d0/)
    DO i = 1, 4
      c(:,i+4) = c(:,i) + (/0d0, 0d0, 1d0/)
    END DO
    v = RESHAPE(c, (/24/))
  END SUBROUTINE cell

END PROGRAM TEST_PCO_COARSE
