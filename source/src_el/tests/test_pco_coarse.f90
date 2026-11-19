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

  USE MPI
  USE TYPES, ONLY: tMultiMesh, TMatrix
  USE var_QuadScalar, ONLY: mg_mesh
  USE UMFPackSolver, ONLY: tUmfCSC, UmfCSC_SetStructure, UmfCSC_SetValues, UmfCSC_Analyse, &
                           UmfCSC_Factorize, UmfCSC_Solve, UmfCSC_Release, &
                           myUmfPack_Factorize, myUmfPack_Solve, myUmfPack_Free
  USE pco_periodic_topology, ONLY: tPcoClasses, PCO_BUILD_CLASSES, PCO_CLASS_MASK_CONFLICTS, &
                                   PCO_ELEMENT_GRAPH
  USE pco_p0_reconstruct, ONLY: PCO_SET_VERTEX_CLASSES, PCO_P0_TO_VERTEX
  USE pco_coarse_solver, ONLY: pco, PCO_INSTALL_LAYOUT, PCO_REFRESH_SOLVER, PCO_SOLVE_CHECK, &
                               PCO_BLOCK_SLOT, PCO_NEIGHBOUR_POS, PCO_SOLVE_FULL, PCO_SOLVE_P0, &
                               PCO_RELEASE
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
    REAL*8, ALLOCATABLE :: dp(:), dpl(:), dax(:), dl0(:), da0(:), num(:), den(:)
    REAL*8 :: c(3), pi, worst_img, worst_val, worst_int, legacy_gap, nn, dd
    REAL*8 :: per(3) = (/1d0, 1d0, 1d0/)
    LOGICAL :: on_image

    WRITE(*,'(A)') '4 P0 -> vertex averaging at periodic planes'
    CALL make_mesh(n, .FALSE., .FALSE.)
    nel = mg_mesh%level(1)%nel
    nvt = mg_mesh%level(1)%nvt
    pi = 4d0*ATAN(1d0)
    ALLOCATE(dp(nel), dpl(nvt), dax(nvt), dl0(nvt), da0(nvt), num(nvt), den(nvt))
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
    CALL PCO_P0_TO_VERTEX(dp, dpl, dax, mg_mesh%level(1)%dvol, mg_mesh%level(1)%kvert, nel, nvt)
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
    DEALLOCATE(dp, dpl, dax, dl0, da0, num, den)
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
  END SUBROUTINE ring_operator

  SUBROUTINE test_gauge_solver()
    INTEGER, PARAMETER :: nel = 6
    INTEGER :: crs, st, i, it
    CHARACTER(len=256) :: msg
    CHARACTER(len=64) :: lab
    REAL*8 :: v(4*nel), resid, devi, b(4*nel), x(4*nel), p0(nel), tol
    LOGICAL :: incompat

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
      CALL PCO_SOLVE_CHECK(v, resid, devi, incompat)
      WRITE(*,'(A,2ES11.3)') '    residual, deviation from v + kernel:', resid, devi
      tol = MERGE(1d-12, 1d-9, crs.EQ.2 .OR. crs.EQ.4)
      CALL EL_TEST_CHECK(resid.LT.tol .AND. devi.LT.100d0*tol .AND. .NOT.incompat, &
        TRIM(lab)//': compatible rhs, ungauged residual and solution up to the kernel', nf, rank)
      CALL EL_TEST_CHECK(ABS(pco%sol(1)).LE.tol*MAXVAL(ABS(pco%sol)), &
        TRIM(lab)//': reference mean value zero', nf, rank)
      ! incompatible: a constant on the means is not in the range
      b = 0d0
      DO i = 1, nel
        b(4*(i-1)+1) = 1d0
      END DO
      it = 500
      i = pco%nincompat
      IF (crs.LE.2) THEN
        CALL PCO_SOLVE_FULL(x, b, it, 1d-14, 1d0)
      ELSE
        CALL PCO_SOLVE_P0(b, p0, it, 1d-14, 1d0)
      END IF
      CALL EL_TEST_CHECK(pco%nincompat.EQ.i+1, TRIM(lab)//': incompatible rhs reported', nf, rank)
      ! open ring: no reference, exact solve
      CALL ring_operator(nel, crs, .TRUE.)
      CALL PCO_REFRESH_SOLVER(st, msg)
      CALL EL_TEST_CHECK(st.EQ.0 .AND. pco%nclosed.EQ.0 .AND. .NOT.ANY(pco%sref), &
        TRIM(lab)//': open ring, no reference', nf, rank)
      CALL PCO_SOLVE_CHECK(v, resid, devi, incompat)
      WRITE(*,'(A,2ES11.3)') '    open: residual, deviation:', resid, devi
      CALL EL_TEST_CHECK(resid.LT.tol .AND. devi.LT.100d0*tol, TRIM(lab)//': open ring solved exactly', nf, rank)
    END DO
    CALL PCO_RELEASE()
  END SUBROUTINE test_gauge_solver

END PROGRAM TEST_PCO_COARSE
