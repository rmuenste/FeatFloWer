!=========================================================================
! test_chi_periodic - Phase 5 periodic-box geometry test (chimera-
! integration-design.md v3, section 9; roadmap row 5).
!
!   1. CHI_PER_DELTA / CHI_PER_WRAP / CHI_PER_DIST: minimum image and
!      wrap identities, and the exact identity when no axis is periodic.
!   2. On a uniform 8x8x8 box [0,1]^3 a sphere body at the box CORNER
!      (0,0,0) - its atmosphere straddles all six periodic faces - is
!      classified (hole/fringe markers) and penalised (nodal lumped
!      operator) with the periodic box; the same body translated to the
!      box centre is classified without periodicity.  The Q2 nodes lie on
!      the 1/16 lattice, so the translation (0.5,0.5,0.5) maps nodes onto
!      nodes: the corner markers must equal the centred markers at the
!      translated node, the penalised volume sum_i D_L(i) must agree to
!      round-off, and the corner case must contain hole nodes on all
!      eight box corners (nothing visible without the periodic images).
!   3. Wrapped outer sample points: for points on the corner body's outer
!      sphere the wrap lands inside [0,1)^3 and preserves the minimum-
!      image distance to the body centre.
!=========================================================================
PROGRAM test_chi_periodic
  USE types, ONLY: tMultiMesh
  USE mesh_structures, ONLY: readTriCoarse, genMeshStructures, &
    getNumberOfEdgesOnVerts
  USE CHI_PERIODIC, ONLY: tChiPeriodic, CHI_PER_ACTIVE, CHI_PER_DELTA, CHI_PER_WRAP, &
    CHI_PER_DIST
  USE CHI_MARKERS, ONLY: tChiBody, CHI_BODY_SPHERE, CHI_CLASSIFY_MARKERS, &
    CHI_MARK_FREE, CHI_MARK_FRINGE, CHI_MARK_HOLE, CHI_POINT_IN_BODY, CHI_BODY_NEAR_BOX
  USE CHI_FEM_EVAL, ONLY: CHI_Q2_DOFMAP, CHI_Q2_REFNODES
  USE CHI_GEOMETRY, ONLY: CHI_Q1_MAP
  USE CHI_PENALTY, ONLY: CHI_PENALTY_NODAL, CHI_PENALTY_BETA
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  INTEGER, PARAMETER :: NX = 8, NL = 2*NX   ! Q2 lattice per axis
  TYPE(tMultiMesh) :: m
  TYPE(tChiPeriodic) :: pb, pb0
  TYPE(tChiBody) :: bcen(1), bcor(1)
  INTEGER :: ierr, nfail, ndof, nel, nvt, net, nat, e, i, q, idx(27), key(3), kk(3)
  INTEGER :: nhole_c, nhole_k, nfr_c, nfr_k, nmis, ncorner, nbad
  INTEGER, ALLOCATABLE :: kind_c(:), pid_c(:), kind_k(:), pid_k(:), tab(:,:,:)
  INTEGER, ALLOCATABLE :: nbo(:)
  LOGICAL, ALLOCATABLE :: inb(:)
  REAL*8, ALLOCATABLE :: q2c(:,:), dl_c(:), dl_k(:), xn(:,:)
  REAL*8 :: refxi(3,27), nodes(3,8), x(3), jac(3,3), detj, d(3), xw(3), hw(1)
  REAL*8 :: r, r2, beta, err, vc, vk, th, ph
  LOGICAL :: inbody
  REAL*8, PARAMETER :: PI = 3.14159265358979323846d0

  nfail = 0
  CALL MPI_INIT(ierr)

  ! ---- 1. algebra of the periodic helpers ------------------------------
  pb%per = .TRUE.; pb%len = (/ 1d0, 2d0, 3d0 /); pb%lo = (/ 0d0, -1d0, 0.5d0 /)
  d = CHI_PER_DELTA(pb, (/ 0.9d0, -0.9d0, 3.4d0 /), (/ 0.1d0, 0.9d0, 0.6d0 /))
  err = MAXVAL(ABS(d - (/ -0.2d0, 0.2d0, -0.2d0 /)))
  WRITE(*,'(A,ES10.2)') 'test_chi_periodic: minimum-image delta error = ', err
  IF (err .GT. 1d-14) CALL fail('CHI_PER_DELTA wrong', nfail)
  xw = CHI_PER_WRAP(pb, (/ 1d0, -1d0, 0.5d0 + 3d0 /))
  err = MAXVAL(ABS(xw - (/ 0d0, -1d0, 0.5d0 /)))
  IF (err .GT. 1d-14) CALL fail('CHI_PER_WRAP of x = lo + len must give lo', nfail)
  xw = CHI_PER_WRAP(pb, (/ -0.25d0, 5.5d0, 0.25d0 /))
  err = MAXVAL(ABS(xw - (/ 0.75d0, -0.5d0, 3.25d0 /)))
  IF (err .GT. 1d-14) CALL fail('CHI_PER_WRAP wrong', nfail)
  IF (CHI_PER_ACTIVE(pb0)) CALL fail('default box must be inactive', nfail)
  d = CHI_PER_DELTA(pb0, (/ 0.9d0, -0.9d0, 3.4d0 /), (/ 0.1d0, 0.9d0, 0.6d0 /))
  IF (ANY(d .NE. (/ 0.9d0, -0.9d0, 3.4d0 /) - (/ 0.1d0, 0.9d0, 0.6d0 /))) &
    CALL fail('inactive delta must be the plain difference (bitwise)', nfail)
  xw = CHI_PER_WRAP(pb0, (/ -0.25d0, 5.5d0, 0.25d0 /))
  IF (ANY(xw .NE. (/ -0.25d0, 5.5d0, 0.25d0 /))) &
    CALL fail('inactive wrap must be the identity (bitwise)', nfail)
  ! partial periodicity: only y
  pb%per = (/ .FALSE., .TRUE., .FALSE. /)
  d = CHI_PER_DELTA(pb, (/ 0.9d0, -0.9d0, 3.4d0 /), (/ 0.1d0, 0.9d0, 0.6d0 /))
  err = MAXVAL(ABS(d - (/ 0.8d0, 0.2d0, 2.8d0 /)))
  IF (err .GT. 1d-14) CALL fail('partial periodicity wrong', nfail)

  ! ---- 2. corner body vs centred body on the box mesh ------------------
  CALL write_box('tcp_box.tri', NX)
  CALL load('tcp_box.tri', m)
  nel = m%level(1)%nel; nvt = m%level(1)%nvt
  net = m%level(1)%net; nat = m%level(1)%nat
  ndof = nvt + net + nat + nel
  ALLOCATE(q2c(3,ndof))
  CALL CHI_Q2_REFNODES(refxi)
  DO e = 1, nel
    DO i = 1, 8
      nodes(:,i) = m%level(1)%dcorvg(:, m%level(1)%kvert(i,e))
    END DO
    CALL CHI_Q2_DOFMAP(e, m%level(1)%kvert, m%level(1)%kedge, m%level(1)%karea, &
      nvt, net, nat, idx)
    DO i = 1, 27
      CALL CHI_Q1_MAP(nodes, refxi(:,i), x, jac, detj)
      q2c(:,idx(i)) = x
    END DO
  END DO

  bcen(1)%shape = CHI_BODY_SPHERE
  bcen(1)%center = (/ 0.5d0, 0.5d0, 0.5d0 /)
  bcen(1)%radius = 0.17d0
  bcor(1) = bcen(1)
  bcor(1)%center = 0d0
  hw(1) = 0.12d0

  pb = tChiPeriodic()
  pb%per = .TRUE.; pb%len = 1d0; pb%lo = 0d0

  ALLOCATE(kind_c(ndof), pid_c(ndof), kind_k(ndof), pid_k(ndof))
  CALL CHI_CLASSIFY_MARKERS(nel, nvt, net, nat, m%level(1)%kvert, m%level(1)%kedge, &
    m%level(1)%karea, q2c, 1, bcen, kind_c, pid_c)
  CALL CHI_CLASSIFY_MARKERS(nel, nvt, net, nat, m%level(1)%kvert, m%level(1)%kedge, &
    m%level(1)%karea, q2c, 1, bcor, kind_k, pid_k, pb)
  nhole_c = COUNT(kind_c .EQ. CHI_MARK_HOLE); nfr_c = COUNT(kind_c .EQ. CHI_MARK_FRINGE)
  nhole_k = COUNT(kind_k .EQ. CHI_MARK_HOLE); nfr_k = COUNT(kind_k .EQ. CHI_MARK_FRINGE)
  WRITE(*,'(A,I0,A,I0,A,I0,A,I0)') 'test_chi_periodic: centred hole/fringe = ', nhole_c, &
    '/', nfr_c, ', corner (periodic) hole/fringe = ', nhole_k, '/', nfr_k
  IF (nhole_c .LE. 0 .OR. nfr_c .LE. 0) CALL fail('centred body classifies nothing', nfail)
  ! the corner body sees the box only through its images: strictly more
  ! marked dofs than the centred body (periodic-face nodes are duplicated)
  IF (nhole_k .LT. nhole_c) CALL fail('corner body has fewer hole nodes than centred', nfail)

  ! lattice table of the centred classification
  ALLOCATE(tab(0:NL, 0:NL, 0:NL))
  tab = -1
  DO i = 1, ndof
    key = NINT(q2c(:,i)*NL)
    IF (MAXVAL(ABS(q2c(:,i)*NL - key)) .GT. 1d-10) CALL fail('node off lattice', nfail)
    tab(key(1), key(2), key(3)) = kind_c(i)
  END DO
  nmis = 0
  ncorner = 0
  DO i = 1, ndof
    key = NINT(q2c(:,i)*NL)
    kk = MOD(key + NX, NL)          ! translate by (0.5,0.5,0.5), wrap
    IF (tab(kk(1), kk(2), kk(3)) .NE. kind_k(i)) nmis = nmis + 1
    IF (ALL(key .EQ. 0 .OR. key .EQ. NL) .AND. kind_k(i) .EQ. CHI_MARK_HOLE) &
      ncorner = ncorner + 1
  END DO
  WRITE(*,'(A,I0,A,I0)') 'test_chi_periodic: translated-marker mismatches = ', nmis, &
    ', hole box corners = ', ncorner
  IF (nmis .NE. 0) CALL fail('corner markers differ from translated centred markers', nfail)
  IF (ncorner .NE. 8) CALL fail('expected all 8 box corners inside the hole', nfail)
  ! duplicates on periodic faces must agree (pure function of the coordinates)
  nbad = 0
  DO i = 1, ndof
    IF (kind_k(i) .NE. CHI_MARK_FREE .AND. pid_k(i) .NE. 1) nbad = nbad + 1
  END DO
  IF (nbad .NE. 0) CALL fail('marked node without body id', nfail)

  ! nodal penalty: total penalised volume agrees
  ALLOCATE(dl_c(ndof), dl_k(ndof), nbo(ndof), inb(ndof), xn(3,ndof))
  CALL CHI_PENALTY_NODAL(nel, nvt, net, nat, m%level(1)%kvert, m%level(1)%kedge, &
    m%level(1)%karea, m%level(1)%dcorvg, 1, bcen, hw, 1d0, ndof, dl_c, nbo, inb, xn, &
    0.25d0, 0.5d0)
  CALL CHI_PENALTY_NODAL(nel, nvt, net, nat, m%level(1)%kvert, m%level(1)%kedge, &
    m%level(1)%karea, m%level(1)%dcorvg, 1, bcor, hw, 1d0, ndof, dl_k, nbo, inb, xn, &
    0.25d0, 0.5d0, pb)
  vc = SUM(dl_c); vk = SUM(dl_k)
  err = ABS(vc - vk)/vc
  WRITE(*,'(A,ES14.7,A,ES14.7,A,ES10.2)') 'test_chi_periodic: penalised volume centred = ', &
    vc, ', corner = ', vk, ', rel. diff = ', err
  IF (err .GT. 1d-11) CALL fail('penalised volume differs between centred and corner', nfail)
  IF (MINVAL(dl_k) .LT. 0d0) CALL fail('negative nodal penalty weight', nfail)
  ! every penalised node of the corner body is inside R + 0.5 H (min image)
  nbad = 0
  DO i = 1, ndof
    IF (dl_k(i) .LE. 0d0) CYCLE
    r = CHI_PER_DIST(pb, q2c(:,i), bcor(1)%center)
    IF (r .GT. bcor(1)%radius + 0.5d0*hw(1) + 1d-12) nbad = nbad + 1
    IF (.NOT. CHI_BODY_NEAR_BOX(bcor(1), q2c(:,i), q2c(:,i), bcor(1)%radius + hw(1), pb)) &
      nbad = nbad + 1
  END DO
  IF (nbad .NE. 0) CALL fail('penalised node outside the minimum-image support', nfail)
  ! point-in-body with the periodic box at the far corner
  IF (.NOT. CHI_POINT_IN_BODY(bcor(1), (/ 1d0, 1d0, 0.9d0 /), pb)) &
    CALL fail('far corner point must be inside the corner body (min image)', nfail)
  IF (CHI_POINT_IN_BODY(bcor(1), (/ 1d0, 1d0, 0.9d0 /))) &
    CALL fail('without periodicity the far corner is outside', nfail)
  CALL CHI_PENALTY_BETA(bcor(1), hw(1), (/ 0.95d0, 1d0, 0.02d0 /), beta, inbody, &
    0.25d0, 0.5d0, pb)
  IF (beta .NE. 1d0 .OR. .NOT. inbody) CALL fail('beta at the far image must be 1', nfail)

  ! ---- 3. wrapped outer sample points ----------------------------------
  nbad = 0
  DO i = 1, 200
    th = PI*DBLE(i)/201d0
    ph = 2d0*PI*DBLE(7*i)/200d0
    r = bcor(1)%radius + hw(1)
    x = bcor(1)%center + r*(/ SIN(th)*COS(ph), SIN(th)*SIN(ph), COS(th) /)
    xw = CHI_PER_WRAP(pb, x)
    IF (ANY(xw .LT. 0d0) .OR. ANY(xw .GE. 1d0)) nbad = nbad + 1
    IF (ABS(CHI_PER_DIST(pb, xw, bcor(1)%center) - r) .GT. 1d-13) nbad = nbad + 1
    ! the unwrapped image in the body frame is the original point
    d = bcor(1)%center + CHI_PER_DELTA(pb, xw, bcor(1)%center)
    IF (MAXVAL(ABS(d - x)) .GT. 1d-13) nbad = nbad + 1
  END DO
  WRITE(*,'(A,I0)') 'test_chi_periodic: wrapped sample-point failures = ', nbad
  IF (nbad .NE. 0) CALL fail('wrap/image round trip failed', nfail)

  CALL MPI_FINALIZE(ierr)
  IF (nfail .GT. 0) THEN
    WRITE(*,'(A,I0,A)') 'test_chi_periodic: ', nfail, ' check(s) FAILED'
    STOP 1
  END IF
  WRITE(*,'(A)') 'test_chi_periodic: all checks passed'

CONTAINS

  SUBROUTINE fail(msg, n)
    CHARACTER(*), INTENT(IN) :: msg
    INTEGER, INTENT(INOUT) :: n
    n = n + 1
    WRITE(*,'(A)') 'FAIL: ' // msg
  END SUBROUTINE fail

  SUBROUTINE write_box(fname, n)
    CHARACTER(*), INTENT(IN) :: fname
    INTEGER, INTENT(IN) :: n
    INTEGER :: iu, ix, iy, iz, i
    iu = 91
    OPEN(iu, FILE=fname, STATUS='REPLACE')
    WRITE(iu,'(A)') 'test_chi_periodic box'
    WRITE(iu,'(A)') 'Parametrisierung PARXC, PARYC, TMAXC'
    WRITE(iu,'(I0,1X,I0,A)') n*n*n, (n+1)**3, ' 1 8 12 6     NEL NVT NBCT NVE NEE NAE'
    WRITE(iu,'(A)') 'DCORVG'
    DO iz = 0, n
      DO iy = 0, n
        DO ix = 0, n
          WRITE(iu,'(3ES24.16)') DBLE(ix)/n, DBLE(iy)/n, DBLE(iz)/n
        END DO
      END DO
    END DO
    WRITE(iu,'(A)') 'KVERT'
    DO iz = 0, n-1
      DO iy = 0, n-1
        DO ix = 0, n-1
          WRITE(iu,'(8(I0,1X))') vid(n,ix,iy,iz), vid(n,ix+1,iy,iz), vid(n,ix+1,iy+1,iz), &
            vid(n,ix,iy+1,iz), vid(n,ix,iy,iz+1), vid(n,ix+1,iy,iz+1), &
            vid(n,ix+1,iy+1,iz+1), vid(n,ix,iy+1,iz+1)
        END DO
      END DO
    END DO
    WRITE(iu,'(A)') 'KNPR'
    DO i = 1, (n+1)**3
      WRITE(iu,'(I0)') 1
    END DO
    CLOSE(iu)
  END SUBROUTINE write_box

  INTEGER FUNCTION vid(n, ix, iy, iz)
    INTEGER, INTENT(IN) :: n, ix, iy, iz
    vid = 1 + ix + (n+1)*(iy + (n+1)*iz)
  END FUNCTION vid

  SUBROUTINE load(fname, mm)
    CHARACTER(*), INTENT(IN) :: fname
    TYPE(tMultiMesh), INTENT(INOUT) :: mm
    INTEGER :: noe
    ALLOCATE(mm%level(1))
    mm%nlmin = 1; mm%nlmax = 1; mm%maxlevel = 1
    mm%level(1)%nel = 0
    CALL readTriCoarse(fname, mm)
    CALL getNumberOfEdgesOnVerts(mm%level(1), noe)
    CALL genMeshStructures(mm, .FALSE., 1, noe)
  END SUBROUTINE load

END PROGRAM test_chi_periodic
