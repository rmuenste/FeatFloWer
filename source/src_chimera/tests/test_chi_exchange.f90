!=========================================================================
! test_chi_exchange - Phase 3 test of the collective background
! evaluation (chimera-integration-design.md v3, sections 3, 4, 9).
! Run with 1, 2 and 3 MPI ranks.
!
! Each rank owns the slab [rank, rank+1] x [0,1] x [0,1] meshed 2x2x2
! (written as a .tri and loaded through the production mesh path), and
! carries the Q2 interpolant of a quadratic velocity field plus the P1
! representation of a linear pressure.  All ranks query the IDENTICAL
! replicated point list (interior points, points exactly on the slab
! interfaces, and one point outside the domain) and check:
!   - values, gradients and pressure reproduce the analytic fields to
!     round-off (Q2 is exact for quadratics, P1 for linears);
!   - the returned arrays are bit-identical on every rank;
!   - interface points are owned by the LOWER rank (MIN rule);
!   - the outside point is reported as missing (and only that one).
!=========================================================================
PROGRAM test_chi_exchange

  USE types, ONLY: tMultiMesh
  USE mesh_structures, ONLY: readTriCoarse, genMeshStructures, &
    getNumberOfEdgesOnVerts
  USE CHI_LOCATOR, ONLY: tChimeraLocator, CHI_LOCATOR_BUILD, CHI_LOCATOR_RELEASE
  USE CHI_EXCHANGE, ONLY: CHI_EXCHANGE_BG_EVAL, CHI_BG_NVAL
  USE CHI_FEM_EVAL, ONLY: CHI_Q2_DOFMAP, CHI_Q2_REFNODES
  USE CHI_GEOMETRY, ONLY: CHI_Q1_MAP

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER :: ierr, myrank, nranks, nfail, ndof, nel, npts, ip, i, e, nmissing
  TYPE(tMultiMesh) :: m
  TYPE(tChimeraLocator) :: loc
  REAL*8, ALLOCATABLE :: u(:), v(:), w(:), p(:), q2c(:,:), pts(:,:), vals(:,:), vals0(:,:)
  INTEGER, ALLOCATABLE :: owner(:)
  REAL*8 :: x(3), ua(3), ga(3,3), pa, err, maxerr, maxdiff
  INTEGER*8 :: seed
  CHARACTER(LEN=64) :: fname

  nfail = 0
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)

  WRITE(fname,'(A,I0,A)') 'tce_slab_', myrank, '.tri'
  CALL write_slab(fname, DBLE(myrank), DBLE(myrank+1))
  ALLOCATE(m%level(1))
  m%nlmin = 1; m%nlmax = 1; m%maxlevel = 1
  m%level(1)%nel = 0
  CALL readTriCoarse(fname, m)
  CALL getNumberOfEdgesOnVerts(m%level(1), i)
  CALL genMeshStructures(m, .FALSE., 1, i)
  nel = m%level(1)%nel
  ndof = m%level(1)%nvt + m%level(1)%net + m%level(1)%nat + nel

  ! nodal Q2 interpolant + P1 pressure
  CALL q2coords(m, q2c)
  ALLOCATE(u(ndof), v(ndof), w(ndof), p(4*nel))
  DO i = 1, ndof
    CALL fields(q2c(:,i), ua, ga, pa)
    u(i) = ua(1); v(i) = ua(2); w(i) = ua(3)
  END DO
  DO e = 1, nel
    CALL centroid(m, e, x)
    CALL fields(x, ua, ga, pa)
    p(4*(e-1)+1) = pa
    p(4*(e-1)+2) = 2d0        ! dp/dx
    p(4*(e-1)+3) = -1d0       ! dp/dy
    p(4*(e-1)+4) = 0.5d0      ! dp/dz
  END DO
  CALL CHI_LOCATOR_BUILD(loc, m%level(1)%dcorvg, m%level(1)%kvert, nel, m%level(1)%nvt)

  ! replicated query list: 60 interior points, nranks-1 interface points
  ! (x = 1, 2, ...), one outside point
  npts = 60 + (nranks-1) + 1
  ALLOCATE(pts(3,npts), vals(CHI_BG_NVAL,npts), vals0(CHI_BG_NVAL,npts), owner(npts))
  seed = 12345_8
  DO ip = 1, 60
    pts(1,ip) = DBLE(nranks)*lcg(seed)
    pts(2,ip) = lcg(seed)
    pts(3,ip) = lcg(seed)
  END DO
  DO ip = 1, nranks-1
    pts(:,60+ip) = (/ DBLE(ip), 0.3d0, 0.7d0 /)
  END DO
  pts(:,npts) = (/ -0.5d0, 0.5d0, 0.5d0 /)

  CALL CHI_EXCHANGE_BG_EVAL(MPI_COMM_WORLD, npts, pts, loc, m%level(1)%dcorvg, &
    m%level(1)%kvert, m%level(1)%kedge, m%level(1)%karea, m%level(1)%nvt, &
    m%level(1)%net, m%level(1)%nat, u, v, w, p, vals, owner, nmissing)

  IF (nmissing .NE. 1 .OR. owner(npts) .NE. nranks) THEN
    WRITE(*,*) 'FAIL: missing-point report', nmissing, owner(npts)
    nfail = nfail + 1
  END IF
  DO ip = 1, nranks-1
    IF (owner(60+ip) .NE. ip-1) THEN
      WRITE(*,*) 'FAIL: interface point owner (MIN rule)', ip, owner(60+ip)
      nfail = nfail + 1
    END IF
  END DO
  maxerr = 0d0
  DO ip = 1, npts-1
    CALL fields(pts(:,ip), ua, ga, pa)
    err = MAXVAL(ABS(vals(1:3,ip) - ua))
    DO i = 1, 3
      err = MAX(err, MAXVAL(ABS(vals(3+3*(i-1)+1:3+3*i,ip) - ga(:,i))))
    END DO
    err = MAX(err, ABS(vals(13,ip) - pa))
    maxerr = MAX(maxerr, err)
  END DO
  IF (maxerr .GT. 1d-10) THEN
    WRITE(*,*) 'FAIL: exactness of the exchanged values, max err =', maxerr
    nfail = nfail + 1
  END IF

  ! bit identity across ranks
  vals0 = vals
  CALL MPI_BCAST(vals0, CHI_BG_NVAL*npts, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  maxdiff = MAXVAL(ABS(vals0 - vals))
  IF (maxdiff .NE. 0d0) THEN
    WRITE(*,*) 'FAIL: values not bit-identical across ranks, rank', myrank, maxdiff
    nfail = nfail + 1
  END IF

  CALL MPI_ALLREDUCE(MPI_IN_PLACE, nfail, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  IF (myrank .EQ. 0) THEN
    WRITE(*,'(A,I0,A,ES10.2,A,I0)') ' ranks = ', nranks, ', max field error = ', &
      maxerr, ', missing = ', nmissing
  END IF
  CALL CHI_LOCATOR_RELEASE(loc)
  OPEN(92, FILE=fname, STATUS='OLD'); CLOSE(92, STATUS='DELETE')
  CALL MPI_FINALIZE(ierr)
  IF (nfail .GT. 0) THEN
    IF (myrank .EQ. 0) WRITE(*,*) 'test_chi_exchange: ', nfail, ' failure(s)'
    STOP 1
  END IF
  IF (myrank .EQ. 0) WRITE(*,*) 'test_chi_exchange: PASS'

CONTAINS

  REAL*8 FUNCTION lcg(s)
    INTEGER*8, INTENT(INOUT) :: s
    s = MOD(s*48271_8, 2147483647_8)
    lcg = DBLE(s)/2147483647d0
  END FUNCTION lcg

  ! quadratic velocity, linear pressure (exactly representable)
  SUBROUTINE fields(x, ua, ga, pa)
    REAL*8, INTENT(IN) :: x(3)
    REAL*8, INTENT(OUT) :: ua(3), ga(3,3), pa
    ua(1) = x(1)*x(1) + 2d0*x(1)*x(2) + x(3)
    ua(2) = x(2)*x(2) - x(1) + 0.5d0*x(2)*x(3)
    ua(3) = x(1)*x(3) + 3d0 - x(2)*x(2)
    ga = 0d0                                 ! ga(a,b) = d u_a / d x_b
    ga(1,1) = 2d0*x(1) + 2d0*x(2); ga(1,2) = 2d0*x(1); ga(1,3) = 1d0
    ga(2,1) = -1d0; ga(2,2) = 2d0*x(2) + 0.5d0*x(3); ga(2,3) = 0.5d0*x(2)
    ga(3,1) = x(3); ga(3,2) = -2d0*x(2); ga(3,3) = x(1)
    pa = 1d0 + 2d0*x(1) - x(2) + 0.5d0*x(3)
  END SUBROUTINE fields

  SUBROUTINE centroid(mm, e, xc)
    TYPE(tMultiMesh), INTENT(IN) :: mm
    INTEGER, INTENT(IN) :: e
    REAL*8, INTENT(OUT) :: xc(3)
    REAL*8 :: nodes(3,8), jac(3,3), detj
    INTEGER :: i
    DO i = 1, 8
      nodes(:,i) = mm%level(1)%dcorvg(:,mm%level(1)%kvert(i,e))
    END DO
    CALL CHI_Q1_MAP(nodes, (/0d0,0d0,0d0/), xc, jac, detj)
  END SUBROUTINE centroid

  SUBROUTINE q2coords(mm, q2c)
    TYPE(tMultiMesh), INTENT(IN) :: mm
    REAL*8, ALLOCATABLE, INTENT(OUT) :: q2c(:,:)
    REAL*8 :: refxi(3,27), nodes(3,8), xx(3), jac(3,3), detj
    INTEGER :: ee, ii, idx(27), nd
    nd = mm%level(1)%nvt + mm%level(1)%net + mm%level(1)%nat + mm%level(1)%nel
    ALLOCATE(q2c(3,nd))
    CALL CHI_Q2_REFNODES(refxi)
    DO ee = 1, mm%level(1)%nel
      DO ii = 1, 8
        nodes(:,ii) = mm%level(1)%dcorvg(:,mm%level(1)%kvert(ii,ee))
      END DO
      CALL CHI_Q2_DOFMAP(ee, mm%level(1)%kvert, mm%level(1)%kedge, mm%level(1)%karea, &
        mm%level(1)%nvt, mm%level(1)%net, mm%level(1)%nat, idx)
      DO ii = 1, 27
        CALL CHI_Q1_MAP(nodes, refxi(:,ii), xx, jac, detj)
        q2c(:,idx(ii)) = xx
      END DO
    END DO
  END SUBROUTINE q2coords

  SUBROUTINE write_slab(fn, x0, x1)
    CHARACTER(*), INTENT(IN) :: fn
    REAL*8, INTENT(IN) :: x0, x1
    INTEGER, PARAMETER :: n = 2
    INTEGER :: iu, ix, iy, iz, ii
    iu = 91
    OPEN(iu, FILE=fn, STATUS='REPLACE')
    WRITE(iu,'(A)') 'test_chi_exchange slab'
    WRITE(iu,'(A)') 'Parametrisierung PARXC, PARYC, TMAXC'
    WRITE(iu,'(I0,1X,I0,A)') n*n*n, (n+1)**3, ' 1 8 12 6     NEL NVT NBCT NVE NEE NAE'
    WRITE(iu,'(A)') 'DCORVG'
    DO iz = 0, n
      DO iy = 0, n
        DO ix = 0, n
          WRITE(iu,'(3ES24.16)') x0 + (x1-x0)*ix/n, DBLE(iy)/n, DBLE(iz)/n
        END DO
      END DO
    END DO
    WRITE(iu,'(A)') 'KVERT'
    DO iz = 0, n-1
      DO iy = 0, n-1
        DO ix = 0, n-1
          WRITE(iu,'(8(I0,1X))') vid(n,ix,iy,iz), vid(n,ix+1,iy,iz), vid(n,ix+1,iy+1,iz), vid(n,ix,iy+1,iz), &
            vid(n,ix,iy,iz+1), vid(n,ix+1,iy,iz+1), vid(n,ix+1,iy+1,iz+1), vid(n,ix,iy+1,iz+1)
        END DO
      END DO
    END DO
    WRITE(iu,'(A)') 'KNPR'
    DO ii = 1, (n+1)**3
      WRITE(iu,'(I0)') 1
    END DO
    CLOSE(iu)
  END SUBROUTINE write_slab

  INTEGER FUNCTION vid(n, ix, iy, iz)
    INTEGER, INTENT(IN) :: n, ix, iy, iz
    vid = 1 + ix + (n+1)*(iy + (n+1)*iz)
  END FUNCTION vid

END PROGRAM test_chi_exchange
