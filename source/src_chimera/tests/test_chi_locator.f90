!=========================================================================
! test_chi_locator - Phase 1 unit test of CHI_LOCATOR
! (chimera-integration-design.md v3, section 9).
!
! Builds a structured, deterministically distorted 4x4x4 hex mesh of the
! unit cube (FEAT kvert vertex ordering) and compares CHI_LOCATE against
! a brute-force scan over all elements for random interior points; also
! checks miss cases (points outside the mesh) and that the returned
! (element, xi) pair reproduces the query point through the forward map.
!
! Serial; prints PASS / STOP 1.
!=========================================================================
PROGRAM test_chi_locator

  USE CHI_GEOMETRY, ONLY: CHI_Q1_MAP, CHI_INVERSE_MAP
  USE CHI_LOCATOR, ONLY: tChimeraLocator, CHI_LOCATOR_BUILD, CHI_LOCATE, &
    CHI_LOCATOR_RELEASE

  IMPLICIT NONE

  INTEGER, PARAMETER :: nc = 4                 ! cells per direction
  INTEGER, PARAMETER :: nvt = (nc+1)**3
  INTEGER, PARAMETER :: nel = nc**3

  REAL*8 :: dcorvg(3,nvt)
  INTEGER :: kvert(8,nel)
  TYPE(tChimeraLocator) :: loc

  REAL*8 :: p(3), xi(3), xib(3), point(3), jac(3,3), detj, nodes(3,8)
  REAL*8 :: h
  INTEGER :: i, j, k, iv, iel, ielb, itest, nfail
  LOGICAL :: found, foundb, conv
  INTEGER*8 :: seed

  nfail = 0
  seed = 424242_8
  h = 1d0/DBLE(nc)

  !--- build the mesh ----------------------------------------------------
  DO k = 0, nc
    DO j = 0, nc
      DO i = 0, nc
        iv = vid(i,j,k)
        dcorvg(:,iv) = (/ DBLE(i)*h, DBLE(j)*h, DBLE(k)*h /)
        ! distort interior vertices deterministically (keeps elements
        ! convex at 0.15h amplitude)
        IF (i.GT.0 .AND. i.LT.nc .AND. j.GT.0 .AND. j.LT.nc .AND. &
            k.GT.0 .AND. k.LT.nc) THEN
          dcorvg(1,iv) = dcorvg(1,iv) + 0.15d0*h*(2d0*lcg(seed)-1d0)
          dcorvg(2,iv) = dcorvg(2,iv) + 0.15d0*h*(2d0*lcg(seed)-1d0)
          dcorvg(3,iv) = dcorvg(3,iv) + 0.15d0*h*(2d0*lcg(seed)-1d0)
        END IF
      END DO
    END DO
  END DO

  iel = 0
  DO k = 0, nc-1
    DO j = 0, nc-1
      DO i = 0, nc-1
        iel = iel + 1
        ! FEAT vertex ordering: bottom (--,+-,++,-+), then top
        kvert(1,iel) = vid(i,  j,  k)
        kvert(2,iel) = vid(i+1,j,  k)
        kvert(3,iel) = vid(i+1,j+1,k)
        kvert(4,iel) = vid(i,  j+1,k)
        kvert(5,iel) = vid(i,  j,  k+1)
        kvert(6,iel) = vid(i+1,j,  k+1)
        kvert(7,iel) = vid(i+1,j+1,k+1)
        kvert(8,iel) = vid(i,  j+1,k+1)
      END DO
    END DO
  END DO

  CALL CHI_LOCATOR_BUILD(loc, dcorvg, kvert, nel, nvt)
  IF (.NOT. loc%initialized) THEN
    WRITE(*,*) 'FAIL: locator build'
    STOP 1
  END IF

  !--- random interior points vs brute force -----------------------------
  DO itest = 1, 200
    p = (/ 0.02d0 + 0.96d0*lcg(seed), 0.02d0 + 0.96d0*lcg(seed), &
           0.02d0 + 0.96d0*lcg(seed) /)
    CALL CHI_LOCATE(loc, dcorvg, kvert, p, iel, xi, found)
    CALL brute_force(p, ielb, xib, foundb)
    IF (found .NEQV. foundb) THEN
      WRITE(*,*) 'FAIL: found mismatch at', p, found, foundb
      nfail = nfail + 1
      CYCLE
    END IF
    IF (.NOT. found) THEN
      WRITE(*,*) 'FAIL: interior point not located', p
      nfail = nfail + 1
      CYCLE
    END IF
    ! the located (iel, xi) must reproduce p through the forward map
    DO iv = 1, 8
      nodes(:,iv) = dcorvg(:,kvert(iv,iel))
    END DO
    CALL CHI_Q1_MAP(nodes, xi, point, jac, detj)
    IF (MAXVAL(ABS(point - p)) .GT. 1d-10) THEN
      WRITE(*,*) 'FAIL: located pair does not reproduce point', p
      nfail = nfail + 1
    END IF
  END DO

  !--- miss cases --------------------------------------------------------
  DO itest = 1, 20
    p = (/ 1.2d0 + lcg(seed), lcg(seed), lcg(seed) /)
    CALL CHI_LOCATE(loc, dcorvg, kvert, p, iel, xi, found)
    IF (found) THEN
      WRITE(*,*) 'FAIL: outside point located', p
      nfail = nfail + 1
    END IF
  END DO

  CALL CHI_LOCATOR_RELEASE(loc)
  CALL CHI_LOCATOR_RELEASE(loc)   ! idempotent

  IF (nfail .GT. 0) THEN
    WRITE(*,*) 'test_chi_locator: ', nfail, ' failure(s)'
    STOP 1
  END IF
  WRITE(*,*) 'test_chi_locator: PASS'

CONTAINS

  INTEGER FUNCTION vid(i, j, k)
    INTEGER, INTENT(IN) :: i, j, k
    vid = 1 + i + (nc+1)*(j + (nc+1)*k)
  END FUNCTION vid

  FUNCTION lcg(s) RESULT(r)
    INTEGER*8, INTENT(INOUT) :: s
    REAL*8 :: r
    ! Lehmer/minstd: product stays < 2^47, no signed-overflow UB
    s = MOD(s*48271_8, 2147483647_8)
    r = DBLE(s)/2147483647d0
  END FUNCTION lcg

  SUBROUTINE brute_force(pt, ie, x, fnd)
    REAL*8, INTENT(IN) :: pt(3)
    INTEGER, INTENT(OUT) :: ie
    REAL*8, INTENT(OUT) :: x(3)
    LOGICAL, INTENT(OUT) :: fnd
    REAL*8 :: nds(3,8)
    INTEGER :: je, jv
    LOGICAL :: cv
    ie = 0
    x = 0d0
    fnd = .FALSE.
    DO je = 1, nel
      DO jv = 1, 8
        nds(:,jv) = dcorvg(:,kvert(jv,je))
      END DO
      CALL CHI_INVERSE_MAP(nds, pt, x, cv)
      IF (cv) THEN
        ie = je
        fnd = .TRUE.
        RETURN
      END IF
    END DO
  END SUBROUTINE brute_force

END PROGRAM test_chi_locator
