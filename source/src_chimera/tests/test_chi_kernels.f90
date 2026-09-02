!=========================================================================
! test_chi_kernels - Phase 2 unit test of CHI_KERNELS
! (chimera-integration-design.md v3, section 9).
!
! Single-element (unit cube) algebraic checks, ff_chimera only:
! 1. Face rule/geometry: each face of the unit cube has area 1 and the
!    expected outward normal (pins the face tables to the FEAT karea
!    convention).
! 2. Mass consistency: with dtinv=1, rho=1, mu=0, u_k=0 the velocity
!    block is the Q2 mass matrix; its total sum is the element volume.
! 3. Deformation-form patch test: A * (linear velocity field) has zero
!    residual in the interior (center-bubble) rows.
! 4. Pressure-gradient consistency: G * (constant pressure) is zero in
!    the interior rows.
! 5. Continuity consistency: D * (divergence-free linear field) = 0
!    exactly; D * (x,0,0) gives element volume in the psi_1 row and 0 in
!    the (symmetric) psi_2..4 rows.
!
! Serial; prints PASS / STOP 1.
!=========================================================================
PROGRAM test_chi_kernels

  USE CHI_KERNELS, ONLY: CHI_BUILD_SADDLE_CSR, CHI_ASM_SADDLE, &
    CHI_FACE_RULE, CHI_FACE_GEOM

  IMPLICIT NONE

  ! one unit-cube element, FEAT conventions; identity edge/face numbering
  INTEGER, PARAMETER :: nvt = 8, net = 12, nat = 6, nel = 1
  INTEGER, PARAMETER :: ndof = nvt + net + nat + nel   ! 27

  REAL*8 :: dcorvg(3,8), nodes(3,8)
  INTEGER :: kvert(8,1), kedge(12,1), karea(6,1)
  INTEGER :: n, i, j, f, q, r, a, nfail
  INTEGER, ALLOCATABLE :: LdA(:), ColA(:)
  REAL*8, ALLOCATABLE :: Avals(:), x(:), y(:)
  REAL*8 :: uk(ndof), vk(ndof), wk(ndof), umesh(3)
  REAL*8 :: xiq(3,9), wq(9), xf(3), nrm(3), sj, area, nexp(3,6)
  REAL*8 :: total, q2c(3,ndof)

  REAL*8, PARAMETER :: signs(3,8) = RESHAPE((/ &
    -1d0,-1d0,-1d0,  1d0,-1d0,-1d0,  1d0, 1d0,-1d0, -1d0, 1d0,-1d0, &
    -1d0,-1d0, 1d0,  1d0,-1d0, 1d0,  1d0, 1d0, 1d0, -1d0, 1d0, 1d0 /), (/3,8/))
  ! FEAT edge -> endpoints, face handled via the kernel tables
  INTEGER, PARAMETER :: edvert(2,12) = RESHAPE((/ &
    1,2, 2,3, 3,4, 4,1, 1,5, 2,6, 3,7, 4,8, 5,6, 6,7, 7,8, 8,5 /), (/2,12/))
  INTEGER, PARAMETER :: fcvert(4,6) = RESHAPE((/ &
    1,2,3,4, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8, 5,6,7,8 /), (/4,6/))

  nfail = 0

  DO i = 1, 8
    dcorvg(:,i) = 0.5d0*(signs(:,i) + 1d0)   ! unit cube [0,1]^3
    nodes(:,i) = dcorvg(:,i)
    kvert(i,1) = i
  END DO
  DO i = 1, 12
    kedge(i,1) = i
  END DO
  DO i = 1, 6
    karea(i,1) = i
  END DO

  ! Q2 node coordinates of the element (vertices, edge mids, face mids,
  ! center) for evaluating test fields at the dofs
  DO i = 1, 8
    q2c(:,i) = dcorvg(:,i)
  END DO
  DO i = 1, 12
    q2c(:,8+i) = 0.5d0*(dcorvg(:,edvert(1,i)) + dcorvg(:,edvert(2,i)))
  END DO
  DO i = 1, 6
    q2c(:,20+i) = 0.25d0*(dcorvg(:,fcvert(1,i)) + dcorvg(:,fcvert(2,i)) + &
                          dcorvg(:,fcvert(3,i)) + dcorvg(:,fcvert(4,i)))
  END DO
  q2c(:,27) = (/ 0.5d0, 0.5d0, 0.5d0 /)

  !--- 1. face rule / geometry -------------------------------------------
  nexp = RESHAPE((/ 0d0,0d0,-1d0,  0d0,-1d0,0d0,  1d0,0d0,0d0, &
                    0d0,1d0,0d0,  -1d0,0d0,0d0,  0d0,0d0,1d0 /), (/3,6/))
  DO f = 1, 6
    CALL CHI_FACE_RULE(f, xiq, wq)
    area = 0d0
    DO q = 1, 9
      CALL CHI_FACE_GEOM(nodes, f, xiq(:,q), xf, nrm, sj)
      area = area + wq(q)*sj
      IF (MAXVAL(ABS(nrm - nexp(:,f))) .GT. 1d-13) THEN
        WRITE(*,*) 'FAIL: face normal, face', f, nrm
        nfail = nfail + 1
        EXIT
      END IF
    END DO
    IF (ABS(area - 1d0) .GT. 1d-13) THEN
      WRITE(*,*) 'FAIL: face area, face', f, area
      nfail = nfail + 1
    END IF
  END DO

  !--- assemble ----------------------------------------------------------
  CALL CHI_BUILD_SADDLE_CSR(nel, nvt, net, nat, kvert, kedge, karea, &
                            n, LdA, ColA)
  ALLOCATE(Avals(LdA(n+1)-1), x(n), y(n))
  uk = 0d0
  vk = 0d0
  wk = 0d0
  umesh = 0d0

  !--- 2. mass consistency ----------------------------------------------
  CALL CHI_ASM_SADDLE(nel, nvt, net, nat, kvert, kedge, karea, dcorvg, &
                      n, LdA, ColA, Avals, 1d0, 0d0, 1d0, uk, vk, wk, umesh)
  total = 0d0
  DO r = 1, ndof                      ! u-u block only
    DO j = LdA(r), LdA(r+1)-1
      IF (ColA(j) .LE. ndof) total = total + Avals(j)
    END DO
  END DO
  IF (ABS(total - 1d0) .GT. 1d-12) THEN
    WRITE(*,*) 'FAIL: mass-matrix sum (volume) =', total
    nfail = nfail + 1
  END IF

  !--- 3./4./5. consistency of viscous, G and D blocks -------------------
  CALL CHI_ASM_SADDLE(nel, nvt, net, nat, kvert, kedge, karea, dcorvg, &
                      n, LdA, ColA, Avals, 1d0, 1d0, 0d0, uk, vk, wk, umesh)

  ! (3) linear velocity u=(x+2y, 3z-y, x+z), p=0: interior residual rows
  ! (the center bubble, dof 27) must vanish.
  x = 0d0
  DO i = 1, ndof
    x(i)        = q2c(1,i) + 2d0*q2c(2,i)
    x(ndof+i)   = 3d0*q2c(3,i) - q2c(2,i)
    x(2*ndof+i) = q2c(1,i) + q2c(3,i)
  END DO
  CALL csr_matvec(n, LdA, ColA, Avals, x, y)
  DO a = 1, 3
    r = (a-1)*ndof + 27
    IF (ABS(y(r)) .GT. 1d-12) THEN
      WRITE(*,*) 'FAIL: viscous patch test, component', a, y(r)
      nfail = nfail + 1
    END IF
  END DO
  ! divergence of that field: 1 - 1 + 1 = 1 -> continuity psi_1 row = vol
  IF (ABS(y(3*ndof+1) - 1d0) .GT. 1d-12) THEN
    WRITE(*,*) 'FAIL: continuity psi_1 row for div=1 field:', y(3*ndof+1)
    nfail = nfail + 1
  END IF
  DO i = 2, 4
    IF (ABS(y(3*ndof+i)) .GT. 1d-12) THEN
      WRITE(*,*) 'FAIL: continuity psi_', i, ' row not zero:', y(3*ndof+i)
      nfail = nfail + 1
    END IF
  END DO

  ! (5) divergence-free linear field u=(y,z,x): all continuity rows 0
  DO i = 1, ndof
    x(i)        = q2c(2,i)
    x(ndof+i)   = q2c(3,i)
    x(2*ndof+i) = q2c(1,i)
  END DO
  x(3*ndof+1:n) = 0d0
  CALL csr_matvec(n, LdA, ColA, Avals, x, y)
  DO i = 1, 4
    IF (ABS(y(3*ndof+i)) .GT. 1d-12) THEN
      WRITE(*,*) 'FAIL: continuity row for div-free field:', i, y(3*ndof+i)
      nfail = nfail + 1
    END IF
  END DO

  ! (4) constant pressure p=1 (dofs (1,0,0,0)), u=0: interior momentum
  ! rows must vanish (G consistency).
  x = 0d0
  x(3*ndof+1) = 1d0
  CALL csr_matvec(n, LdA, ColA, Avals, x, y)
  DO a = 1, 3
    r = (a-1)*ndof + 27
    IF (ABS(y(r)) .GT. 1d-12) THEN
      WRITE(*,*) 'FAIL: pressure-gradient consistency, component', a, y(r)
      nfail = nfail + 1
    END IF
  END DO

  IF (nfail .GT. 0) THEN
    WRITE(*,*) 'test_chi_kernels: ', nfail, ' failure(s)'
    STOP 1
  END IF
  WRITE(*,*) 'test_chi_kernels: PASS'

CONTAINS

  SUBROUTINE csr_matvec(nn, Ld, Co, Av, xx, yy)
    INTEGER, INTENT(IN) :: nn, Ld(*), Co(*)
    REAL*8, INTENT(IN) :: Av(*), xx(*)
    REAL*8, INTENT(OUT) :: yy(*)
    INTEGER :: ii, jj
    DO ii = 1, nn
      yy(ii) = 0d0
      DO jj = Ld(ii), Ld(ii+1)-1
        yy(ii) = yy(ii) + Av(jj)*xx(Co(jj))
      END DO
    END DO
  END SUBROUTINE csr_matvec

END PROGRAM test_chi_kernels
