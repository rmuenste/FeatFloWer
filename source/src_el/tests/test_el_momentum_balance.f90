!=========================================================================
! test_el_momentum_balance - helpers of the momentum-balance audit (EL_MOMENTUM_BALANCE;
! SimPar@ELMomentumBalanceAudit). Serial, on the synthetic structured unit-cube mesh.
!
!   1 masked sums: all rows, Dirichlet rows per component, scale
!   2 plane integrals at z = 0 and z = 1 for u = (k z, 0, a + b z + g x) with an element
!     coefficient rho_eps: area 1, momentum flux int rho_eps u (u.n) dS and viscous traction
!     int c (grad u) n dS against the closed forms (integrands of degree <= 2: exact)
!   3 an interior plane (z = 1/2) is found from both sides: the two outward normals give fluxes
!     and tractions of opposite sign (each side's sum equals the closed form)
!   4 a plane between element layers that no element face lies in is not picked up
!=========================================================================
PROGRAM TEST_EL_MOMENTUM_BALANCE

  USE MPI
  USE TYPES, ONLY: tMultiMesh
  USE EL_MOMENTUM_BALANCE, ONLY: EL_MOMBAL_SUMS, EL_MOMBAL_PLANE_LOCAL
  USE EL_TEST_HELPERS, ONLY: EL_TEST_INIT, EL_TEST_CHECK, EL_TEST_REPORT, &
                             EL_TEST_RELEASE_MESH, BUILD_STRUCTURED_MESH, &
                             BUILD_DOF_COORDINATES, nel, nvt, net, nat

  IMPLICIT NONE

  INTEGER :: rank, nproc, ierr, failures

  CALL MPI_Init(ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  CALL EL_TEST_INIT()

  failures = 0
  CALL test_sums(failures)
  CALL test_planes(failures)
  CALL EL_TEST_REPORT(failures, rank, nproc, 'EL momentum balance')
  CALL MPI_Finalize(ierr)
  IF (failures.GT.0) STOP 1

CONTAINS

  SUBROUTINE test_sums(nf)
    INTEGER, INTENT(INOUT) :: nf
    REAL*8 :: f(3,4), s_all(3), s_dir(3)
    LOGICAL :: mask(3,4)
    INTEGER :: i

    DO i = 1, 4
      f(:,i) = (/ DBLE(i), 10.0d0*DBLE(i), -DBLE(i*i) /)
    END DO
    mask = .FALSE.
    mask(:,2) = .TRUE.          ! a fully constrained row
    mask(3,4) = .TRUE.          ! only the z component of row 4
    CALL EL_MOMBAL_SUMS(f, mask, 0.5d0, s_all, s_dir)
    CALL EL_TEST_CHECK(ALL(ABS(s_all - 0.5d0*(/ 10.0d0, 100.0d0, -30.0d0 /)).LT.1.0d-14), &
      'sums: all rows, scaled', nf, rank)
    CALL EL_TEST_CHECK(ALL(ABS(s_dir - 0.5d0*(/ 2.0d0, 20.0d0, -4.0d0 - 16.0d0 /)).LT.1.0d-14), &
      'sums: Dirichlet rows per component', nf, rank)
  END SUBROUTINE test_sums

  SUBROUTINE test_planes(nf)
    INTEGER, INTENT(INOUT) :: nf
    TYPE(tMultiMesh) :: mesh
    REAL*8, ALLOCATABLE :: coor(:,:), u(:), v(:), w(:), coef(:)
    REAL*8, PARAMETER :: a = 0.3d0, b = -1.7d0, g = 0.45d0, k = 2.1d0, rho_eps = 0.62d0, &
                         c = 1.8d-2
    REAL*8 :: area, flux(3), trac(3), area2, flux2(3), trac2(3), z, wint, w2int
    INTEGER :: n, ndof, i

    n = 4
    CALL BUILD_STRUCTURED_MESH(mesh, n)
    ndof = nvt + net + nat + nel
    ALLOCATE(coor(3,ndof), u(ndof), v(ndof), w(ndof), coef(nel))
    CALL BUILD_DOF_COORDINATES(mesh, coor)
    DO i = 1, ndof
      u(i) = k*coor(3,i)
      v(i) = 0.0d0
      w(i) = a + b*coor(3,i) + g*coor(1,i)
    END DO
    coef = rho_eps

    ! bottom plane, outward normal -e_z: w(x) = a + g x on z = 0
    area = 0.0d0
    flux = 0.0d0
    trac = 0.0d0
    CALL EL_MOMBAL_PLANE_LOCAL(mesh, 1, coor, u, v, w, coef, c, 0.0d0, 1.0d-9, area, flux, trac)
    w2int = a*a + a*g + g*g/3.0d0
    CALL EL_TEST_CHECK(ABS(area - 1.0d0).LT.1.0d-13, 'bottom plane: area', nf, rank)
    CALL EL_TEST_CHECK(ABS(flux(3) + rho_eps*w2int).LT.1.0d-13 .AND. ABS(flux(1)).LT.1.0d-13 &
      .AND. ABS(flux(2)).LT.1.0d-13, 'bottom plane: momentum flux -rho_eps int w^2', nf, rank)
    CALL EL_TEST_CHECK(ABS(trac(3) + c*b).LT.1.0d-13 .AND. ABS(trac(1) + c*k).LT.1.0d-13 .AND. &
      ABS(trac(2)).LT.1.0d-13, 'bottom plane: traction -c du/dz', nf, rank)

    ! top plane, outward normal +e_z: w(x) = a + b + g x, u_x = k
    area = 0.0d0
    flux = 0.0d0
    trac = 0.0d0
    CALL EL_MOMBAL_PLANE_LOCAL(mesh, 1, coor, u, v, w, coef, c, 1.0d0, 1.0d-9, area, flux, trac)
    wint = (a + b) + 0.5d0*g
    w2int = (a + b)**2 + (a + b)*g + g*g/3.0d0
    CALL EL_TEST_CHECK(ABS(area - 1.0d0).LT.1.0d-13 .AND. ABS(flux(3) - rho_eps*w2int).LT.1.0d-13 &
      .AND. ABS(flux(1) - rho_eps*k*wint).LT.1.0d-13, 'top plane: area and momentum flux', nf, rank)
    CALL EL_TEST_CHECK(ABS(trac(3) - c*b).LT.1.0d-13 .AND. ABS(trac(1) - c*k).LT.1.0d-13, &
      'top plane: traction', nf, rank)

    ! interior plane z = 1/2: the layer below (outward +e_z) and the layer above (-e_z) both
    ! carry a face there; the sums of the two sides cancel
    z = 0.5d0
    area = 0.0d0
    flux = 0.0d0
    trac = 0.0d0
    CALL EL_MOMBAL_PLANE_LOCAL(mesh, 1, coor, u, v, w, coef, c, z, 1.0d-9, area, flux, trac)
    CALL EL_TEST_CHECK(ABS(area - 2.0d0).LT.1.0d-13 .AND. MAXVAL(ABS(flux)).LT.1.0d-13 .AND. &
      MAXVAL(ABS(trac)).LT.1.0d-13, 'interior plane: both sides found, opposite normals', nf, rank)

    ! z = 0.3 cuts the second layer: no face lies there
    area2 = 0.0d0
    flux2 = 0.0d0
    trac2 = 0.0d0
    CALL EL_MOMBAL_PLANE_LOCAL(mesh, 1, coor, u, v, w, coef, c, 0.3d0, 1.0d-9, area2, flux2, trac2)
    CALL EL_TEST_CHECK(area2.EQ.0.0d0 .AND. ALL(flux2.EQ.0.0d0) .AND. ALL(trac2.EQ.0.0d0), &
      'plane without element faces: nothing picked up', nf, rank)

    CALL EL_TEST_RELEASE_MESH(mesh)
    DEALLOCATE(coor, u, v, w, coef)
  END SUBROUTINE test_planes

END PROGRAM TEST_EL_MOMENTUM_BALANCE
