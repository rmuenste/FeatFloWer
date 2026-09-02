!=========================================================================
! CHI_SOLVER - per-submesh monolithic saddle-point solver, Layer M
! (design: chimera-integration-design.md v3, sections 3 and 7).
!
! Solves the (steady or theta-stepped) Navier-Stokes problem on ONE
! atmosphere submesh: Q2/P1 monolithic system, Picard linearization of
! convection (and of the Robin alpha-term), direct solve per iteration
! through the instance-based tSparseDirectSolver.  Boundary conditions:
!   - inner surface + z-faces: Dirichlet (values from the bc callback,
!     evaluated at the PROJECTED Q2 node coordinates),
!   - outer surface: Robin (paper eq. (5d), data from the robin callback)
!     or Dirichlet (diagnostic mode).
! In the all-Dirichlet mode the constant-pressure nullspace is removed
! by pinning the first pressure dof (gauge); with Robin active the
! traction datum fixes the pressure level and no gauge is applied.
!=========================================================================
MODULE CHI_SOLVER

  USE CHI_SUBMESH, ONLY: tChimeraSubmesh, CHI_SURF_INNER, CHI_SURF_OUTER, &
    CHI_SURF_ZLO, CHI_SURF_ZHI
  USE CHI_KERNELS, ONLY: CHI_BUILD_SADDLE_CSR, CHI_ASM_SADDLE, &
    CHI_ASM_ROBIN_TAB, CHI_ROBIN_POINTS, CHI_APPLY_DIRICHLET_ROW, &
    CHI_ASM_MASS_RHS, &
    chi_bc_velocity, chi_robin_data
  USE CHI_SPARSE_DIRECT, ONLY: tSparseDirectSolver, SD_INIT, SD_FACTORIZE, &
    SD_SOLVE, SD_FREE

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: tChiSubSolver
  PUBLIC :: CHI_SOLVER_INIT
  PUBLIC :: CHI_SOLVE_STEADY
  PUBLIC :: CHI_SOLVE_STEADY_TAB
  PUBLIC :: CHI_DIR_U, CHI_DIR_V, CHI_DIR_W, CHI_DIR_ALL
  PUBLIC :: CHI_SOLVER_ADVANCE
  PUBLIC :: CHI_SOLVER_RELEASE

  ! Component bits of the per-dof Dirichlet mask of the tabulated solve.
  INTEGER, PARAMETER :: CHI_DIR_U = 1, CHI_DIR_V = 2, CHI_DIR_W = 4
  INTEGER, PARAMETER :: CHI_DIR_ALL = 7

  TYPE tChiSubSolver
    LOGICAL :: initialized = .FALSE.
    INTEGER :: n = 0, ndof = 0, nel = 0
    INTEGER, ALLOCATABLE :: LdA(:), ColA(:)
    REAL*8,  ALLOCATABLE :: Avals(:), rhs(:), sol(:)
    REAL*8,  ALLOCATABLE :: u(:), v(:), w(:)     ! (ndof)
    REAL*8,  ALLOCATABLE :: uold(:), vold(:), wold(:)  ! previous time level
    REAL*8,  ALLOCATABLE :: p(:)                 ! (4*nel)
    TYPE(tSparseDirectSolver) :: sd
  END TYPE tChiSubSolver

CONTAINS

  !-----------------------------------------------------------------------
  ! Build the sparsity pattern (once per submesh) and prepare storage.
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_SOLVER_INIT(slv, sub, ok)
    TYPE(tChiSubSolver), INTENT(INOUT) :: slv
    TYPE(tChimeraSubmesh), INTENT(IN) :: sub
    LOGICAL, INTENT(OUT) :: ok

    INTEGER :: nvt, net, nat, nel

    ok = .FALSE.
    CALL CHI_SOLVER_RELEASE(slv)
    IF (.NOT. sub%initialized) RETURN

    nvt = sub%mesh%level(sub%nlmax)%nvt
    net = sub%mesh%level(sub%nlmax)%net
    nat = sub%mesh%level(sub%nlmax)%nat
    nel = sub%mesh%level(sub%nlmax)%nel

    slv%ndof = nvt + net + nat + nel
    slv%nel = nel
    CALL CHI_BUILD_SADDLE_CSR(nel, nvt, net, nat, &
      sub%mesh%level(sub%nlmax)%kvert, sub%mesh%level(sub%nlmax)%kedge, &
      sub%mesh%level(sub%nlmax)%karea, slv%n, slv%LdA, slv%ColA)

    ALLOCATE(slv%Avals(slv%LdA(slv%n+1)-1))
    ALLOCATE(slv%rhs(slv%n), slv%sol(slv%n))
    ALLOCATE(slv%u(slv%ndof), slv%v(slv%ndof), slv%w(slv%ndof))
    ALLOCATE(slv%uold(slv%ndof), slv%vold(slv%ndof), slv%wold(slv%ndof))
    slv%uold = 0d0
    slv%vold = 0d0
    slv%wold = 0d0
    ALLOCATE(slv%p(4*nel))
    slv%u = 0d0
    slv%v = 0d0
    slv%w = 0d0
    slv%p = 0d0

    CALL SD_INIT(slv%sd, slv%n, slv%LdA, slv%ColA, ok)
    IF (.NOT. ok) RETURN

    slv%initialized = .TRUE.
  END SUBROUTINE CHI_SOLVER_INIT

  !-----------------------------------------------------------------------
  ! Steady Navier-Stokes Picard solve, callback form (Phase-2 interface,
  ! used by the analytic tests): tabulates the Dirichlet values at the
  ! projected boundary nodes (all three components on inner + z faces,
  ! plus the outer surface in the diagnostic all-Dirichlet mode) and the
  ! Robin data at the outer quadrature points, then runs the tabulated
  ! solve.
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_SOLVE_STEADY(slv, sub, rho, mu, alpha, outer_dirichlet, &
                              bc_vel, robin_h, npicard, resid, ok)
    TYPE(tChiSubSolver), INTENT(INOUT) :: slv
    TYPE(tChimeraSubmesh), INTENT(IN) :: sub
    REAL*8, INTENT(IN) :: rho, mu, alpha
    LOGICAL, INTENT(IN) :: outer_dirichlet
    PROCEDURE(chi_bc_velocity) :: bc_vel
    PROCEDURE(chi_robin_data) :: robin_h
    INTEGER, INTENT(IN) :: npicard
    REAL*8, INTENT(OUT) :: resid
    LOGICAL, INTENT(OUT) :: ok

    INTEGER, ALLOCATABLE :: dirmask(:)
    REAL*8,  ALLOCATABLE :: ubc(:,:), xq(:,:,:), nq(:,:,:), hq(:,:,:)
    INTEGER :: i, dmask, nouter, ifc, q
    LOGICAL :: isdir

    ok = .FALSE.
    resid = HUGE(1d0)
    IF (.NOT. slv%initialized) RETURN

    ALLOCATE(dirmask(slv%ndof), ubc(3,slv%ndof))
    dirmask = 0
    ubc = 0d0
    DO i = 1, slv%ndof
      dmask = sub%dofmask(i)
      IF (dmask .EQ. 0) CYCLE
      isdir = (IAND(dmask, CHI_SURF_INNER) .NE. 0) .OR. &
              (IAND(dmask, CHI_SURF_ZLO)   .NE. 0) .OR. &
              (IAND(dmask, CHI_SURF_ZHI)   .NE. 0)
      IF (outer_dirichlet .AND. IAND(dmask, CHI_SURF_OUTER) .NE. 0) &
        isdir = .TRUE.
      IF (.NOT. isdir) CYCLE
      dirmask(i) = CHI_DIR_ALL
      CALL bc_vel(sub%q2coor(:,i), ubc(:,i))
    END DO

    nouter = SIZE(sub%outerFaces,2)
    ALLOCATE(xq(3,9,MAX(nouter,1)), nq(3,9,MAX(nouter,1)), hq(3,9,MAX(nouter,1)))
    hq = 0d0
    IF (.NOT. outer_dirichlet) THEN
      CALL CHI_ROBIN_POINTS(sub%outerFaces, nouter, &
        sub%mesh%level(sub%nlmax)%kvert, sub%mesh%level(sub%nlmax)%dcorvg, xq, nq)
      DO ifc = 1, nouter
        DO q = 1, 9
          CALL robin_h(xq(:,q,ifc), nq(:,q,ifc), hq(:,q,ifc))
        END DO
      END DO
    END IF

    CALL CHI_SOLVE_STEADY_TAB(slv, sub, rho, mu, 0d0, alpha, outer_dirichlet, &
                              dirmask, ubc, hq, npicard, resid, ok)
    DEALLOCATE(dirmask, ubc, xq, nq, hq)
  END SUBROUTINE CHI_SOLVE_STEADY

  !-----------------------------------------------------------------------
  ! Steady Navier-Stokes Picard solve with TABULATED boundary data (the
  ! Phase-3 coupling path): dirmask(ndof) holds CHI_DIR_* component bits,
  ! ubc(3,ndof) the Dirichlet values, hq(3,9,nouter) the Robin data at
  ! the CHI_ROBIN_POINTS quadrature points of sub%outerFaces (ignored in
  ! the all-Dirichlet mode, where the pressure gauge is applied instead).
  ! dtinv = 0: steady problem; dtinv > 0: one backward-Euler step of
  ! size 1/dtinv from the stored previous level slv%uold (advance it with
  ! CHI_SOLVER_ADVANCE after the step).  Warm-started from slv%u/v/w;
  ! resid is the max-norm velocity update of the LAST Picard iteration.
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_SOLVE_STEADY_TAB(slv, sub, rho, mu, dtinv, alpha, outer_dirichlet, &
                                  dirmask, ubc, hq, npicard, resid, ok)
    TYPE(tChiSubSolver), INTENT(INOUT) :: slv
    TYPE(tChimeraSubmesh), INTENT(IN) :: sub
    REAL*8, INTENT(IN) :: rho, mu, dtinv, alpha
    LOGICAL, INTENT(IN) :: outer_dirichlet
    INTEGER, INTENT(IN) :: dirmask(*)
    REAL*8,  INTENT(IN) :: ubc(3,*), hq(3,9,*)
    INTEGER, INTENT(IN) :: npicard
    REAL*8, INTENT(OUT) :: resid
    LOGICAL, INTENT(OUT) :: ok

    INTEGER :: nvt, net, nat, nel, ndof, ip, i, a, r
    REAL*8 :: umesh(3)
    LOGICAL :: sdok

    ok = .FALSE.
    resid = HUGE(1d0)
    IF (.NOT. slv%initialized) RETURN

    nvt = sub%mesh%level(sub%nlmax)%nvt
    net = sub%mesh%level(sub%nlmax)%net
    nat = sub%mesh%level(sub%nlmax)%nat
    nel = sub%mesh%level(sub%nlmax)%nel
    ndof = slv%ndof
    umesh = 0d0

    DO ip = 1, npicard
      CALL CHI_ASM_SADDLE(nel, nvt, net, nat, &
        sub%mesh%level(sub%nlmax)%kvert, sub%mesh%level(sub%nlmax)%kedge, &
        sub%mesh%level(sub%nlmax)%karea, sub%mesh%level(sub%nlmax)%dcorvg, &
        slv%n, slv%LdA, slv%ColA, slv%Avals, &
        rho, mu, dtinv, slv%u, slv%v, slv%w, umesh)
      slv%rhs = 0d0
      IF (dtinv .GT. 0d0) THEN
        CALL CHI_ASM_MASS_RHS(nel, nvt, net, nat, &
          sub%mesh%level(sub%nlmax)%kvert, sub%mesh%level(sub%nlmax)%kedge, &
          sub%mesh%level(sub%nlmax)%karea, sub%mesh%level(sub%nlmax)%dcorvg, &
          rho*dtinv, slv%uold, slv%vold, slv%wold, slv%rhs)
      END IF

      IF (.NOT. outer_dirichlet) THEN
        CALL CHI_ASM_ROBIN_TAB(sub%outerFaces, SIZE(sub%outerFaces,2), &
          nvt, net, nat, nel, &
          sub%mesh%level(sub%nlmax)%kvert, sub%mesh%level(sub%nlmax)%kedge, &
          sub%mesh%level(sub%nlmax)%karea, sub%mesh%level(sub%nlmax)%dcorvg, &
          slv%n, slv%LdA, slv%ColA, slv%Avals, slv%rhs, &
          alpha, slv%u, slv%v, slv%w, hq)
      END IF

      ! Dirichlet rows (after Robin, so rim dofs shared between the
      ! outer surface and a Dirichlet surface end up Dirichlet).
      DO i = 1, ndof
        IF (dirmask(i) .EQ. 0) CYCLE
        DO a = 1, 3
          IF (IAND(dirmask(i), 2**(a-1)) .EQ. 0) CYCLE
          r = (a-1)*ndof + i
          CALL CHI_APPLY_DIRICHLET_ROW(r, ubc(a,i), slv%LdA, slv%ColA, &
                                       slv%Avals, slv%rhs)
        END DO
      END DO

      ! Pressure gauge only in the all-Dirichlet mode.
      IF (outer_dirichlet) THEN
        CALL CHI_APPLY_DIRICHLET_ROW(3*ndof + 1, 0d0, slv%LdA, slv%ColA, &
                                     slv%Avals, slv%rhs)
      END IF

      CALL SD_FACTORIZE(slv%sd, slv%Avals, sdok)
      IF (.NOT. sdok) THEN
        WRITE(*,*) 'CHI_SOLVE_STEADY_TAB: factorization failed, iteration', ip
        RETURN
      END IF
      CALL SD_SOLVE(slv%sd, slv%sol, slv%rhs, sdok)
      IF (.NOT. sdok) THEN
        WRITE(*,*) 'CHI_SOLVE_STEADY_TAB: solve failed, iteration', ip
        RETURN
      END IF

      resid = 0d0
      DO i = 1, ndof
        resid = MAX(resid, ABS(slv%sol(i)        - slv%u(i)), &
                           ABS(slv%sol(ndof+i)   - slv%v(i)), &
                           ABS(slv%sol(2*ndof+i) - slv%w(i)))
        slv%u(i) = slv%sol(i)
        slv%v(i) = slv%sol(ndof+i)
        slv%w(i) = slv%sol(2*ndof+i)
      END DO
      slv%p(1:4*nel) = slv%sol(3*ndof+1:3*ndof+4*nel)
      IF (resid .LT. 1d-11) EXIT   ! Picard converged to solver precision
    END DO

    ok = .TRUE.
  END SUBROUTINE CHI_SOLVE_STEADY_TAB

  !-----------------------------------------------------------------------
  ! Accept the current solution as the previous time level.
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_SOLVER_ADVANCE(slv)
    TYPE(tChiSubSolver), INTENT(INOUT) :: slv
    IF (.NOT. slv%initialized) RETURN
    slv%uold = slv%u
    slv%vold = slv%v
    slv%wold = slv%w
  END SUBROUTINE CHI_SOLVER_ADVANCE

  SUBROUTINE CHI_SOLVER_RELEASE(slv)
    TYPE(tChiSubSolver), INTENT(INOUT) :: slv
    CALL SD_FREE(slv%sd)
    IF (ALLOCATED(slv%uold)) DEALLOCATE(slv%uold, slv%vold, slv%wold)
    IF (ALLOCATED(slv%LdA))  DEALLOCATE(slv%LdA)
    IF (ALLOCATED(slv%ColA)) DEALLOCATE(slv%ColA)
    IF (ALLOCATED(slv%Avals)) DEALLOCATE(slv%Avals)
    IF (ALLOCATED(slv%rhs))  DEALLOCATE(slv%rhs)
    IF (ALLOCATED(slv%sol))  DEALLOCATE(slv%sol)
    IF (ALLOCATED(slv%u))    DEALLOCATE(slv%u)
    IF (ALLOCATED(slv%v))    DEALLOCATE(slv%v)
    IF (ALLOCATED(slv%w))    DEALLOCATE(slv%w)
    IF (ALLOCATED(slv%p))    DEALLOCATE(slv%p)
    slv%initialized = .FALSE.
    slv%n = 0
    slv%ndof = 0
    slv%nel = 0
  END SUBROUTINE CHI_SOLVER_RELEASE

END MODULE CHI_SOLVER
