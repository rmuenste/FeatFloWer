!=========================================================================
! QuadSc_el_mombal.f90
!
! Momentum-balance audit (SimPar@ELMomentumBalanceAudit): fluid-side adapter. Included in
! Transport_Q2P1 (host association to QuadSc, LinSc, the matrices, mg_mesh, mgDensity, ...).
! EL_MOMENTUM_BALANCE holds the accumulators, the meaning and sign of every term, the closures and
! the line format; this file extracts the terms from the solver arrays.
!
! Four calls, one line each, every one a no-op unless the deck key is set, and on the master,
! which holds no fields:
!   EL_MomBal_Begin    fluid core, predictor, before the old-velocity defect (Matdef idef = 1):
!                      Dirichlet mask, u^n, p^n, sum M u^n
!   EL_MomBal_Rhs      fluid core, the complete momentum right-hand side before its Dirichlet rows
!                      are zeroed: split into M u^n, B p^n, -fp(p^n), c fv(u^n), the particle source,
!                      the hrz partner and the remainder
!   EL_MomBal_Defect   fluid core, every nonlinear iteration after the defect of the iterate is
!                      assembled and before its Dirichlet rows are zeroed (the last call holds the
!                      converged predictor u~): defect, M (u~ - u^n), -B u~, -D u~, -K u~
!   EL_MomBal_End      Euler-Lagrange transport, right after the fluid core (projection and the
!                      corrector passes done): u^{n+1}, p^{n+1}, q = p^{n+1} - p^n, the boundary
!                      planes, then the collective write over the workers
! Read-only: the solver arrays are operands only; products are formed in local scratch arrays with
! the solver's own kernels (B_Mul_U, LAX17, EL_MM_PHI_OPERATORS). No global state (ILEV, matrix
! pointers, communication) is touched. The scheme the split assumes (backward Euler, Newtonian,
! matrix M + dt (D + K + B)) is enforced at start-up (app_init).
!=========================================================================
SUBROUTINE EL_MomBal_Begin()

  USE EL_MOMENTUM_BALANCE, ONLY: EL_MOMBAL_ACTIVE, EL_MOMBAL_BEGIN_STEP, EL_MOMBAL_RECORD, &
                                 MB_MOM_START
  USE PP3D_MPI, ONLY: myid, master, MPI_COMM_SUBS

  LOGICAL, ALLOCATABLE :: mask(:,:)
  REAL*8, ALLOCATABLE :: f(:,:)
  INTEGER :: i, n, ierr

  IF (.NOT.EL_MOMBAL_ACTIVE()) RETURN
  IF (myid.EQ.master) RETURN

  n = QuadSc%ndof
  ALLOCATE(mask(3,n), f(3,n))
  ! the rows Boundary_QuadScalar_Def zeroes; slip rows (a normal projection) are not covered
  DO i = 1, n
    mask(1,i) = (QuadSc%knprU(NLMAX)%x(i).EQ.1)
    mask(2,i) = (QuadSc%knprV(NLMAX)%x(i).EQ.1)
    mask(3,i) = (QuadSc%knprW(NLMAX)%x(i).EQ.1)
    IF (ALLOCATED(FictKNPR)) THEN
      IF (FictKNPR(i).NE.0) mask(:,i) = .TRUE.
    END IF
    IF (ALLOCATED(MixerKNPR)) THEN
      IF (MixerKNPR(i).NE.0) mask(:,i) = .TRUE.
    END IF
    IF (ALLOCATED(myBoundary%bSlip)) THEN
      IF (myBoundary%bSlip(i) .AND. .NOT.(myBoundary%bWall(i) .OR. myBoundary%iInflow(i).GT.0)) THEN
        WRITE(*,'(A,I0)') 'Fatal: ELMomentumBalanceAudit does not cover slip boundaries; rank ', myid
        CALL MPI_Abort(MPI_COMM_SUBS, 915, ierr)
      END IF
    END IF
  END DO
  CALL EL_MOMBAL_BEGIN_STEP(mask, QuadSc%valU, QuadSc%valV, QuadSc%valW, LinSc%valP(NLMAX)%x)

  f(1,:) = mg_MlRhoMat(NLMAX)%a(1:n)*QuadSc%valU(1:n)
  f(2,:) = mg_MlRhoMat(NLMAX)%a(1:n)*QuadSc%valV(1:n)
  f(3,:) = mg_MlRhoMat(NLMAX)%a(1:n)*QuadSc%valW(1:n)
  CALL EL_MOMBAL_RECORD(MB_MOM_START, f, 1.0d0)
  DEALLOCATE(mask, f)

END SUBROUTINE EL_MomBal_Begin

!=========================================================================
SUBROUTINE EL_MomBal_Rhs()

  USE EL_MOMENTUM_BALANCE, ONLY: EL_MOMBAL_ACTIVE, EL_MOMBAL_RECORD, mb_u_old, mb_p_old, &
                                 MB_RHS_FULL, MB_PGRAD_OLD, MB_FEEDBACK, MB_PARTNER, MB_MMP_OLD, &
                                 MB_MMV_OLD, MB_RHS_OTHER
  USE EL_MM_OPERATOR, ONLY: EL_MM_PHI_OPERATORS
  USE EL_FORCES, ONLY: EL_FLUID_OPERATOR_VISCOSITY
  USE PP3D_MPI, ONLY: myid, master

  REAL*8, ALLOCATABLE :: f(:,:), other(:,:), fp(:,:), fv(:,:), bx(:), by(:), bz(:)
  REAL*8 :: dt, c
  INTEGER :: n

  IF (.NOT.EL_MOMBAL_ACTIVE()) RETURN
  IF (myid.EQ.master) RETURN

  n = QuadSc%ndof
  dt = TSTEP
  ALLOCATE(f(3,n), other(3,n), bx(n), by(n), bz(n))

  ! the right-hand side as assembled, and its part beyond the old-time mass term
  f(1,:) = QuadSc%defU(1:n)
  f(2,:) = QuadSc%defV(1:n)
  f(3,:) = QuadSc%defW(1:n)
  CALL EL_MOMBAL_RECORD(MB_RHS_FULL, f, 1.0d0/dt)
  other(1,:) = (QuadSc%defU(1:n) - mg_MlRhoMat(NLMAX)%a(1:n)*mb_u_old(1,:))/dt
  other(2,:) = (QuadSc%defV(1:n) - mg_MlRhoMat(NLMAX)%a(1:n)*mb_u_old(2,:))/dt
  other(3,:) = (QuadSc%defW(1:n) - mg_MlRhoMat(NLMAX)%a(1:n)*mb_u_old(3,:))/dt

  ! unweighted pressure term B p^n (AddPressureGradient)
  CALL B_Mul_U(mg_qlMat(NLMAX)%ColA, mg_qlMat(NLMAX)%LdA, mg_BXMat(NLMAX)%a, mg_BYMat(NLMAX)%a, &
               mg_BZMat(NLMAX)%a, mb_p_old, bx, by, bz, n, 1.0d0, 0.0d0)
  f(1,:) = bx
  f(2,:) = by
  f(3,:) = bz
  CALL EL_MOMBAL_RECORD(MB_PGRAD_OLD, f, 1.0d0)
  other = other - f

  ! explicit particle source and the hrz partner (EL_APPLY_FLUID_FEEDBACK_SOURCE,
  ! EL_APPLY_DRAG_SPLIT_PARTNER, under the guards of the fluid core)
  f = 0.0d0
  IF (el_apply_fluid_feedback .AND. ALLOCATED(el_field_data%fluid_feedback_source)) &
    f = el_field_data%fluid_feedback_source(:,1:n)
  CALL EL_MOMBAL_RECORD(MB_FEEDBACK, f, 1.0d0)
  other = other - f
  f = 0.0d0
  IF (el_apply_fluid_feedback .AND. el_drag_hrz .AND. ALLOCATED(el_field_data%drag_B_source)) THEN
    f(1,:) = el_field_data%drag_B_source(1:n)*mb_u_old(1,:)
    f(2,:) = el_field_data%drag_B_source(1:n)*mb_u_old(2,:)
    f(3,:) = el_field_data%drag_B_source(1:n)*mb_u_old(3,:)
  END IF
  CALL EL_MOMBAL_RECORD(MB_PARTNER, f, 1.0d0)
  other = other - f

  ! momentum model A: the solid-fraction corrections -fp(p^n) and c fv(u^n)
  ! (EL_MM_Add_Operator_Terms, same operands)
  f = 0.0d0
  IF (EL_MM_OPERATOR_FORM()) THEN
    ALLOCATE(fp(3,n), fv(3,n))
    CALL EL_MM_PHI_OPERATORS(mg_mesh, NLMAX, mg_mesh%level(NLMAX)%dcorvg, &
      el_field_data%epsilon_transport_q2, mb_p_old, mb_u_old(1,:), mb_u_old(2,:), mb_u_old(3,:), &
      .TRUE., fp, fv)
    c = EL_FLUID_OPERATOR_VISCOSITY(Properties%Density(1), &
                                    EL_CLOSURE_VISCOSITY(Properties%Viscosity(1)))
    CALL EL_MOMBAL_RECORD(MB_MMP_OLD, fp, -1.0d0)
    CALL EL_MOMBAL_RECORD(MB_MMV_OLD, fv, c)
    other = other - (c*fv - fp)
    DEALLOCATE(fp, fv)
  ELSE
    CALL EL_MOMBAL_RECORD(MB_MMP_OLD, f, 1.0d0)
    CALL EL_MOMBAL_RECORD(MB_MMV_OLD, f, 1.0d0)
  END IF
  CALL EL_MOMBAL_RECORD(MB_RHS_OTHER, other, 1.0d0)

  DEALLOCATE(f, other, bx, by, bz)

END SUBROUTINE EL_MomBal_Rhs

!=========================================================================
SUBROUTINE EL_MomBal_Defect()

  USE EL_MOMENTUM_BALANCE, ONLY: EL_MOMBAL_ACTIVE, EL_MOMBAL_RECORD, EL_MOMBAL_SET_PREDICTOR, &
                                 mb_u_old, MB_DEFECT, MB_INERTIA_PRED, MB_MOM_PRED, &
                                 MB_DRAGIMP_PRED, MB_VISC_PRED, MB_CONV_PRED
  USE EL_CONFIG, ONLY: el_drag_semi_implicit
  USE PP3D_MPI, ONLY: myid, master

  REAL*8, ALLOCATABLE :: f(:,:), ax(:)
  REAL*8 :: dt
  INTEGER :: n

  IF (.NOT.EL_MOMBAL_ACTIVE()) RETURN
  IF (myid.EQ.master) RETURN

  n = QuadSc%ndof
  dt = TSTEP
  ALLOCATE(f(3,n), ax(n))
  CALL EL_MOMBAL_SET_PREDICTOR(QuadSc%valU, QuadSc%valV, QuadSc%valW)

  f(1,:) = QuadSc%defU(1:n)
  f(2,:) = QuadSc%defV(1:n)
  f(3,:) = QuadSc%defW(1:n)
  CALL EL_MOMBAL_RECORD(MB_DEFECT, f, 1.0d0/dt)

  f(1,:) = mg_MlRhoMat(NLMAX)%a(1:n)*(QuadSc%valU(1:n) - mb_u_old(1,:))
  f(2,:) = mg_MlRhoMat(NLMAX)%a(1:n)*(QuadSc%valV(1:n) - mb_u_old(2,:))
  f(3,:) = mg_MlRhoMat(NLMAX)%a(1:n)*(QuadSc%valW(1:n) - mb_u_old(3,:))
  CALL EL_MOMBAL_RECORD(MB_INERTIA_PRED, f, 1.0d0/dt)
  f(1,:) = mg_MlRhoMat(NLMAX)%a(1:n)*QuadSc%valU(1:n)
  f(2,:) = mg_MlRhoMat(NLMAX)%a(1:n)*QuadSc%valV(1:n)
  f(3,:) = mg_MlRhoMat(NLMAX)%a(1:n)*QuadSc%valW(1:n)
  CALL EL_MOMBAL_RECORD(MB_MOM_PRED, f, 1.0d0)

  ! the drag diagonal of the matrix (Matdef idef = -1, same guards)
  f = 0.0d0
  IF (el_apply_fluid_feedback .AND. el_drag_semi_implicit .AND. &
      ALLOCATED(el_field_data%drag_B_source)) THEN
    f(1,:) = el_field_data%drag_B_source(1:n)*QuadSc%valU(1:n)
    f(2,:) = el_field_data%drag_B_source(1:n)*QuadSc%valV(1:n)
    f(3,:) = el_field_data%drag_B_source(1:n)*QuadSc%valW(1:n)
  END IF
  CALL EL_MOMBAL_RECORD(MB_DRAGIMP_PRED, f, -1.0d0)

  ! diffusion and convection matrices of this iterate
  CALL MATVEC(mg_DMat(NLMAX)%a, MB_VISC_PRED)
  CALL MATVEC(mg_KMat(NLMAX)%a, MB_CONV_PRED)

  DEALLOCATE(f, ax)

CONTAINS

  ! record -A u~ for a scalar matrix A acting on each component
  SUBROUTINE MATVEC(a, term)
    REAL*8, INTENT(IN) :: a(:)
    INTEGER, INTENT(IN) :: term
    CALL LAX17(a, mg_qMat(NLMAX)%ColA, mg_qMat(NLMAX)%LdA, n, QuadSc%valU, ax, 1.0d0, 0.0d0)
    f(1,:) = ax
    CALL LAX17(a, mg_qMat(NLMAX)%ColA, mg_qMat(NLMAX)%LdA, n, QuadSc%valV, ax, 1.0d0, 0.0d0)
    f(2,:) = ax
    CALL LAX17(a, mg_qMat(NLMAX)%ColA, mg_qMat(NLMAX)%LdA, n, QuadSc%valW, ax, 1.0d0, 0.0d0)
    f(3,:) = ax
    CALL EL_MOMBAL_RECORD(term, f, -1.0d0)
  END SUBROUTINE MATVEC

END SUBROUTINE EL_MomBal_Defect

!=========================================================================
SUBROUTINE EL_MomBal_End(mfile, istep)

  USE EL_MOMENTUM_BALANCE, ONLY: EL_MOMBAL_ACTIVE, EL_MOMBAL_RECORD, EL_MOMBAL_RECORD_PLANES, &
                                 EL_MOMBAL_PLANE_LOCAL, EL_MOMBAL_WRITE, mb_u_pred, mb_p_old, &
                                 MB_INERTIA_PROJ, MB_MOM_END, MB_DRAGIMP_END, MB_PGRAD_END, &
                                 MB_PGRAD_Q, MB_MMP_END, MB_MMP_Q, MB_NPLANE
  USE EL_MM_OPERATOR, ONLY: EL_MM_PHI_OPERATORS
  USE EL_FORCES, ONLY: EL_FLUID_OPERATOR_VISCOSITY
  USE EL_BED_AUDIT, ONLY: EL_PLANE_PRESSURE_LOCAL
  USE EL_GEOMETRY, ONLY: el_domain_box
  USE EL_CONFIG, ONLY: el_drag_semi_implicit
  USE PP3D_MPI, ONLY: myid, master

  INTEGER, INTENT(IN) :: mfile, istep

  REAL*8, ALLOCATABLE :: f(:,:), fp(:,:), fv(:,:), q(:), bx(:), by(:), bz(:)
  REAL*8 :: dt, c, zmin, zmax, tol, pl(MB_NPLANE), face_area
  INTEGER :: n

  IF (.NOT.EL_MOMBAL_ACTIVE()) RETURN
  IF (myid.EQ.master) RETURN

  n = QuadSc%ndof
  dt = TSTEP
  ALLOCATE(f(3,n), bx(n), by(n), bz(n), q(SIZE(mb_p_old)))

  f(1,:) = mg_MlRhoMat(NLMAX)%a(1:n)*(QuadSc%valU(1:n) - mb_u_pred(1,:))
  f(2,:) = mg_MlRhoMat(NLMAX)%a(1:n)*(QuadSc%valV(1:n) - mb_u_pred(2,:))
  f(3,:) = mg_MlRhoMat(NLMAX)%a(1:n)*(QuadSc%valW(1:n) - mb_u_pred(3,:))
  CALL EL_MOMBAL_RECORD(MB_INERTIA_PROJ, f, 1.0d0/dt)
  f(1,:) = mg_MlRhoMat(NLMAX)%a(1:n)*QuadSc%valU(1:n)
  f(2,:) = mg_MlRhoMat(NLMAX)%a(1:n)*QuadSc%valV(1:n)
  f(3,:) = mg_MlRhoMat(NLMAX)%a(1:n)*QuadSc%valW(1:n)
  CALL EL_MOMBAL_RECORD(MB_MOM_END, f, 1.0d0)

  f = 0.0d0
  IF (el_apply_fluid_feedback .AND. el_drag_semi_implicit .AND. &
      ALLOCATED(el_field_data%drag_B_source)) THEN
    f(1,:) = el_field_data%drag_B_source(1:n)*QuadSc%valU(1:n)
    f(2,:) = el_field_data%drag_B_source(1:n)*QuadSc%valV(1:n)
    f(3,:) = el_field_data%drag_B_source(1:n)*QuadSc%valW(1:n)
  END IF
  CALL EL_MOMBAL_RECORD(MB_DRAGIMP_END, f, -1.0d0)

  ! pressure terms of the whole step: p^{n+1} and the projection increment q = p^{n+1} - p^n
  q = LinSc%valP(NLMAX)%x - mb_p_old
  CALL PGRAD(LinSc%valP(NLMAX)%x, MB_PGRAD_END)
  CALL PGRAD(q, MB_PGRAD_Q)
  f = 0.0d0
  IF (EL_MM_OPERATOR_FORM()) THEN
    ALLOCATE(fp(3,n), fv(3,n))
    CALL EL_MM_PHI_OPERATORS(mg_mesh, NLMAX, mg_mesh%level(NLMAX)%dcorvg, &
      el_field_data%epsilon_transport_q2, LinSc%valP(NLMAX)%x, QuadSc%valU, QuadSc%valV, &
      QuadSc%valW, .FALSE., fp, fv)
    CALL EL_MOMBAL_RECORD(MB_MMP_END, fp, -1.0d0)
    CALL EL_MM_PHI_OPERATORS(mg_mesh, NLMAX, mg_mesh%level(NLMAX)%dcorvg, &
      el_field_data%epsilon_transport_q2, q, QuadSc%valU, QuadSc%valV, QuadSc%valW, &
      .FALSE., fp, fv)
    CALL EL_MOMBAL_RECORD(MB_MMP_Q, fp, -1.0d0)
    DEALLOCATE(fp, fv)
  ELSE
    CALL EL_MOMBAL_RECORD(MB_MMP_END, f, 1.0d0)
    CALL EL_MOMBAL_RECORD(MB_MMP_Q, f, 1.0d0)
  END IF

  ! inlet (bottom) and outlet (top) planes of the domain box: pressure as the bed audit,
  ! momentum flux with the convection coefficient, viscous traction with the operator viscosity
  zmin = el_domain_box(5)
  zmax = el_domain_box(6)
  tol = 1.0d-6*(zmax - zmin)
  c = EL_FLUID_OPERATOR_VISCOSITY(Properties%Density(1), &
                                  EL_CLOSURE_VISCOSITY(Properties%Viscosity(1)))
  pl = 0.0d0
  CALL EL_PLANE_PRESSURE_LOCAL(mg_mesh, NLMAX, LinSc%valP(NLMAX)%x, zmin, tol, pl(1), pl(2))
  CALL EL_PLANE_PRESSURE_LOCAL(mg_mesh, NLMAX, LinSc%valP(NLMAX)%x, zmax, tol, pl(3), pl(4))
  ! (the face areas of the flux kernel repeat pl(1) / pl(3) and are not kept)
  face_area = 0.0d0
  CALL EL_MOMBAL_PLANE_LOCAL(mg_mesh, NLMAX, mg_mesh%level(NLMAX)%dcorvg, QuadSc%valU, &
    QuadSc%valV, QuadSc%valW, mgDensity(NLMAX)%x, c, zmin, tol, face_area, pl(5:7), pl(11:13))
  face_area = 0.0d0
  CALL EL_MOMBAL_PLANE_LOCAL(mg_mesh, NLMAX, mg_mesh%level(NLMAX)%dcorvg, QuadSc%valU, &
    QuadSc%valV, QuadSc%valW, mgDensity(NLMAX)%x, c, zmax, tol, face_area, pl(8:10), pl(14:16))
  CALL EL_MOMBAL_RECORD_PLANES(pl)

  CALL EL_MOMBAL_WRITE(timens, istep, dt, mfile)
  DEALLOCATE(f, bx, by, bz, q)

CONTAINS

  ! record sum B x for a P1 vector x (B_Mul_U, the kernel of AddPressureGradient)
  SUBROUTINE PGRAD(x, term)
    REAL*8, INTENT(IN) :: x(:)
    INTEGER, INTENT(IN) :: term
    CALL B_Mul_U(mg_qlMat(NLMAX)%ColA, mg_qlMat(NLMAX)%LdA, mg_BXMat(NLMAX)%a, mg_BYMat(NLMAX)%a, &
                 mg_BZMat(NLMAX)%a, x, bx, by, bz, n, 1.0d0, 0.0d0)
    f(1,:) = bx
    f(2,:) = by
    f(3,:) = bz
    CALL EL_MOMBAL_RECORD(term, f, 1.0d0)
  END SUBROUTINE PGRAD

END SUBROUTINE EL_MomBal_End
