MODULE EL_VA_CONTINUITY

  ! Volume-averaged continuity for the unresolved Euler-Lagrange coupling.
  !
  ! The single-phase solver enforces div u = 0. In a particle-laden fluid the
  ! volume-averaged continuity equation is
  !
  !     d eps/dt + div(eps u) = 0        (eps = fluid volume fraction)
  !
  ! i.e.  div u = g,   g = -( d eps/dt + u . grad eps ) / eps.
  !
  ! Stage 'continuity' (SimPar@ELVolumeAveraged = continuity) keeps the
  ! single-phase momentum operator and imposes the source g through the
  ! pressure Poisson right-hand side: the discrete constraint is B^T u = G
  ! with G_q = int_K q g dV, and the projection solves C dp = -(1/k) B^T u~ + G/k,
  ! so G/k is added to defP after Matdef_General_LinScalar(...,1) and before
  ! it is copied to rhsP (see Transport_q2p1_UxyzP_fluid_core). Only the
  ! constant P1 test function receives the source (the three gradient entries
  ! vanish for the axis-aligned brick cells of every EL mesh).
  !
  ! Inputs are the element fields of EL_FIELDS on level NLMAX, fresh from the
  ! post-advance transfer pass: eps, d eps/dt = (eps^{n+1}-eps^n)/dt and the
  ! cell-averaged kernel gradient grad(alpha_p) (grad eps = -grad alpha_p).
  ! The element velocity is the mean of its eight vertex values.
  !
  ! Compatibility: in a closed or periodic box the pressure matrix has a
  ! constant nullspace, so the source must have zero volume integral. The
  ! continuous form div(eps u) has that property exactly (constant particle
  ! volume), the divided form g loses it at O(phi); the volume mean of g is
  ! therefore removed (one Allreduce over the worker communicator; the CFD
  ! master never calls this routine).

  USE MPI
  USE PP3D_MPI, ONLY: myid, master, showid, MPI_COMM_SUBS
  USE TYPES, ONLY: tMultiMesh
  USE EL_CONFIG, ONLY: el_volume_averaged, el_momentum_audit_freq
  USE EL_FIELDS, ONLY: el_field_data

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: EL_ADD_CONTINUITY_SOURCE, EL_VOLUME_AVERAGED_ACTIVE

CONTAINS

  LOGICAL FUNCTION EL_VOLUME_AVERAGED_ACTIVE()
    EL_VOLUME_AVERAGED_ACTIVE = (TRIM(el_volume_averaged).EQ.'continuity') .AND. &
      ALLOCATED(el_field_data%grad_alpha_p) .AND. ALLOCATED(el_field_data%deps_f_dt)
  END FUNCTION EL_VOLUME_AVERAGED_ACTIVE

  SUBROUTINE EL_ADD_CONTINUITY_SOURCE(defp, val_u, val_v, val_w, mesh, ilev, &
                                      dt, time, istep, mfile)

    REAL*8, INTENT(INOUT) :: defp(:)
    REAL*8, INTENT(IN) :: val_u(:), val_v(:), val_w(:)
    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev, istep, mfile
    REAL*8, INTENT(IN) :: dt, time

    INTEGER :: iel, nel, i, ig, ierr
    REAL*8 :: u(3), eps, deps_dt, adv, g, vol
    REAL*8 :: local(4), global(4), gbar, vol_total
    REAL*8, ALLOCATABLE :: gk(:)

    IF (myid.EQ.master) RETURN
    IF (.NOT.EL_VOLUME_AVERAGED_ACTIVE()) RETURN
    IF (dt.LE.0.0d0) RETURN

    nel = mesh%level(ilev)%nel
    IF (nel.LE.0) RETURN
    IF (SIZE(el_field_data%deps_f_dt).LT.nel) RETURN
    IF (SIZE(defp).LT.4*nel) RETURN

    ALLOCATE(gk(nel))
    local = 0.0d0

    DO iel = 1, nel
      u = 0.0d0
      DO i = 1, 8
        ig = mesh%level(ilev)%kvert(i,iel)
        u(1) = u(1) + val_u(ig)
        u(2) = u(2) + val_v(ig)
        u(3) = u(3) + val_w(ig)
      END DO
      u = 0.125d0*u
      eps = MAX(el_field_data%epsilon_f(iel), 1.0d-3)
      deps_dt = el_field_data%deps_f_dt(iel)
      ! grad eps = - grad alpha_p  (cell-averaged kernel gradient)
      adv = -DOT_PRODUCT(u, el_field_data%grad_alpha_p(:,iel))
      g = -(deps_dt + adv)/eps
      gk(iel) = g
      vol = mesh%level(ilev)%dvol(iel)
      local(1) = local(1) + g*vol
      local(2) = local(2) + vol
      local(3) = local(3) + ABS(deps_dt/eps)*vol
      local(4) = local(4) + ABS(adv/eps)*vol
    END DO

    CALL MPI_Allreduce(local, global, 4, MPI_DOUBLE_PRECISION, MPI_SUM, &
      MPI_COMM_SUBS, ierr)
    vol_total = MAX(global(2), 1.0d-300)
    gbar = global(1)/vol_total

    DO iel = 1, nel
      vol = mesh%level(ilev)%dvol(iel)
      defp(4*(iel-1)+1) = defp(4*(iel-1)+1) + (gk(iel)-gbar)*vol/dt
    END DO

    IF (el_momentum_audit_freq.GT.0) THEN
      IF (MOD(ABS(istep), el_momentum_audit_freq).EQ.0 .AND. myid.EQ.showid) THEN
        WRITE(mfile,'(A,ES14.6,A,I0,A,ES13.6,A,ES13.6,A,ES13.6,A,ES13.6)') &
          'EL_VA_SOURCE t= ', time, ' step= ', istep, &
          ' int_g= ', global(1), ' gbar= ', gbar, &
          ' L1_dedt= ', global(3), ' L1_adv= ', global(4)
        WRITE(*,'(A,ES14.6,A,I0,A,ES13.6,A,ES13.6,A,ES13.6,A,ES13.6)') &
          'EL_VA_SOURCE t= ', time, ' step= ', istep, &
          ' int_g= ', global(1), ' gbar= ', gbar, &
          ' L1_dedt= ', global(3), ' L1_adv= ', global(4)
      END IF
    END IF

    DEALLOCATE(gk)

  END SUBROUTINE EL_ADD_CONTINUITY_SOURCE

END MODULE EL_VA_CONTINUITY
