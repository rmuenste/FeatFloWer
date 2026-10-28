!=========================================================================
! QuadSc_el_va.f90
!
! Euler-Lagrange volume-averaged continuity: source/projection correctors.
! Included in Transport_Q2P1 (host association to LinSc, QuadSc, PLinSc,
! mg_mesh, myStat, Velocity_Correction, Pressure_Correction, ...).
!
! Switched by the deck alone: SimPar@ELVASourceCorrectors = n (default 0 =
! this file is never entered) with SimPar@ELVolumeAveraged = continuity.
!=========================================================================

!=========================================================================
! EL_VA_Projection_Correctors - n source/projection corrector passes after
! the pressure projection of Transport_q2p1_UxyzP_fluid_core.
!
! The predictor solves C dp = -(1/k) B^T u~ + G(u~)/k with the source G
! evaluated from the uncorrected velocity u~. Each corrector pass:
!   1. keeps the current pressure as the base of the next increment,
!   2. reassembles (1/k) B^T u with the corrected velocity u and adds the
!      source G(u)/k recomputed from it (EL_ADD_CONTINUITY_SOURCE, pass icorr),
!   3. logs the constraint residual dt*defP on the audit steps,
!   4. solves the pressure Poisson problem for the increment and applies
!      Velocity_Correction / Pressure_Correction (same solver, same routines).
! The particle advance and the momentum predictor are not repeated. The
! pressure of the old time level P^n is kept in a local copy so that after
! the last pass valP_old holds P^n again and the extrapolated P_new
! (1.5 P^{n+1} - 0.5 P^n, used by the force/torque evaluation) is rebuilt
! from the final P^{n+1}. On the audit steps one more assembly (no solve)
! logs the final residual with pass label -1.
!
! State on exit: valU/valV/valW and valP(NLMAX) corrected n times, valP_old
! = P^n, P_new rebuilt; defP(NLMAX) and rhsP(NLMAX) hold the last pass's
! assembly (both are work arrays of the projection, no reader after the
! step). Nothing else is touched.
!
! Collective: the pressure solve is joined by every rank, so the caller
! must guard the call on the deck keys only (EL_VA_SOURCE_CORRECTORS_ENABLED),
! never on the field allocation, which the master rank does not have.
!=========================================================================
SUBROUTINE EL_VA_Projection_Correctors(mfile, itns)
  use EL_CONFIG, only: el_va_source_correctors
  use EL_VA_CONTINUITY, only: EL_ADD_CONTINUITY_SOURCE, EL_VA_LOG_CONSTRAINT_RESIDUAL, &
                              EL_VA_AUDIT_STEP
  implicit none
  INTEGER, INTENT(IN) :: mfile, itns

  INTEGER :: icorr, i
  REAL*8 :: int_g
  REAL*8, ALLOCATABLE :: p_n(:)

  IF (el_va_source_correctors.LE.0) RETURN

  IF (myid.ne.0) THEN
    ALLOCATE(p_n(SIZE(LinSc%valP_old)))
    p_n = LinSc%valP_old
  END IF

  DO icorr = 1, el_va_source_correctors
    IF (myid.ne.0) THEN
      CALL ZTIME(tttt0)
      LinSc%valP_old = LinSc%valP(NLMAX)%x
      LinSc%valP(NLMAX)%x = 0d0
      CALL Matdef_General_LinScalar(LinSc,QuadSc,PLinSc,1)
      CALL EL_ADD_CONTINUITY_SOURCE(LinSc%defP(NLMAX)%x, QuadSc%valU, QuadSc%valV, &
        QuadSc%valW, mg_mesh, NLMAX, tstep, timens, itns, mfile, icorr, int_g)
      CALL EL_VA_LOG_CONSTRAINT_RESIDUAL(LinSc%defP(NLMAX)%x, mg_mesh%level(NLMAX)%nel, &
        tstep, timens, itns, icorr, int_g, mfile)
      LinSc%rhsP(NLMAX)%x = LinSc%defP(NLMAX)%x
      CALL ZTIME(tttt1)
      myStat%tDefP = myStat%tDefP + (tttt1-tttt0)
    END IF

    CALL Solve_General_LinScalar(LinSc,PLinSc,QuadSc,Boundary_LinScalar_Mat,Boundary_LinScalar_Def,mfile)
    CALL Protocol_LinScalar(mfile,LinSc," Pressure-Poisson corrector")

    IF (myid.ne.0) THEN
      CALL ZTIME(tttt0)
      CALL Velocity_Correction()
      CALL Pressure_Correction()
      CALL ZTIME(tttt1)
      myStat%tCorrUVWP = myStat%tCorrUVWP + (tttt1-tttt0)
    END IF
  END DO

  IF (myid.ne.0) THEN
    ! final constraint residual of the corrected velocity (assembly only, audit steps)
    IF (EL_VA_AUDIT_STEP(itns)) THEN
      CALL Matdef_General_LinScalar(LinSc,QuadSc,PLinSc,1)
      CALL EL_ADD_CONTINUITY_SOURCE(LinSc%defP(NLMAX)%x, QuadSc%valU, QuadSc%valV, &
        QuadSc%valW, mg_mesh, NLMAX, tstep, timens, itns, mfile, -1, int_g)
      CALL EL_VA_LOG_CONSTRAINT_RESIDUAL(LinSc%defP(NLMAX)%x, mg_mesh%level(NLMAX)%nel, &
        tstep, timens, itns, -1, int_g, mfile)
    END IF
    ! restore P^n and rebuild the extrapolated pressure from the final P^{n+1}
    LinSc%valP_old = p_n
    DO i = 1, SIZE(LinSc%P_new)
      LinSc%P_new(i) = 1.5d0*LinSc%valP(NLMAX)%x(i) - 0.5d0*LinSc%valP_old(i)
    END DO
    DEALLOCATE(p_n)
  END IF

END SUBROUTINE EL_VA_Projection_Correctors
