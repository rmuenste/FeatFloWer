!=========================================================================
! QuadSc_pco_operator.f90 - pressure operator renewal and the periodic coarse operator
!
! Defect. The master (rank 0) solves the coarsest pressure level of the multigrid with the
! Schur matrix C = B^T M^-1 B that it assembles on its unpartitioned coarse mesh. There the
! two images of a periodic velocity dof are separate dofs, so the C of the master lacks every
! coupling across a periodic plane, while the workers apply the periodic operator (local C
! plus the parallel part CP). The coarse correction then belongs to another operator and the
! pressure multigrid stagnates on periodic meshes (case s10, xy-periodic boxes).
!
! Entry point. Renew_Pressure_Operator is the one assembly of the pressure Schur operator,
! used by the ordinary matrix renewal (OperatorRegenaration, C branch), by the weighted EL
! operator (EL_MM_Pressure_Operator, every step) and by q2p1_cc:
!     Create_CMat (every rank)  ->  Create_ParCMat (workers)  ->  coarse operator of the master
!
! Decision (first renewal, collective; resolve_mode), in this order:
!   1 no finite period (dPeriodicity >= 1d8 on every axis): legacy; nothing below is entered.
!   2 key = legacy: legacy, logged.
!   3 periodicity: census of the shifted worker links (pco_periodic_topology: element pairs of
!     CP that share a vertex only through a periodic shift). No shifted link, or no pressure
!     structure on the coarse level: periodicity is not established -> legacy, logged, for
!     every key (a walled run stays bitwise whatever the deck says; q2p1_dns_drag declares
!     unit periods on every box). Couplings without a shared vertex in any image: failure.
!   4 default scope (key absent): gathered only inside default_scope (one table below: the
!     configurations covered by gates); outside it legacy with the reason logged. An explicit
!     key is honoured after validation.
!   5 validation: coarse solver type 1-4, MinLev = MedLev <= 3, master solving; no pressure
!     Dirichlet element (knprP) on the coarse level; the master holds the coordinates of the
!     Q2 dofs of the coarse level; complete, unambiguous classes with equal constraint flags
!     that realise every active axis; element graph and layout; every worker coupling inside
!     the layout, one owner per row (gather plan, also in wrapped mode); wrap structure.
!   A failure in 4-5 refuses the run on every rank (reason, counts, first dof and its
!   coordinates) when the key was given, and falls back to legacy (logged) when it is absent.
!
! Every renewal of a periodic mode (renew_periodic, collective): the pattern of the workers
! must be the one of the plan (PCO_GATHER_UNCHANGED) and the coarse level free of pressure
! Dirichlet elements, otherwise every rank stops (a changed pattern needs a restart); values
! (gathered: collected; wrapped: assembled on the master and compared with the gathered values,
! refused above 1e-10 or on a non-finite entry - the master geometry is not the worker geometry
! on curved boundaries); solver refresh on the master (gauge decided at the first refresh,
! factorisation; pco_coarse_solver), status collective; the check when due.
! Legacy keeps its arithmetic: Setup_UMFPACK_CoarseSolver as before; the release of the
! previous legacy factorisation stays at its call sites (EL adapter for types 2/4, MG_Solver
! for MatrixRenewal C >= 2).
!
! Collective operations per periodic renewal: Create_ParCMat (as always), PCO_GATHER_UNCHANGED
! (one pressure exchange + one allreduce), the knprP allreduce, the value gather (gathered
! mode and the wrapped comparison: one MPI_Gatherv), the status broadcast and, when the
! check is due, the check (one gather, three distribute/exchange/collect cycles).
!
! Runtime reload (PCO_PRESSURE_PARAMETERS_RELOADED, from ProcessControl Reload_Pres, every
! rank): the key and the check value are read once and kept; under a periodic operator a new
! coarse solver type 1-4 rebuilds the solver copy (same layout and gauge), any other change of
! the coarse solver type or the master levels is not taken (old values restored, logged).
!
! Output: PCO lines (mode, reason, census, classes, gauge) on the protocol and terminal
! (showid); PCO_TIMING lines at renewals 1, 2 and every 20th (master); PCO_CHECK lines of
! the check (pco_diagnostics).
!
! Extension points: a further operator = a key value + a branch in refresh_raw; a further
! default configuration = an entry of default_scope with its gate; a deeper master hierarchy
! (MinLev < MedLev) = an operator per master level and a tested master cycle (see
! docs/md_docs/pressure_periodic_coarse_operator_2026-11-19.md).
!=========================================================================
MODULE pco_coarse_operator

  USE PP3D_MPI, ONLY: myid, master, showid, MPI_COMM_WORLD, MPI_INTEGER, MPI_LOGICAL, MPI_LOR, &
                      MPI_DOUBLE_PRECISION, MPI_CHARACTER, dPeriodicity
  USE def_FEAT, ONLY: MTERM, MFILE1
  USE var_QuadScalar, ONLY: LinSc, QuadSc, mg_mesh, mg_MlRhoPmat, mg_MlPMat, bNoOutflow, &
                            bMasterTurnedOn
  USE QuadSc_assembly, ONLY: Create_CMat
  USE def_QuadScalar, ONLY: Create_ParCMat
  USE QuadSc_solver_coarse, ONLY: Setup_UMFPACK_CoarseSolver
  USE pressure_coarse_key, ONLY: PCO_KEY_UNSET, PCO_KEY_LEGACY, PCO_KEY_GATHERED, PCO_KEY_WRAPPED, &
                                 PCO_KEY_NAME
  USE pco_periodic_topology, ONLY: tPcoClasses, tPcoLinkCensus, PCO_LINK_CENSUS, PCO_BUILD_CLASSES, &
                                   PCO_CLASS_MASK_CONFLICTS, PCO_ELEMENT_GRAPH, &
                                   PCO_MASTER_HAS_DOF_COORDINATES
  USE pco_gather, ONLY: tPcoGatherPlan, PCO_GATHER_BUILD, PCO_GATHER_VALUES, PCO_GATHER_UNCHANGED, &
                        PCO_GATHER_FREE
  USE pco_wrap, ONLY: tPcoWrap, PCO_WRAP_BUILD, PCO_WRAP_ASSEMBLE, PCO_MASS_ORDINARY, &
                      PCO_MASS_EL_WEIGHTED, PCO_MASS_EL_PER_IMAGE
  USE pco_coarse_solver, ONLY: pco, PCO_INSTALL_LAYOUT, PCO_REFRESH_SOLVER, PCO_RELEASE, &
                               PCO_CHANGE_SOLVER_TYPE
  USE pco_p0_reconstruct, ONLY: PCO_SET_VERTEX_CLASSES
  USE pco_diagnostics, ONLY: PCO_COMPARE_OPERATORS, PCO_RUN_CHECK

  IMPLICIT NONE
  PRIVATE

  TYPE tPcoState
    LOGICAL :: resolved = .FALSE.
    INTEGER :: requested = PCO_KEY_UNSET, mode = PCO_KEY_LEGACY
    INTEGER :: lev = 0, crs = 0, check = 0, minlev = 0, medlev = 0
    INTEGER :: nrenew = 0
    CHARACTER(len=64) :: application = ' '
    TYPE(tPcoLinkCensus) :: census
    TYPE(tPcoClasses) :: cls
    TYPE(tPcoGatherPlan) :: plan
    TYPE(tPcoWrap) :: wrap
    REAL*8 :: t_values = 0d0, t_refresh = 0d0, t_check = 0d0, t_setup = 0d0, t_wrapcheck = 0d0
  END TYPE tPcoState

  TYPE(tPcoState), SAVE :: st

  REAL*8, PARAMETER :: WRAP_TOL = 1d-10

  PUBLIC :: Renew_Pressure_Operator, PCO_REGISTER_APPLICATION, PCO_PRESSURE_PARAMETERS_RELOADED

CONTAINS

  !=======================================================================
  ! An application names itself for the default scope (default_scope); call before the first
  ! renewal. Applications that do not register get legacy when the key is absent.
  !=======================================================================
  SUBROUTINE PCO_REGISTER_APPLICATION(name)
    CHARACTER(len=*), INTENT(IN) :: name
    st%application = name
  END SUBROUTINE PCO_REGISTER_APPLICATION

  !=======================================================================
  ! Default scope: the configurations in which the absent key selects the gathered operator,
  ! each covered by gates (design note, table "configuration -> default -> gate"):
  !   application q2p1_el_pipeflow, coarse pressure level 1 (MinLev = MedLev = 1), coarse
  !   solver type 1-4.
  ! Everything else stays legacy when the key is absent (reason logged) and gets a periodic
  ! operator only on explicit request. reason: why a configuration is outside.
  !=======================================================================
  LOGICAL FUNCTION default_scope(reason) RESULT(inside)
    CHARACTER(len=*), INTENT(OUT) :: reason
    inside = .FALSE.
    IF (TRIM(st%application).NE.'q2p1_el_pipeflow') THEN
      IF (LEN_TRIM(st%application).EQ.0) THEN
        reason = 'this application (unregistered) is not in the validated default scope'
      ELSE
        reason = 'application '//TRIM(st%application)//' is not in the validated default scope'
      END IF
    ELSE IF (st%minlev.NE.1 .OR. st%medlev.NE.1) THEN
      reason = 'a coarse pressure level other than 1 is not in the validated default scope'
    ELSE IF (st%crs.LT.1 .OR. st%crs.GT.4) THEN
      reason = 'coarse solver type outside 1-4'
    ELSE
      inside = .TRUE.
      reason = ' '
    END IF
  END FUNCTION default_scope

  !=======================================================================
  ! The pressure operator renewal (see module header). Collective. w_coarse / rho: the weighted
  ! lumped mass W of the coarse level and the density under the weighted EL operator, whose C
  ! was assembled with the mass U^2/W, U = rho Ml (EL_MM_Pressure_Operator); absent for the
  ! ordinary mass.
  !=======================================================================
  SUBROUTINE Renew_Pressure_Operator(w_coarse, rho)
    REAL*8, INTENT(IN), OPTIONAL :: w_coarse(:), rho

    CALL Create_CMat(QuadSc%knprU, QuadSc%knprV, QuadSc%knprW)
    IF (myid.NE.master) CALL Create_ParCMat(QuadSc%knprU, QuadSc%knprV, QuadSc%knprW)

    IF (.NOT.st%resolved) CALL resolve_mode()

    IF (st%mode.EQ.PCO_KEY_LEGACY) THEN
      CALL Setup_UMFPACK_CoarseSolver(LinSc%knprP, LinSc%prm%MGprmIn%CrsSolverType, &
                                      LinSc%prm%MGprmIn%MinLev, bNoOutflow)
      RETURN
    END IF

    CALL renew_periodic(w_coarse, rho)

  END SUBROUTINE Renew_Pressure_Operator

  !=======================================================================
  ! Mode decision and validation (first renewal, collective). See the module header.
  !=======================================================================
  SUBROUTINE resolve_mode()
    INTEGER :: target, ierr, info(10), lev, nel, nconf, ifirst
    REAL*8 :: xbad(3)
    INTEGER, ALLOCATABLE :: eptr(:), enbr(:), idum(:)
    CHARACTER(len=256) :: reason
    REAL*8 :: t0
    DOUBLE PRECISION MPI_WTIME
    EXTERNAL MPI_WTIME

    st%resolved = .TRUE.
    st%requested = LinSc%prm%PeriodicCoarseOperator
    st%check = LinSc%prm%PeriodicCoarseCheck
    st%minlev = LinSc%prm%MGprmIn%MinLev
    st%medlev = LinSc%prm%MGprmIn%MedLev
    st%lev = st%minlev
    st%crs = LinSc%prm%MGprmIn%CrsSolverType
    st%mode = PCO_KEY_LEGACY
    lev = st%lev

    ! 1 not periodic at all: nothing to decide, nothing to enter
    IF (.NOT.ANY(dPeriodicity.LT.1d8)) THEN
      IF (st%requested.EQ.PCO_KEY_GATHERED .OR. st%requested.EQ.PCO_KEY_WRAPPED) &
        CALL say('PCO periodic coarse operator: '//TRIM(PCO_KEY_NAME(st%requested))// &
                 ' requested, but no axis is periodic: legacy operator')
      RETURN
    END IF
    ! 2 legacy requested
    IF (st%requested.EQ.PCO_KEY_LEGACY) THEN
      CALL say('PCO periodic coarse operator = legacy (requested)')
      RETURN
    END IF
    t0 = MPI_WTIME()

    ! 3 is the mesh periodic? (before any periodic-only restriction)
    CALL PCO_LINK_CENSUS(lev, st%census)
    IF (.NOT.st%census%possible) THEN
      WRITE(reason,'(A,I0,A)') 'PCO periodicity not established: no pressure structure on level ', lev, &
        ' on every worker; legacy operator'
      CALL say(reason)
      RETURN
    END IF
    IF (.NOT.ANY(st%census%axis_active) .AND. st%census%nbad.EQ.0) THEN
      WRITE(reason,'(A,3ES11.3,A,I0,A)') 'PCO periods declared (', dPeriodicity, &
        ') but no worker coupling crosses a periodic plane on level ', lev, &
        ': periodicity inactive, legacy operator'
      CALL say(reason)
      RETURN
    END IF

    ! 4 default scope
    IF (st%requested.EQ.PCO_KEY_UNSET) THEN
      IF (.NOT.default_scope(reason)) THEN
        CALL say('PCO periodic coarse operator: legacy (key absent) - '//TRIM(reason)// &
                 '; Pres@PeriodicCoarseOperator = gathered requests it')
        RETURN
      END IF
      target = PCO_KEY_GATHERED
    ELSE
      target = st%requested
    END IF

    IF (st%census%nbad.GT.0) THEN
      WRITE(reason,'(I0,A,I0,A)') st%census%nbad, ' worker couplings on level ', lev, &
        ' share no vertex in any periodic image: dPeriodicity inconsistent with the mesh'
      CALL fail(target, reason)
      RETURN
    END IF

    ! 5 validation: configuration (deck values, the same on every rank)
    reason = ' '
    IF (.NOT.bMasterTurnedOn) THEN
      reason = 'the master does not solve the coarse pressure level (bMasterTurnedOn off)'
    ELSE IF (st%crs.LT.1 .OR. st%crs.GT.4) THEN
      WRITE(reason,'(A,I0,A)') 'Pres@MGCrsSolverType ', st%crs, &
        ' does not solve with the master matrix (supported: 1, 2, 3, 4)'
    ELSE IF (st%minlev.NE.st%medlev) THEN
      WRITE(reason,'(A,I0,A,I0,A)') 'Pres@MGMinLev ', st%minlev, ' < MGMedLev ', st%medlev, &
        ': a master multigrid hierarchy is not supported'
    ELSE IF (lev.LT.1 .OR. lev.GT.3) THEN
      WRITE(reason,'(A,I0,A)') 'coarse level ', lev, ' outside 1..3 (coarse vector transfers)'
    END IF
    IF (reason.NE.' ') THEN
      CALL fail(target, reason)
      RETURN
    END IF

    ! no pressure Dirichlet elements on the coarse level
    IF (knprp_on_level(lev)) THEN
      CALL fail(target, 'pressure Dirichlet elements (knprP) on the coarse level are not supported')
      RETURN
    END IF

    ! master: coordinates, classes, constraint flags, layout
    info = 0
    xbad = 0d0
    ifirst = 0
    IF (myid.EQ.master) THEN
      IF (.NOT.PCO_MASTER_HAS_DOF_COORDINATES(lev)) THEN
        info(9) = 1
      ELSE
        nel = mg_mesh%level(lev)%nel
        CALL PCO_BUILD_CLASSES(lev, st%census%axis_active, dPeriodicity, st%cls)
        nconf = PCO_CLASS_MASK_CONFLICTS(st%cls, QuadSc%knprU(lev)%x, QuadSc%knprV(lev)%x, &
                                         QuadSc%knprW(lev)%x, ifirst)
        info(1) = st%cls%nambiguous
        info(2) = st%cls%nincomplete
        info(3) = st%cls%nmixed
        info(4) = nconf
        info(5) = st%cls%nmulti
        info(6) = st%cls%maxsize
        info(7) = st%cls%first_bad
        IF (info(7).EQ.0) info(7) = ifirst
        IF (info(7).GT.0) xbad = mg_mesh%level(lev+1)%dcorvg(:, info(7))
        info(10) = COUNT(st%census%axis_active .AND. .NOT.st%cls%axis_found)
        IF (ALL(info(1:4).EQ.0) .AND. info(10).EQ.0) THEN
          CALL PCO_ELEMENT_GRAPH(lev, st%cls, eptr, enbr)
          CALL PCO_INSTALL_LAYOUT(lev, nel, eptr, enbr, st%crs, target)
          info(8) = SIZE(pco%colA)
        END IF
      END IF
    END IF
    CALL MPI_Bcast(info, 10, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_Bcast(xbad, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    IF (info(9).NE.0) THEN
      WRITE(reason,'(A,I0,A)') 'the master holds no coordinates of the Q2 dofs of level ', lev, &
        ' (its mesh is not refined to level MedLev+1 by this application)'
      CALL fail(target, reason)
      RETURN
    END IF
    IF (ANY(info(1:4).NE.0)) THEN
      WRITE(reason,'(A,I0,A,4(I0,A),I0,A,3ES11.3,A)') 'periodic classes on level ', lev, ': ambiguous ', info(1), &
        ', incomplete ', info(2), ', mixed ', info(3), ', constraint conflicts ', info(4), &
        ' (first: dof ', info(7), ' at', xbad, ')'
      reason = TRIM(reason)//': a missing or ambiguous partner, or unequal conditions on the two faces of a periodic pair'
      CALL fail(target, reason)
      RETURN
    END IF
    IF (info(10).NE.0) THEN
      WRITE(reason,'(I0,A)') info(10), ' axes with shifted worker links have no periodic images in the master classes'
      CALL fail(target, reason)
      RETURN
    END IF

    ! gather plan: also validates the master classes against the worker couplings
    IF (myid.EQ.master) THEN
      CALL PCO_GATHER_BUILD(st%plan, lev, pco%ldA, pco%colA)
    ELSE
      ALLOCATE(idum(2))
      idum = 1
      CALL PCO_GATHER_BUILD(st%plan, lev, idum, idum)
      DEALLOCATE(idum)
    END IF
    IF (st%plan%nmissing.GT.0 .OR. st%plan%nownership.GT.0 .OR. st%plan%nunowned.GT.0) THEN
      WRITE(reason,'(A,I0,A,3(I0,A))') 'worker operator on level ', lev, ' vs master layout: ', &
        st%plan%nmissing, ' couplings outside the periodic classes, ', st%plan%nownership, &
        ' rows of two owners, ', st%plan%nunowned, ' rows without owner'
      CALL fail(target, reason)
      RETURN
    END IF

    ! wrap structure
    IF (target.EQ.PCO_KEY_WRAPPED .OR. st%check.GT.0) THEN
      info(1) = 0
      IF (myid.EQ.master) THEN
        CALL PCO_WRAP_BUILD(st%wrap, lev, st%cls)
        info(1) = st%wrap%nbad
      END IF
      CALL MPI_Bcast(info, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      IF (info(1).GT.0) THEN
        WRITE(reason,'(I0,A)') info(1), ' element pairs of the wrapped assembly outside the layout'
        CALL fail(target, reason)
        RETURN
      END IF
    END IF

    st%mode = target
    IF (myid.EQ.master) THEN
      CALL PCO_SET_VERTEX_CLASSES(st%cls%nvt, st%cls%ncls, st%cls%ptr, st%cls%mem)
      pco%active = .TRUE.
    END IF
    st%t_setup = MPI_WTIME() - t0

    WRITE(reason,'(A,A,A,A,I0,A,3L1,A,3(I0,1X))') 'PCO periodic coarse operator = ', &
      TRIM(PCO_KEY_NAME(target)), MERGE(' (default)  ', ' (requested)', st%requested.EQ.PCO_KEY_UNSET), &
      ' level ', lev, ' axes ', st%census%axis_active, ' shifted links ', st%census%nshift
    CALL say(reason)
    WRITE(reason,'(A,I0,A,I0,A,I0,A,I0,A,I0)') 'PCO classes ', info(5), ' (largest ', info(6), &
      '), layout entries ', info(8), ', uncovered by workers ', st%plan%nuncovered, &
      ', coarse solver type ', st%crs
    CALL say(reason)

  END SUBROUTINE resolve_mode

  ! some rank has a pressure Dirichlet element on level lev (collective)
  LOGICAL FUNCTION knprp_on_level(lev) RESULT(bad_any)
    INTEGER, INTENT(IN) :: lev
    LOGICAL :: bad_local
    INTEGER :: ierr
    bad_local = .FALSE.
    IF (ALLOCATED(LinSc%knprP)) bad_local = ANY(LinSc%knprP(lev)%x(1:mg_mesh%level(lev)%nel).NE.0)
    CALL MPI_Allreduce(bad_local, bad_any, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
  END FUNCTION knprp_on_level

  ! validation failed: refusal (explicit key) or logged fallback to legacy (key absent)
  SUBROUTINE fail(target, reason)
    INTEGER, INTENT(IN) :: target
    CHARACTER(len=*), INTENT(IN) :: reason
    IF (st%requested.EQ.PCO_KEY_UNSET) THEN
      CALL say('PCO periodic coarse operator: fallback to legacy - '//TRIM(reason))
      st%mode = PCO_KEY_LEGACY
      IF (myid.EQ.master) CALL PCO_RELEASE()
      CALL PCO_GATHER_FREE(st%plan)
    ELSE
      CALL stop_all('PCO REFUSED: Pres@PeriodicCoarseOperator = '//TRIM(PCO_KEY_NAME(target))//': '//TRIM(reason))
    END IF
  END SUBROUTINE fail

  ! every rank reaches this with the same decision: report (master and protocol rank), stop
  SUBROUTINE stop_all(text)
    CHARACTER(len=*), INTENT(IN) :: text
    IF (myid.EQ.master .OR. myid.EQ.showid) WRITE(*,'(A)') TRIM(text)
    IF (myid.EQ.showid) WRITE(MFILE1,'(A)') TRIM(text)
    FLUSH(6)
    STOP 1
  END SUBROUTINE stop_all

  SUBROUTINE say(text)
    CHARACTER(len=*), INTENT(IN) :: text
    IF (myid.EQ.showid) THEN
      WRITE(MTERM,'(A)') TRIM(text)
      WRITE(MFILE1,'(A)') TRIM(text)
    END IF
  END SUBROUTINE say

  !=======================================================================
  ! One renewal of a periodic mode (collective).
  !=======================================================================
  SUBROUTINE renew_periodic(w_coarse, rho)
    REAL*8, INTENT(IN), OPTIONAL :: w_coarse(:), rho
    INTEGER :: status, ierr
    CHARACTER(len=256) :: msg
    LOGICAL :: do_check
    REAL*8 :: t0, t1, t2
    DOUBLE PRECISION MPI_WTIME
    EXTERNAL MPI_WTIME

    st%nrenew = st%nrenew + 1
    do_check = (st%check.GT.0) .AND. (st%nrenew.LE.2 .OR. MOD(st%nrenew, MAX(st%check,1)).EQ.0)

    ! the plan, the layout, the classes and the gauge belong to one pattern
    IF (.NOT.PCO_GATHER_UNCHANGED(st%plan)) CALL stop_all('PCO STOP: the coarse pressure pattern, '// &
      'ghost numbering or ownership of the workers changed since the first renewal; not supported (restart)')
    IF (knprp_on_level(st%lev)) CALL stop_all('PCO STOP: pressure Dirichlet elements (knprP) '// &
      'appeared on the coarse level; not supported with a periodic coarse operator')

    t0 = MPI_WTIME()
    CALL refresh_raw(st%mode, w_coarse, rho)
    t1 = MPI_WTIME()
    IF (st%mode.EQ.PCO_KEY_WRAPPED) CALL validate_wrapped()
    t2 = MPI_WTIME()
    st%t_values = st%t_values + (t1 - t0)
    st%t_wrapcheck = st%t_wrapcheck + (t2 - t1)

    status = 0
    msg = ' '
    t0 = MPI_WTIME()
    IF (myid.EQ.master) CALL PCO_REFRESH_SOLVER(status, msg)
    st%t_refresh = st%t_refresh + (MPI_WTIME() - t0)
    CALL MPI_Bcast(status, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    IF (status.NE.0) THEN
      CALL MPI_Bcast(msg, LEN(msg), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      CALL stop_all('PCO STOP (renewal): '//TRIM(msg))
    END IF

    IF (st%nrenew.EQ.1 .AND. myid.EQ.master) CALL report_gauge()
    IF (do_check) THEN
      t0 = MPI_WTIME()
      CALL run_check(w_coarse, rho)
      st%t_check = st%t_check + (MPI_WTIME() - t0)
    END IF
    IF (myid.EQ.master .AND. (st%nrenew.LE.2 .OR. MOD(st%nrenew,20).EQ.0)) &
      WRITE(MTERM,'(A,I0,A,A,5(A,ES10.3),2(A,I0),A,ES10.3,A,ES10.3,A,I0)') 'PCO_TIMING renewal ', st%nrenew, &
        ' mode ', TRIM(PCO_KEY_NAME(st%mode)), ' setup_s ', st%t_setup, ' values_s ', st%t_values, &
        ' wrapcheck_s ', st%t_wrapcheck, ' refresh_s ', st%t_refresh, ' check_s ', st%t_check, &
        ' coarse solves ', pco%nsolve, ' reported imbalances ', pco%nincompat, &
        ' largest relative imbalance ', pco%max_incompat, ' net ', pco%max_incompat_net, &
        ' component ', pco%max_incompat_comp

  END SUBROUTINE renew_periodic

  !-----------------------------------------------------------------------
  ! Raw values of an operator into pco%a (or a_out). Wrapped forms the unweighted mass
  ! U = rho Ml on the master only, here, when the weighted EL operator is renewed.
  !-----------------------------------------------------------------------
  SUBROUTINE refresh_raw(mode, w_coarse, rho, a_out, mass_variant)
    INTEGER, INTENT(IN) :: mode
    REAL*8, INTENT(IN), OPTIONAL :: w_coarse(:), rho
    REAL*8, INTENT(INOUT), OPTIONAL, TARGET :: a_out(:)
    INTEGER, INTENT(IN), OPTIONAL :: mass_variant
    REAL*8 :: dum(1)
    REAL*8, ALLOCATABLE :: u(:)
    INTEGER :: kind, lev

    lev = st%lev
    SELECT CASE (mode)
    CASE (PCO_KEY_GATHERED)
      IF (myid.EQ.master) THEN
        IF (PRESENT(a_out)) THEN
          CALL PCO_GATHER_VALUES(st%plan, a_out)
        ELSE
          CALL PCO_GATHER_VALUES(st%plan, pco%a)
        END IF
      ELSE
        CALL PCO_GATHER_VALUES(st%plan, dum)
      END IF
    CASE (PCO_KEY_WRAPPED)
      IF (myid.NE.master) RETURN
      kind = PCO_MASS_ORDINARY
      IF (PRESENT(w_coarse)) kind = PCO_MASS_EL_WEIGHTED
      IF (PRESENT(mass_variant)) kind = mass_variant
      IF (kind.EQ.PCO_MASS_ORDINARY) THEN
        IF (PRESENT(a_out)) THEN
          CALL PCO_WRAP_ASSEMBLE(st%wrap, st%cls, kind, mg_MlRhoPmat(lev)%a, dum, dum, &
            QuadSc%knprU(lev)%x, QuadSc%knprV(lev)%x, QuadSc%knprW(lev)%x, a_out)
        ELSE
          CALL PCO_WRAP_ASSEMBLE(st%wrap, st%cls, kind, mg_MlRhoPmat(lev)%a, dum, dum, &
            QuadSc%knprU(lev)%x, QuadSc%knprV(lev)%x, QuadSc%knprW(lev)%x, pco%a)
        END IF
      ELSE
        ALLOCATE(u(SIZE(mg_MlPMat(lev)%a)))
        u = rho*mg_MlPMat(lev)%a
        IF (PRESENT(a_out)) THEN
          CALL PCO_WRAP_ASSEMBLE(st%wrap, st%cls, kind, dum, u, w_coarse, &
            QuadSc%knprU(lev)%x, QuadSc%knprV(lev)%x, QuadSc%knprW(lev)%x, a_out)
        ELSE
          CALL PCO_WRAP_ASSEMBLE(st%wrap, st%cls, kind, dum, u, w_coarse, &
            QuadSc%knprU(lev)%x, QuadSc%knprV(lev)%x, QuadSc%knprW(lev)%x, pco%a)
        END IF
        DEALLOCATE(u)
      END IF
    END SELECT
  END SUBROUTINE refresh_raw

  !-----------------------------------------------------------------------
  ! The wrapped operator is the master assembly of what the workers apply. It is exact only
  ! where the master has the coefficients and the geometry of the workers: on a curved
  ! boundary the master does not project the coordinates of its finest level (app_init:
  ! ParametrizeBndr of level MedLev+1 on the workers only), and on the s05 pipe its entries
  ! differ by up to 18 % next to the curved wall, where the legacy master matrix has the
  ! wrapped value (consistent with the missing projection; not proven to be the only cause).
  ! Every renewal the wrapped values are compared with the gathered ones entry by entry
  ! (PCO_COMPARE_OPERATORS: rows scaled by their largest modulus, rows that are zero in the
  ! gathered operator by the largest modulus overall, non-finite entries rejected); a
  ! difference above WRAP_TOL stops the run. Collective.
  !-----------------------------------------------------------------------
  SUBROUTINE validate_wrapped()
    REAL*8, ALLOCATABLE :: ag(:)
    REAL*8 :: d
    INTEGER :: ierr, kmax
    IF (myid.EQ.master) THEN
      ALLOCATE(ag(SIZE(pco%a)))
    ELSE
      ALLOCATE(ag(1))
    END IF
    CALL refresh_raw(PCO_KEY_GATHERED, a_out=ag)
    d = 0d0
    IF (myid.EQ.master) THEN
      d = PCO_COMPARE_OPERATORS(pco%a, ag, kmax)
      IF (st%nrenew.LE.2 .OR. MOD(st%nrenew,20).EQ.0 .OR. d.GT.WRAP_TOL) &
        WRITE(MTERM,'(A,I0,A,ES10.3,A,ES8.1,A)') 'PCO wrapped operator vs gathered, renewal ', st%nrenew, &
        ': ', d, ' (limit ', WRAP_TOL, ', entry difference / largest modulus of the row)'
    END IF
    DEALLOCATE(ag)
    CALL MPI_Bcast(d, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    IF (d.GT.WRAP_TOL) THEN
      IF (st%nrenew.EQ.1) THEN
        CALL stop_all('PCO REFUSED: Pres@PeriodicCoarseOperator = wrapped: the master assembly differs '// &
          'from the operator of the workers (curved boundary? the master does not project its '// &
          'finest-level coordinates); use gathered')
      ELSE
        CALL stop_all('PCO STOP: the wrapped operator left the gathered one at a later renewal; use gathered')
      END IF
    END IF
  END SUBROUTINE validate_wrapped

  ! gauge decision of the first renewal against SimPar@NoOutflow (master)
  SUBROUTINE report_gauge()
    CHARACTER(len=256) :: line
    WRITE(line,'(A,I0,A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,L1)') 'PCO gauge: components ', pco%ncomp, &
      ', closed (C e = 0, one zero reference each) ', pco%nclosed, ', closure measure closed ', &
      pco%closure_closed, ' open ', MERGE(pco%closure_open, 0d0, pco%nclosed.LT.pco%ncomp), &
      ' left closed ', pco%left_closed, ', SimPar@NoOutflow ', bNoOutflow
    WRITE(MTERM,'(A)') TRIM(line)
    IF ((pco%nclosed.GT.0) .NEQV. bNoOutflow) WRITE(MTERM,'(A)') &
      'PCO gauge note: the operator decides the reference (closed components), SimPar@NoOutflow disagrees'
    IF (pco%nclosed.GT.0 .AND. pco%left_closed.GT.1d-8) WRITE(MTERM,'(A)') &
      'PCO gauge note: e^T A /= 0 on a closed component (non-symmetric operator): the imbalance '// &
      'measure of the right-hand side is an indicator only'
  END SUBROUTINE report_gauge

  ! the check (pco_diagnostics): both operators on the master, collective
  SUBROUTINE run_check(w_coarse, rho)
    REAL*8, INTENT(IN), OPTIONAL :: w_coarse(:), rho
    REAL*8, ALLOCATABLE :: ag(:), aw(:), ai(:)
    INTEGER :: nnz
    LOGICAL :: weighted
    weighted = PRESENT(w_coarse)
    nnz = 1
    IF (myid.EQ.master) nnz = SIZE(pco%colA)
    ALLOCATE(ag(nnz), aw(nnz), ai(nnz))
    CALL refresh_raw(PCO_KEY_GATHERED, a_out=ag)
    IF (myid.EQ.master) THEN
      CALL refresh_raw(PCO_KEY_WRAPPED, w_coarse, rho, a_out=aw)
      IF (weighted) THEN
        CALL refresh_raw(PCO_KEY_WRAPPED, w_coarse, rho, a_out=ai, mass_variant=PCO_MASS_EL_PER_IMAGE)
      ELSE
        ai = aw
      END IF
    END IF
    CALL PCO_RUN_CHECK(st%lev, st%crs, st%nrenew, st%cls, ag, aw, ai, weighted)
    DEALLOCATE(ag, aw, ai)
  END SUBROUTINE run_check

  !=======================================================================
  ! After a runtime reload of the pressure parameters (ProcessControl, every rank).
  !=======================================================================
  SUBROUTINE PCO_PRESSURE_PARAMETERS_RELOADED()
    INTEGER :: status, ierr
    CHARACTER(len=256) :: msg
    INTEGER :: crs_new
    IF (.NOT.st%resolved) RETURN
    IF (LinSc%prm%PeriodicCoarseOperator.NE.st%requested .OR. LinSc%prm%PeriodicCoarseCheck.NE.st%check) THEN
      LinSc%prm%PeriodicCoarseOperator = st%requested
      LinSc%prm%PeriodicCoarseCheck = st%check
      CALL say('PCO: Pres@PeriodicCoarseOperator / PeriodicCoarseCheck are read once; the reloaded values are not taken')
    END IF
    IF (st%mode.EQ.PCO_KEY_LEGACY) RETURN
    crs_new = LinSc%prm%MGprmIn%CrsSolverType
    IF (LinSc%prm%MGprmIn%MinLev.NE.st%minlev .OR. LinSc%prm%MGprmIn%MedLev.NE.st%medlev .OR. &
        crs_new.LT.1 .OR. crs_new.GT.4) THEN
      LinSc%prm%MGprmIn%MinLev = st%minlev
      LinSc%prm%MGprmIn%MedLev = st%medlev
      LinSc%prm%MGprmIn%CrsSolverType = st%crs
      CALL say('PCO: a reloaded coarse level or a coarse solver type outside 1-4 is not taken under the '// &
               'periodic coarse operator (restart for that); previous values kept')
      RETURN
    END IF
    IF (crs_new.EQ.st%crs) RETURN
    status = 0
    msg = ' '
    IF (myid.EQ.master) CALL PCO_CHANGE_SOLVER_TYPE(crs_new, status, msg)
    CALL MPI_Bcast(status, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    IF (status.NE.0) THEN
      CALL MPI_Bcast(msg, LEN(msg), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      CALL stop_all('PCO STOP (reload of the coarse solver type): '//TRIM(msg))
    END IF
    WRITE(msg,'(A,I0,A,I0)') 'PCO: coarse solver type changed by reload from ', st%crs, ' to ', crs_new
    CALL say(msg)
    st%crs = crs_new
  END SUBROUTINE PCO_PRESSURE_PARAMETERS_RELOADED

END MODULE pco_coarse_operator
