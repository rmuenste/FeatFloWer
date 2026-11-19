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
! used by the ordinary matrix renewal (OperatorRegenaration, C branch) and by the weighted EL
! operator (EL_MM_Pressure_Operator, every step):
!     Create_CMat (every rank)  ->  Create_ParCMat (workers)  ->  coarse operator of the master
! The coarse operator is chosen once (first renewal) by Pres@PeriodicCoarseOperator
! (pressure_coarse_key):
!
!   no finite period (dPeriodicity >= 1d8 on every axis)
!        -> legacy, whatever the key says; nothing below is entered (bitwise as before)
!   legacy requested
!        -> legacy, logged; nothing else entered
!   otherwise (gathered: key absent or gathered; wrapped: key wrapped), validation:
!     deck   coarse solver type 1-4, Pres@MGMinLev = MGMedLev (<= 3), master solving
!     census shifted worker links (pco_periodic_topology); none on any axis -> periodicity
!            inactive -> legacy, logged (not an error: q2p1_dns_drag declares unit periods
!            on non-periodic boxes); couplings inconsistent with the periods -> failure
!     flags  no pressure Dirichlet element (knprP) on the coarse level
!     master periodic dof classes complete and unambiguous, equal velocity constraint flags
!            over every class, element graph and layout (pco_coarse_solver)
!     plan   every worker coupling inside the master layout, every master row owned by
!            exactly one worker (pco_gather; built in both modes: it cross-checks the master
!            classes against the links the workers really use)
!     wrap   structure inside the layout (wrapped mode, or the check)
!   A failure refuses the run on every rank (reason, level, counts) when the value was given
!   in the deck, and falls back to legacy with a logged reason when the key is absent.
!   Wrapped is in addition compared with the gathered operator at the first renewal and
!   refused when they differ (validate_wrapped: the master geometry is not the worker
!   geometry on curved boundaries).
!
! Every renewal of a periodic mode: values (gathered: collected from the workers; wrapped:
! assembled on the master from B, B^T and the class masses), then the solver refresh on the
! master (gauge, factorisation; pco_coarse_solver), status collective. Cached: classes,
! layout, gather plan, wrap structure, symbolic analysis; refreshed: values and numeric
! factors. Legacy keeps its arithmetic: the release of the previous UMFPACK factorisation
! (formerly in EL_MM_Pressure_Operator, flag-protected) and Setup_UMFPACK_CoarseSolver.
!
! Check (Pres@PeriodicCoarseCheck = n > 0, renewals 1, 2 and every n-th; any periodic mode):
! both operators are built; entry-by-entry difference gathered - wrapped; for three test
! vectors (mean modes = 1; a four-mode impulse on an element at the largest vertex class;
! a reproducible vector exciting all modes) the worker product (E012_DAX arithmetic: C x +
! CP xp on the workers, distributed and collected with the coarse vector transfers)
! against the master product of the gathered and of the wrapped matrix, every row, scaled by
! sum_j |A_ij| max|v|; under the weighted EL mass also the per-image transformed variant (to show
! the check discriminates); the closure measure; a compatible solve with the ungauged
! residual, the reference value and the left null vector (e^T A = 0 on closed components).
! Output: PCO_CHECK lines on the terminal of the master.
!
! Output: PCO lines (mode, reason, census, classes, gauge) on the protocol and terminal
! (showid); PCO_TIMING lines at renewals 1, 2 and every 20th (master).
!
! Extension points: a further operator = a key value + a branch in refresh_raw; a deeper
! master hierarchy (MinLev < MedLev) = an operator per master level and a tested master
! cycle (see docs/md_docs/pressure_periodic_coarse_operator_2026-11-19.md).
!=========================================================================
MODULE pco_coarse_operator

  USE PP3D_MPI, ONLY: myid, master, showid, numnodes, MPI_COMM_WORLD, MPI_INTEGER, MPI_LOGICAL, &
                      MPI_LOR, MPI_SUM, MPI_MAX, MPI_DOUBLE_PRECISION, dPeriodicity
  USE def_FEAT, ONLY: MTERM, MFILE1, ILEV
  USE var_QuadScalar, ONLY: LinSc, QuadSc, mg_mesh, mg_MlRhoPmat, mg_CMat, mg_CPMat, mg_lMat, &
                            mg_lPMat, bNoOutflow, bMasterTurnedOn
  USE QuadSc_assembly, ONLY: Create_CMat
  USE def_QuadScalar, ONLY: Create_ParCMat
  USE QuadSc_solver_coarse, ONLY: Setup_UMFPACK_CoarseSolver
  USE UMFPackSolver, ONLY: myUmfPack_Free
  USE pressure_coarse_key, ONLY: PCO_KEY_UNSET, PCO_KEY_LEGACY, PCO_KEY_GATHERED, PCO_KEY_WRAPPED, &
                                 PCO_KEY_NAME
  USE pco_periodic_topology, ONLY: tPcoClasses, tPcoLinkCensus, PCO_LINK_CENSUS, PCO_BUILD_CLASSES, &
                                   PCO_CLASS_MASK_CONFLICTS, PCO_ELEMENT_GRAPH, PCO_GHOST_BUFFER_SIZE
  USE pco_gather, ONLY: tPcoGatherPlan, PCO_GATHER_BUILD, PCO_GATHER_VALUES, PCO_GATHER_STALE, &
                        PCO_GATHER_FREE
  USE pco_wrap, ONLY: tPcoWrap, PCO_WRAP_BUILD, PCO_WRAP_ASSEMBLE, PCO_MASS_ORDINARY, &
                      PCO_MASS_EL_WEIGHTED, PCO_MASS_EL_PER_IMAGE
  USE pco_coarse_solver, ONLY: pco, PCO_INSTALL_LAYOUT, PCO_REFRESH_SOLVER, PCO_RAW_MATVEC, &
                               PCO_RELEASE, PCO_SOLVE_CHECK
  USE pco_p0_reconstruct, ONLY: PCO_SET_VERTEX_CLASSES

  IMPLICIT NONE
  PRIVATE

  TYPE tPcoState
    LOGICAL :: resolved = .FALSE.
    INTEGER :: requested = PCO_KEY_UNSET, mode = PCO_KEY_LEGACY
    INTEGER :: lev = 0, crs = 0, check = 0
    INTEGER :: nrenew = 0
    TYPE(tPcoLinkCensus) :: census
    TYPE(tPcoClasses) :: cls
    TYPE(tPcoGatherPlan) :: plan
    TYPE(tPcoWrap) :: wrap
    REAL*8 :: t_values = 0d0, t_refresh = 0d0, t_check = 0d0, t_setup = 0d0
  END TYPE tPcoState

  TYPE(tPcoState), SAVE :: st

  PUBLIC :: Renew_Pressure_Operator

CONTAINS

  !=======================================================================
  ! The pressure operator renewal (see module header). Collective. u_coarse / w_coarse:
  ! the unweighted (rho Ml) and weighted lumped masses of the coarse level under the weighted
  ! EL operator, whose C was assembled with the mass U^2/W (EL_MM_Pressure_Operator); absent
  ! for the ordinary mass.
  !=======================================================================
  SUBROUTINE Renew_Pressure_Operator(u_coarse, w_coarse)
    REAL*8, INTENT(IN), OPTIONAL :: u_coarse(:), w_coarse(:)

    CALL Create_CMat(QuadSc%knprU, QuadSc%knprV, QuadSc%knprW)
    IF (myid.NE.master) CALL Create_ParCMat(QuadSc%knprU, QuadSc%knprV, QuadSc%knprW)

    IF (.NOT.st%resolved) CALL resolve_mode()

    IF (st%mode.EQ.PCO_KEY_LEGACY) THEN
      ! the coarse solver setup of before: release the previous UMFPACK factorisation (a no-op
      ! when there is none) and factorise the own matrix of the master
      IF (myid.EQ.master) CALL myUmfPack_Free()
      CALL Setup_UMFPACK_CoarseSolver(LinSc%knprP, LinSc%prm%MGprmIn%CrsSolverType, &
                                      LinSc%prm%MGprmIn%MinLev, bNoOutflow)
      RETURN
    END IF

    CALL renew_periodic(u_coarse, w_coarse)

  END SUBROUTINE Renew_Pressure_Operator

  !=======================================================================
  ! Mode decision and validation (first renewal, collective).
  !=======================================================================
  SUBROUTINE resolve_mode()
    INTEGER :: target, ierr, info(8), lev, nel, nconf, ifirst
    REAL*8 :: xbad(3)
    LOGICAL :: bad_local, bad_any
    INTEGER, ALLOCATABLE :: eptr(:), enbr(:), idum(:)
    CHARACTER(len=256) :: reason
    REAL*8 :: t0
    DOUBLE PRECISION MPI_WTIME
    EXTERNAL MPI_WTIME

    st%resolved = .TRUE.
    st%requested = LinSc%prm%PeriodicCoarseOperator
    st%check = LinSc%prm%PeriodicCoarseCheck
    st%lev = LinSc%prm%MGprmIn%MinLev
    st%crs = LinSc%prm%MGprmIn%CrsSolverType
    st%mode = PCO_KEY_LEGACY

    ! not periodic at all: nothing to decide, nothing to enter
    IF (.NOT.ANY(dPeriodicity.LT.1d8)) THEN
      IF (st%requested.EQ.PCO_KEY_GATHERED .OR. st%requested.EQ.PCO_KEY_WRAPPED) &
        CALL say('PCO periodic coarse operator: '//TRIM(PCO_KEY_NAME(st%requested))// &
                 ' requested, but no axis is periodic: legacy operator')
      RETURN
    END IF
    IF (st%requested.EQ.PCO_KEY_LEGACY) THEN
      CALL say('PCO periodic coarse operator = legacy (requested)')
      RETURN
    END IF
    target = MERGE(PCO_KEY_GATHERED, st%requested, st%requested.EQ.PCO_KEY_UNSET)
    t0 = MPI_WTIME()
    lev = st%lev

    ! configuration supported? (deck values, the same on every rank)
    reason = ' '
    IF (.NOT.bMasterTurnedOn) THEN
      reason = 'the master does not solve the coarse pressure level (bMasterTurnedOn off)'
    ELSE IF (st%crs.LT.1 .OR. st%crs.GT.4) THEN
      WRITE(reason,'(A,I0,A)') 'Pres@MGCrsSolverType ', st%crs, &
        ' does not solve with the master matrix (supported: 1, 2, 3, 4)'
    ELSE IF (LinSc%prm%MGprmIn%MinLev.NE.LinSc%prm%MGprmIn%MedLev) THEN
      WRITE(reason,'(A,I0,A,I0,A)') 'Pres@MGMinLev ', LinSc%prm%MGprmIn%MinLev, ' < MGMedLev ', &
        LinSc%prm%MGprmIn%MedLev, ': a master multigrid hierarchy is not supported'
    ELSE IF (lev.LT.1 .OR. lev.GT.3) THEN
      WRITE(reason,'(A,I0,A)') 'coarse level ', lev, ' outside 1..3 (coarse vector transfers)'
    END IF
    IF (reason.NE.' ') THEN
      CALL fail(target, reason)
      RETURN
    END IF

    ! shifted worker links
    CALL PCO_LINK_CENSUS(lev, st%census)
    IF (st%census%nbad.GT.0) THEN
      WRITE(reason,'(I0,A,I0,A,ES10.3,A)') st%census%nbad, ' worker couplings on level ', lev, &
        ' match no period (largest offset ', st%census%worst, ' element sizes): dPeriodicity inconsistent with the mesh'
      CALL fail(target, reason)
      RETURN
    END IF
    IF (.NOT.ANY(st%census%axis_active)) THEN
      WRITE(reason,'(A,3ES11.3,A,I0,A)') 'PCO periods declared (', dPeriodicity, &
        ') but no worker coupling crosses a periodic plane on level ', lev, &
        ': periodicity inactive, legacy operator'
      CALL say(reason)
      RETURN
    END IF

    ! no pressure Dirichlet elements on the coarse level (Boundary_LinScalar_Mat would change
    ! the matrices of master and workers behind the back of this operator)
    bad_local = .FALSE.
    IF (ALLOCATED(LinSc%knprP)) bad_local = ANY(LinSc%knprP(lev)%x(1:mg_mesh%level(lev)%nel).NE.0)
    CALL MPI_Allreduce(bad_local, bad_any, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (bad_any) THEN
      CALL fail(target, 'pressure Dirichlet elements (knprP) on the coarse level are not supported')
      RETURN
    END IF

    ! master: classes, constraint flags, layout
    info = 0
    xbad = 0d0
    ifirst = 0
    IF (myid.EQ.master) THEN
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
      IF (ALL(info(1:4).EQ.0)) THEN
        CALL PCO_ELEMENT_GRAPH(lev, st%cls, eptr, enbr)
        CALL PCO_INSTALL_LAYOUT(lev, nel, eptr, enbr, st%crs, target)
        info(8) = SIZE(pco%colA)
      END IF
    END IF
    CALL MPI_Bcast(info, 8, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_Bcast(xbad, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    IF (ANY(info(1:4).NE.0)) THEN
      WRITE(reason,'(A,I0,A,4(I0,A),I0,A,3ES11.3,A)') 'periodic classes on level ', lev, ': ambiguous ', info(1), &
        ', incomplete ', info(2), ', mixed ', info(3), ', constraint conflicts ', info(4), &
        ' (first: dof ', info(7), ' at', xbad, ')'
      reason = TRIM(reason)//': a missing or ambiguous partner, or unequal conditions on the two faces of a periodic pair'
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

    WRITE(reason,'(A,A,A,A,I0,A,3L1,A,3(I0,1X),A,ES9.2)') 'PCO periodic coarse operator = ', &
      TRIM(PCO_KEY_NAME(target)), MERGE(' (default)  ', ' (requested)', st%requested.EQ.PCO_KEY_UNSET), &
      ' level ', lev, ' axes ', st%census%axis_active, ' shifted links ', st%census%nshift, &
      'largest offset/size', st%census%worst
    CALL say(reason)
    WRITE(reason,'(A,I0,A,I0,A,I0,A,I0,A,I0)') 'PCO classes ', info(5), ' (largest ', info(6), &
      '), layout entries ', info(8), ', uncovered by workers ', st%plan%nuncovered, &
      ', coarse solver type ', st%crs
    CALL say(reason)

  END SUBROUTINE resolve_mode

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
      IF (myid.EQ.master .OR. myid.EQ.showid) WRITE(*,'(A)') 'PCO REFUSED: Pres@PeriodicCoarseOperator = '// &
        TRIM(PCO_KEY_NAME(target))//': '//TRIM(reason)
      IF (myid.EQ.showid) WRITE(MFILE1,'(A)') 'PCO REFUSED: Pres@PeriodicCoarseOperator = '// &
        TRIM(PCO_KEY_NAME(target))//': '//TRIM(reason)
      FLUSH(6)
      STOP 1
    END IF
  END SUBROUTINE fail

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
  SUBROUTINE renew_periodic(u_coarse, w_coarse)
    REAL*8, INTENT(IN), OPTIONAL :: u_coarse(:), w_coarse(:)
    INTEGER :: status, ierr, lev, idum(2)
    CHARACTER(len=256) :: msg
    LOGICAL :: do_check
    REAL*8 :: t0, t1, t2
    DOUBLE PRECISION MPI_WTIME
    EXTERNAL MPI_WTIME

    lev = st%lev
    st%nrenew = st%nrenew + 1
    do_check = (st%check.GT.0) .AND. (st%nrenew.LE.2 .OR. MOD(st%nrenew, MAX(st%check,1)).EQ.0)

    ! a rebuilt worker pattern needs a new plan (validated again; a failure now is a refusal)
    IF (PCO_GATHER_STALE(st%plan)) THEN
      IF (myid.EQ.master) THEN
        CALL PCO_GATHER_BUILD(st%plan, lev, pco%ldA, pco%colA)
      ELSE
        idum = 1
        CALL PCO_GATHER_BUILD(st%plan, lev, idum, idum)
      END IF
      IF (st%plan%nmissing.GT.0 .OR. st%plan%nownership.GT.0 .OR. st%plan%nunowned.GT.0) THEN
        IF (myid.EQ.master .OR. myid.EQ.showid) WRITE(*,'(A,I0)') &
          'PCO REFUSED: the rebuilt worker pressure pattern no longer fits the master layout, renewal ', st%nrenew
        FLUSH(6)
        STOP 1
      END IF
    END IF

    t0 = MPI_WTIME()
    CALL refresh_raw(st%mode, u_coarse, w_coarse)
    IF (st%mode.EQ.PCO_KEY_WRAPPED .AND. st%nrenew.EQ.1) CALL validate_wrapped()
    t1 = MPI_WTIME()
    status = 0
    msg = ' '
    IF (myid.EQ.master) CALL PCO_REFRESH_SOLVER(status, msg)
    t2 = MPI_WTIME()
    CALL MPI_Bcast(status, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    IF (status.NE.0) THEN
      IF (myid.EQ.master) WRITE(*,'(A,I0,A,A)') 'PCO STOP (renewal ', st%nrenew, '): ', TRIM(msg)
      FLUSH(6)
      STOP 1
    END IF
    st%t_values = st%t_values + (t1 - t0)
    st%t_refresh = st%t_refresh + (t2 - t1)

    IF (st%nrenew.EQ.1 .AND. myid.EQ.master) CALL report_gauge()
    IF (do_check) THEN
      t0 = MPI_WTIME()
      CALL run_check(u_coarse, w_coarse)
      st%t_check = st%t_check + (MPI_WTIME() - t0)
    END IF
    IF (myid.EQ.master .AND. (st%nrenew.LE.2 .OR. MOD(st%nrenew,20).EQ.0)) &
      WRITE(MTERM,'(A,I0,A,A,4(A,ES10.3),A,I0,A,ES10.3)') 'PCO_TIMING renewal ', st%nrenew, ' mode ', &
        TRIM(PCO_KEY_NAME(st%mode)), ' setup_s ', st%t_setup, ' values_s ', st%t_values, &
        ' refresh_s ', st%t_refresh, ' check_s ', st%t_check, ' coarse solves ', pco%nsolve, &
        ' largest rhs imbalance ', pco%max_incompat

  END SUBROUTINE renew_periodic

  !-----------------------------------------------------------------------
  ! The wrapped operator is the master assembly of what the workers apply. It is exact only
  ! where the master has the coefficients and the geometry of the workers: on a curved
  ! boundary the master does not project the coordinates of its finest level (the Q2 dofs of
  ! the coarse level, app_init: ParametrizeBndr of level MedLev+1 on the workers only), so its
  ! B and lumped mass differ there (s05 pipe: up to 18 % of an entry). At the first renewal
  ! the wrapped values are therefore compared with the gathered operator entry by entry
  ! (scaled by the largest modulus of the row); a difference above WRAP_TOL refuses the run.
  ! Collective.
  !-----------------------------------------------------------------------
  SUBROUTINE validate_wrapped()
    REAL*8, PARAMETER :: WRAP_TOL = 1d-10
    REAL*8, ALLOCATABLE :: ag(:)
    REAL*8 :: d, rmax
    INTEGER :: row, k, ierr
    IF (myid.EQ.master) THEN
      ALLOCATE(ag(SIZE(pco%a)))
    ELSE
      ALLOCATE(ag(1))
    END IF
    CALL refresh_raw(PCO_KEY_GATHERED, a_out=ag)
    d = 0d0
    IF (myid.EQ.master) THEN
      DO row = 1, pco%n
        rmax = MAXVAL(ABS(ag(pco%ldA(row):pco%ldA(row+1)-1)))
        IF (rmax.LE.0d0) CYCLE
        DO k = pco%ldA(row), pco%ldA(row+1)-1
          d = MAX(d, ABS(ag(k)-pco%a(k))/rmax)
        END DO
      END DO
      WRITE(MTERM,'(A,ES10.3,A,ES8.1,A)') 'PCO wrapped operator vs gathered at the first renewal: ', d, &
        ' (limit ', WRAP_TOL, ', entry difference / largest modulus of the row)'
    END IF
    DEALLOCATE(ag)
    CALL MPI_Bcast(d, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    IF (d.GT.WRAP_TOL) THEN
      IF (myid.EQ.master .OR. myid.EQ.showid) WRITE(*,'(A,ES10.3,A)') 'PCO REFUSED: Pres@PeriodicCoarseOperator'// &
        ' = wrapped: the master assembly differs from the operator of the workers by ', d, &
        ' (curved boundary? the master does not project its finest-level coordinates); use gathered'
      FLUSH(6)
      STOP 1
    END IF
  END SUBROUTINE validate_wrapped

  ! raw values of the active (or another, for the check) operator into pco%a / a
  SUBROUTINE refresh_raw(mode, u_coarse, w_coarse, a_out, mass_variant)
    INTEGER, INTENT(IN) :: mode
    REAL*8, INTENT(IN), OPTIONAL :: u_coarse(:), w_coarse(:)
    REAL*8, INTENT(INOUT), OPTIONAL, TARGET :: a_out(:)
    INTEGER, INTENT(IN), OPTIONAL :: mass_variant
    REAL*8 :: dum(1)
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
      IF (PRESENT(u_coarse)) kind = PCO_MASS_EL_WEIGHTED
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
        IF (PRESENT(a_out)) THEN
          CALL PCO_WRAP_ASSEMBLE(st%wrap, st%cls, kind, dum, u_coarse, w_coarse, &
            QuadSc%knprU(lev)%x, QuadSc%knprV(lev)%x, QuadSc%knprW(lev)%x, a_out)
        ELSE
          CALL PCO_WRAP_ASSEMBLE(st%wrap, st%cls, kind, dum, u_coarse, w_coarse, &
            QuadSc%knprU(lev)%x, QuadSc%knprV(lev)%x, QuadSc%knprW(lev)%x, pco%a)
        END IF
      END IF
    END SELECT
  END SUBROUTINE refresh_raw

  ! gauge decision of the first renewal against SimPar@NoOutflow (master)
  SUBROUTINE report_gauge()
    CHARACTER(len=256) :: line
    WRITE(line,'(A,I0,A,I0,A,ES10.3,A,ES10.3,A,L1)') 'PCO gauge: components ', pco%ncomp, &
      ', closed (C e = 0, one zero reference each) ', pco%nclosed, ', closure measure closed ', &
      pco%closure_closed, ' open ', MERGE(pco%closure_open, 0d0, pco%nclosed.LT.pco%ncomp), &
      ', SimPar@NoOutflow ', bNoOutflow
    WRITE(MTERM,'(A)') TRIM(line)
    IF ((pco%nclosed.GT.0) .NEQV. bNoOutflow) WRITE(MTERM,'(A)') &
      'PCO gauge note: the operator decides the reference (closed components), SimPar@NoOutflow disagrees'
  END SUBROUTINE report_gauge

  !=======================================================================
  ! The check (module header). Collective.
  !=======================================================================
  SUBROUTINE run_check(u_coarse, w_coarse)
    REAL*8, INTENT(IN), OPTIONAL :: u_coarse(:), w_coarse(:)
    REAL*8, ALLOCATABLE :: ag(:), aw(:), ai(:), v(:), yd(:), yg(:), yw(:), yi(:), sc(:)
    INTEGER :: lev, n, nnz, k, row, iv, e
    REAL*8 :: dgw, dgi, rmax, eg, ew, ei, resid, devi, cnorm, tol_res, tol_dev, cleft, xref
    INTEGER :: kmax, rowk
    REAL*8, ALLOCATABLE :: cs(:), ca(:)
    LOGICAL :: weighted, incompat
    CHARACTER(len=8) :: vname(3) = (/'mean    ', 'impulse ', 'random  '/)

    lev = st%lev
    weighted = PRESENT(u_coarse)
    IF (myid.EQ.master) THEN
      n = pco%n
      nnz = SIZE(pco%colA)
      ALLOCATE(ag(nnz), aw(nnz), ai(nnz), v(n), yd(n), yg(n), yw(n), yi(n), sc(n))
    ELSE
      n = 4*mg_mesh%level(lev)%nel
      ALLOCATE(ag(1), aw(1), ai(1), v(n), yd(n), yg(1), yw(1), yi(1), sc(1))
    END IF
    CALL refresh_raw(PCO_KEY_GATHERED, a_out=ag)
    IF (myid.EQ.master) THEN
      CALL refresh_raw(PCO_KEY_WRAPPED, u_coarse, w_coarse, a_out=aw)
      IF (weighted) THEN
        CALL refresh_raw(PCO_KEY_WRAPPED, u_coarse, w_coarse, a_out=ai, mass_variant=PCO_MASS_EL_PER_IMAGE)
      ELSE
        ai = aw
      END IF
      ! entry-by-entry, scaled by the largest modulus of the row (gathered)
      dgw = 0d0
      dgi = 0d0
      kmax = 0
      rowk = 0
      DO row = 1, n
        rmax = MAXVAL(ABS(ag(pco%ldA(row):pco%ldA(row+1)-1)))
        IF (rmax.LE.0d0) CYCLE
        DO k = pco%ldA(row), pco%ldA(row+1)-1
          IF (ABS(ag(k)-aw(k))/rmax.GT.dgw) THEN
            dgw = ABS(ag(k)-aw(k))/rmax
            kmax = k
            rowk = row
          END IF
          dgi = MAX(dgi, ABS(ag(k)-ai(k))/rmax)
        END DO
      END DO
      WRITE(MTERM,'(A,I0,A,I0,A,I0,A,ES10.3,A,ES10.3,A,A)') 'PCO_CHECK renewal ', st%nrenew, ' level ', lev, &
        ' entries ', nnz, ' |gathered-wrapped|/rowmax ', dgw, ' per-image-mass variant ', dgi, &
        ' ', MERGE('PASS', 'FAIL', dgw.LE.1d-12)
      ! where: row and column unknown (element, mode), element centroids, both values
      IF (kmax.GT.0) WRITE(MTERM,'(A,I0,A,2(I0,A,I0,A,3ES11.3,A),3(A,ES12.4))') 'PCO_CHECK renewal ', st%nrenew, &
        ' largest difference at ', (rowk-1)/4+1, ' mode ', MOD(rowk-1,4)+1, ' (', centroid((rowk-1)/4+1), ')', &
        (pco%colA(kmax)-1)/4+1, ' mode ', MOD(pco%colA(kmax)-1,4)+1, ' (', centroid((pco%colA(kmax)-1)/4+1), ')', &
        ' gathered ', ag(kmax), ' wrapped ', aw(kmax), ' legacy master C ', legacy_entry(rowk, pco%colA(kmax))
    END IF

    DO iv = 1, 3
      IF (myid.EQ.master) CALL test_vector(iv, v)
      CALL worker_product(lev, v, yd)
      IF (myid.EQ.master) THEN
        CALL matvec(ag, v, yg)
        CALL matvec(aw, v, yw)
        CALL matvec(ai, v, yi)
        ! normwise row scale sum_j |A_ij| max|v| (a product scale would divide the rounding
        ! of entries that cancel to nearly zero by nearly zero)
        DO row = 1, n
          sc(row) = 0d0
          DO k = pco%ldA(row), pco%ldA(row+1)-1
            sc(row) = sc(row) + ABS(ag(k))
          END DO
          sc(row) = sc(row)*MAXVAL(ABS(v))
        END DO
        eg = 0d0
        ew = 0d0
        ei = 0d0
        DO row = 1, n
          IF (sc(row).LE.0d0) CYCLE
          eg = MAX(eg, ABS(yd(row)-yg(row))/sc(row))
          ew = MAX(ew, ABS(yd(row)-yw(row))/sc(row))
          ei = MAX(ei, ABS(yd(row)-yi(row))/sc(row))
        END DO
        cnorm = MAXVAL(ABS(yd)/MAX(sc, 1d-300))
        WRITE(MTERM,'(A,I0,A,A,A,I0,3(A,ES10.3),A,ES10.3,A,A)') 'PCO_CHECK renewal ', st%nrenew, ' vector ', &
          TRIM(vname(iv)), ' rows ', n, ' dax-gathered ', eg, ' dax-wrapped ', ew, ' dax-per-image ', ei, &
          ' max |Av|/(|A| |v|) ', cnorm, ' ', MERGE('PASS', 'FAIL', eg.LE.1d-12 .AND. ew.LE.1d-12)
      END IF
    END DO

    ! solver: compatible right-hand side, ungauged residual; and an incompatible one
    IF (myid.EQ.master) THEN
      ! direct solvers: rounding level; BiCGStab (types 1, 3): its tolerance in PCO_SOLVE_CHECK
      tol_res = MERGE(1d-10, 1d-8, st%crs.EQ.2 .OR. st%crs.EQ.4)
      ! solution error = conditioning x residual: direct solvers rounding level, BiCGStab only
      ! to the iteration tolerance amplified by the condition number of the coarse matrix
      tol_dev = MERGE(1d-8, 1d-3, st%crs.EQ.2 .OR. st%crs.EQ.4)
      ! left null vector of a closed component: e^T A = 0 (column sums of the mean rows),
      ! the compatibility condition the gauge relies on
      ALLOCATE(cs(n), ca(n))
      cs = 0d0
      ca = 0d0
      DO row = 1, n, 4
        DO k = pco%ldA(row), pco%ldA(row+1)-1
          cs(pco%colA(k)) = cs(pco%colA(k)) + ag(k)
          ca(pco%colA(k)) = ca(pco%colA(k)) + ABS(ag(k))
        END DO
      END DO
      cleft = 0d0
      DO k = 1, n
        IF (pco%comp_closed(pco%comp_of((k-1)/4+1)) .AND. ca(k).GT.0d0) cleft = MAX(cleft, ABS(cs(k))/ca(k))
      END DO
      DEALLOCATE(cs, ca)
      CALL test_vector(3, v)
      CALL PCO_SOLVE_CHECK(v, resid, devi, incompat)
      xref = 0d0
      IF (ANY(pco%sref)) xref = MAXVAL(ABS(pco%sol), MASK=pco%sref)/MAX(MAXVAL(ABS(pco%sol)), 1d-300)
      WRITE(MTERM,'(A,I0,A,I0,4(A,ES10.3),A,L1,A,A)') 'PCO_CHECK renewal ', st%nrenew, ' solver type ', &
        st%crs, ' left |e^T A|/(e^T |A|) ', cleft, ' compatible rhs: ungauged residual ', resid, &
        ' solution - v off the kernel ', devi, ' reference value ', xref, &
        ' incompatibility reported ', incompat, ' ', MERGE('PASS', 'FAIL', resid.LE.tol_res .AND. &
        devi.LE.tol_dev .AND. .NOT.incompat .AND. cleft.LE.1d-12 .AND. xref.LE.tol_res)
      IF (pco%nclosed.GT.0) THEN
        ! v = mean modes of the elements of the first closed component, b = raw v + that constant:
        ! the constant is not in the range of a closed operator, the solve must report it
        v = 0d0
        DO e = 1, pco%nel
          v(4*(e-1)+1) = 1d0
        END DO
        CALL incompatible_probe(v, incompat)
        WRITE(MTERM,'(A,I0,A,L1,A,A)') 'PCO_CHECK renewal ', st%nrenew, ' incompatible rhs reported ', &
          incompat, ' ', MERGE('PASS', 'FAIL', incompat)
      END IF
    END IF
    DEALLOCATE(ag, aw, ai, v, yd, yg, yw, yi, sc)

  CONTAINS

    ! value of the own (legacy, non-periodic) matrix of the master at (r, c), 0 if not in its pattern
    REAL*8 FUNCTION legacy_entry(r, c) RESULT(val)
      INTEGER, INTENT(IN) :: r, c
      INTEGER :: kk
      val = 0d0
      DO kk = mg_lMat(lev)%LdA(r), mg_lMat(lev)%LdA(r+1)-1
        IF (mg_lMat(lev)%ColA(kk).EQ.c) val = mg_CMat(lev)%a(kk)
      END DO
    END FUNCTION legacy_entry

    FUNCTION centroid(e0) RESULT(c)
      INTEGER, INTENT(IN) :: e0
      REAL*8 :: c(3)
      INTEGER :: j
      c = 0d0
      DO j = 1, 8
        c = c + 0.125d0*mg_mesh%level(lev)%dcorvg(:, mg_mesh%level(lev)%kvert(j,e0))
      END DO
    END FUNCTION centroid

    SUBROUTINE matvec(a, x, y)
      REAL*8, INTENT(IN) :: a(:), x(:)
      REAL*8, INTENT(OUT) :: y(:)
      INTEGER :: i, kk
      DO i = 1, pco%n
        y(i) = 0d0
        DO kk = pco%ldA(i), pco%ldA(i+1)-1
          y(i) = y(i) + a(kk)*x(pco%colA(kk))
        END DO
      END DO
    END SUBROUTINE matvec

  END SUBROUTINE run_check

  ! the three test vectors in the master numbering
  SUBROUTINE test_vector(iv, v)
    INTEGER, INTENT(IN) :: iv
    REAL*8, INTENT(OUT) :: v(:)
    INTEGER :: e, k, best, msz, i
    INTEGER(8) :: s
    v = 0d0
    SELECT CASE (iv)
    CASE (1)
      DO e = 1, pco%nel
        v(4*(e-1)+1) = 1d0
      END DO
    CASE (2)
      ! element of the first member of the largest vertex class
      best = 1
      msz = 0
      DO k = 1, st%cls%ncls
        IF (st%cls%mem(st%cls%ptr(k)).GT.st%cls%nvt) CYCLE
        IF (st%cls%ptr(k+1)-st%cls%ptr(k).GT.msz) THEN
          msz = st%cls%ptr(k+1)-st%cls%ptr(k)
          best = st%cls%mem(st%cls%ptr(k))
        END IF
      END DO
      e = 1
      DO i = 1, mg_mesh%level(st%lev)%nel
        IF (ANY(mg_mesh%level(st%lev)%kvert(:,i).EQ.best)) THEN
          e = i
          EXIT
        END IF
      END DO
      v(4*(e-1)+1:4*e) = (/1d0, 0.5d0, -0.25d0, 0.75d0/)
    CASE DEFAULT
      s = 20261119_8
      DO i = 1, SIZE(v)
        s = MOD(1103515245_8*s + 12345_8, 2147483648_8)
        v(i) = 2d0*DBLE(s)/2147483648d0 - 1d0
      END DO
    END SELECT
  END SUBROUTINE test_vector

  ! y = worker operator applied to v (master numbering on the master; collective):
  ! distribute v to the workers, E012_DAX arithmetic there, collect the result
  SUBROUTINE worker_product(lev, v, y)
    INTEGER, INTENT(IN) :: lev
    REAL*8, INTENT(INOUT) :: v(:), y(:)
    INTEGER :: nel, i, k, ilev_save
    REAL*8, ALLOCATABLE :: xp(:)

    ilev_save = ILEV
    ILEV = lev
    nel = mg_mesh%level(lev)%nel
    CALL coarse_to_workers(lev, v, nel)
    IF (myid.NE.master) THEN
      ALLOCATE(xp(PCO_GHOST_BUFFER_SIZE(lev)))
      xp = 0d0
    ELSE
      ALLOCATE(xp(4))
    END IF
    CALL GetParPressure(v, xp)
    IF (myid.NE.master) THEN
      DO i = 1, 4*nel
        y(i) = 0d0
        DO k = mg_lMat(lev)%LdA(i), mg_lMat(lev)%LdA(i+1)-1
          y(i) = y(i) + mg_CMat(lev)%a(k)*v(mg_lMat(lev)%ColA(k))
        END DO
        DO k = mg_lPMat(lev)%LdA(i), mg_lPMat(lev)%LdA(i+1)-1
          y(i) = y(i) + mg_CPMat(lev)%a(k)*xp(mg_lPMat(lev)%ColA(k))
        END DO
      END DO
    END IF
    DEALLOCATE(xp)
    CALL workers_to_coarse(lev, y, nel)
    ILEV = ilev_save
  END SUBROUTINE worker_product

  ! the coarse vector transfers of the multigrid (QuadSc_mpi.f90), by level
  SUBROUTINE coarse_to_workers(lev, d, nel)
    INTEGER, INTENT(IN) :: lev, nel
    REAL*8, INTENT(INOUT) :: d(:)
    EXTERNAL E012GATHR_L1, E012GATHR_L2, E012GATHR_L3
    IF (lev.EQ.1) CALL E012GATHR_L1(d, nel)
    IF (lev.EQ.2) CALL E012GATHR_L2(d, nel)
    IF (lev.EQ.3) CALL E012GATHR_L3(d, nel)
  END SUBROUTINE coarse_to_workers

  SUBROUTINE workers_to_coarse(lev, d, nel)
    INTEGER, INTENT(IN) :: lev, nel
    REAL*8, INTENT(INOUT) :: d(:)
    EXTERNAL E012DISTR_L1, E012DISTR_L2, E012DISTR_L3
    IF (lev.EQ.1) CALL E012DISTR_L1(d, nel)
    IF (lev.EQ.2) CALL E012DISTR_L2(d, nel)
    IF (lev.EQ.3) CALL E012DISTR_L3(d, nel)
  END SUBROUTINE workers_to_coarse

  ! a right-hand side with a constant added on the closed components must be reported
  SUBROUTINE incompatible_probe(c, seen)
    REAL*8, INTENT(IN) :: c(:)
    LOGICAL, INTENT(OUT) :: seen
    REAL*8, ALLOCATABLE :: b(:), x(:)
    INTEGER :: nin, it
    nin = pco%nincompat
    ALLOCATE(b(pco%n), x(pco%n))
    CALL test_vector(3, x)
    CALL PCO_RAW_MATVEC(x, b)
    b = b + c*MAXVAL(ABS(b))
    it = 2000
    IF (pco%reduced) THEN
      CALL solve_p0_probe(b)
    ELSE
      CALL solve_full_probe(b)
    END IF
    seen = (pco%nincompat.GT.nin)
    DEALLOCATE(b, x)
  CONTAINS
    SUBROUTINE solve_full_probe(bb)
      USE pco_coarse_solver, ONLY: PCO_SOLVE_FULL
      REAL*8, INTENT(IN) :: bb(:)
      REAL*8, ALLOCATABLE :: xx(:)
      ALLOCATE(xx(pco%n))
      CALL PCO_SOLVE_FULL(xx, bb, it, 1d-14, 1d0)
      DEALLOCATE(xx)
    END SUBROUTINE solve_full_probe
    SUBROUTINE solve_p0_probe(bb)
      USE pco_coarse_solver, ONLY: PCO_SOLVE_P0
      REAL*8, INTENT(IN) :: bb(:)
      REAL*8, ALLOCATABLE :: xx(:)
      ALLOCATE(xx(pco%nel))
      CALL PCO_SOLVE_P0(bb, xx, it, 1d-14, 1d0)
      DEALLOCATE(xx)
    END SUBROUTINE solve_p0_probe
  END SUBROUTINE incompatible_probe

END MODULE pco_coarse_operator
