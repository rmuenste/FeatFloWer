!=========================================================================
! CHIMERA_API - the single facade of the Chimera overlapping-mesh
! component (design: chimera-integration-design.md, v3, section 2).
!
! This is the ONLY Chimera module that existing solver code is allowed
! to USE.  All component state stays private behind it (CHI_COUPLING).
! Every operation is internally rank-safe: callers never need myid
! guards (the one documented exception is the H1 call-site placement,
! which is dictated by argument association, not by this module).
!
! Contract for the hooks in existing code:
!   - the first statement of every operation is the master-switch test,
!     so a disabled run (the default) takes exactly one extra branch and
!     no floating-point path changes (off-regression bit identity);
!   - enabled-but-uninitialized: the boundary hooks are no-ops (they are
!     legitimately reached during application initialization, before
!     Chimera_Initialize), whereas Chimera_BeginStep aborts - it is the
!     first per-step call and catches an application that enabled the
!     component without initializing it.
!=========================================================================
MODULE CHIMERA_API

  USE CHIMERA_CONFIG, ONLY: chimera_enable, bChimeraW, &
    CHIMERA_VALIDATE_CONFIG
  USE CHI_COUPLING, ONLY: CHI_COUPLING_INIT, CHI_COUPLING_BEGIN_STEP, &
    CHI_COUPLING_APPLY_DEF, CHI_COUPLING_APPLY_VAL, CHI_COUPLING_FILTER_MAT, &
    CHI_COUPLING_FILTER_MAT_9, CHI_COUPLING_FINALIZE

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: Chimera_IsEnabled
  PUBLIC :: Chimera_VariantIsWeak
  PUBLIC :: Chimera_Initialize
  PUBLIC :: Chimera_BeginStep
  PUBLIC :: Chimera_ApplyBoundaryDef
  PUBLIC :: Chimera_ApplyBoundaryValues
  PUBLIC :: Chimera_FilterMatrixRows
  PUBLIC :: Chimera_FilterMatrixRows9
  PUBLIC :: Chimera_Finalize

  ! Lifecycle state of the component (private; set by Chimera_Initialize,
  ! cleared by Chimera_Finalize).
  LOGICAL :: chi_initialized = .FALSE.

CONTAINS

  !-----------------------------------------------------------------------
  ! Master switch, read by every hook in existing solver code.
  !-----------------------------------------------------------------------
  LOGICAL FUNCTION Chimera_IsEnabled()
    Chimera_IsEnabled = chimera_enable
  END FUNCTION Chimera_IsEnabled

  !-----------------------------------------------------------------------
  ! True when the weak (interior-penalty) variant is active.
  !-----------------------------------------------------------------------
  LOGICAL FUNCTION Chimera_VariantIsWeak()
    Chimera_VariantIsWeak = chimera_enable .AND. bChimeraW
  END FUNCTION Chimera_VariantIsWeak

  !-----------------------------------------------------------------------
  ! Application-local initialization (hook H5, called by q2p1_chimera
  ! after init_q2p1_app; paired with Chimera_Finalize).  No-op when the
  ! component is disabled.  Milestone-1 restrictions are enforced here:
  ! the strong variant only (the weak variant arrives with Phase 4), and
  ! no coexistence with the FBM particle mode.
  !-----------------------------------------------------------------------
  SUBROUTINE Chimera_Initialize(mfile)
    USE var_QuadScalar, ONLY: myFBM
    INTEGER, INTENT(IN) :: mfile

    IF (.NOT. chimera_enable) RETURN

    CALL CHIMERA_VALIDATE_CONFIG()

    IF (bChimeraW) THEN
      WRITE(*,'(A)') 'CHIMERA_API error: ChimeraVariant = weak is not ' // &
        'implemented yet (arrives with design phase 4).'
      STOP 1
    END IF
    IF (myFBM%nParticles .GT. 0) THEN
      WRITE(*,'(A)') 'CHIMERA_API error: ChimeraEnable = Yes cannot be ' // &
        'combined with FBM particles (milestone-1 restriction).'
      STOP 1
    END IF

    CALL CHI_COUPLING_INIT(mfile)
    chi_initialized = .TRUE.
  END SUBROUTINE Chimera_Initialize

  !-----------------------------------------------------------------------
  ! Per-time-step coupling update (hook H1 in
  ! Transport_q2p1_UxyzP_fluid_core).  No-op when disabled; aborts when
  ! enabled without initialization - i.e. when an application other than
  ! one that calls Chimera_Initialize runs with ChimeraEnable = Yes.
  !-----------------------------------------------------------------------
  SUBROUTINE Chimera_BeginStep(valU, valV, valW, valP)
    REAL*8, INTENT(IN) :: valU(*), valV(*), valW(*), valP(*)

    IF (.NOT. chimera_enable) RETURN

    IF (.NOT. chi_initialized) THEN
      WRITE(*,'(A)') 'CHIMERA_API error: SimPar@ChimeraEnable = Yes, but ' // &
        'Chimera was never initialized.'
      WRITE(*,'(A)') '  ChimeraEnable requires an application that calls ' // &
        'Chimera_Initialize (e.g. q2p1_chimera).'
      STOP 1
    END IF

    CALL CHI_COUPLING_BEGIN_STEP(valU, valV, valW, valP)
  END SUBROUTINE Chimera_BeginStep

  !-----------------------------------------------------------------------
  ! Hook H2: zero the momentum defect at hole/fringe dofs (sibling of the
  ! FictKNPR branch in Boundary_QuadScalar_Def).
  !-----------------------------------------------------------------------
  SUBROUTINE Chimera_ApplyBoundaryDef(defU, defV, defW, ndof)
    REAL*8, INTENT(INOUT) :: defU(*), defV(*), defW(*)
    INTEGER, INTENT(IN) :: ndof
    IF (.NOT. chimera_enable) RETURN
    IF (.NOT. chi_initialized) RETURN
    CALL CHI_COUPLING_APPLY_DEF(defU, defV, defW, ndof)
  END SUBROUTINE Chimera_ApplyBoundaryDef

  !-----------------------------------------------------------------------
  ! Hook H3: impose hole (rigid-body) and fringe (submesh) velocities.
  !-----------------------------------------------------------------------
  SUBROUTINE Chimera_ApplyBoundaryValues(valU, valV, valW, ndof)
    REAL*8, INTENT(INOUT) :: valU(*), valV(*), valW(*)
    INTEGER, INTENT(IN) :: ndof
    IF (.NOT. chimera_enable) RETURN
    IF (.NOT. chi_initialized) RETURN
    CALL CHI_COUPLING_APPLY_VAL(valU, valV, valW, ndof)
  END SUBROUTINE Chimera_ApplyBoundaryValues

  !-----------------------------------------------------------------------
  ! Hook H4: Dirichlet row filter of the (block-diagonal) momentum matrix.
  !-----------------------------------------------------------------------
  SUBROUTINE Chimera_FilterMatrixRows(DA11, DA22, DA33, KLD, ndof)
    REAL*8, INTENT(INOUT) :: DA11(*), DA22(*), DA33(*)
    INTEGER, INTENT(IN) :: KLD(*), ndof
    IF (.NOT. chimera_enable) RETURN
    IF (.NOT. chi_initialized) RETURN
    CALL CHI_COUPLING_FILTER_MAT(DA11, DA22, DA33, KLD, ndof)
  END SUBROUTINE Chimera_FilterMatrixRows

  !-----------------------------------------------------------------------
  ! Hook H4 (9-block variant).
  !-----------------------------------------------------------------------
  SUBROUTINE Chimera_FilterMatrixRows9(DA11, DA22, DA33, DA12, DA13, DA23, &
                                       DA21, DA31, DA32, KLD, ndof)
    REAL*8, INTENT(INOUT) :: DA11(*), DA22(*), DA33(*), DA12(*), DA13(*), &
                             DA23(*), DA21(*), DA31(*), DA32(*)
    INTEGER, INTENT(IN) :: KLD(*), ndof
    IF (.NOT. chimera_enable) RETURN
    IF (.NOT. chi_initialized) RETURN
    CALL CHI_COUPLING_FILTER_MAT_9(DA11, DA22, DA33, DA12, DA13, DA23, &
                                   DA21, DA31, DA32, KLD, ndof)
  END SUBROUTINE Chimera_FilterMatrixRows9

  !-----------------------------------------------------------------------
  ! Application-local finalization (hook H6).  Idempotent and safe on
  ! partial initialization.
  !-----------------------------------------------------------------------
  SUBROUTINE Chimera_Finalize()
    IF (.NOT. chi_initialized) RETURN
    CALL CHI_COUPLING_FINALIZE()
    chi_initialized = .FALSE.
  END SUBROUTINE Chimera_Finalize

END MODULE CHIMERA_API
