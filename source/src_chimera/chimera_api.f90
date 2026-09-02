!=========================================================================
! CHIMERA_API - the single facade of the Chimera overlapping-mesh
! component (design: chimera-integration-design.md, v3, section 2).
!
! This is the ONLY Chimera module that existing solver code is allowed
! to USE.  All component state stays private behind it.  Every operation
! is internally rank-safe: callers never need myid guards (the one
! documented exception is the H1 call-site placement, which is dictated
! by argument association, not by this module).
!
! Phase 0/1 status: the facade carries the master switch and the
! lifecycle contract.  Operations whose implementation arrives in later
! phases abort with a clear message instead of silently doing nothing -
! an enabled-but-unimplemented Chimera run must never masquerade as a
! plain flow solve (no print-and-continue stubs).
!=========================================================================
MODULE CHIMERA_API

  USE CHIMERA_CONFIG, ONLY: chimera_enable, bChimeraW, &
    CHIMERA_VALIDATE_CONFIG

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: Chimera_IsEnabled
  PUBLIC :: Chimera_VariantIsWeak
  PUBLIC :: Chimera_Initialize
  PUBLIC :: Chimera_BeginStep
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
  ! component is disabled.
  !
  ! Phase 3 will load submeshes, build the locator and donor caches and
  ! classify markers here.  Until then an enabled run aborts loudly.
  !-----------------------------------------------------------------------
  SUBROUTINE Chimera_Initialize(mfile)
    INTEGER, INTENT(IN) :: mfile

    IF (.NOT. chimera_enable) RETURN

    CALL CHIMERA_VALIDATE_CONFIG()

    WRITE(*,'(A)') 'CHIMERA_API error: Chimera_Initialize is not ' // &
      'implemented yet (arrives with design phase 3).'
    WRITE(*,'(A)') '  Set SimPar@ChimeraEnable = No.'
    STOP 1

    ! Unreachable until phase 3:
    ! chi_initialized = .TRUE.
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

    ! Phase 3: gather Robin data, solve submeshes, broadcast, update
    ! fringe values, integrate forces.
  END SUBROUTINE Chimera_BeginStep

  !-----------------------------------------------------------------------
  ! Application-local finalization (hook H6).  Idempotent and safe on
  ! partial initialization.
  !-----------------------------------------------------------------------
  SUBROUTINE Chimera_Finalize()
    IF (.NOT. chi_initialized) RETURN
    ! Phase 3: release submeshes, locator, caches, solver handles.
    chi_initialized = .FALSE.
  END SUBROUTINE Chimera_Finalize

END MODULE CHIMERA_API
