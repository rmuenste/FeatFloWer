!=========================================================================
! CHIMERA_CONFIG - runtime configuration for the Chimera overlapping-mesh
! component (design: chimera-integration-design.md, v3).
!
! Layer L: this module lives in the ff_util library so that the central
! parameter parser (source/src_util/param_parser.f90) can write into it
! directly - the same pattern el_config.f90 uses.  It has NO dependencies
! on solver state; every key defaults to the value that keeps the
! standard operational mode untouched (chimera_enable = .FALSE.).
!
! The component is always compiled (no preprocessor gating); it is
! activated purely at runtime via SimPar@ChimeraEnable.  Unknown keys in
! a deck are ignored by the parser, so decks remain portable.
!=========================================================================
MODULE CHIMERA_CONFIG

  IMPLICIT NONE

  PRIVATE

  !-----------------------------------------------------------------------
  ! SimPar@ChimeraEnable (Yes/No, default No)
  ! Master switch. Every hook in existing solver code guards on this flag
  ! (via CHIMERA_API); with the default value the binary is bit-identical
  ! to a build without the component.
  !-----------------------------------------------------------------------
  LOGICAL, PUBLIC :: chimera_enable = .FALSE.

  !-----------------------------------------------------------------------
  ! SimPar@ChimeraVariant ('strong' | 'weak', default 'strong')
  ! 'strong' = Chimera-S: hole/fringe Dirichlet constraints on the
  !            background mesh (paper Section 6).
  ! 'weak'   = Chimera-W: distributed interior-penalty operator D and
  !            vector g (paper eqs. (7)-(12)).
  !-----------------------------------------------------------------------
  CHARACTER(LEN=8),  PUBLIC :: chimera_variant = 'strong'

  !-----------------------------------------------------------------------
  ! SimPar@ChimeraOuterBC ('robin' | 'dirichlet', default 'robin')
  ! Outer (atmosphere) boundary condition of the submesh problems.
  ! 'robin' is the paper's eq. (5d) and the acceptance path; 'dirichlet'
  ! (interpolated background velocity, with a pressure gauge) is kept as
  ! a diagnostic comparison mode only.
  !-----------------------------------------------------------------------
  CHARACTER(LEN=16), PUBLIC :: chimera_outer_bc = 'robin'

  !-----------------------------------------------------------------------
  ! SimPar@ChimeraParticleFile (path, default empty)
  ! Table of particle centers / radii / atmosphere widths H_k.
  !-----------------------------------------------------------------------
  CHARACTER(LEN=256), PUBLIC :: chimera_particle_file = ''

  !-----------------------------------------------------------------------
  ! SimPar@ChimeraSubmeshFile (path, default empty)
  ! Coarse body-fitted shell mesh (.tri) instantiated per particle by
  ! scale/translate (tools/chimera_meshgen).
  !-----------------------------------------------------------------------
  CHARACTER(LEN=256), PUBLIC :: chimera_submesh_file = ''

  !-----------------------------------------------------------------------
  ! SimPar@ChimeraSubmeshLev (integer >= 1, default 3)
  ! Number of refinement levels of each submesh hierarchy.
  !-----------------------------------------------------------------------
  INTEGER, PUBLIC :: chimera_submesh_nlmax = 3

  !-----------------------------------------------------------------------
  ! SimPar@ChimeraRobinAlpha (real >= 0, default 1.0)
  ! Robin coefficient alpha of paper eq. (5d).
  !-----------------------------------------------------------------------
  REAL*8, PUBLIC :: chimera_robin_alpha = 1d0

  !-----------------------------------------------------------------------
  ! SimPar@ChimeraGammaMax (real >= 0, default 0.0)
  ! Interior-penalty parameter gamma_max of paper eq. (7).  Required > 0
  ! for the weak variant.
  !-----------------------------------------------------------------------
  REAL*8, PUBLIC :: chimera_gamma_max = 0d0

  !-----------------------------------------------------------------------
  ! SimPar@ChimeraOuterIters (integer >= 1, default 1)
  ! In-step outer coupling iterations (hook H13, milestone M4).  The
  ! default of 1 keeps the solver control flow byte-identical.
  !-----------------------------------------------------------------------
  INTEGER, PUBLIC :: chimera_outer_iters = 1

  !-----------------------------------------------------------------------
  ! SimPar@ChimeraSubNL (integer >= 1, default 3)
  ! Picard iterations of each submesh saddle-point solve.
  !-----------------------------------------------------------------------
  INTEGER, PUBLIC :: chimera_sub_nl = 3

  !-----------------------------------------------------------------------
  ! SimPar@ChimeraWriteVTK (Yes/No, default No)
  ! Dump submesh solutions for visualization.
  !-----------------------------------------------------------------------
  LOGICAL, PUBLIC :: chimera_write_vtk = .FALSE.

  ! Derived flags, set by CHIMERA_VALIDATE_CONFIG.
  LOGICAL, PUBLIC :: bChimeraS = .FALSE.
  LOGICAL, PUBLIC :: bChimeraW = .FALSE.

  PUBLIC :: CHIMERA_VALIDATE_CONFIG

CONTAINS

  !-----------------------------------------------------------------------
  ! Validate the parsed configuration.  Called by the parameter parser
  ! once the SimPar section has been read and chimera_enable is set, and
  ! again defensively from Chimera_Initialize.  Aborts with a clear
  ! message on any violation (never print-and-continue).
  !-----------------------------------------------------------------------
  SUBROUTINE CHIMERA_VALIDATE_CONFIG()

    IF (.NOT. chimera_enable) THEN
      bChimeraS = .FALSE.
      bChimeraW = .FALSE.
      RETURN
    END IF

    CALL chimera_lowercase(chimera_variant)
    CALL chimera_lowercase(chimera_outer_bc)

    SELECT CASE (TRIM(chimera_variant))
    CASE ('strong')
      bChimeraS = .TRUE.
      bChimeraW = .FALSE.
    CASE ('weak')
      bChimeraS = .FALSE.
      bChimeraW = .TRUE.
    CASE DEFAULT
      WRITE(*,'(A,A)') 'CHIMERA_CONFIG error: invalid ChimeraVariant: ', &
        TRIM(chimera_variant)
      WRITE(*,'(A)') '  valid values: strong | weak'
      STOP 1
    END SELECT

    SELECT CASE (TRIM(chimera_outer_bc))
    CASE ('robin', 'dirichlet')
      ! valid
    CASE DEFAULT
      WRITE(*,'(A,A)') 'CHIMERA_CONFIG error: invalid ChimeraOuterBC: ', &
        TRIM(chimera_outer_bc)
      WRITE(*,'(A)') '  valid values: robin | dirichlet'
      STOP 1
    END SELECT

    IF (bChimeraW .AND. chimera_gamma_max .LE. 0d0) THEN
      WRITE(*,'(A)') 'CHIMERA_CONFIG error: ChimeraVariant = weak requires ' // &
        'ChimeraGammaMax > 0'
      STOP 1
    END IF

    IF (chimera_robin_alpha .LT. 0d0) THEN
      WRITE(*,'(A)') 'CHIMERA_CONFIG error: ChimeraRobinAlpha must be >= 0'
      STOP 1
    END IF

    IF (chimera_submesh_nlmax .LT. 1) THEN
      WRITE(*,'(A)') 'CHIMERA_CONFIG error: ChimeraSubmeshLev must be >= 1'
      STOP 1
    END IF

    IF (chimera_outer_iters .LT. 1) THEN
      WRITE(*,'(A)') 'CHIMERA_CONFIG error: ChimeraOuterIters must be >= 1'
      STOP 1
    END IF

    IF (chimera_sub_nl .LT. 1) THEN
      WRITE(*,'(A)') 'CHIMERA_CONFIG error: ChimeraSubNL must be >= 1'
      STOP 1
    END IF

    ! NOTE: the milestone-1 cross-check that Chimera and the FBM particle
    ! mode are mutually exclusive needs runtime context (particle counts,
    ! calculateFBM state) and therefore runs in Chimera_Initialize, not
    ! here.

  END SUBROUTINE CHIMERA_VALIDATE_CONFIG

  SUBROUTINE chimera_lowercase(s)
    CHARACTER(*), INTENT(INOUT) :: s
    INTEGER :: i
    CHARACTER :: ch
    DO i = 1, LEN_TRIM(s)
      ch = s(i:i)
      IF (ch .GE. 'A' .AND. ch .LE. 'Z') s(i:i) = CHAR(IACHAR(ch) + 32)
    END DO
  END SUBROUTINE chimera_lowercase

END MODULE CHIMERA_CONFIG
