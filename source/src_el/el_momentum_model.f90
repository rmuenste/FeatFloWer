!=========================================================================
! EL_MOMENTUM_MODEL - momentum formulation of the Euler-Lagrange fluid-particle set
!
! SimPar@ELMomentumModel selects one of the CFDEM model types (EL_CONFIG: key, table of names,
! el_mm_implemented). This module owns everything a model changes and answers it through a few
! questions, so that no other code tests the model name:
!
!   EL_MM_ACTIVE()     .TRUE. when a model other than 'none' is selected (deck-only, valid on
!                      every rank including particle-free workers and the master).
!
!   EL_MM_WEIGHTED_INERTIA()
!                      .TRUE. when the fluid inertia is eps rho Du/Dt (all three models). The
!                      mesh and operator work is done by the fluid adapter in Transport_Q2P1
!                      (QuadSc_el_mm.f90), which asks this function only.
!   EL_MM_DRAG_SCALE(eps)
!                      factor on the drag force and its coefficient: 1/eps for model B (the
!                      whole fluid-particle interaction is carried by the drag there), 1
!                      otherwise. eps is the bounded kernel sample the closure itself uses.
!   EL_MM_RESIDUAL_DIVISOR(eps_e)
!                      divisor that turns the conservative continuity residual r_m of an
!                      element into the projection right-hand side -r_m/(divisor dt). It is the
!                      preconditioning of the residual iteration by the existing pressure
!                      matrix C0 = B^T (rho Ml)^-1 B: a pressure increment q changes the flux
!                      eps u by about eps dt (m Ml)^-1 G q with the inertia coefficient m and
!                      the pressure-force operator G of the model,
!                        none   m = rho,     G = B      -> eps C0 q   divisor eps_e
!                        bfull  m = eps rho, G = B      -> C0 q       divisor 1
!                        b      m = eps rho, G = B      -> C0 q       divisor 1
!                        a      m = eps rho, G = eps B  -> eps C0 q   divisor eps_e
!                      The cancellation is exact only for uniform porosity (lumped P0-weighted
!                      mass against the Q1 transport porosity of the flux): C0 is a
!                      preconditioner, the fixed point r = 0 does not depend on it, and the
!                      contraction is measured (EL_VA_FLUX res1_L1 per corrector), not assumed.
!
!   EL_MM_VISCOUS_FORCE() / EL_MM_VISCOUS_FEEDBACK()
!                      the particle receives V div tau (Bfull, A) and the reaction is deposited
!                      on the fluid (Bfull). The sampling is done by EL_VISCOUS.
!
!   EL_MM_OPERATOR_FORM()
!                      model A: -eps grad p and div(eps tau) in the fluid operator
!                      (EL_MM_OPERATOR; momentum right-hand side and velocity correction).
!
! The pressure part of the particle force split (pressure force, its feedback) is expressed
! through the existing keys, which EL_CONFIG:EL_RESOLVE_MOMENTUM_MODEL sets for a model; the
! force module needs no model test for them.
!
! Non-Euler-Lagrange solvers never reach this module: it is called from the EL transport
! routine and the EL force/transfer modules only.
!=========================================================================
MODULE EL_MOMENTUM_MODEL

  USE EL_CONFIG, ONLY: EL_MM_ID

  IMPLICIT NONE
  PRIVATE

  ! model ids = indices of el_mm_names (EL_CONFIG)
  INTEGER, PARAMETER, PUBLIC :: EL_MM_NONE = 0, EL_MM_BFULL = 1, EL_MM_B = 2, EL_MM_A = 3

  PUBLIC :: EL_MM_ACTIVE, EL_MM_WEIGHTED_INERTIA, EL_MM_DRAG_SCALE, EL_MM_RESIDUAL_DIVISOR
  PUBLIC :: EL_MM_VISCOUS_FORCE, EL_MM_VISCOUS_FEEDBACK, EL_MM_OPERATOR_FORM

CONTAINS

  LOGICAL FUNCTION EL_MM_ACTIVE()
    EL_MM_ACTIVE = (EL_MM_ID().NE.EL_MM_NONE)
  END FUNCTION EL_MM_ACTIVE

  LOGICAL FUNCTION EL_MM_WEIGHTED_INERTIA()
    EL_MM_WEIGHTED_INERTIA = EL_MM_ACTIVE()
  END FUNCTION EL_MM_WEIGHTED_INERTIA

  ! The particle receives the viscous force V div tau (models Bfull and A) ...
  LOGICAL FUNCTION EL_MM_VISCOUS_FORCE()
    EL_MM_VISCOUS_FORCE = (EL_MM_ID().EQ.EL_MM_BFULL .OR. EL_MM_ID().EQ.EL_MM_A)
  END FUNCTION EL_MM_VISCOUS_FORCE

  ! ... and its reaction is deposited on the fluid (Bfull only; in model A the fluid operator
  ! carries div(eps tau) itself).
  LOGICAL FUNCTION EL_MM_VISCOUS_FEEDBACK()
    EL_MM_VISCOUS_FEEDBACK = (EL_MM_ID().EQ.EL_MM_BFULL)
  END FUNCTION EL_MM_VISCOUS_FEEDBACK

  ! Model A: the fluid operator itself carries -eps grad p and div(eps tau) (EL_MM_OPERATOR,
  ! applied by the fluid adapter as solid-fraction corrections of the unweighted operators).
  LOGICAL FUNCTION EL_MM_OPERATOR_FORM()
    EL_MM_OPERATOR_FORM = (EL_MM_ID().EQ.EL_MM_A)
  END FUNCTION EL_MM_OPERATOR_FORM

  REAL*8 FUNCTION EL_MM_DRAG_SCALE(eps)
    REAL*8, INTENT(IN) :: eps
    EL_MM_DRAG_SCALE = 1.0d0
    IF (EL_MM_ID().EQ.EL_MM_B) EL_MM_DRAG_SCALE = 1.0d0/eps
  END FUNCTION EL_MM_DRAG_SCALE

  REAL*8 FUNCTION EL_MM_RESIDUAL_DIVISOR(eps_e)
    REAL*8, INTENT(IN) :: eps_e
    SELECT CASE (EL_MM_ID())
    CASE (EL_MM_BFULL, EL_MM_B)
      EL_MM_RESIDUAL_DIVISOR = 1.0d0
    CASE DEFAULT
      EL_MM_RESIDUAL_DIVISOR = eps_e
    END SELECT
  END FUNCTION EL_MM_RESIDUAL_DIVISOR

END MODULE EL_MOMENTUM_MODEL
