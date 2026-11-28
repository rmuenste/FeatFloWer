MODULE EL_CONFIG

  USE MPI, ONLY: MPI_COMM_WORLD
  USE EL_INLET_TABLE, ONLY: tElInletTable, EL_INLET_TABLE_LOAD_SHARED, EL_INLET_TABLE_VALUE, &
                            EL_INLET_TABLE_PATH

  IMPLICIT NONE

  CHARACTER(LEN=32) :: el_kernel = 'deen_poly'
  CHARACTER(LEN=32) :: el_drag_model = 'difelice'
  ! ---- The drag closure table (SimPar@ELDragModel): the one place that names the closures --
  ! Index = closure id. EL_FORCES:EL_DRAG_CLOSURE dispatches on the id (one function per closure
  ! family, docs/md_docs/el_min_fluidization_prep_2026-11-19.md section 1), EL_VALIDATE_CONFIG
  ! accepts exactly these names, and EL_RESOLVE_MOMENTUM_MODEL reads the porosity column. A
  ! further closure is one entry here, one function and one CASE line in EL_DRAG_CLOSURE.
  !   difelice          Di Felice, campaign voidage eps^(-chi) (no eps^2 on the squared slip)
  !   difelice_std      Di Felice, standard form eps^(2-chi) (Kafui 2002, CFDEM DiFeliceDrag)
  !   stokes            single sphere 3 pi mu d (no porosity)
  !   schiller_naumann  single sphere with the Schiller-Naumann correction (no porosity)
  !   gidaspow          Ergun for eps <= 0.8, Wen-Yu above (CFDEM GidaspowDrag, MFiX GIDASPOW)
  !   wenyu             Wen-Yu at every porosity
  !   ergun             Ergun at every porosity (no dilute limit: vanishes for eps -> 1 at Re 0)
  ! Every closure returns the model-A per-particle coefficient V_p beta / (1 - eps) (EL_FORCES);
  ! momentum model B divides it by eps (EL_MOMENTUM_MODEL:EL_MM_DRAG_SCALE).
  INTEGER, PARAMETER :: EL_DRAG_DIFELICE = 0, EL_DRAG_DIFELICE_STD = 1, EL_DRAG_STOKES = 2, &
                        EL_DRAG_SCHILLER_NAUMANN = 3, EL_DRAG_GIDASPOW = 4, &
                        EL_DRAG_WENYU = 5, EL_DRAG_ERGUN = 6
  CHARACTER(LEN=16), PARAMETER :: el_drag_names(0:6) = (/ 'difelice        ', &
    'difelice_std    ', 'stokes          ', 'schiller_naumann', 'gidaspow        ', &
    'wenyu           ', 'ergun           ' /)
  ! the closure is a drag law in the standard porosity convention that the momentum models
  ! assume (refused under a model otherwise); the single-sphere laws carry no porosity at all
  LOGICAL, PARAMETER :: el_drag_standard_form(0:6) = (/ .FALSE., .TRUE., .TRUE., .TRUE., &
    .TRUE., .TRUE., .TRUE. /)
  ! Id of the selected closure, resolved once from el_drag_model by EL_RESOLVE_DRAG_CLOSURE
  ! (EL_VALIDATE_CONFIG; a unit test that sets el_drag_model directly calls it as well). -1:
  ! not resolved; EL_FORCES:EL_DRAG_CLOSURE stops the run on an id it does not dispatch.
  INTEGER :: el_drag_closure_id = -1
  CHARACTER(LEN=32) :: el_drag_coupling = 'explicit'
  ! SimPar@ELDragImplicitWeights = rowsum | hrz
  !   Nodal weights of the implicit part of the semi-implicit drag (the coefficient B_i on the
  !   diagonal of the fluid matrix and its split partner +B u_f on the explicit side; the drag
  !   impulse itself and the other reactions keep the row-sum weights in both modes):
  !   'rowsum' (default) w_i = int K phi_i / norm, the weights of every other deposit (campaign
  !            behaviour, bitwise). With Q2 basis functions w_i is negative at the nodes that a
  !            particle kernel reaches only with the negative lobes of phi_i: a node at the edge
  !            of a sparse particle cloud can get B_i < 0, i.e. an anti-damping term. Harmless
  !            while the fluid mass dominates (liquids); in a gas (case 11, rho_p/rho_f 860)
  !            thstep |B_i| is comparable with the lumped fluid mass and the velocity at such a
  !            node grows (spikes of 2-20 m/s in the first protocol run, job 163872).
  !   'hrz'    diagonal-scaling (Hinton-Rock-Zienkiewicz) lumping per element:
  !            w_i = int K phi_i^2 * (sum_j int K phi_j) / (sum_j int K phi_j^2) / norm >= 0, the
  !            same total per particle and element as 'rowsum'. Needs ELDragCoupling =
  !            semi_implicit (refused otherwise). docs/md_docs/el_min_fluidization_prep_2026-11-19.md
  !            section 7. Also needs a momentum defect that carries the drag diagonal
  !            (EL_DRAG_HRZ_DEFECT_SUPPORTED; refused at start-up by the application otherwise).
  CHARACTER(LEN=32) :: el_drag_implicit_weights = 'rowsum'
  CHARACTER(LEN=6), PARAMETER :: el_drag_iw_names(0:1) = (/ 'rowsum', 'hrz   ' /)
  ! Resolved by EL_VALIDATE_CONFIG (deck-only, the same on every rank): 'hrz' selected.
  LOGICAL :: el_drag_hrz = .FALSE.
  CHARACTER(LEN=32) :: el_lift_model = 'none'
  CHARACTER(LEN=32) :: el_inertial_lift = 'none'
  CHARACTER(LEN=32) :: el_prescribed_field = 'none'
  ! One-shot initial velocity field for COUPLED runs (istart=0 only):
  ! 'linear_shear' seeds u_x = el_shear_rate*(z - z_c) at t=0 (Couette
  ! startup skip); unlike el_prescribed_field it does NOT freeze the solve.
  CHARACTER(LEN=32) :: el_initial_field = 'none'
  CHARACTER(LEN=32) :: el_domain_type = 'box'
  CHARACTER(LEN=8) :: el_cylinder_axis = 'x'
  REAL*8 :: el_kernel_width_factor = 2.5d0
  REAL*8 :: el_eps_f_min = 0.4d0
  REAL*8 :: el_eps_f_relax = 0.0d0
  REAL*8 :: el_shear_rate = 0.0d0
  REAL*8 :: el_prescribed_umax = 0.0d0
  REAL*8 :: el_inertial_lift_umax = 0.0d0
  REAL*8 :: el_cylinder_center(3) = (/0.0d0,0.0d0,0.0d0/)
  REAL*8 :: el_cylinder_radius = -1.0d0
  ! ELApplyParticleForces: hand the evaluated hydrodynamic forces to PE (No =
  ! the forces are still evaluated and, with ELApplyFluidFeedback = Yes, fed
  ! back to the fluid, but PE advances the particles without them: a particle
  ! at rest without contacts stays at rest, a moving one keeps its velocity;
  ! a bed is held still by its PE configuration, not by this key. With
  ! semi-implicit drag the fluid split keeps +B u_f on the explicit side so
  ! that such a bed feels the drag once).
  LOGICAL :: el_apply_particle_forces = .TRUE.
  LOGICAL :: el_apply_fluid_feedback = .FALSE.
  ! EL_WALL_STRESS diagnostic (Couette wall force balance): fluid wall shear
  ! stress + PE wall lubrication/contact channels on the two z-walls.
  LOGICAL :: el_wall_stress_diag = .FALSE.
  ! Apply Prop@Gravity as a fluid body force (AddGravForce). Set to No for
  ! fully periodic suspensions: without boundaries there is no hydrostatic
  ! pressure to absorb rho_f*g, so fluid gravity must be dropped (absorbed
  ! into the modified pressure) while particles keep the analytic
  ! grav_buoy = (rho_p-rho_f)*V*g. Default Yes preserves bounded-domain
  ! behavior (tier2 terminal, ten Cate).
  LOGICAL :: el_fluid_gravity = .TRUE.
  ! SimPar@ELFluidViscosity = kinematic | dynamic
  !   Coefficient of the fluid's viscous term. The solver assembles its diffusion matrix with
  !   Properties%Viscosity and no density factor, against mass and convection matrices that
  !   carry the density: it expects the DYNAMIC viscosity there (the DNS decks give it so).
  !   The Euler-Lagrange decks give the KINEMATIC viscosity (Prop@Viscosity, or Prop@DynVisc
  !   which the parser divides by the density), because the drag closures form mu = rho nu
  !   themselves.
  !   'kinematic' the campaign behaviour up to 2026-11-15: the diffusion matrix carries nu,
  !               i.e. the fluid's viscous term is too small by the factor rho_f (no effect
  !               for rho_f = 1). Set it explicitly to reproduce a certified log of a case
  !               with rho_f /= 1 bit for bit.
  !   'dynamic'   (default since 2026-11-15, owner decision) the application stores mu = rho nu in Properties%Viscosity after the solver
  !               initialisation and rebuilds the diffusion matrix (EL_APPLY_FLUID_VISCOSITY in
  !               the application); the closures keep receiving nu (EL_CLOSURE_VISCOSITY), the
  !               particle viscous force and model A follow the operator
  !               (EL_FORCES: EL_FLUID_OPERATOR_VISCOSITY). Supported: Newtonian fluid without
  !               viscoelastic stress, Prop@PowerLawExp = 1, MatrixRenewal D >= 1 and S = 0
  !               (refused otherwise). A property reload at run time (Reload_Prop of the
  !               process control) would put nu back: the run then stops at the next
  !               Euler-Lagrange pass. The deck convention does not change (still nu or
  !               Prop@DynVisc); fluidViscosity_ of the particle json is an independent input
  !               and must hold mu as before.
  CHARACTER(LEN=32) :: el_fluid_viscosity = 'dynamic'
  CHARACTER(LEN=9), PARAMETER :: el_fluid_visc_names(0:1) = (/ 'kinematic', 'dynamic  ' /)
  ! Resolved by EL_VALIDATE_CONFIG (deck-only):
  LOGICAL :: el_fluid_visc_dynamic = .FALSE.
  ! Recorded by EL_SET_FLUID_VISCOSITY when the application converts the property: the
  ! kinematic value for the closures and the dynamic one now stored in Properties%Viscosity.
  ! Read through EL_CLOSURE_VISCOSITY / EL_FLUID_VISCOSITY_CONVERTED only.
  REAL*8, PRIVATE :: el_closure_nu = -1.0d0, el_operator_mu = -1.0d0
  LOGICAL :: el_drag_semi_implicit = .FALSE.
  LOGICAL :: el_pressure_force = .TRUE.
  ! Feed the particle pressure force -V grad p back to the fluid (kernel
  ! deposited with drag + lift). Without it the particle receives the
  ! suspension-weight share of grad p but the fluid keeps its full -grad p
  ! over the particle volume: an O(phi) momentum imbalance (liggghts-public
  ! finding 2026-10-15). With it the cell-averaged fluid pressure force is
  ! -eps grad p, i.e. the Model-A form. Fed back in both transfer passes (the
  ! pre-advance deposit is the one the fluid receives; until 2026-11-04 the
  ! pre-advance rebuild with particles advanced dropped it, P5). Needs
  ! ELPressureForce = Yes and ELApplyFluidFeedback = Yes (refused otherwise).
  ! Default No keeps every existing numerical output bitwise; the
  ! pressure-force rungs are to be re-run with Yes.
  LOGICAL :: el_pressure_feedback = .FALSE.
  LOGICAL :: el_magnus = .FALSE.
  LOGICAL :: el_write_diagnostics = .TRUE.
  ! Measured-leak momentum compensator (SimPar@ELMomentumFix): each step
  ! the EL_FLUID_PAIR audit measures the fluid-internal momentum error of
  ! the previous fluid solve (dominated by the Galerkin convective form,
  ! -rho*int(u_i div u)); when enabled, its negative is applied as a
  ! uniform body force through the Grav_QuadSc pathway one step lagged
  ! (deadbeat: cumulative drift stays bounded by one step's leak).
  ! Use with the legacy convective form; the momentum-exact divergence
  ! form (ELConvectionForm=divergence) is not energy-robust in developed
  ! pseudo-turbulence (NaN blowup, see v2_rz_settling RUNBOOK part 4).
  LOGICAL :: el_momentum_fix = .FALSE.
  INTEGER :: el_momentum_audit_freq = 100
  ! Per-particle trajectory log (SimPar@ELParticleLogFreq, 0 = off): every N
  ! steps each worker appends its locally owned particles to
  ! el_particles_rank<myid>.log in the DNS particle_force.log column layout
  ! (time ip fx fy fz tx ty tz px py pz vx vy vz); merge with
  ! tools/el_particle_log_merge.py. Collective-free by construction.
  INTEGER :: el_particle_log_freq = 0
  ! Volume-averaged continuity (SimPar@ELVolumeAveraged): 'none' (default,
  ! single-phase div u = 0) or 'continuity' (div u = -(d eps/dt + u.grad eps)/eps
  ! imposed through the pressure Poisson right-hand side; EL_VOLUME_AVERAGED).
  CHARACTER(LEN=32) :: el_volume_averaged = 'none'
  ! --- Continuity-source time integration (only read when ELVolumeAveraged = continuity;
  !     both defaults reproduce the backward-Euler campaign code bitwise) -----------------
  ! SimPar@ELVASourceTimeScheme = BE | BDF2
  !   Time scheme of the rate d eps/dt in the source g = -(d eps/dt + u.grad eps)/eps:
  !   'be'   backward difference (eps^{n+1}-eps^n)/dt (default)
  !   'bdf2' endpoint BDF2 1.5*(d eps/dt)^{n+1} - 0.5*(d eps/dt)^n, BE on the first step
  !          and after a restart without rate history (s04 RUNBOOK section 8: removes the
  !          O(dt) return-flow deficit of the VB twin at level 2).
  CHARACTER(LEN=32) :: el_va_source_time_scheme = 'be'
  ! SimPar@ELVASourceCorrectors = n (0..20, default 0)
  !   Source/projection corrector passes after the pressure projection: reassemble the
  !   divergence of the corrected velocity, recompute g from it, solve the pressure Poisson
  !   problem for the residual and correct u and p again (Transport_Q2P1,
  !   EL_VA_Projection_Correctors). Each pass costs one pressure solve.
  INTEGER :: el_va_source_correctors = 0
  ! SimPar@ELVACompatibility = uniform | alpha
  !   How the volume integral of the raw source is removed so that the pressure problem is
  !   compatible (closed/periodic boxes): 'uniform' (default) subtracts the volume mean
  !   everywhere, 'alpha' subtracts int_g * alpha_p(e)/W (W = sum alpha_p|e|), i.e. the
  !   compensating sink sits where the particles are (EL_VA_SINK; uniform fallback without
  !   particles). 'alpha' closes the level-3 return-flow deficit of the VB twin (s04 RUNBOOK
  !   section 10).
  CHARACTER(LEN=32) :: el_va_compatibility = 'uniform'
  !   'none' applies the raw source unchanged: required when the domain has an outflow
  !   boundary (SimPar@NoOutflow = No; the pressure problem is solvable for any source and a
  !   removal would act as a spurious source, case 05 E2), refused in a closed box (checked
  !   in the application, where NoOutflow is known).
  ! The valid values, blank-padded; index = sink mode id of EL_VA_CONTINUITY (0 uniform, 1 alpha,
  ! 2 none).
  CHARACTER(LEN=7), PARAMETER :: el_va_compat_names(0:2) = (/ 'uniform', 'alpha  ', 'none   ' /)
  ! SimPar@ELVAContinuityForm = advective | conservative
  !   Spatial form of the continuity constraint (EL_VA_CONTINUITY):
  !   'advective'    (default) div u = -(d eps/dt + u.grad eps)/eps per element, the advective
  !                  term from the vertex-mean velocity and the kernel gradient (campaign code).
  !                  The rate comes from the stored porosity, the gradient from the raw kernel
  !                  field -grad alpha_p: the two agree only without relaxation (ELEpsFRelax =
  !                  0, enforced) and outside clipped cells (ELEpsFMin; clipping is reported
  !                  per pass, the advective term is then inconsistent there).
  !   'conservative' div(eps u) = -d eps/dt in the weak form of all four P1 pressure test
  !                  functions, with shared face fluxes of a continuous (Q1-reconstructed)
  !                  transport porosity (EL_VA_FLUX): the net source telescopes to the boundary
  !                  flux exactly, so the compatibility sink vanishes at convergence. Iterated
  !                  through the existing projection, hence needs ELVASourceCorrectors >= 1;
  !                  certified for closed impermeable boxes (docs/md_docs/
  !                  el_va_conservative_form_2026-11-05.md, s04 RUNBOOK section 12).
  CHARACTER(LEN=32) :: el_va_continuity_form = 'advective'
  ! The valid values, blank-padded; index = form id of EL_VA_CONTINUITY (0 advective,
  ! 1 conservative).
  CHARACTER(LEN=12), PARAMETER :: el_va_form_names(0:1) = (/ 'advective   ', 'conservative' /)
  ! SimPar@ELPeriodicKernel = off | wrap
  !   Periodic image treatment of the particle-mesh coupling (EL_PERIODIC, EL_IMAGE_TRANSFER;
  !   docs/md_docs/el_periodic_kernel_2026-11-17.md):
  !   'off'  (default) the kernel of a particle ends at the periodic planes of the fluid
  !          domain like at a wall: a particle within one kernel radius of a plane deposits and
  !          samples on its own side only, and its footprint jumps by the period when the
  !          particle engine wraps it (the campaign behaviour, every certified log bitwise).
  !   'wrap' on every periodic axis of the fluid (SimPar@Periodic / PeriodicAxis /
  !          PeriodicLength) a particle acts through its periodic images: deposit, sampling,
  !          kernel normalisation and the halo exchange see the particle on both sides of a
  !          plane; the faces of a periodic axis are no walls for the wall-corrected lift.
  !          Refused (list in the design note): without a periodic axis, with a period that
  !          is not the extent of the fluid domain, with ELWallStressDiag and z periodic
  !          (the application checks these, where the periodicity is known); when the
  !          kernel radius reaches a period or a subdomain spans one (EL_HALO).
  !          The periodicity of the particle engine (periodicX_/Y_/Z_ of its json) is a
  !          separate input, must describe the same axes and is not cross-checked.
  !   The echo line appears only for 'wrap'.
  CHARACTER(LEN=32) :: el_periodic_kernel = 'off'
  ! The valid values, blank-padded; index 1 = 'wrap' is the one the resolution tests.
  CHARACTER(LEN=4), PARAMETER :: el_periodic_kernel_names(0:1) = (/ 'off ', 'wrap' /)
  ! Resolved by EL_VALIDATE_CONFIG (deck-only, the same on every rank): the deck asks for
  ! 'wrap'. The application turns it into the lattice of EL_PERIODIC; nothing else reads it.
  LOGICAL :: el_periodic_wrap = .FALSE.
  ! SimPar@ELMomentumModel = none | Bfull | B | A
  !   Momentum formulation of the fluid-particle set, named after the CFDEM model types
  !   (cfdemSolverPiso; docs/md_docs/el_momentum_models_plan_2026-11-15.md):
  !   'none'  (default) unweighted fluid momentum rho Du/Dt = -grad p + ..., the particle force
  !           split is whatever ELPressureForce / ELPressureFeedback say (campaign behaviour,
  !           every certified log bitwise).
  !   'bfull' eps rho Du/Dt = -grad p + div tau - n (F_d + F_p + F_v): particle pressure and
  !           viscous force, both fed back to the fluid.
  !   'b'     eps rho Du/Dt = -grad p + div tau - n F_d/eps: drag divided by eps, no particle
  !           pressure or viscous force.
  !   'a'     eps rho Du/Dt = -eps grad p + div(eps tau) - n F_d: eps in the operator, particle
  !           pressure and viscous force without feedback.
  !   A model implies its consistent key set; el_mm_implemented says which names this build
  !   accepts (a recognised but unimplemented name is refused, never silently run as 'none').
  !   Behaviour lives in EL_MOMENTUM_MODEL; nothing else tests the name.
  CHARACTER(LEN=32) :: el_momentum_model = 'none'
  ! SimPar@ELMMPressureOperator = startup | weighted
  !   Pressure matrix of the projection under a momentum model (the preconditioner of the
  !   conservative continuity iteration; the converged state does not depend on it):
  !   'startup'  (default) the matrix of start-up, C0 = B^T (rho Ml)^-1 B, with the residual
  !              divisor of the model (EL_MOMENTUM_MODEL: EL_MM_RESIDUAL_DIVISOR). For model A
  !              that is eps C0 against the true operator div(eps grad): in a dense bed the
  !              flux through the bed is over-corrected and the iteration converges slowly
  !              with alternating sign (case 05, moving bed: about 4 correctors needed).
  !   'weighted' the matrix is rebuilt every step from the current porosity,
  !              C = B^T (rho Ml / eps)^-1 B with the nodal porosity (a nodal approximation
  !              of the operator), the divisor is 1. Available for the models whose table
  !              column el_mmt_pressure_rebuild is set (A); accepted without effect for
  !              bfull and b, whose start-up matrix with divisor 1 already is that
  !              approximation (exact cancellation for uniform porosity); refused without a
  !              model. Needs Pres@MGCrsSolverType 1-4 or 9 and Pres@MGMinLev = MGMedLev =
  !              coarsest level (checked by the fluid adapter). The echo line appears only
  !              for 'weighted'.
  CHARACTER(LEN=32) :: el_mm_pressure_operator = 'startup'
  ! The valid values; index 1 = 'weighted' is the one the resolution tests.
  CHARACTER(LEN=8), PARAMETER :: el_mm_pop_names(0:1) = (/ 'startup ', 'weighted' /)
  ! Resolved by EL_RESOLVE_MOMENTUM_MODEL: the deck asks for 'weighted' / the selected model
  ! rebuilds its pressure matrix for it (read through EL_MOMENTUM_MODEL:
  ! EL_MM_PRESSURE_REBUILD; deck-only, the same on every rank).
  LOGICAL :: el_mm_pop_weighted = .FALSE.
  LOGICAL :: el_mm_pressure_rebuild_on = .FALSE.
  ! ---- The model table: the one place that says what a model is -------------------------
  ! Index = model id. EL_MOMENTUM_MODEL answers the run-time questions from these columns and
  ! EL_RESOLVE_MOMENTUM_MODEL the key set; a further model is one more column.
  INTEGER, PARAMETER :: EL_MM_NONE = 0, EL_MM_BFULL = 1, EL_MM_B = 2, EL_MM_A = 3
  CHARACTER(LEN=5), PARAMETER :: el_mm_names(0:3) = (/ 'none ', 'bfull', 'b    ', 'a    ' /)
  LOGICAL, PARAMETER :: el_mm_implemented(0:3)       = (/ .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.  /)
  ! the particle receives -V grad p / its reaction is deposited on the fluid
  LOGICAL, PARAMETER :: el_mmt_particle_pressure(0:3) = (/ .FALSE., .TRUE.,  .FALSE., .TRUE.  /)
  LOGICAL, PARAMETER :: el_mmt_pressure_feedback(0:3) = (/ .FALSE., .TRUE.,  .FALSE., .FALSE. /)
  ! the particle receives V c lap u / its reaction is deposited on the fluid
  LOGICAL, PARAMETER :: el_mmt_viscous_force(0:3)     = (/ .FALSE., .TRUE.,  .FALSE., .TRUE.  /)
  LOGICAL, PARAMETER :: el_mmt_viscous_feedback(0:3)  = (/ .FALSE., .TRUE.,  .FALSE., .FALSE. /)
  ! drag force and coefficient divided by eps
  LOGICAL, PARAMETER :: el_mmt_drag_over_eps(0:3)     = (/ .FALSE., .FALSE., .TRUE.,  .FALSE. /)
  ! eps-weighted fluid inertia / -eps grad p and div(eps tau) in the fluid operator
  LOGICAL, PARAMETER :: el_mmt_weighted_inertia(0:3)  = (/ .FALSE., .TRUE.,  .TRUE.,  .TRUE.  /)
  LOGICAL, PARAMETER :: el_mmt_operator_form(0:3)     = (/ .FALSE., .FALSE., .FALSE., .TRUE.  /)
  ! the conservative continuity residual is divided by the element porosity (else by 1) before
  ! it enters the pressure right-hand side (EL_MOMENTUM_MODEL: EL_MM_RESIDUAL_DIVISOR)
  LOGICAL, PARAMETER :: el_mmt_divisor_eps(0:3)       = (/ .TRUE.,  .FALSE., .FALSE., .TRUE.  /)
  ! the start-up pressure matrix does not approximate the projection operator of the model
  ! and SimPar@ELMMPressureOperator = weighted rebuilds it (bfull, b: the start-up matrix
  ! already does; 'none': the rebuild is not implemented)
  LOGICAL, PARAMETER :: el_mmt_pressure_rebuild(0:3)  = (/ .FALSE., .FALSE., .FALSE., .TRUE.  /)
  ! (row 'none' of the two pressure columns is not used: the campaign keys decide there)
  ! Model id of the validated deck value, set by EL_VALIDATE_CONFIG (0 = none before that and
  ! by default); read through EL_MM_ID. Code that changes el_momentum_model must validate again.
  INTEGER :: el_mm_id_validated = 0
  ! Deck keys of the SimPar@EL* family that the parameter file set explicitly (recorded by the
  ! application parser through EL_NOTE_KEY_SET). A momentum model implies values for other
  ! keys; the registry lets the resolution tell "left at its default" from "set to a
  ! contradicting value", which is refused.
  INTEGER, PARAMETER :: el_max_keys_set = 128
  CHARACTER(LEN=32) :: el_keys_set(el_max_keys_set) = ' '
  INTEGER :: el_n_keys_set = 0
  ! Echo lines of the key values a momentum model implied (EL_RESOLVE_MOMENTUM_MODEL),
  ! printed by EL_PRINT_CONFIG after the model name.
  INTEGER, PARAMETER :: el_max_implied = 16
  CHARACTER(LEN=96) :: el_mm_implied(el_max_implied) = ' '
  INTEGER :: el_n_implied = 0
  INTEGER :: el_mm_conflicts = 0
  REAL*8 :: el_tavg_window = 0.5d0
  ! Uniform ramped inlet (boundary type Inflow975, GetVeloBCVal): axial velocity normal to a
  ! z-plane inlet, ramped linearly from ELInletVelocityStart to ELInletVelocityEnd over
  ! ELInletRampTime and constant afterwards (EL_INLET_VELOCITY). The deck value is the
  ! SUPERFICIAL velocity: under ELVolumeAveraged = continuity the solver velocity is
  ! interstitial, so the boundary value is divided by the nodal porosity on the inlet plane
  ! (EL_INLET_SUPERFICIAL, EL_FIELDS:EL_INLET_INTERSTITIAL_FACTOR; a bed resting on the inlet
  ! has eps < 1 there). Without the continuity stage the value is applied as it is. All
  ! three default to an inert inlet (0, 0, ramp 1 s); a mesh that carries Inflow975 with both
  ! velocities zero is refused, because it would silently be a no-slip floor.
  REAL*8 :: el_inlet_u0 = 0.0d0
  REAL*8 :: el_inlet_u1 = 0.0d0
  REAL*8 :: el_inlet_tramp = 1.0d0
  ! SimPar@ELInletVelocityTable = "<file>": piecewise-linear U(t) instead of the ramp (rows
  ! "t U", superficial velocity; format and refusals in EL_INLET_TABLE). Mutually exclusive
  ! with the three ramp keys (a deck that sets both is refused), an empty path is refused.
  ! Loaded once by EL_VALIDATE_CONFIG (rank 0 reads, every rank receives the rows; a start-up
  ! setting, not touched by Reload_SimPar); echoed by EL_PRINT_CONFIG. Unset (default): the
  ! ramp, bitwise as before.
  CHARACTER(LEN=256) :: el_inlet_table_file = ' '
  TYPE(tElInletTable) :: el_inlet_table
  ! The inlet definition (ramp keys or table) is a start-up setting: Reload_SimPar runs the
  ! generic parser, which reads no SimPar@EL* key. At a reload EL_INLET_RELOAD_NOTICE compares
  ! the definition the reload file asks for with the one of the run and reports a difference
  ! (changed table path, table removed, table <-> ramp, changed ramp values); nothing changes.
  TYPE tElInletDefinition
    LOGICAL :: table = .FALSE.         ! ELInletVelocityTable given
    CHARACTER(LEN=256) :: path = ' '   ! its file name (EL_INLET_TABLE_PATH of the value)
    LOGICAL :: ramp_keys = .FALSE.     ! one of the three ramp keys given
    REAL*8 :: u0 = 0.0d0, u1 = 0.0d0, tramp = 1.0d0   ! ramp values (the defaults when unset)
    LOGICAL :: readable = .TRUE.       ! every ramp value of the file could be read
  END TYPE tElInletDefinition
  ! EL_VALIDATE_CONFIG has run (the application parsed and checked the EL keys); the reload
  ! notice is silent in applications that never validate an EL configuration
  LOGICAL :: el_config_validated = .FALSE.
  ! SimPar@ELBedAuditFreq = N (default 0 = off): every N steps the EL_BED_AUDIT line (EL_BED_AUDIT
  ! module: area-mean pressure on the inlet and outlet planes, their difference with and
  ! without the hydrostatic part, buoyant bed weight, floor contact force, bed heights, porosity
  ! of the bed window). Diagnostic only; no effect on the solution.
  INTEGER :: el_bed_audit_freq = 0
  ! SimPar@ELVelocityAuditFreq = N (default 0 = off): every N steps the EL_UMAX_AUDIT line
  ! (EL_VELOCITY_AUDIT: location, components and local porosity of the fastest fluid dof,
  ! smallest element porosity, inlet-plane porosity). Diagnostic only.
  INTEGER :: el_velocity_audit_freq = 0
  ! SimPar@ELMomentumBalanceAudit = Yes (default No): every step the EL_MOM_BALANCE_* lines
  ! (EL_MOMENTUM_BALANCE: the impulse the particles receive, the coupling terms of the gas
  ! momentum equation as assembled, boundary planes and storage, closure residuals). Reads the
  ! solver state only; no effect on the solution. Needs ELDomainType = box (inlet / outlet
  ! planes); the application refuses solver settings the balance does not cover (app_init).
  LOGICAL :: el_momentum_balance_audit = .FALSE.
  ! Mesoscale lane-mode filter for periodic settling boxes: each step,
  ! subtract from u_z its zero-mean x-binned and y-binned plane averages
  ! (the gravest "lane" modes of the two-way-coupled concentration
  ! instability, which lateral container walls would suppress in the
  ! wall-bounded experiments RZ describes). The subtracted field is
  ! z-independent (divergence-free on the axis-aligned box) and has its
  ! global mean removed (momentum-neutral; any lumping residue lands in
  ! the EL_FLUID_PAIR mismatch and is absorbed by ELMomentumFix).
  ! Periodic-box suspension runs only.
  LOGICAL :: el_meso_filter = .FALSE.
  INTEGER :: el_meso_filter_bins = 16

CONTAINS

  SUBROUTINE EL_VALIDATE_CONFIG()

    CALL EL_LOWERCASE(el_kernel)
    CALL EL_LOWERCASE(el_drag_model)
    CALL EL_LOWERCASE(el_drag_coupling)
    CALL EL_LOWERCASE(el_drag_implicit_weights)
    CALL EL_LOWERCASE(el_lift_model)
    CALL EL_LOWERCASE(el_inertial_lift)
    CALL EL_LOWERCASE(el_prescribed_field)
    CALL EL_LOWERCASE(el_initial_field)
    CALL EL_LOWERCASE(el_domain_type)
    CALL EL_LOWERCASE(el_volume_averaged)
    CALL EL_LOWERCASE(el_va_source_time_scheme)
    CALL EL_LOWERCASE(el_va_compatibility)
    CALL EL_LOWERCASE(el_va_continuity_form)
    CALL EL_LOWERCASE(el_cylinder_axis)
    CALL EL_LOWERCASE(el_momentum_model)
    CALL EL_LOWERCASE(el_mm_pressure_operator)
    CALL EL_LOWERCASE(el_fluid_viscosity)
    CALL EL_LOWERCASE(el_periodic_kernel)
    IF (.NOT.ANY(TRIM(el_periodic_kernel).EQ.el_periodic_kernel_names)) THEN
      WRITE(*,'(A,A)') 'Invalid ELPeriodicKernel (off|wrap): ', TRIM(el_periodic_kernel)
      STOP 1
    END IF
    el_periodic_wrap = (TRIM(el_periodic_kernel).EQ.TRIM(el_periodic_kernel_names(1)))
    IF (.NOT.ANY(TRIM(el_fluid_viscosity).EQ.el_fluid_visc_names)) THEN
      WRITE(*,'(A,A)') 'Invalid ELFluidViscosity (kinematic|dynamic): ', TRIM(el_fluid_viscosity)
      STOP 1
    END IF
    el_fluid_visc_dynamic = (TRIM(el_fluid_viscosity).EQ.TRIM(el_fluid_visc_names(1)))

    el_mm_id_validated = EL_MM_LOOKUP(el_momentum_model)
    IF (el_mm_id_validated.LT.0) THEN
      WRITE(*,'(A,A)') 'Invalid ELMomentumModel (none|Bfull|B|A): ', TRIM(el_momentum_model)
      STOP 1
    END IF
    IF (.NOT.el_mm_implemented(EL_MM_ID())) THEN
      WRITE(*,'(A,A,A)') 'ELMomentumModel = ', TRIM(el_momentum_model), &
        ' is not implemented in this build.'
      STOP 1
    END IF
    ! a model fixes other keys; resolved before the ordinary validation below reads them
    CALL EL_RESOLVE_MOMENTUM_MODEL()

    SELECT CASE (TRIM(el_kernel))
    CASE ('deen_poly', 'gaussian')
      CONTINUE
    CASE DEFAULT
      WRITE(*,'(A,A)') 'Invalid ELKernel: ', TRIM(el_kernel)
      STOP 1
    END SELECT

    CALL EL_RESOLVE_DRAG_CLOSURE()
    IF (el_drag_closure_id.LT.0) THEN
      WRITE(*,'(A,A)') 'Invalid ELDragModel: ', TRIM(el_drag_model)
      STOP 1
    END IF

    SELECT CASE (TRIM(el_drag_coupling))
    CASE ('explicit')
      el_drag_semi_implicit = .FALSE.
    CASE ('semi_implicit')
      el_drag_semi_implicit = .TRUE.
    CASE DEFAULT
      WRITE(*,'(A,A)') 'Invalid ELDragCoupling: ', TRIM(el_drag_coupling)
      STOP 1
    END SELECT

    IF (.NOT.ANY(TRIM(el_drag_implicit_weights).EQ.el_drag_iw_names)) THEN
      WRITE(*,'(A,A)') 'Invalid ELDragImplicitWeights (rowsum|hrz): ', TRIM(el_drag_implicit_weights)
      STOP 1
    END IF
    el_drag_hrz = (TRIM(el_drag_implicit_weights).EQ.TRIM(el_drag_iw_names(1)))
    IF (el_drag_hrz .AND. .NOT.el_drag_semi_implicit) THEN
      WRITE(*,*) 'ELDragImplicitWeights = hrz needs ELDragCoupling = semi_implicit.'
      STOP 1
    END IF

    IF (el_drag_semi_implicit .AND. .NOT.el_apply_fluid_feedback) THEN
      WRITE(*,*) 'ELDragCoupling=semi_implicit requires ELApplyFluidFeedback=Yes.'
      STOP 1
    END IF
    IF (el_pressure_feedback .AND. .NOT.(el_apply_fluid_feedback .AND. el_pressure_force)) THEN
      WRITE(*,*) 'ELPressureFeedback=Yes requires ELApplyFluidFeedback=Yes and ELPressureForce=Yes.'
      STOP 1
    END IF
    IF (el_momentum_fix .AND. el_drag_semi_implicit) THEN
      ! EL_FLUID_PAIR counts the captured explicit source only; the drag operator
      ! acting through the matrix would be read as a mismatch and compensated.
      WRITE(*,*) 'ELMomentumFix=Yes is not supported with ELDragCoupling=semi_implicit', &
                 ' (the momentum audit does not include the implicit drag operator).'
      STOP 1
    END IF

    SELECT CASE (TRIM(el_lift_model))
    CASE ('none', 'saffman', 'saffman_mei', 'saffman_mei_wall')
      CONTINUE
    CASE DEFAULT
      WRITE(*,'(A,A)') 'Invalid ELLiftModel: ', TRIM(el_lift_model)
      STOP 1
    END SELECT

    SELECT CASE (TRIM(el_inertial_lift))
    CASE ('none', 'matas_asmolov')
      CONTINUE
    CASE DEFAULT
      WRITE(*,'(A,A)') 'Invalid ELInertialLift: ', TRIM(el_inertial_lift)
      STOP 1
    END SELECT

    IF (TRIM(el_inertial_lift).NE.'none') THEN
      IF (TRIM(el_domain_type).NE.'cylinder' .OR. &
          el_cylinder_radius.LE.0.0d0) THEN
        WRITE(*,*) 'ELInertialLift requires ELDomainType=cylinder with ', &
          'ELCylinderRadius > 0.'
        STOP 1
      END IF
      IF (el_inertial_lift_umax.LE.0.0d0 .AND. &
          el_prescribed_umax.LE.0.0d0) THEN
        WRITE(*,*) 'ELInertialLift requires ELInertialLiftUmax > 0 ', &
          '(or ELPrescribedUmax > 0 in frozen runs).'
        STOP 1
      END IF
    END IF

    SELECT CASE (TRIM(el_prescribed_field))
    CASE ('none', 'linear_shear', 'poiseuille')
      CONTINUE
    CASE DEFAULT
      WRITE(*,'(A,A)') 'Invalid ELPrescribedField: ', TRIM(el_prescribed_field)
      STOP 1
    END SELECT

    SELECT CASE (TRIM(el_initial_field))
    CASE ('none', 'linear_shear')
      CONTINUE
    CASE DEFAULT
      WRITE(*,'(A,A)') 'Invalid ELInitialField: ', TRIM(el_initial_field)
      STOP 1
    END SELECT
    IF (TRIM(el_initial_field).EQ.'linear_shear' .AND. &
        ABS(el_shear_rate).LE.0.0d0) THEN
      WRITE(*,*) 'ELInitialField=linear_shear requires ELShearRate /= 0.'
      STOP 1
    END IF

    SELECT CASE (TRIM(el_domain_type))
    CASE ('box', 'cylinder')
      CONTINUE
    CASE DEFAULT
      WRITE(*,'(A,A)') 'Invalid ELDomainType: ', TRIM(el_domain_type)
      STOP 1
    END SELECT

    SELECT CASE (TRIM(el_cylinder_axis))
    CASE ('x', 'y', 'z')
      CONTINUE
    CASE DEFAULT
      WRITE(*,'(A,A)') 'Invalid ELCylinderAxis: ', TRIM(el_cylinder_axis)
      STOP 1
    END SELECT

    IF (TRIM(el_domain_type).EQ.'cylinder' .AND. &
        el_cylinder_radius.LE.0.0d0) THEN
      WRITE(*,*) 'ELDomainType=cylinder requires ELCylinderRadius > 0.'
      STOP 1
    END IF

    IF (el_magnus) THEN
      WRITE(*,*) 'ELMagnus is not implemented in Phase 2.'
      STOP 1
    END IF

    IF (el_kernel_width_factor.LE.0.0d0) THEN
      WRITE(*,*) 'ELKernelWidthFactor must be positive.'
      STOP 1
    END IF

    IF (el_eps_f_min.LE.0.0d0 .OR. el_eps_f_min.GT.1.0d0) THEN
      WRITE(*,*) 'ELEpsFMin must be in (0,1].'
      STOP 1
    END IF

    IF (el_eps_f_relax.LT.0.0d0 .OR. el_eps_f_relax.GE.1.0d0) THEN
      WRITE(*,*) 'ELEpsFRelax must be in [0,1).'
      STOP 1
    END IF

    IF (el_momentum_audit_freq.LE.0) THEN
      WRITE(*,*) 'ELMomentumAuditFreq must be positive.'
      STOP 1
    END IF

    IF (el_particle_log_freq.LT.0) THEN
      WRITE(*,*) 'ELParticleLogFreq must be >= 0 (0 = off).'
      STOP 1
    END IF

    SELECT CASE (TRIM(el_volume_averaged))
    CASE ('none', 'continuity')
    CASE DEFAULT
      WRITE(*,'(A,A)') 'Invalid ELVolumeAveraged: ', TRIM(el_volume_averaged)
      STOP 1
    END SELECT
    SELECT CASE (TRIM(el_va_source_time_scheme))
    CASE ('be', 'bdf2')
    CASE DEFAULT
      WRITE(*,'(A,A)') 'Invalid ELVASourceTimeScheme (BE|BDF2): ', TRIM(el_va_source_time_scheme)
      STOP 1
    END SELECT
    IF (el_va_source_correctors.LT.0 .OR. el_va_source_correctors.GT.20) THEN
      WRITE(*,*) 'ELVASourceCorrectors must be in 0..20.'
      STOP 1
    END IF
    IF (.NOT.ANY(TRIM(el_va_compatibility).EQ.el_va_compat_names)) THEN
      WRITE(*,'(A,A)') 'Invalid ELVACompatibility (uniform|alpha|none): ', TRIM(el_va_compatibility)
      STOP 1
    END IF
    IF (.NOT.ANY(TRIM(el_va_continuity_form).EQ.el_va_form_names)) THEN
      WRITE(*,'(A,A)') 'Invalid ELVAContinuityForm (advective|conservative): ', &
        TRIM(el_va_continuity_form)
      STOP 1
    END IF
    IF (TRIM(el_volume_averaged).NE.'continuity' .AND. &
        (TRIM(el_va_source_time_scheme).NE.'be' .OR. el_va_source_correctors.GT.0 .OR. &
         TRIM(el_va_compatibility).NE.'uniform' .OR. &
         TRIM(el_va_continuity_form).NE.'advective')) THEN
      WRITE(*,*) 'ELVASourceTimeScheme / ELVASourceCorrectors / ELVACompatibility / ', &
                 'ELVAContinuityForm need ELVolumeAveraged = continuity.'
      STOP 1
    END IF
    IF (TRIM(el_va_continuity_form).EQ.'conservative') THEN
      ! the conservative constraint is a fixed point through the projection: without a
      ! corrector pass only the predictor estimate would be applied
      IF (el_va_source_correctors.LT.1) THEN
        WRITE(*,*) 'ELVAContinuityForm = conservative needs ELVASourceCorrectors >= 1.'
        STOP 1
      END IF
      ! the lane-mode filter acts on the velocity after the projection and would break the
      ! constraint it just established
      IF (el_meso_filter) THEN
        WRITE(*,*) 'ELVAContinuityForm = conservative is not supported with ELMesoFilter = Yes.'
        STOP 1
      END IF
    END IF
    ! The advective form takes the rate from the stored (relaxed, clipped) porosity and the
    ! gradient from the raw kernel field; with relaxation they describe different fields.
    ! The conservative form reconstructs its transport porosity from the stored field and is
    ! consistent with both relaxation and clipping.
    IF (TRIM(el_volume_averaged).EQ.'continuity' .AND. &
        TRIM(el_va_continuity_form).EQ.'advective' .AND. el_eps_f_relax.GT.0.0d0) THEN
      WRITE(*,*) 'ELVolumeAveraged = continuity with ELVAContinuityForm = advective needs ', &
                 'ELEpsFRelax = 0 (the advective term uses the unrelaxed kernel gradient).'
      STOP 1
    END IF

    IF (el_inlet_tramp.LE.0.0d0) THEN
      WRITE(*,*) 'ELInletRampTime must be positive.'
      STOP 1
    END IF
    CALL EL_LOAD_INLET_TABLE()

    IF (el_bed_audit_freq.LT.0) THEN
      WRITE(*,*) 'ELBedAuditFreq must be >= 0 (0 = off).'
      STOP 1
    END IF
    ! the audit planes and the bed window are the z planes of the domain box (EL_BED_AUDIT)
    IF (el_bed_audit_freq.GT.0 .AND. TRIM(el_domain_type).NE.'box') THEN
      WRITE(*,'(A,A,A)') 'ELBedAuditFreq > 0 needs ELDomainType = box (inlet / outlet planes and ', &
        'the bed window are z planes of the domain box); ELDomainType = ', TRIM(el_domain_type)
      STOP 1
    END IF
    IF (el_velocity_audit_freq.LT.0) THEN
      WRITE(*,*) 'ELVelocityAuditFreq must be >= 0 (0 = off).'
      STOP 1
    END IF
    ! the boundary terms of the momentum balance are integrals over the z planes of the box
    IF (el_momentum_balance_audit .AND. TRIM(el_domain_type).NE.'box') THEN
      WRITE(*,'(A,A)') 'ELMomentumBalanceAudit = Yes needs ELDomainType = box (inlet / outlet ', &
        'are the z planes of the domain box); ELDomainType = '//TRIM(el_domain_type)
      STOP 1
    END IF
    ! the particle lines record the forces of the transfer pass as the forces armed in pe
    IF (el_momentum_balance_audit .AND. .NOT.el_apply_particle_forces) THEN
      WRITE(*,*) 'ELMomentumBalanceAudit = Yes needs ELApplyParticleForces = Yes (the particle ', &
        'balance takes the evaluated forces as the forces applied in pe).'
      STOP 1
    END IF

    IF (el_tavg_window.LE.0.0d0 .OR. el_tavg_window.GT.1.0d0) THEN
      WRITE(*,*) 'ELTAvgWindow must be in (0,1].'
      STOP 1
    END IF

    el_config_validated = .TRUE.

  END SUBROUTINE EL_VALIDATE_CONFIG

  SUBROUTINE EL_PRINT_CONFIG(mfile, mterm)

    INTEGER, INTENT(IN) :: mfile, mterm
    INTEGER :: i_implied

    WRITE(mfile,'(A,A)') 'EL kernel                 = ', TRIM(el_kernel)
    WRITE(mfile,'(A,ES14.6)') 'EL kernel width factor    = ', el_kernel_width_factor
    WRITE(mfile,'(A,ES14.6)') 'EL epsilon minimum        = ', el_eps_f_min
    WRITE(mfile,'(A,ES14.6)') 'EL epsilon relaxation     = ', el_eps_f_relax
    WRITE(mfile,'(A,A)') 'EL drag model             = ', TRIM(el_drag_model)
    WRITE(mfile,'(A,A)') 'EL drag coupling          = ', TRIM(el_drag_coupling)
    IF (el_drag_hrz) WRITE(mfile,'(A,A)') 'EL drag implicit weights  = ', TRIM(el_drag_implicit_weights)
    WRITE(mfile,'(A,L1)') 'EL pressure force         = ', el_pressure_force
    WRITE(mfile,'(A,L1)') 'EL pressure feedback      = ', el_pressure_feedback
    WRITE(mfile,'(A,A)') 'EL lift model             = ', TRIM(el_lift_model)
    WRITE(mfile,'(A,A)') 'EL inertial lift          = ', TRIM(el_inertial_lift)
    WRITE(mfile,'(A,ES14.6)') 'EL inertial lift umax     = ', el_inertial_lift_umax
    WRITE(mfile,'(A,A)') 'EL prescribed field       = ', TRIM(el_prescribed_field)
    WRITE(mfile,'(A,ES14.6)') 'EL shear rate             = ', el_shear_rate
    WRITE(mfile,'(A,ES14.6)') 'EL prescribed umax        = ', el_prescribed_umax
    WRITE(mfile,'(A,A)') 'EL domain type            = ', TRIM(el_domain_type)
    WRITE(mfile,'(A,3ES14.6)') 'EL cylinder center        = ', el_cylinder_center
    WRITE(mfile,'(A,ES14.6)') 'EL cylinder radius        = ', el_cylinder_radius
    WRITE(mfile,'(A,A)') 'EL cylinder axis          = ', TRIM(el_cylinder_axis)
    WRITE(mfile,'(A,L1)') 'EL Magnus                 = ', el_magnus
    WRITE(mfile,'(A,L1)') 'EL apply particle forces  = ', el_apply_particle_forces
    WRITE(mfile,'(A,L1)') 'EL apply fluid feedback   = ', el_apply_fluid_feedback
    WRITE(mfile,'(A,L1)') 'EL fluid gravity          = ', el_fluid_gravity
    WRITE(mfile,'(A,L1)') 'EL momentum fix           = ', el_momentum_fix
    WRITE(mfile,'(A,L1,A,I0,A)') 'EL meso filter            = ', el_meso_filter, &
      ' (bins ', el_meso_filter_bins, ')'
    WRITE(mfile,'(A,I0)') 'EL momentum audit freq    = ', el_momentum_audit_freq
    WRITE(mfile,'(A,I0)') 'EL particle log freq      = ', el_particle_log_freq
    WRITE(mfile,'(A,A)') 'EL volume averaged        = ', TRIM(el_volume_averaged)
    WRITE(mfile,'(A,A)') 'EL VA source time scheme  = ', TRIM(el_va_source_time_scheme)
    WRITE(mfile,'(A,I0)') 'EL VA source correctors   = ', el_va_source_correctors
    WRITE(mfile,'(A,A)') 'EL VA compatibility       = ', TRIM(el_va_compatibility)
    WRITE(mfile,'(A,A)') 'EL VA continuity form     = ', TRIM(el_va_continuity_form)
    WRITE(mfile,'(A,A)') 'EL momentum model         = ', TRIM(el_momentum_model)
    IF (el_mm_pop_weighted) WRITE(mfile,'(A,A)') 'EL MM pressure operator   = ', &
      TRIM(EL_MM_PRESSURE_OPERATOR_TEXT())
    WRITE(mfile,'(A,A)') 'EL fluid viscosity        = ', TRIM(el_fluid_viscosity)
    IF (el_periodic_wrap) WRITE(mfile,'(A,A)') 'EL periodic kernel        = ', &
      TRIM(el_periodic_kernel)
    DO i_implied=1,el_n_implied
      WRITE(mfile,'(A,A)') '  implied by the model    : ', TRIM(el_mm_implied(i_implied))
    END DO
    WRITE(mfile,'(A,ES14.6)') 'EL tavg window            = ', el_tavg_window
    WRITE(mfile,'(A,3ES14.6)') 'EL inlet u0 / u1 / ramp   = ', el_inlet_u0, el_inlet_u1, el_inlet_tramp
    CALL EL_PRINT_INLET_DEFINITION(mfile)
    IF (el_bed_audit_freq.GT.0) WRITE(mfile,'(A,I0)') 'EL bed audit freq         = ', &
      el_bed_audit_freq
    IF (el_velocity_audit_freq.GT.0) WRITE(mfile,'(A,I0)') 'EL velocity audit freq    = ', &
      el_velocity_audit_freq
    IF (el_momentum_balance_audit) WRITE(mfile,'(A,L1)') 'EL momentum balance audit = ', &
      el_momentum_balance_audit

    IF (mterm.NE.mfile) THEN
      WRITE(mterm,'(A,A)') 'EL kernel                 = ', TRIM(el_kernel)
      WRITE(mterm,'(A,ES14.6)') 'EL kernel width factor    = ', el_kernel_width_factor
      WRITE(mterm,'(A,ES14.6)') 'EL epsilon minimum        = ', el_eps_f_min
      WRITE(mterm,'(A,ES14.6)') 'EL epsilon relaxation     = ', el_eps_f_relax
      WRITE(mterm,'(A,A)') 'EL drag model             = ', TRIM(el_drag_model)
      WRITE(mterm,'(A,A)') 'EL drag coupling          = ', TRIM(el_drag_coupling)
      IF (el_drag_hrz) WRITE(mterm,'(A,A)') 'EL drag implicit weights  = ', TRIM(el_drag_implicit_weights)
      WRITE(mterm,'(A,L1)') 'EL pressure force         = ', el_pressure_force
      WRITE(mterm,'(A,L1)') 'EL pressure feedback      = ', el_pressure_feedback
      WRITE(mterm,'(A,A)') 'EL lift model             = ', TRIM(el_lift_model)
      WRITE(mterm,'(A,A)') 'EL inertial lift          = ', TRIM(el_inertial_lift)
      WRITE(mterm,'(A,ES14.6)') 'EL inertial lift umax     = ', el_inertial_lift_umax
      WRITE(mterm,'(A,A)') 'EL prescribed field       = ', TRIM(el_prescribed_field)
      WRITE(mterm,'(A,ES14.6)') 'EL shear rate             = ', el_shear_rate
      WRITE(mterm,'(A,ES14.6)') 'EL prescribed umax        = ', el_prescribed_umax
      WRITE(mterm,'(A,A)') 'EL domain type            = ', TRIM(el_domain_type)
      WRITE(mterm,'(A,3ES14.6)') 'EL cylinder center        = ', el_cylinder_center
      WRITE(mterm,'(A,ES14.6)') 'EL cylinder radius        = ', el_cylinder_radius
      WRITE(mterm,'(A,A)') 'EL cylinder axis          = ', TRIM(el_cylinder_axis)
      WRITE(mterm,'(A,L1)') 'EL Magnus                 = ', el_magnus
      WRITE(mterm,'(A,L1)') 'EL apply particle forces  = ', el_apply_particle_forces
      WRITE(mterm,'(A,L1)') 'EL apply fluid feedback   = ', el_apply_fluid_feedback
      WRITE(mterm,'(A,L1)') 'EL fluid gravity          = ', el_fluid_gravity
      WRITE(mterm,'(A,L1)') 'EL momentum fix           = ', el_momentum_fix
      WRITE(mterm,'(A,L1,A,I0,A)') 'EL meso filter            = ', el_meso_filter, &
        ' (bins ', el_meso_filter_bins, ')'
      WRITE(mterm,'(A,I0)') 'EL momentum audit freq    = ', el_momentum_audit_freq
      WRITE(mterm,'(A,I0)') 'EL particle log freq      = ', el_particle_log_freq
      WRITE(mterm,'(A,A)') 'EL volume averaged        = ', TRIM(el_volume_averaged)
      WRITE(mterm,'(A,A)') 'EL VA source time scheme  = ', TRIM(el_va_source_time_scheme)
      WRITE(mterm,'(A,I0)') 'EL VA source correctors   = ', el_va_source_correctors
      WRITE(mterm,'(A,A)') 'EL VA compatibility       = ', TRIM(el_va_compatibility)
      WRITE(mterm,'(A,A)') 'EL VA continuity form     = ', TRIM(el_va_continuity_form)
      WRITE(mterm,'(A,A)') 'EL momentum model         = ', TRIM(el_momentum_model)
      IF (el_mm_pop_weighted) WRITE(mterm,'(A,A)') 'EL MM pressure operator   = ', &
        TRIM(EL_MM_PRESSURE_OPERATOR_TEXT())
      WRITE(mterm,'(A,A)') 'EL fluid viscosity        = ', TRIM(el_fluid_viscosity)
      IF (el_periodic_wrap) WRITE(mterm,'(A,A)') 'EL periodic kernel        = ', &
        TRIM(el_periodic_kernel)
      DO i_implied=1,el_n_implied
        WRITE(mterm,'(A,A)') '  implied by the model    : ', TRIM(el_mm_implied(i_implied))
      END DO
      WRITE(mterm,'(A,ES14.6)') 'EL tavg window            = ', el_tavg_window
      WRITE(mterm,'(A,3ES14.6)') 'EL inlet u0 / u1 / ramp   = ', el_inlet_u0, el_inlet_u1, el_inlet_tramp
      CALL EL_PRINT_INLET_DEFINITION(mterm)
      IF (el_bed_audit_freq.GT.0) WRITE(mterm,'(A,I0)') 'EL bed audit freq         = ', &
        el_bed_audit_freq
      IF (el_velocity_audit_freq.GT.0) WRITE(mterm,'(A,I0)') 'EL velocity audit freq    = ', &
        el_velocity_audit_freq
      IF (el_momentum_balance_audit) WRITE(mterm,'(A,L1)') 'EL momentum balance audit = ', &
        el_momentum_balance_audit
    END IF

  END SUBROUTINE EL_PRINT_CONFIG

  SUBROUTINE EL_LOWERCASE(text)

    CHARACTER(LEN=*), INTENT(INOUT) :: text
    INTEGER :: i, code

    DO i=1,LEN_TRIM(text)
      code = IACHAR(text(i:i))
      IF (code.GE.IACHAR('A') .AND. code.LE.IACHAR('Z')) THEN
        text(i:i) = ACHAR(code + IACHAR('a') - IACHAR('A'))
      END IF
    END DO

  END SUBROUTINE EL_LOWERCASE


  ! Deck-only switch for the source/projection corrector loop. Evaluated identically on every
  ! rank, including the master that holds no EL fields: the loop contains collective pressure
  ! solves, so the fluid core must not guard it on field allocation (EL_VOLUME_AVERAGED_ACTIVE).
  LOGICAL FUNCTION EL_VA_SOURCE_CORRECTORS_ENABLED()
    EL_VA_SOURCE_CORRECTORS_ENABLED = (TRIM(el_volume_averaged).EQ.'continuity') .AND. &
                                      (el_va_source_correctors.GT.0)
  END FUNCTION EL_VA_SOURCE_CORRECTORS_ENABLED

  ! Continuity form id = index of the deck value in el_va_form_names (0 advective,
  ! 1 conservative; the deck is validated against the same table). Deck-only, identical on
  ! every rank; the single lookup behind EL_VA_CONSERVATIVE_FORM and the form dispatch of
  ! EL_VA_CONTINUITY.
  ! Radius of the kernel support of a particle: the one definition (kernel integration and
  ! deposit, halo graph and send list, periodic images). Callers must not spell the
  ! expression out: sender and receiver of a record decide on the same number.
  PURE REAL*8 FUNCTION EL_KERNEL_SUPPORT_RADIUS(particle_radius)
    REAL*8, INTENT(IN) :: particle_radius
    EL_KERNEL_SUPPORT_RADIUS = 2.0d0*el_kernel_width_factor*particle_radius
  END FUNCTION EL_KERNEL_SUPPORT_RADIUS

  ! ELDragImplicitWeights = hrz puts the explicit partner dt B_i u_i^n of the drag diagonal on
  ! the right-hand side of the momentum equation (EL_FIELDS:EL_APPLY_DRAG_SPLIT_PARTNER); the
  ! nonlinear defect of the momentum solve must then subtract the matching theta dt B_i u_i, or
  ! the partner stays in the residual of an unchanged field. The defects built from the
  ! assembled matrix (Newtonian path; non-Newtonian path with SimPar@MatrixRenewal S >= 1)
  ! carry it through the matrix diagonal (QuadSc_def: drag_B_source). The hand-built
  ! non-Newtonian defect of S = 0 (mass, convection, STRESS; QuadSc_def Matdef idef = -1) has no
  ! drag term: that combination is not supported (review r2 finding 4; the missing defect term
  ! is a follow-up, no test case exists).
  ! The check concerns runs that execute the momentum solve: a prescribed field
  ! (EL_FLUID_SOLVE_ACTIVE false) never reaches the defect, and the application accepts it.
  PURE LOGICAL FUNCTION EL_DRAG_HRZ_DEFECT_SUPPORTED(non_newtonian, matrix_renewal_s)
    LOGICAL, INTENT(IN) :: non_newtonian
    INTEGER, INTENT(IN) :: matrix_renewal_s
    EL_DRAG_HRZ_DEFECT_SUPPORTED = .NOT.(non_newtonian .AND. matrix_renewal_s.EQ.0)
  END FUNCTION EL_DRAG_HRZ_DEFECT_SUPPORTED

  ! Backstop in the momentum assembly of an hrz run (Transport_Q2P1, workers only): a FlowType or
  ! MatrixRenewal changed by Reload_SimPar (the generic parser reads both) after the start-up check.
  ! An unsupported combination aborts MPI_COMM_WORLD, so every rank, the master included, stops at
  ! once whichever worker detects it (a STOP inside the workers-only block would leave the master
  ! waiting in its next collective).
  SUBROUTINE EL_DRAG_HRZ_RUNTIME_CHECK(non_newtonian, matrix_renewal_s)
    LOGICAL, INTENT(IN) :: non_newtonian
    INTEGER, INTENT(IN) :: matrix_renewal_s
    INTEGER :: ierr
    IF (EL_DRAG_HRZ_DEFECT_SUPPORTED(non_newtonian, matrix_renewal_s)) RETURN
    WRITE(*,'(A)') 'Fatal: ELDragImplicitWeights = hrz with FlowType = non-Newtonian and MatrixRenewal S = 0 '// &
      '(changed by Reload_SimPar?): the momentum defect has no drag diagonal; aborting all ranks.'
    CALL MPI_Abort(MPI_COMM_WORLD, 917, ierr)
  END SUBROUTINE EL_DRAG_HRZ_RUNTIME_CHECK

  ! Solver settings the momentum-balance audit covers (EL_MOMENTUM_BALANCE splits the momentum
  ! defect into lumped mass, diffusion matrix D, convection matrix K and drag diagonal): the
  ! backward-Euler Newtonian path whose matrix is M_l + dt (D + K + B), every operator rebuilt
  ! each step. The one definition for the start-up refusal of the application and the run-time
  ! backstop below (Reload_SimPar reads TimeScheme, FlowType and MatrixRenewal again).
  PURE LOGICAL FUNCTION EL_MOMBAL_SOLVER_SUPPORTED(theta, non_newtonian, viscoelastic, &
                                                   ns_stabilization, gamma, renewal_m, renewal_d, renewal_k)
    REAL*8, INTENT(IN) :: theta, gamma
    LOGICAL, INTENT(IN) :: non_newtonian, viscoelastic, ns_stabilization
    INTEGER, INTENT(IN) :: renewal_m, renewal_d, renewal_k
    EL_MOMBAL_SOLVER_SUPPORTED = theta.EQ.1.0d0 .AND. gamma.EQ.0.0d0 .AND. &
      .NOT.(non_newtonian .OR. viscoelastic .OR. ns_stabilization) .AND. &
      renewal_m.GE.1 .AND. renewal_d.GE.1 .AND. renewal_k.GE.1
  END FUNCTION EL_MOMBAL_SOLVER_SUPPORTED

  ! Backstop at the start of every audited fluid step (Transport_Q2P1, workers only): a setting
  ! changed by Reload_SimPar after the start-up check. Aborts MPI_COMM_WORLD for the reason given
  ! at EL_DRAG_HRZ_RUNTIME_CHECK.
  SUBROUTINE EL_MOMBAL_RUNTIME_CHECK(theta, non_newtonian, viscoelastic, ns_stabilization, &
                                     gamma, renewal_m, renewal_d, renewal_k)
    REAL*8, INTENT(IN) :: theta, gamma
    LOGICAL, INTENT(IN) :: non_newtonian, viscoelastic, ns_stabilization
    INTEGER, INTENT(IN) :: renewal_m, renewal_d, renewal_k
    INTEGER :: ierr
    IF (EL_MOMBAL_SOLVER_SUPPORTED(theta, non_newtonian, viscoelastic, ns_stabilization, &
                                   gamma, renewal_m, renewal_d, renewal_k)) RETURN
    WRITE(*,'(A)') 'Fatal: ELMomentumBalanceAudit = Yes with a time scheme, flow type, stabilisation or '// &
      'MatrixRenewal the balance does not cover (changed by Reload_SimPar?); aborting all ranks.'
    CALL MPI_Abort(MPI_COMM_WORLD, 918, ierr)
  END SUBROUTINE EL_MOMBAL_RUNTIME_CHECK

  ! The fluid momentum solve runs: Transport_q2p1_UxyzP_el skips it for a prescribed field
  ! (SimPar@ELPrescribedField /= none). The one definition for the transport, the start-up checks
  ! that concern the solve and the solver terms of the velocity audit.
  PURE LOGICAL FUNCTION EL_FLUID_SOLVE_ACTIVE()
    EL_FLUID_SOLVE_ACTIVE = (TRIM(el_prescribed_field).EQ.'none')
  END FUNCTION EL_FLUID_SOLVE_ACTIVE

  ! Inlet velocity of boundary type Inflow975 at time t (end-of-step time of the caller): the
  ! table when the deck names one, else the linear ramp (arithmetic unchanged).
  PURE REAL*8 FUNCTION EL_INLET_VELOCITY(t)
    REAL*8, INTENT(IN) :: t
    IF (el_inlet_table%n.GT.0) THEN
      EL_INLET_VELOCITY = EL_INLET_TABLE_VALUE(el_inlet_table, t)
    ELSE
      EL_INLET_VELOCITY = el_inlet_u0 + (el_inlet_u1 - el_inlet_u0)* &
                          MIN(1.0d0, MAX(0.0d0, t)/el_inlet_tramp)
    END IF
  END FUNCTION EL_INLET_VELOCITY

  ! True when the deck configures the Inflow975 inlet (either ramp velocity nonzero, or a
  ! table). Deck-only, the same on every rank.
  PURE LOGICAL FUNCTION EL_INLET_CONFIGURED()
    EL_INLET_CONFIGURED = (el_inlet_u0.NE.0.0d0 .OR. el_inlet_u1.NE.0.0d0) .OR. &
                          (LEN_TRIM(el_inlet_table_file).GT.0)
  END FUNCTION EL_INLET_CONFIGURED

  ! SimPar@ELInletVelocityTable: refuse an empty path and the table together with a ramp key
  ! (deck-only, every rank alike), then load the file: rank 0 of MPI_COMM_WORLD reads one
  ! snapshot, status and rows are broadcast (EL_INLET_TABLE_LOAD_SHARED; every rank parses the
  ! deck, so every rank is here). A missing, unreadable, malformed, non-finite or unsorted table
  ! stops every rank with the same message. Without the key nothing happens (no collective).
  SUBROUTINE EL_LOAD_INLET_TABLE()
    INTEGER :: ierr
    CHARACTER(LEN=320) :: message

    IF (LEN_TRIM(el_inlet_table_file).EQ.0) THEN
      IF (EL_KEY_SET('ELInletVelocityTable')) THEN
        WRITE(*,'(A)') 'ELInletVelocityTable: empty file name (remove the key for the ramp).'
        STOP 1
      END IF
      RETURN
    END IF
    IF (EL_KEY_SET('ELInletVelocityStart') .OR. EL_KEY_SET('ELInletVelocityEnd') .OR. &
        EL_KEY_SET('ELInletRampTime')) THEN
      WRITE(*,'(A)') 'ELInletVelocityTable cannot be combined with ELInletVelocityStart / '// &
        'ELInletVelocityEnd / ELInletRampTime (one inlet definition only).'
      STOP 1
    END IF
    CALL EL_INLET_TABLE_LOAD_SHARED(el_inlet_table_file, el_inlet_table, ierr, message, &
                                    MPI_COMM_WORLD)
    IF (ierr.NE.0) THEN
      WRITE(*,'(A)') TRIM(message)
      STOP 1
    END IF
  END SUBROUTINE EL_LOAD_INLET_TABLE

  ! Echo of the inlet definition after the ramp line of EL_PRINT_CONFIG: for a table the file,
  ! the row count, the start-up notice and every row; for a configured ramp its start-up
  ! notice; nothing without an inlet.
  SUBROUTINE EL_PRINT_INLET_DEFINITION(unit)
    INTEGER, INTENT(IN) :: unit
    INTEGER :: i
    IF (el_inlet_table%n.LE.0) THEN
      IF (EL_INLET_CONFIGURED()) WRITE(unit,'(A)') 'EL inlet velocity ramp: start-up setting; '// &
        'a reload of the parameter file (Reload_SimPar) does not change it'
      RETURN
    END IF
    WRITE(unit,'(A,A,A,I0,A)') 'EL inlet velocity table   = ', TRIM(el_inlet_table_file), &
      ' (', el_inlet_table%n, ' rows t / U; linear, constant beyond the ends)'
    WRITE(unit,'(A)') 'EL inlet velocity table: start-up setting, read once here; a reload of '// &
      'the parameter file (Reload_SimPar) neither re-reads nor replaces it'
    DO i=1,el_inlet_table%n
      WRITE(unit,'(A,2ES16.8)') '  EL_INLET_TABLE ', el_inlet_table%t(i), el_inlet_table%u(i)
    END DO
  END SUBROUTINE EL_PRINT_INLET_DEFINITION

  ! The inlet definition the run uses (start-up values of the EL keys).
  FUNCTION EL_INLET_STARTUP_DEFINITION() RESULT(def)
    TYPE(tElInletDefinition) :: def
    def%table = (LEN_TRIM(el_inlet_table_file).GT.0)
    def%path = el_inlet_table_file
    def%ramp_keys = EL_KEY_SET('ELInletVelocityStart') .OR. EL_KEY_SET('ELInletVelocityEnd') .OR. &
                    EL_KEY_SET('ELInletRampTime')
    def%u0 = el_inlet_u0
    def%u1 = el_inlet_u1
    def%tramp = el_inlet_tramp
  END FUNCTION EL_INLET_STARTUP_DEFINITION

  ! The inlet definition a parameter file asks for, read the way the application's deck parser
  ! reads it (myGDATNEW: a line 'SimPar@<key> = <value>' with the last '@' and the last '=';
  ! a later line overrides an earlier one; the table value through EL_INLET_TABLE_PATH, the
  ! ramp values list-directed). ierr /= 0: the file cannot be opened (def is the default).
  SUBROUTINE EL_INLET_DEFINITION_OF_FILE(filename, def, ierr)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    TYPE(tElInletDefinition), INTENT(OUT) :: def
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=500) :: string
    CHARACTER(LEN=7) :: cvar
    CHARACTER(LEN=32) :: cpar
    INTEGER :: iunit, ios, iat, ieq, i
    REAL*8 :: val

    OPEN(NEWUNIT=iunit, FILE=TRIM(filename), STATUS='OLD', ACTION='READ', IOSTAT=ierr)
    IF (ierr.NE.0) RETURN
    DO
      READ(iunit, '(A500)', IOSTAT=ios) string
      IF (ios.NE.0) EXIT
      iat = 0
      ieq = 0
      DO i = 1, LEN(string)
        IF (string(i:i).EQ.'@') iat = i
        IF (string(i:i).EQ.'=') ieq = i
      END DO
      IF (iat.LE.1 .OR. ieq.LE.iat + 1) CYCLE
      cvar = ' '
      cpar = ' '
      READ(string(1:iat-1), *, IOSTAT=ios) cvar
      IF (ios.NE.0 .OR. TRIM(ADJUSTL(cvar)).NE.'SimPar') CYCLE
      READ(string(iat+1:ieq-1), *, IOSTAT=ios) cpar
      IF (ios.NE.0) CYCLE
      SELECT CASE (TRIM(ADJUSTL(cpar)))
      CASE ('ELInletVelocityTable')
        def%table = .TRUE.
        def%path = EL_INLET_TABLE_PATH(string(ieq+1:))
      CASE ('ELInletVelocityStart', 'ELInletVelocityEnd', 'ELInletRampTime')
        def%ramp_keys = .TRUE.
        READ(string(ieq+1:), *, IOSTAT=ios) val
        IF (ios.NE.0) THEN
          def%readable = .FALSE.
        ELSE IF (TRIM(ADJUSTL(cpar)).EQ.'ELInletVelocityStart') THEN
          def%u0 = val
        ELSE IF (TRIM(ADJUSTL(cpar)).EQ.'ELInletVelocityEnd') THEN
          def%u1 = val
        ELSE
          def%tramp = val
        END IF
      END SELECT
    END DO
    CLOSE(iunit)
    ierr = 0
  END SUBROUTINE EL_INLET_DEFINITION_OF_FILE

  ! How the requested inlet definition differs from the one of the run; blank when it does not.
  FUNCTION EL_INLET_DEFINITION_CHANGE(run, requested) RESULT(text)
    TYPE(tElInletDefinition), INTENT(IN) :: run, requested
    CHARACTER(LEN=400) :: text
    text = ' '
    IF (requested%table .AND. requested%ramp_keys) THEN
      text = 'the reload file sets ELInletVelocityTable together with ramp keys (refused at start-up)'
    ELSE IF (run%table .AND. requested%table) THEN
      IF (TRIM(run%path).NE.TRIM(requested%path)) text = 'table '//TRIM(run%path)// &
        ' -> table '//TRIM(requested%path)
    ELSE IF (run%table) THEN
      IF (requested%ramp_keys) THEN
        text = 'table '//TRIM(run%path)//' -> ramp keys'
      ELSE
        text = 'table '//TRIM(run%path)//' removed'
      END IF
    ELSE IF (requested%table) THEN
      text = 'ramp -> table '//TRIM(requested%path)
    ELSE IF (run%ramp_keys .AND. .NOT.requested%ramp_keys) THEN
      text = 'ramp keys removed'
    ELSE IF (.NOT.requested%readable) THEN
      text = 'ramp keys with an unreadable value'
    ELSE IF (run%u0.NE.requested%u0 .OR. run%u1.NE.requested%u1 .OR. &
             run%tramp.NE.requested%tramp) THEN
      text = 'ramp values changed'
    END IF
  END FUNCTION EL_INLET_DEFINITION_CHANGE

  ! Reload_SimPar of filename (ProcCtrl): report on unit mfile and mterm (when tell) an inlet
  ! definition of the file that differs from the one of the run, which is kept. Rank-local
  ! (only the telling rank reads the file), no collective; silent without a validated EL
  ! configuration (other applications) and when the file cannot be opened.
  SUBROUTINE EL_INLET_RELOAD_NOTICE(filename, mfile, mterm, tell)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: mfile, mterm
    LOGICAL, INTENT(IN) :: tell
    TYPE(tElInletDefinition) :: run, requested
    CHARACTER(LEN=400) :: change
    CHARACTER(LEN=300) :: kept
    INTEGER :: ierr

    IF (.NOT.(tell .AND. el_config_validated)) RETURN
    CALL EL_INLET_DEFINITION_OF_FILE(filename, requested, ierr)
    IF (ierr.NE.0) RETURN
    run = EL_INLET_STARTUP_DEFINITION()
    change = EL_INLET_DEFINITION_CHANGE(run, requested)
    IF (LEN_TRIM(change).EQ.0) RETURN
    IF (run%table) THEN
      kept = 'the table '//TRIM(run%path)//' as read at start-up'
    ELSE
      WRITE(kept,'(A,3ES14.6)') 'the ramp u0 / u1 / ramp time ', run%u0, run%u1, run%tramp
    END IF
    CALL TELL_UNIT(mfile)
    IF (mterm.NE.mfile) CALL TELL_UNIT(mterm)
  CONTAINS
    SUBROUTINE TELL_UNIT(unit)
      INTEGER, INTENT(IN) :: unit
      WRITE(unit,'(A,A,A,A)') 'WARNING: EL inlet definition of the reload file ', TRIM(filename), &
        ' ignored: ', TRIM(change)
      WRITE(unit,'(A,A)') '  the inlet is a start-up setting; the run keeps ', TRIM(kept)
    END SUBROUTINE TELL_UNIT
  END SUBROUTINE EL_INLET_RELOAD_NOTICE

  ! True when the inlet value has to be converted from superficial to interstitial velocity
  ! (deck-only: the same on every rank; decides the nodal-porosity collective of the transfer
  ! pass).
  PURE LOGICAL FUNCTION EL_INLET_SUPERFICIAL()
    EL_INLET_SUPERFICIAL = EL_INLET_CONFIGURED() .AND. (TRIM(el_volume_averaged).EQ.'continuity')
  END FUNCTION EL_INLET_SUPERFICIAL

  INTEGER FUNCTION EL_VA_FORM_ID()
    INTEGER :: m
    EL_VA_FORM_ID = 0
    IF (TRIM(el_volume_averaged).NE.'continuity') RETURN
    DO m = 0, UBOUND(el_va_form_names,1)
      IF (TRIM(el_va_continuity_form).EQ.TRIM(el_va_form_names(m))) EL_VA_FORM_ID = m
    END DO
  END FUNCTION EL_VA_FORM_ID

  ! Deck-only switch for the conservative continuity form (EL_VA_FLUX): decides the nodal
  ! reconstruction and its collective in the transfer pass and the flux-moment path of the
  ! source; identical on every rank.
  ! Key set of a momentum model. Called by EL_VALIDATE_CONFIG after the name is known to be
  ! valid and implemented, before any other key is validated. For 'none' nothing happens (the
  ! campaign keys keep their deck values). For a model every key of its consistent set is
  ! assigned; a key the deck set explicitly to a different value is refused (el_keys_set),
  ! and the settings a model cannot be combined with are refused (all conflicts are listed
  ! before the one stop). The values come from the model table at the top of the module. Checks that need solver
  ! state (matrix renewal, GAMMA, material law) are made by the fluid adapter at its first
  ! call (Transport_Q2P1: EL_MM_Check_Solver_Settings).
  !
  !   key                    bfull          b              a
  !   ELVolumeAveraged       continuity     continuity     continuity
  !   ELVAContinuityForm     conservative   conservative   conservative
  !   ELVASourceCorrectors   >= 1 (1 if the deck does not set it)
  !   ELApplyFluidFeedback   Yes            Yes            Yes
  !   ELDragCoupling         semi_implicit  semi_implicit  semi_implicit
  !   ELPressureForce        Yes            No             Yes
  !   ELPressureFeedback     Yes            No             No
  !   ELFluidGravity         No             No             No
  ! ELFluidGravity = No is the dynamic-pressure formulation all three share here: the fluid
  ! carries no gravity term, the particle the buoyant weight (rho_p - rho_f) V g, and the
  ! pressure force uses the gradient of that pressure; with a hydrostatic pressure in the
  ! field the buoyancy would be counted twice. SimPar@ConstantForcing stays available (a fluid
  ! acceleration, weighted with eps like the inertia).
  SUBROUTINE EL_RESOLVE_MOMENTUM_MODEL()

    INTEGER :: id

    id = EL_MM_ID()
    IF (.NOT.ANY(TRIM(el_mm_pressure_operator).EQ.el_mm_pop_names)) THEN
      WRITE(*,'(A,A)') 'Invalid ELMMPressureOperator (startup|weighted): ', &
        TRIM(el_mm_pressure_operator)
      STOP 1
    END IF
    el_mm_pop_weighted = (TRIM(el_mm_pressure_operator).EQ.TRIM(el_mm_pop_names(1)))
    el_mm_pressure_rebuild_on = el_mm_pop_weighted .AND. el_mmt_pressure_rebuild(id)
    IF (id.EQ.EL_MM_NONE) THEN
      IF (el_mm_pop_weighted) THEN
        WRITE(*,'(A)') 'ELMMPressureOperator = weighted needs a momentum model (ELMomentumModel).'
        STOP 1
      END IF
      RETURN
    END IF
    el_n_implied = 0
    el_mm_conflicts = 0

    CALL EL_MM_IMPLY_TEXT('ELVolumeAveraged', el_volume_averaged, 'continuity')
    CALL EL_MM_IMPLY_TEXT('ELVAContinuityForm', el_va_continuity_form, 'conservative')
    CALL EL_MM_IMPLY_TEXT('ELDragCoupling', el_drag_coupling, 'semi_implicit')
    CALL EL_MM_IMPLY_LOGICAL('ELApplyFluidFeedback', el_apply_fluid_feedback, .TRUE.)
    CALL EL_MM_IMPLY_LOGICAL('ELPressureForce', el_pressure_force, el_mmt_particle_pressure(id))
    CALL EL_MM_IMPLY_LOGICAL('ELPressureFeedback', el_pressure_feedback, &
                             el_mmt_pressure_feedback(id))
    CALL EL_MM_IMPLY_LOGICAL('ELFluidGravity', el_fluid_gravity, .FALSE.)
    IF (.NOT.EL_KEY_SET('ELVASourceCorrectors')) THEN
      el_va_source_correctors = 1
      CALL EL_MM_NOTE_IMPLIED('ELVASourceCorrectors = 1')
    END IF

    ! settings a model cannot be combined with; every conflict is reported before the stop
    ! (an unknown drag name is refused by the drag validation after this resolution)
    IF (EL_DRAG_LOOKUP(el_drag_model).GE.0) THEN
      IF (.NOT.el_drag_standard_form(EL_DRAG_LOOKUP(el_drag_model))) CALL EL_MM_CONFLICT( &
        'needs a drag law in the standard porosity convention; ELDragModel = '// &
        TRIM(el_drag_model)//' (campaign form without eps^2) is refused, use difelice_std.')
    END IF
    IF (el_meso_filter) CALL EL_MM_CONFLICT('is not supported with ELMesoFilter = Yes.')
    IF (el_momentum_fix) CALL EL_MM_CONFLICT('is not supported with ELMomentumFix = Yes '// &
      '(the compensator acts on the unweighted fluid momentum).')
    IF (el_eps_f_relax.GT.0.0d0) CALL EL_MM_CONFLICT('needs ELEpsFRelax = 0.')
    IF (TRIM(el_prescribed_field).NE.'none') CALL EL_MM_CONFLICT( &
      'is not supported with a prescribed fluid field.')
    ! the weak kernel average of the viscous force assumes a kernel that vanishes at the edge
    ! of its support; the gaussian is cut at a finite value there (an 8 % low Laplacian)
    IF (el_mmt_viscous_force(id) .AND. TRIM(el_kernel).EQ.'gaussian') CALL EL_MM_CONFLICT( &
      'with a particle viscous force needs ELKernel = deen_poly (the truncated gaussian '// &
      'kernel is not supported by the viscous average).')
    IF (el_mm_conflicts.GT.0) STOP 1

  END SUBROUTINE EL_RESOLVE_MOMENTUM_MODEL

  ! Echo text of SimPar@ELMMPressureOperator = weighted for the selected model.
  FUNCTION EL_MM_PRESSURE_OPERATOR_TEXT() RESULT(text)
    CHARACTER(LEN=96) :: text
    text = 'weighted (pressure matrix rebuilt from the porosity every step)'
    IF (.NOT.el_mm_pressure_rebuild_on) &
      text = 'weighted accepted; start-up preconditioner retained for this model'
  END FUNCTION EL_MM_PRESSURE_OPERATOR_TEXT

  ! Kinematic viscosity for the closures and diagnostics of the Euler-Lagrange path. The
  ! argument is Properties%Viscosity(1): the kinematic value itself under 'kinematic' (returned
  ! unchanged), the dynamic one under ELFluidViscosity = dynamic, where the kinematic value
  ! recorded at the conversion is returned (not mu/rho: the closures see the same bits in
  ! both conventions). Under 'dynamic' the argument must be the recorded mu: a property that
  ! was never converted, or changed afterwards (Reload_Prop), stops the run here instead of
  ! running with mixed conventions.
  REAL*8 FUNCTION EL_CLOSURE_VISCOSITY(properties_viscosity)
    REAL*8, INTENT(IN) :: properties_viscosity
    IF (.NOT.el_fluid_visc_dynamic) THEN
      EL_CLOSURE_VISCOSITY = properties_viscosity
    ELSE
      IF (el_closure_nu.LE.0.0d0) THEN
        WRITE(*,*) 'EL_CLOSURE_VISCOSITY: ELFluidViscosity = dynamic, but the property was ', &
                   'not converted (EL_APPLY_FLUID_VISCOSITY not called).'
        STOP 1
      END IF
      IF (properties_viscosity.NE.el_operator_mu) THEN
        WRITE(*,*) 'EL_CLOSURE_VISCOSITY: ELFluidViscosity = dynamic, but the fluid viscosity ', &
                   'changed after the conversion (property reload is not supported): ', &
                   properties_viscosity, ' /= ', el_operator_mu
        STOP 1
      END IF
      EL_CLOSURE_VISCOSITY = el_closure_nu
    END IF
  END FUNCTION EL_CLOSURE_VISCOSITY

  ! Record the conversion (the application's EL_APPLY_FLUID_VISCOSITY; unit tests): nu for the
  ! closures, mu = rho nu for the operator. Returns mu, the value to store in the property.
  REAL*8 FUNCTION EL_SET_FLUID_VISCOSITY(nu, rho)
    REAL*8, INTENT(IN) :: nu, rho
    el_closure_nu = nu
    el_operator_mu = rho*nu
    EL_SET_FLUID_VISCOSITY = el_operator_mu
  END FUNCTION EL_SET_FLUID_VISCOSITY

  LOGICAL FUNCTION EL_FLUID_VISCOSITY_CONVERTED()
    EL_FLUID_VISCOSITY_CONVERTED = (el_closure_nu.GT.0.0d0)
  END FUNCTION EL_FLUID_VISCOSITY_CONVERTED

  SUBROUTINE EL_MM_CONFLICT(text)
    CHARACTER(LEN=*), INTENT(IN) :: text
    WRITE(*,'(4A)') 'ELMomentumModel = ', TRIM(el_momentum_model), ' ', text
    el_mm_conflicts = el_mm_conflicts + 1
  END SUBROUTINE EL_MM_CONFLICT

  SUBROUTINE EL_MM_NOTE_IMPLIED(line)
    CHARACTER(LEN=*), INTENT(IN) :: line
    IF (el_n_implied.GE.el_max_implied) RETURN
    el_n_implied = el_n_implied + 1
    el_mm_implied(el_n_implied) = line
  END SUBROUTINE EL_MM_NOTE_IMPLIED

  ! Assign the model's value to a text key; refuse a different explicit deck value.
  SUBROUTINE EL_MM_IMPLY_TEXT(key, var, value)
    CHARACTER(LEN=*), INTENT(IN) :: key, value
    CHARACTER(LEN=*), INTENT(INOUT) :: var
    IF (EL_KEY_SET(key) .AND. TRIM(var).NE.TRIM(value)) THEN
      WRITE(*,'(6A)') 'ELMomentumModel = ', TRIM(el_momentum_model), ' implies ', key, &
        ' = ', value
      WRITE(*,'(4A)') '  the parameter file sets ', key, ' = ', TRIM(var)
      el_mm_conflicts = el_mm_conflicts + 1
      RETURN
    END IF
    IF (TRIM(var).NE.TRIM(value) .OR. .NOT.EL_KEY_SET(key)) &
      CALL EL_MM_NOTE_IMPLIED(key//' = '//value)
    var = value
  END SUBROUTINE EL_MM_IMPLY_TEXT

  ! Assign the model's value to a Yes/No key; refuse a different explicit deck value.
  SUBROUTINE EL_MM_IMPLY_LOGICAL(key, var, value)
    CHARACTER(LEN=*), INTENT(IN) :: key
    LOGICAL, INTENT(INOUT) :: var
    LOGICAL, INTENT(IN) :: value
    CHARACTER(LEN=3) :: text
    text = 'No'
    IF (value) text = 'Yes'
    IF (EL_KEY_SET(key) .AND. (var.NEQV.value)) THEN
      WRITE(*,'(6A)') 'ELMomentumModel = ', TRIM(el_momentum_model), ' implies ', key, &
        ' = ', TRIM(text)
      WRITE(*,'(3A)') '  the parameter file sets ', key, ' to the opposite value.'
      el_mm_conflicts = el_mm_conflicts + 1
      RETURN
    END IF
    IF ((var.NEQV.value) .OR. .NOT.EL_KEY_SET(key)) &
      CALL EL_MM_NOTE_IMPLIED(key//' = '//TRIM(text))
    var = value
  END SUBROUTINE EL_MM_IMPLY_LOGICAL

  ! Id of the momentum model of the validated configuration (index of the model table;
  ! EL_MM_NONE before EL_VALIDATE_CONFIG has run). Deck-only, the same on every rank.
  INTEGER FUNCTION EL_MM_ID()
    EL_MM_ID = el_mm_id_validated
  END FUNCTION EL_MM_ID

  ! Index of a (lower-case) model name in the table, -1 if it is not a model name.
  INTEGER FUNCTION EL_MM_LOOKUP(name)
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: i
    EL_MM_LOOKUP = -1
    DO i=LBOUND(el_mm_names,1),UBOUND(el_mm_names,1)
      IF (TRIM(name).EQ.TRIM(el_mm_names(i))) EL_MM_LOOKUP = i
    END DO
  END FUNCTION EL_MM_LOOKUP

  ! Id of a (lower-case) drag closure name in the drag table, -1 if it is not a closure name.
  PURE INTEGER FUNCTION EL_DRAG_LOOKUP(name)
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: i
    EL_DRAG_LOOKUP = -1
    DO i=LBOUND(el_drag_names,1),UBOUND(el_drag_names,1)
      IF (TRIM(name).EQ.TRIM(el_drag_names(i))) THEN
        EL_DRAG_LOOKUP = i
        RETURN
      END IF
    END DO
  END FUNCTION EL_DRAG_LOOKUP

  ! Resolve the closure id of el_drag_model once (EL_VALIDATE_CONFIG refuses -1); the per-particle
  ! dispatch reads el_drag_closure_id, never the name.
  SUBROUTINE EL_RESOLVE_DRAG_CLOSURE()
    el_drag_closure_id = EL_DRAG_LOOKUP(el_drag_model)
  END SUBROUTINE EL_RESOLVE_DRAG_CLOSURE

  ! Record that the parameter file set the SimPar@EL* key `name` (called by the parser for
  ! every such key, before EL_VALIDATE_CONFIG).
  SUBROUTINE EL_NOTE_KEY_SET(name)
    CHARACTER(LEN=*), INTENT(IN) :: name
    IF (EL_KEY_SET(name)) RETURN
    IF (el_n_keys_set.GE.el_max_keys_set) THEN
      WRITE(*,*) 'EL_NOTE_KEY_SET: more than el_max_keys_set distinct SimPar@EL* keys.'
      STOP 1
    END IF
    el_n_keys_set = el_n_keys_set + 1
    el_keys_set(el_n_keys_set) = name
  END SUBROUTINE EL_NOTE_KEY_SET

  LOGICAL FUNCTION EL_KEY_SET(name)
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: i
    EL_KEY_SET = .FALSE.
    DO i=1,el_n_keys_set
      IF (TRIM(el_keys_set(i)).EQ.TRIM(name)) EL_KEY_SET = .TRUE.
    END DO
  END FUNCTION EL_KEY_SET

  LOGICAL FUNCTION EL_VA_CONSERVATIVE_FORM()
    EL_VA_CONSERVATIVE_FORM = (EL_VA_FORM_ID().EQ.1)
  END FUNCTION EL_VA_CONSERVATIVE_FORM

END MODULE EL_CONFIG
