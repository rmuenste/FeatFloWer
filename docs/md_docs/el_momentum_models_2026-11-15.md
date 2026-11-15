# Euler-Lagrange momentum models A / B / Bfull (2026-11-15)

Branch `feature/el-momentum-models` (worktree `FeatFloWer-mommodels`, from
`feature/el-ergun-column`). Deck key `SimPar@ELMomentumModel = none | Bfull | B | A`, default
`none` (every certified configuration bitwise). The three names are the model types of CFDEM's
`cfdemSolverPiso`; "model A" here means the same equations as CFDEM's model A. The plan this
started from is `el_momentum_models_plan_2026-11-15.md`; where the two differ, this note holds.

## 1. Equations

eps porosity, phi = 1 - eps, u interstitial fluid velocity, F_d drag on a particle, V its
volume, n the number density, c the coefficient of the fluid's viscous operator (section 6).

| | Fluid momentum | Force on a particle |
|---|---|---|
| `none` | rho Du/Dt = -grad p + c lap u - n F_d [+ phi grad p with `ELPressureFeedback`] | F_d [- V grad p with `ELPressureForce`] |
| `Bfull` | eps rho Du/Dt = -grad p + c lap u - n (F_d + F_p + F_v) | F_d + F_p + F_v |
| `B` | eps rho Du/Dt = -grad p + c lap u - n F_d / eps | F_d / eps |
| `A` | eps rho Du/Dt = -eps grad p + div(eps c grad u) - n F_d | F_d + F_p + F_v |

F_p = -V grad p, F_v = V c lap u, both taken as kernel averages at the particle. Continuity in
all three: d eps/dt + div(eps u) = 0 in the conservative form
(`el_va_conservative_form_2026-11-05.md`). The pressure is the dynamic pressure: the fluid has no
gravity term, the particle carries (rho_p - rho_f) V g; with a uniform fluid density this is
CFDEM's full pressure plus Archimedes.

## 2. Where it lives

| Piece | File | Notes |
|---|---|---|
| key, the model table (one column per model: id, name, property columns), registry of explicit keys, implied keys and refusals | `source/src_el/el_config.f90` | `EL_RESOLVE_MOMENTUM_MODEL`, `EL_NOTE_KEY_SET`, `EL_MM_ID` |
| policy: the run-time questions, answered from the table | `source/src_el/el_momentum_model.f90` | `EL_MM_ACTIVE`, `_WEIGHTED_INERTIA`, `_DRAG_SCALE`, `_RESIDUAL_DIVISOR`, `_VISCOUS_FORCE`, `_VISCOUS_FEEDBACK`, `_OPERATOR_FORM` |
| kernel average of lap u at a particle | `source/src_el/el_viscous.f90` | weak form with wall term and consistency correction |
| solid-fraction operators of model A | `source/src_el/el_mm_operator.f90` | fp = int p div(phi_h N), fv = int phi_h grad u : grad N |
| mesh / MPI / operator work | `source/src_quadLS/QuadSc_el_mm.f90` (included in `Transport_Q2P1`) | begin-of-step (setting check, boundary-face mask, porosity priming), coefficient field, model A terms |
| force split | `source/src_el/el_forces.f90`, `el_transfer.f90` | `particle_nondrag`, `feedback_nondrag`, `viscous` |
| residual scaling | `source/src_el/el_volume_averaged.f90`, `el_va_flux.f90` | divisor w_e in the source, the right-hand side and the sink identity |
| audit | `source/src_el/el_diagnostics.f90` | `EL_FLUID_PAIR` on the eps-weighted momentum under a model; new line `EL_MM_FORCE` |

Shared (non-Euler-Lagrange) code: `QuadSc_main.f90` gains one include and four guarded calls, two
in `Transport_q2p1_UxyzP_el` and two in the fluid core inside `IF (.NOT.enable_fbm)` (the model
predicate is not evaluated on the FBM path), next to the Euler-Lagrange calls that were already there. No assembly
routine, `Create_CMat`, the multigrid, `AddPressureGradient` or `Velocity_Correction` is changed.
Applications that never call the Euler-Lagrange transport routine cannot reach the code (DNS twins).

## 3. Mechanisms

**Keys.** A model implies ELVolumeAveraged = continuity, ELVAContinuityForm = conservative,
ELVASourceCorrectors >= 1, ELApplyFluidFeedback = Yes, ELDragCoupling = semi_implicit,
ELFluidGravity = No, and ELPressureForce / ELPressureFeedback = Yes/Yes (Bfull), No/No (B), Yes/No
(A). A key the deck sets to a different value is refused; the implied values are echoed. Refused
with any model: ELDragModel = difelice (campaign form without eps^2), the meso filter, the
momentum-fix compensator, porosity relaxation, a prescribed field; with a viscous-force model
(Bfull, A) also ELKernel = gaussian. All conflicts of a deck are listed before the one stop.
Checked at the start of the first step, before any pass: MatrixRenewal M >= 1, K = 2 or 3, C = 1;
GAMMA = 0; convective form; Newtonian path without viscoelastic stress.

**First step.** A fresh run has the porosity field at 1 when the pre-advance pass samples the
forces. Under a model one deposit-only pass (label 0: no PE arming, no history advance, nothing
captured for the fluid) builds the porosity of the initial positions first, so the first particle
advance already sees the right porosity in the drag, in the 1/eps scale of model B and in the
inertia weight. `none` keeps its start-up sequence.

**ConstantForcing under a model** is a fluid acceleration and is weighted like the inertia (eps
rho f on the fluid, nothing on the particles). It is not a driving pressure gradient of any of
the three equation sets; use it as what it is.

**eps-weighted inertia (all models).** `EL_MM_Update_Fluid_Coefficients` writes the element-wise
coefficient field `mgDensity = rho * eps_f` once per step after the last particle-mesh pass,
restricts it volume-weighted to the coarser levels, hands the coarsest level to the master and
rebuilds the lumped mass. Mass matrix, convection, body forces and the velocity correction then
carry eps rho. Time level: eps^{n+1} (u^{n+1} - u^n)/dt + eps^{n+1} u.grad u. CFDEM's
backward-Euler assembly gives eps^n (u^{n+1} - u^n)/dt; the difference is O(dt d eps/dt du/dt) and
was kept deliberately (one coefficient field).

**Projection.** The start-up pressure matrix C0 = B^T (rho Ml)^-1 B is not rebuilt. It is the
preconditioner of the conservative residual iteration: the residual r_m of an element enters the
pressure right-hand side as -r_m/(w_e dt), with w_e = 1 for B and Bfull (the eps of the flux and
the 1/eps of the weighted mass cancel) and w_e = eps_e for `none` and A. The cancellation is exact
for uniform porosity only; the fixed point r = 0 does not depend on w_e and the contraction is
observed per run in `EL_VA_FLUX` (res1_L1 per corrector), not assumed.

**Particle viscous force (Bfull, A).** <lap u_c> = (A_c - B.G_c)/N with A_c = -int grad W.grad u_c
+ int_wall W grad u_c.n, B = -int grad W + int_wall W n, G_c the kernel-averaged gradient, N the
kernel integral. The B.G term removes the quadrature defect, so a uniform velocity gradient gives
exactly zero. The boundary integrals run over all faces on the physical boundary of the domain
(a face with one adjacent element over all ranks, whatever its boundary condition); the kernel
must vanish at the edge of its support, hence deen_poly only.
Additive across ranks (summed at the owner). Enters the semi-implicit particle update, the
Newton-pair expectation and, under Bfull, the deposited feedback.

**Model A, operator form.** -eps grad p and div(eps c grad u) are the unweighted operators minus
their solid-fraction parts fp(p) and c fv(u): `EL_MM_Add_Operator_Terms` adds dt (c fv(u^n) -
fp(p)) to the momentum right-hand side after `AddPressureGradient`, `EL_MM_Correct_Velocity_Operator`
subtracts dt M^-1 fp(q) after `Velocity_Correction`. fp is the transpose of the solid-fraction part
of the conservative divergence (to rounding on straight cells), so pressure gradient and continuity
form one saddle point, and a constant pressure exerts no force for any continuous porosity. The
solid-fraction share of the viscous term is explicit; with an implicit unweighted term this is
unconditionally stable.

## 4. Verification

Unit tests `test_el_momentum_model` (policy table, key resolution, force split, viscous average,
start-up balance, model A operators): pass. Viscous average of a quadratic field: 0.03-0.07 % in
the interior and at a wall, brick and warped mesh; without the wall term 14 %. Transpose identity
of fp: 5e-14 (brick), 7e-9 (warped).

Gates at the default (ten certified configurations against `ergun_gates2`, strict comparison):
stage 1 (key only) jobs 160443-160452, combined code jobs 160466-160475: bitwise. DNS twins e4_l3
and DKT: stage 1 jobs 160455/160456 bitwise; combined code: see the case RUNBOOK. Final round on
the commit: RUNBOOK of `s06_startup_inertia`.

Validation: `applications/q2p1_el_pipeflow/validation_cases/s06_startup_inertia/RUNBOOK.md`
(inertia weighting, time constant of a suspension starting from rest) and section "Momentum
models" of the s05 RUNBOOK (Ergun column: the three models against each other and CFDEM).

## 5. Design consultation (GPT-6 Astra, high, on the plan) and what was done with it

| Finding | Decision |
|---|---|
| new-time eps in the mass term is not CFDEM's discrete form | kept, documented (section 3) |
| unweighted C is a preconditioner, not the exact operator | adopted in the wording; convergence observed per run |
| divisor must change in the source bookkeeping and the right-hand side together | done (w_e in both and in the sink identity) |
| viscous force: operator, units, wall traction, additive assembly | weak kernel form with wall term; coefficient of the fluid operator (section 6) |
| force hooks must reach the pre-advance reconstruction and the Newton pair | done (`particle_nondrag`, `feedback_nondrag`, expected impulse) |
| model A cannot be a scaling of the pressure dofs | done as transpose of the conservative divergence (fp) |
| mass lifecycle, renewal digits, GAMMA, convective form | checked at the first step |
| volume-weighted restriction | done |
| fluid-pair audit on weighted momentum | done for the fluid substep; the implicit drag share is outside the audit as before |
| ELFluidGravity with a hydrostatic pressure double-counts buoyancy | ELFluidGravity = No implied |
| mass projection with stiff implicit drag is a splitting | unchanged from the campaign scheme; open (section 6) |
| restart equivalence, output weights | open (section 6) |

## 5b. Code reviews of 0fb06e5a (GPT-6 Astra high; an independent Claude agent) and the fix round

Both: no change of a floating-point result or an MPI sequence on the default path; signs and
scalings of the model A terms, the residual divisor, the force split, the drag scale and the
parallel restriction correct. Findings and what was done:

| Finding (who) | Done |
|---|---|
| truncated gaussian kernel: viscous average 8 % low (both) | gaussian refused with Bfull / A |
| boundary term only on all-Dirichlet faces: wrong at outlets and slip planes (both) | geometric boundary-face mask (one adjacent element over all ranks) |
| model B without its 1/eps on the first particle advance, porosity still 1 (Astra) | porosity priming pass under a model |
| Newtonian check lets a viscoelastic deck through (Astra) | refused |
| `.AND.` need not short-circuit: FBM path may evaluate the predicate (Astra) | nested guard |
| "advective" audit term changes meaning with divisor 1 (agent) | formed with the porosity for every model; default expression untouched |
| model ids as literals in EL_CONFIG, table in two modules (both) | one model table in EL_CONFIG, EL_MOMENTUM_MODEL reads it; id cached at validation |
| duplicated 3 x 3 inverse; misplaced header comment (agent) | `EL_VA_INVERT3`; comment restored. The face-normal code of `EL_VA_ELEMENT_FACE_FLUX` stays inline on purpose (certified arithmetic) |
| tests never select A, loose tolerances (both) | A in policy / keys / force split; tolerances at the measured accuracy |
| setting check late and printed by every rank; refusals one at a time; mask lifecycle; readiness of the transport porosity (agent) | check at the first step start on one rank; conflicts collected; mask released and rebuilt with the fields; readiness stop |
| not unit-tested: adapter (restriction, mass rebuild, model A corrections as applied), multi-rank viscous record (both) | stated in the test header; covered by the s06 and s05 runs and the rank-count gates only |
| per-step mass rebuild assembles the consistent mass and prints a line (agent) | left; cost in the s06 RUNBOOK |

## 6. Open points

1. **Viscosity convention of the fluid operator.** The diffusion matrix carries
   `Properties%Viscosity(1)` without a density factor against a mass matrix that carries the
   density. With `Prop@DynVisc` the parser stores mu/rho there (the drag closures form mu =
   rho nu themselves), so for rho /= 1 the fluid's viscous term is too small by the factor rho in
   every Euler-Lagrange run; the DNS decks give the dynamic viscosity in `Prop@Viscosity` and are
   consistent. Pre-existing, owner decision. `EL_FLUID_OPERATOR_VISCOSITY` is the one place the
   particle viscous force and model A follow whatever is decided.
2. Stress form: the Newtonian path uses the component Laplacian; CFDEM's deviatoric form differs
   where div u /= 0. The particle force matches this solver's operator.
3. Momentum-aware pressure mobility (drag in the projection), outer momentum/pressure iterations.
4. Restart with a model: the coefficient field is rebuilt from the porosity; not certified.
5. P1 -> Q2 pressure output and the profile post-processing use the lumped mass as nodal weight,
   which is eps-weighted under a model.
6. Periodic domains with a model: not run (the conservative form is certified for closed boxes
   and, since case 05, open outlets).
7. Cost of the viscous sampling: one more loop over the supported elements per particle.
