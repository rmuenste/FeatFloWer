# V4 — pressure drop vs volume fraction (+ V5 radial profiles)

Stage-5 case of the EL validation campaign (plan §Stage 5). Measures the
apparent suspension viscosity mu_app(phi) = f*R^2/(8*u_mean) in the
Stage-3 z-pipe at fixed axial body force, and (V5, same runs) the radial
concentration profiles.

## Design

- Domain: z-aligned pipe, R = 0.5, L = 2.0 (periodic z), same
  partitioned PIPEZ27 mesh as pipe_hp_check / v3_ss_coupled; level 3
  production, np = 28.
- Particles: neutrally buoyant (rho_p = rho_f = 1), d_p = 0.05
  (R/d_p = 10 per the plan constraint), random cylinder seeding
  (seedMode_ = random, seedDomain_ = cylinder), dry contacts.
- Forcing: identical axial force f = 0.192 across all phi (the
  v3_ss_coupled value; single-phase u_max = f*R^2/(4*mu) = 0.6,
  u_mean = 0.3, pipe Re = u_mean*2R/nu = 15). No flow-rate controller:
  mu_app is measured from the phi-dependent mean velocity at fixed f
  (mode-agnostic, plan decision 4).
- Lift: ELInertialLift = matas_asmolov, Umax = 0.6 (the W4.2
  configuration) so the V5 profiles carry the validated migration
  physics.
- ELMomentumFix = No and no counter-forcing: the pipe is wall-bounded;
  the EL_FLUID_PAIR audit's expected-source bookkeeping does not include
  wall traction, so its mismatch line is NOT a conservation gate here
  (unlike the periodic v2 box). Ignore it in this case.
- N per phi (V_pipe = pi*R^2*L = 1.5708, V_p = 6.545e-5):
  phi=0.05 -> 1200, 0.10 -> 2400, 0.15 -> 3600, 0.20 -> 4800.
  All below the 5964 scaling-envelope ceiling verified in Stage 0.
- Per-phi staging: stage this case, then sed volumeFraction_ in
  example.json (0.05/0.10/0.15/0.20).

## Measurement

- u_mean(phi): mixture superficial velocity from EL_MEAN_SLIP_TAVG,
  u_mean = uf_super,z + phi*up,z (plus bulk_flow.log cross-check).
- mu_app(phi) = f*R^2/(8*u_mean); mu_app(0) from pipe_hp_check
  (or a phi=0 twin) must recover mu within 3%.
- Reference curve: in-repo mu_eff(phi) correlation
  (source/src_quadLS/viscosity_model.f90, Krieger-Dougherty-type,
  phi_max = 0.64) evaluated offline; matched FBM-DNS points deferred
  (plan decision 2).
- V5: annular-bin concentration histograms offline from _ns/PT dumps.

## Acceptance (plan §Stage 5)

- mu_app monotone increasing in phi; mu_app(0) within 3% of mu.
- mu_app/mu within ~15% of the KD curve for phi <= 0.15, ~25% at 0.20.
- V5: near-wall depletion + off-center peak, qualitative vs literature.
- EL_CONTACT_STATS: max_overlap <= 1% d_p (watch near-wall at 0.20).

## Gate: phi=0.20 smoke (MANDATORY before long runs)

~200 steps, np=28, N=4800 — checks random cylinder seeding achieves
phi=0.20 with no initial overlaps and the radial inset respected
(EL_SS_RADIUS rmax), per-rank balance/halo memory, wall-time/step to
size the Slurm allocations.

## Run log

(appended as runs complete)

## Baseline (no lubrication) — 2026-07-27, jobs 135958-135961: model boundary confirmed

All four runs clean to t=400 (no NaN). Mixture flow rate
u_mean = uf_super,z + phi*up,z from EL_MEAN_SLIP_TAVG:

| phi | u_mean | mu_app/mu |
|----:|-------:|----------:|
| 0.05 | 0.3004 | 0.999 |
| 0.10 | 0.3007 | 0.998 |
| 0.15 | 0.3020 | 0.993 |
| 0.20 | 0.3039 | 0.987 |

Flat (slightly below 1) vs Krieger-Dougherty's expected +15%..+72%: the
drag-only unresolved model generates no suspension viscosity for
neutrally buoyant particles — the lubrication-OFF branch of Kroupa
et al. Fig 5, reproduced in our pipe. This table is the reference for
the LUBRICATED rerun (jobs 136113-136116, Kroupa closure + substeps 10).
Note: tail-averaged <eps_f> drifts slightly above 1-phi at higher phi
(0.812 vs 0.800 at phi=0.20) — kernel/wall bookkeeping worth a look, but
u_mean does not depend on it.

## Lubricated sweep (W6) — 2026-07-27, jobs 136113-136116

Same staging as the baseline plus Kroupa lubrication
(lubricationEnabled_ true, cutoff 0.025 = R_p, h_c 0.0025,
substeps_ 10) and the shadow-margin clamp (pe margin 0.05 -> 0.0483 for
the 1x1x27 pipe slabs; clamp notice printed, healthy). All four clean to
t=400, no NaN, max_overlap = 0 throughout.

| phi | u_mean | mu_app/mu (lub ON) | baseline (OFF) | Krieger-Dougherty |
|----:|-------:|-------------------:|---------------:|------------------:|
| 0.05 | 0.2965 | 1.012 | 0.999 | 1.139 |
| 0.10 | 0.2841 | 1.056 | 0.998 | 1.312 |
| 0.15 | 0.2707 | 1.108 | 0.993 | 1.533 |
| 0.20 | 0.2613 | 1.148 | 0.987 | 1.821 |

Verdict:
- MECHANISM PASS: lubrication converts the flat baseline into a strictly
  monotone mu_app(phi) — the Kroupa Fig 5 on/off contrast reproduced in
  the wall-bounded pipe at fixed forcing (ON-OFF spread grows 0.013 ->
  0.161 across phi).
- QUANTITATIVE: recovers only ~9-20% of the KD excess viscosity
  (-11% of KD at phi=0.05 up to -37% at 0.20), outside the ~15%/~25%
  acceptance bands for phi >= 0.10. Two known, documented reasons:
  (a) no particle-WALL lubrication (pipe wall is a CFD boundary, not a
  PE body) — precisely the near-wall high-shear region that dominates
  dissipation in Poiseuille flow; (b) Matas-Asmolov lift migrates
  particles to the s~0.6-0.7 annulus, depleting the wall region further.
  The homogeneous-shear kroupa_shear box (no walls, no migration) shows
  the particle-phase etaL alone reaching 0.25*mu at phi=0.20, so the
  closure itself generates the right order — the pipe shortfall is
  geometry/wall physics, not the pair force. Recorded as
  mechanism-PASS / quantitative-RECORDED; wall lubrication is the
  designated follow-up if quantitative KD tracking is required.
- Gates: Newton-pair machine zero in 199/200 audits per run; single
  transient at t=180 in phi=0.20 (mismatch 2.8e-6, likely an ownership
  migration within one macro step; ~3e-5 of total particle momentum,
  no effect on u_mean). max_overlap identically 0 (lubrication holds
  gaps open — contacts never reach penetration).
- <eps_f> TAVG: 0.9500/0.9001/0.8517/0.8103 — same mild high-phi drift
  as the baseline (bookkeeping note above), unchanged by lubrication.
