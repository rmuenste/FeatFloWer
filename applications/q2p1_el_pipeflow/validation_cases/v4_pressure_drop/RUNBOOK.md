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
