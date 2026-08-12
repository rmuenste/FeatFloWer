# kroupa_shear — suspension viscosity from the lubrication virial

Measures the particle-phase suspension viscosity eta_L(phi) of the
pairwise-lubrication closure (Kroupa/Vonka/Soos/Kosek, Langmuir 32
(2016) 8451; `kroupa_soos_2016.pdf` at repo root) in the walls-free
frozen-linear-shear unit box (fourway_shear geometry, G = 1), via the
Irving-Kirkwood impulse virial accumulated in PE and printed as
EL_SUSP_STRESS (sigma = -virial/(dt*V); eta_L = sig_xz_tavg/G).

Equivalence to the paper's wall force balance: in statistically steady
homogeneous shear the momentum flux through any plane equals the bulk
stress, so the wall sum (their eqs 27-32) and the volume virial measure
the same eta_L; the virial needs no walls. The hydrodynamic component
eta_H is NOT resolved in a frozen field — the total-viscosity test is
the V4 pipe rerun (mu_app vs Krieger-Dougherty).

## Configuration
- Neutrally buoyant d_p = 0.05 spheres (particleDensity_ = 1.0), unit
  triply periodic box, frozen linear_shear G = 1, dt = 0.002,
  substeps_ = 10 (lubrication damping time tau = m/c ~ dt/4 at these
  parameters — the force is re-evaluated per substep inside PE).
- Lubrication: cutoff 0.025 = R_p... note cutoff is d_p/2 (surface gap
  trigger), slip length h_c = 0.0025 (eps_c = h_c/R_p = 0.1, the value
  Kroupa fit to Krieger/de Kruif data), minEpsLub_/alphaImpulseCap_
  defaults. Dynamic viscosity mu = rho_f*nu = 1*0.02 = 0.02.
- 5000 steps to t = 10 (t*G*phi >= 1-2 = steady per Kroupa Fig 2);
  eta_L from the EL_SUSP_STRESS tail average (ELTAvgWindow default 0.5).
- Runs: phi in {0.05, 0.10, 0.20, 0.30} lubrication-ON + phi = 0.20
  OFF-twin (lubricationEnabled_ = false).

## Gates
- OFF-twin: eta_L = 0 identically (virial only accumulates lubrication).
- ON: eta_L > 0, monotone increasing in phi; phi-trend consistent with
  Kroupa Fig 2b / Fig 3 shape (rapid growth with phi). Quantitative
  Krieger-Dougherty comparison happens at total-viscosity level in V4.
- EL_NEWTON_PAIR machine zero (lubrication is internal, antisymmetric).
- max_overlap <= 1% d_p; contact count stationary; no NaN.

## Run log
(appended as runs complete)

## Smoke ladder (2026-07-26) — two real bugs caught by the Newton-pair audit

1. Job 135995 (cutoff-only shadow margin): etaL sane (5.15e-3) and
   OFF-twin (135996) exactly zero, but EL_NEWTON_PAIR worst 4.3e-6 —
   one-sided pairs. Root cause: margin must be radius + cutoff (a
   partner can sit up to R+cutoff beyond the boundary when the owned
   body's center is at the domain edge). Fix pe 11895f3.
2. Job 136001 (margin fixed): 3.2e-7 residual — Gauss-Seidel fold: pair
   forces applied into v_ mid-sweep, so later pairs on the owner rank
   read updates the partner's owner rank never saw. Fix: Jacobi
   accumulation against the frozen synced state (pe commit after
   11895f3).
3. Job 136017 (both fixes): EL_NEWTON_PAIR worst 4.0e-20 at phi=0.20
   under 83k pair-substep evaluations; etaL_tavg 5.11e-3; max_overlap
   0.0; no NaN. ALL SMOKE GATES PASS.

Full sweep submitted: jobs 136021-136024 (phi 0.05/0.10/0.20/0.30 ON)
+ 136025 (phi=0.20 OFF twin), 5000 steps to t=10, med partition.

## Sweep continuation — two more one-sidedness bugs, then the final scheme

4. Full-length phi>=0.20 runs (136023/136024) reintroduced intermittent
   EL_NEWTON_PAIR bursts (13/200 audits at phi=0.20 up to 4.3e-6; 55/200
   at 0.30) absent in the 500-step smoke and dilute runs. An ownership
   fix (bodystorage membership -> !isRemote(), pe commit) was falsified
   by a byte-identical rerun (136096/136097).
5. EL_LUB_IMPULSE instrumentation (job 136098) matched the bursts
   digit-for-digit to the net momentum folded by the sweep: pair
   one-sidedness confirmed — the fold-once scheme requires perfectly
   symmetric cross-rank pair visibility, which nothing guarantees.
6. Final scheme (pe designated-treater commit): each pair evaluated by
   exactly ONE rank (owner of the smaller system ID), corrections for
   BOTH members written to dv_/dw_ pre-divided by relaxationParam_ and
   delivered by the synchronizeVelocities relay — Newton-exact by
   construction, like the PGS contact path; visibility asymmetry
   degrades to a skipped substep, not a momentum leak.

## RESULT (2026-07-27, jobs 136099/136100/136101/136102 + OFF-twin 136025)

| phi | eta_L | eta_L/mu | worst NP | violated audits |
|----:|------:|--------:|---------:|----:|
| 0.05 | 1.581e-4 | 0.0079 | 1.3e-20 | 0/200 |
| 0.10 | 7.252e-4 | 0.0363 | 6.1e-20 | 0/200 |
| 0.20 | 4.991e-3 | 0.2496 | 1.7e-19 | 0/200 |
| 0.30 | 1.522e-2 | 0.7610 | 2.1e-19 | 0/200 |
| 0.20 OFF | 0.0 (exact) | 0 | 0.0 (exact) | 0/200 |

- Monotone, strongly superlinear (effective eta_L ~ phi^2.2-2.8 between
  consecutive points) — the rapid-growth shape of Kroupa Fig 2b/Fig 3.
- OFF-twin identically zero: the virial isolates the lubrication stress.
- max_overlap <= 7e-13 at all phi; contact counts stationary; no NaN.
- ALL GATES PASS. Total-viscosity (Krieger-Dougherty) comparison happens
  in the V4 lubricated pipe rerun (jobs 136103-136106).

## Contact-channel stress + substep ladder (2026-07-27, jobs 136206/136207/136208)

Step 1 of the publication-grade Kroupa reproduction (user decision
2026-07-27): before building the walled Couette twin, (a) accumulate the
converged PGS contact impulses into their own virial (EL_CONT_STRESS,
pe 11b518a + FF b9e6499f) and (b) check temporal convergence of the
stiff lubrication forces (Kroupa used ~250 particle substeps per fluid
step; we ran 10). Ladder at phi=0.30, substeps_ = 10/50/100, 5000 steps
to t=10, new binary:

| substeps | etaL/mu | etaC/mu | etaP/mu | 1+etaP/mu | worst NP | max_overlap |
|---------:|--------:|--------:|--------:|----------:|---------:|------------:|
| 10  | 0.761 | 0.831 | 1.592 | 2.592 | 2.1e-19 | 4.4e-10 |
| 50  | 0.758 | 0.847 | 1.606 | 2.606 | 3.5e-19 | 1.5e-11 |
| 100 | 0.753 | 0.862 | 1.614 | 2.614 | 5.2e-19 | 3.5e-11 |

Findings:

1. SUBSTEPS CONVERGED: total particle-phase viscosity drifts +0.8%
   from substeps 10 -> 100 (etaL -1.1%, etaC +3.8%, partially
   cancelling). substeps_ = 10 is adequate for these parameters;
   Kroupa's 250 is not needed. No submission-cost blocker for the
   Couette build.
2. THE CONTACT CHANNEL IS THE LARGER PARTICLE-PHASE CHANNEL at
   phi=0.30: etaC/mu = 0.83-0.86 vs etaL/mu = 0.75-0.76. The
   hard-contact impulses fire at essentially zero overlap
   (max_overlap ~ 1e-10, velocity-level constraint), which is why the
   max_overlap gate never saw them: contact stress without penetration.
   The earlier "etaL" sweep numbers are the lubrication channel only,
   NOT the particle-phase total.
3. GAP TO KRIEGER-DOUGHERTY NEARLY CLOSED: 1 + etaP/mu = 2.59-2.61 vs
   KD(0.30) = 2.75 -- within ~5-6% with the hydrodynamic channel eta_H
   still structurally invisible to the frozen field. (For the pipe V4
   this also reframes the shortfall: mu_app there already includes
   contact effects implicitly via the flow rate; the missing pieces
   remain eta_H-type feedback, wall lubrication, and migration
   depletion.)
4. Gates: NP machine zero (<= 5.2e-19), no NaN, OFF-twin logic
   untouched.

Contact-channel phi sweep (kc05/kc10/kc20, jobs 136222-136224,
substeps 10, new binary) completes the etaC(phi) curve; results below
when finished.

## Contact-channel phi sweep (2026-07-27, jobs 136222/136223/136224)

Reruns of phi=0.05/0.10/0.20 with the EL_CONT_STRESS binary (substeps
10); etaL reproduces the original sweep to all printed digits (bitwise
repro of the lubrication channel). Full two-channel picture:

| phi | etaL/mu | etaC/mu | etaP/mu | 1+etaP/mu | KD | KD-(1+etaP/mu) | 2.5*phi |
|----:|--------:|--------:|--------:|----------:|-----:|------:|------:|
| 0.05 | 0.0079 | 0.0017 | 0.0096 | 1.010 | 1.139 | 0.129 | 0.125 |
| 0.10 | 0.0363 | 0.0125 | 0.0488 | 1.049 | 1.312 | 0.264 | 0.250 |
| 0.20 | 0.2496 | 0.2106 | 0.4602 | 1.460 | 1.821 | 0.361 | 0.500 |
| 0.30 | 0.7610 | 0.8306 | 1.5916 | 2.592 | 2.751 | 0.160 | 0.750 |

- Contact fraction of the particle-phase stress grows 18% -> 26% -> 46%
  -> 52% with phi: the dry-contact channel is co-dominant from
  phi ~ 0.2 up.
- CONSISTENCY CHECK, dilute limit: the component MISSING from the
  frozen box (KD minus our total) matches the Einstein hydrodynamic
  term 2.5*phi almost exactly at phi = 0.05 (0.129 vs 0.125) and 0.10
  (0.264 vs 0.250). Exactly what the construction predicts: the frozen
  field cannot produce the hydrodynamic (eta_H) contribution, which IS
  the entire suspension viscosity in the dilute limit. This is strong
  evidence the pair channels are measured correctly and the walled
  Couette twin (which adds eta_H by solving the fluid) is the right
  remaining step.
- Gates: no NaN, NP <= 1.7e-19 at all phi.
