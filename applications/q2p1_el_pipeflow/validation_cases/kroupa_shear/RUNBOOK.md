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
