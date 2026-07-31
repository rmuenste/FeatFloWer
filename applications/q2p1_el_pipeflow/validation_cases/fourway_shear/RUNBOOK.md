# Runbook: fourway_shear — four-way sheared-box stress test (Stage 2)

Goal: stress the contact system at the campaign's densest loading. φ = 0.20
(N = 3056) in the fully periodic unit box under a frozen `linear_shear`
prescribed field (rate 1.0), gravity off, contacts active. The frozen field
skips the CFD solve, so every step is particle/contact work — the heaviest
sustained contact load anywhere in the campaign, including the Newton-pair
audit exercising cross-rank and periodic-image contacts continuously.

## Configuration

- Mesh: QBOX9 unit box (symlink to tier2 momentum_conservation), level 3,
  np = 28. CFD `Periodic = Yes`, PE `periodicX_/Y_/Z_: true`
  (decomposePeriodic3D).
- Particles: d_p = 0.05, N = 3056 (`seedMode_ = random`,
  `volumeFraction_ = 0.20`, seed_ 12345), dry contacts, restitution 0.
- Field: `ELPrescribedField = linear_shear`, `ELShearRate = 1.0`,
  `ELFluidGravity = No`, gravity off.
- Time: dt = 0.002, 10 000 steps → t = 20. Particle VTK: `vtk_ true`,
  `visspacing_ 200` → paraview/collector.pvd.
- Diagnostics: `EL_CONTACT_STATS` + `EL_NEWTON_PAIR` every 25 steps.

## Gates (campaign plan Stage 2)

1. max_overlap ≤ 1% d_p at every sampled step
2. Tgran bounded — last window ≤ early window × (1+tol), no spontaneous
   heating
3. contact count stationary
4. no NaN
5. (added) EL_NEWTON_PAIR at machine zero under heavy contact load

## Results (2026-07-08, Slurm job 132694, nx-06, np = 28, 10k steps)

ALL GATES PASS:

- max_overlap = 0.0 at all 400 samples (PGS resolves every contact within
  the step; gate was ≤ 1% d_p).
- Tgran flat at 2.6727e-2 → 2.6729e-2 (constant to 4 significant digits
  over t = 0.05…20 — no spontaneous heating; the value is the equilibrium
  shear-induced fluctuation level).
- ncontacts stationary: mean 5315, range [3883, 5820], no trend
  (5516 at the first sample, 5473 at the last).
- No NaN anywhere in the protocol.
- EL_NEWTON_PAIR: exactly 0.0 at all 400 audits — the pe 8b037ae
  fold-once/exact-zero/migration-handover fix holds under ~5.3k
  simultaneous contacts with periodic images on all three axes.

Wall time 3h25 for 10k steps (≈ 1.2 s/step at np = 28 on nx-06, frozen
field). Particle VTK (29 frames) in the run dir's `paraview/`.

Together with unit_periodic_contact (contact THROUGH a periodic plane
fires and resolves; see that case's RUNBOOK), this closes Stage 2 of the
campaign plan and removes the periodic-image asterisk from the Stage-1
φ = 0.20 results.
