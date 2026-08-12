# Runbook: v3_ss_coupled — two-way Segré–Silberberg pipe (W4.2)

Goal: the coupled lift arbiter. Same physics as the frozen convergence
demonstration (v3_ss_frozen RUNBOOK, Stage 4a) but with the full two-way
solver: body-force-driven pipe flow, kernel feedback, drag + Matas–Asmolov
inertial lift. Lift-ON must focus the suspension onto the annulus; the
lift-OFF twin must not.

## Configuration

- Mesh: PIPEZ27 (symlink to pipe_hp_check), R = 0.5, L = 2, z-periodic
  (CFD `PeriodicAxis=z` + PE `periodicZ_`), np = 28, level 3.
- Flow: `ConstantForcing` z = 0.192 → laminar u_max = f·R²/(4μ) = 0.6,
  Rc = u_max·D/ν = 30 (matches the digitized lift table; HP validation of
  this exact setup: pipe_hp_check). Startup transient ≈ 7·τ ≈ 15 t.u.
- Particles: N ≈ 30 (volumeFraction_ 0.01, random cylinder seeding with
  radial inset), d_p = 0.1 (a/R = 0.1, as in the frozen demonstration),
  neutrally buoyant, gravity off, dry contacts.
- Lift: `ELInertialLift = matas_asmolov` with `ELInertialLiftUmax = 0.6`
  (coupled runs REQUIRE the explicit key — no prescribed-field fallback).
  `ELLiftModel = none` (slip lift ≈ 0 at zero steady slip anyway; keep the
  arbiter clean).
- Time: dt = 0.02 (CFL ≈ 0.46 at u_max, dt/τ_p = 0.72 — PE drag update is
  semi-implicit), 50 000 steps → t = 1000. Frozen-field reference reached
  the annulus band by t ≈ 1500 from a wider fan; random mid-radius seeds
  plus two-way noise should show clear focusing by t = 1000 (fallback per
  plan: direction/rate demonstration if full focusing needs more time).
- Metric: `EL_SS_RADIUS` every 100 steps (rmean/rmin/rmax of r/R).
  ON acceptance: band contracts toward ≈ 0.675 (frozen-run reference;
  ±0.05); OFF acceptance: no systematic contraction (spread statistics
  fluctuate about the seed distribution).
- Momentum audit runs (coupled, ρ=1): EL_MOMENTUM_ELEMINT drift bounded as
  in tier2; EL_VOLUME_CONSERVATION at machine zero through wraps.

## Run (Slurm)

```bash
tools/el_stage_rundir.sh v3_ss_coupled <rundir_on>    # staged default = lift ON
tools/el_stage_rundir.sh v3_ss_coupled <rundir_off>
sed -i 's/^SimPar@ELInertialLift = .*/SimPar@ELInertialLift = none/' <rundir_off>/_data/q2p1_param.dat
tools/el_slurm_submit.sh <rundir_on>  w42_lift_on  1-23:00:00 long
tools/el_slurm_submit.sh <rundir_off> w42_lift_off 1-23:00:00 long
```

## Results (2026-07-06, Slurm jobs 131967/131968, nx nodes, 50k steps)

| t | lift ON rmean [rmin, rmax] | lift OFF rmean [rmin, rmax] |
|---|---|---|
| 2 | 0.532 [0.209, 0.869] | 0.532 [0.208, 0.871] |
| 400 | 0.575 [0.262, 0.852] | — |
| 1000 | **0.620 [0.338, 0.839]** | **0.536 [0.092, 0.871]** |

PASS (direction/rate criterion, the plan's fallback): lift ON focuses the
30-particle suspension steadily toward the annulus (Δrmean = +0.088 over
t = 1000, inner tail migrating out at the frozen-run rate for a/R = 0.1);
lift OFF shows no systematic contraction (Δrmean = +0.004, spread
random-walks). Full annulus convergence would need ≈ 2× the steps
(frozen-field reference reached the band at t ≈ 1500–2000 from a similar
inner tail); a 100k-step continuation is optional follow-up, the ON/OFF
discrimination is already unambiguous (22:1 in Δrmean).

## Results — report-grade 100k-step run (2026-07-09, Slurm job 132699 `w42fix_on`, nx-01, 21h23)

Clean rerun with the fixed binaries (pe 8b037ae Newton-pair fix +
src_el 2c6ee44d audit ILEV pin), lift ON, 100 000 steps → t = 2000,
full VTK (29 particle frames in `paraview/`):

| t | rmean | [rmin, rmax] |
|---|---|---|
| 2 | 0.532 | [0.209, 0.869] |
| 100 | 0.542 | [0.231, 0.841] |
| 500 | 0.582 | [0.262, 0.853] |
| 1000 | 0.622 | [0.348, 0.820] |
| 1500 | 0.651 | [0.379, 0.798] |
| 2000 | **0.664** | [0.454, 0.778] |

PASS (full criterion, no fallback needed): rmean converges to 0.6643 vs
the Matas–Asmolov equilibrium 0.675 (−1.6%, well inside ±0.05); the
radial band contracts monotonically from width 0.66 to 0.32 and is still
tightening at t = 2000 (the slowest particles are the inner tail, whose
lift force ∝ ĝ(s) is smallest near the axis — expected). This is the W4.2
lift-arbiter result to quote: the coupled two-way solver with the
Matas–Asmolov closure reproduces Segré–Silberberg annulus focusing at
Rc = 30 with a/R = 0.1, N = 30.
