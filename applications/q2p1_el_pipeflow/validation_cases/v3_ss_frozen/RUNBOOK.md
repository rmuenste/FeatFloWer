# Runbook: v3_ss_frozen — Segré–Silberberg frozen-Poiseuille plumbing test (W4.1)

Goal: exercise the full lift path (Saffman–Mei/Zeng blend, wall arm, cylinder
wall distance) in a pipe geometry with a prescribed (frozen) Poiseuille field,
BEFORE committing compute to the coupled Segré–Silberberg run (W4.2). This is
a plumbing/sign test: it must show (a) the `EL_SS_RADIUS` diagnostic working,
(b) radial particle motion with lift ON, (c) radially frozen trajectories with
lift OFF.

## Configuration

- Mesh: PIPEZ27 reused from `pipe_hp_check` (the `mesh` entry here is a
  symlink to `../pipe_hp_check/mesh`; staging resolves it). R = 0.5, L = 2,
  axis z, np = 28.
- Frozen field: `ELPrescribedField = poiseuille`, `ELPrescribedUmax = 1.2`
  → nominal pipe Re = u_mean·D/ν = 0.6·1/0.02 = 30 (ν enters the closures
  only; the CFD solve is skipped in prescribed mode).
- Particles: 16 neutrally buoyant spheres (ρ_p = ρ_f = 1), d_p = 0.05
  (a/R = 0.05), seeded at z = 0.2 on a radial fan r/R = 0.1 … 0.9 (uniform
  spacing 0.8/15) with golden-angle azimuthal stagger; min pairwise distance
  0.107 ≫ d_p. Gravity off. Feedback off. Contacts irrelevant at these
  separations.
- Time: dt = 0.002 (PE `stepsize_` matched), 650 steps → t = 1.3. Particle
  response time τ_p = ρ_p d²/(18μ) = 6.9e-3 → dt/τ_p ≈ 0.29 (explicit drag
  stable). The committed baseline (650 steps) predates particle wrap; runs
  are no longer track-capped when `"periodicZ_": true` is set in
  `example.json` (PE `decomposePeriodicZ3D`, requires `processesZ_ ≥ 3`).
  Wrap verification (lift OFF, 2000 steps → fastest particle wraps z = 0/2
  twice): count = 16 throughout, radii frozen to all printed digits, worst
  EL_VOLUME_CONSERVATION rel_error 2.5e-15 across the periodic plane.
- Diagnostic: `EL_SS_RADIUS t= step= count= rmean= rmin= rmax=` every 25
  steps (values are r/R, reduced over all particles).

## Physics expectation (read before judging the numbers)

The implemented lift closures are slip-based (Saffman/Mei factor, Zeng wall
and shear arms — all scale with |u_p − u_f|). A neutrally buoyant particle in
a frozen field relaxes to the local fluid velocity on τ_p ≈ 7e-3, after which
the steady slip — and hence the closure lift — is near zero. The expected
signal is therefore:

1. An initial radial kick during the start-up transient (particles start at
   rest → axial slip ≈ u(r) → lagging-particle lift points toward the
   centerline), decaying on a few τ_p.
2. Little further migration. TRUE Segré–Silberberg equilibrium (r/R ≈ 0.6)
   is NOT expected from slip-based closures at zero steady slip; capturing it
   in the coupled run may require a dedicated inertial-lift closure
   (Ho & Leal / Asmolov type). This is a campaign-level finding to record,
   not a defect of the implementation.

The pass criteria below are therefore about plumbing, not about reaching the
Segré–Silberberg annulus.

## Run

```bash
BUILD=/home/user/rmuenste/nobackup/code/FF-EL/FeatFloWer/build-el-phase2-pe-gcc14
BIN=$BUILD/applications/q2p1_el_pipeflow/q2p1_el_pipeflow
MPIPREFIX=$(ldd "$BIN" | grep -m1 'libmpi\.so' | awk '{print $3}' | sed 's#/lib/libmpi.*##')
export PATH="$MPIPREFIX/bin:$PATH"; export LD_LIBRARY_PATH="$MPIPREFIX/lib:$LD_LIBRARY_PATH"
export OMPI_MCA_rmaps_base_oversubscribe=1

cmake --build "$BUILD" --target q2p1_el_pipeflow_val_v3_ss_frozen_stage
RUN=$BUILD/applications/q2p1_el_pipeflow

# Lift ON (staged default: ELLiftModel = saffman_mei_wall)
mpirun --oversubscribe --wdir "$RUN" -np 28 "$BIN" > "$RUN/run_ss_on.log" 2>&1
grep EL_SS_RADIUS "$RUN/simulation_output_level_3.log" > ss_on.dat

# Lift OFF
sed -i 's/^SimPar@ELLiftModel = .*/SimPar@ELLiftModel = none/' "$RUN/_data/q2p1_param.dat"
mpirun --oversubscribe --wdir "$RUN" -np 28 "$BIN" > "$RUN/run_ss_off.log" 2>&1
grep EL_SS_RADIUS "$RUN/simulation_output_level_3.log" > ss_off.dat
```

## Variant: slip-driven lift (sustained-migration proof)

The neutral pair only exercises lift during the start-up transient. To prove
the lift path can drive SUSTAINED radial migration, a third run gives the
particles a persistent slip: `particleDensity_ = 2.0` (json) and
`Prop@Gravity = 0d0,0d0,-0.981d0` (param; PE json `gravity_` stays zero —
gravity enters only through the EL grav_buoy force, tier2-terminal
convention). Settling slip v_t = Δρ·g·d²/(18μ) ≈ 6.8e-3 (slip Re_p ≈ 0.017,
Saffman–Mei arm of the blend), τ_p = 1.39e-2 (dt/τ_p = 0.14). The particle
lags the fluid axially, so the shear lift points toward the CENTERLINE:
expect a steady rmean decrease, of order Δ(r/R) ≈ 0.03 over t = 1.3
(Saffman estimate at r/R = 0.5), versus a pure settle-through with radii
frozen when `ELLiftModel = none`.

```bash
sed -i 's/^SimPar@ELLiftModel = .*/SimPar@ELLiftModel = saffman_mei_wall/' "$RUN/_data/q2p1_param.dat"
sed -i 's/"particleDensity_": 1.0/"particleDensity_": 2.0/' "$RUN/example.json"
sed -i 's/^Prop@Gravity = .*/Prop@Gravity = 0d0,0d0,-0.981d0/' "$RUN/_data/q2p1_param.dat"
mpirun --oversubscribe --wdir "$RUN" -np 28 "$BIN" > "$RUN/run_ss_dense.log" 2>&1
grep EL_SS_RADIUS "$RUN/run_ss_dense.log" > ss_dense.dat
```

## Variant: Matas–Asmolov inertial lift (W4.2 closure, option b)

`ELInertialLift = matas_asmolov` (with `ELLiftModel = none`) activates the
neutrally-buoyant inertial-lift closure: the Rc = 30 lift profile digitized
from Matas, Morris & Guazzelli (2004) JFM 515, figure 14 (their
Asmolov-method computation), F_r = ĝ(r/R)·ρ·Umax²·a⁴/(8√2·R²), table zero
crossing at s_eq = 0.675. Set `ELPrescribedUmax = 0.6` so Rc = Umax·D/ν =
30 matches the table's Reynolds number, and `"periodicZ_": true` (runs are
longer than the axial track).

Results (2026-07-04, 450 steps / t = 0.9 of a 2000-step run, np=28):
- Migration directions correct: rmin (s=0.1) drifts outward, rmax (s=0.9)
  inward, both toward s_eq; trajectories linear in t.
- Profile shape: measured slope ratio |rmax/rmin| = 7.34 vs table
  ĝ(0.9)/ĝ(0.1) = 7.12 (3%).
- Absolute rates: ds/dt = +1.56e-5 (rmin) / −1.15e-4 (rmax) vs ε=1
  quasi-steady predictions +1.87e-5 / −1.33e-4 (ratio 0.84–0.86). The
  deficit is the Di Felice ε^(−χ−1) voidage factor acting on the particle's
  own kernel deposit (ε ≈ 0.97 at the particle ⇒ ~15% extra lateral drag) —
  an unresolved-EL self-influence effect, not a closure defect.
- Full convergence to the annulus at a/R = 0.05 takes O(10⁴) time units
  (physical: v_migr ∝ (a/D)³; SS experiments needed meters of tube). The
  convergence demonstration run uses d_p = 0.1 (64× faster, ~250 t.u.,
  ~50k steps at dt = 0.005) — pending.

## Acceptance (plumbing)

- Both runs complete cleanly: no fatal/OOB/clipping lines, dt assertion
  passes, `count= 16` on every `EL_SS_RADIUS` line (no particle lost).
- Lift OFF: radial positions frozen — |Δ(r/R)| < 0.02 for rmin, rmean, rmax
  over the whole run (drag in an axisymmetric axial field has no radial
  component; any drift would indicate a transfer/geometry bug).
- Lift ON: measurable radial response distinct from OFF (transient kick per
  the expectation above); sign consistent with lagging-particle lift
  (toward centerline, i.e. rmax decreasing) while slip is axialy negative.
- Record rmean/rmin/rmax trajectories for both runs in Results.

## Results (2026-07-04, np=28, build-el-phase2-pe-gcc14)

All three runs: exit 0, `count= 16` on every sample, no fatal/OOB/clipping,
EL_VOLUME_CONSERVATION at machine zero (deposited = 16·V_p = 1.047198e-3).

| run | rmean t=1.3 | rmin t=1.3 | rmax t=1.3 | verdict |
|---|---|---|---|---|
| neutral, lift ON  | 0.4995905 | 0.0998234 | 0.8996307 | transient-only kick, PASS |
| neutral, lift OFF | 0.5000000 | 0.1000000 | 0.9000000 | radially frozen, PASS |
| dense (ρ_p=2, g_z=−0.981), lift ON | 0.4980670 | 0.0991647 | 0.8977161 | sustained inward drift, PASS |

- **Neutral ON**: inward kick Δ(r/R) = −(2…4)e-4 completed by t ≈ 0.1
  (≈14 τ_p, the start-from-rest slip relaxation); afterwards static. Force
  budget confirms the mechanism: lift ≈ 1e-4 at step 1 → O(1e-11) at step
  650 (slip → 0 ⇒ closure lift → 0, as predicted above).
- **Neutral OFF**: rmean/rmin/rmax frozen at 0.5/0.1/0.9 to all printed
  digits (gate < 0.02) — drag in the axisymmetric field has no radial
  component; transfer/geometry clean.
- **Dense variant**: steady linear inward migration. Measured d(r/R)/dt =
  1.30e-3 (outermost, seeded 0.9) and 4.18e-4 (innermost, 0.1) vs
  pure-Saffman estimates 1.72e-3 / 5.72e-4 — a consistent 0.73–0.76×, the
  expected trim from the Mei finite-Re correction (shear Re ≈ 0.3). Sign
  (lagging particle → toward centerline) and radius-dependence both correct.
  Step-650 budget: drag_z +1.024e-3 vs grav_buoy_z −1.027e-3 (steady
  settling slip), sustained transverse lift ≈ 2.5e-6.

**Campaign finding (for W4.2 planning):** with slip-based closures, a
neutrally buoyant particle produces no sustained migration — the coupled
Segré–Silberberg run as planned would show a flat concentration profile for
physical reasons, not bugs. Options: (a) accept W4.2 as a null test of the
closure set and document; (b) add a dedicated neutrally-buoyant inertial-lift
closure (Ho & Leal 1974 / Asmolov 1999 f_L(Re, r/R) tables) as a small
Stage-4 work item and use it for W4.2. Decision pending.
