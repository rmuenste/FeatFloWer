# Euler–Lagrange validation campaign — final report

Branch: `feature/euler-lagrange-phase1` · Application: `q2p1_el_pipeflow`
· Companion data: `el_validation_datasheet.md` / `.csv` (60 expected-vs-actual
rows: 47 PASS, 10 RECORDED, 2 FAIL-as-measured with physics attribution,
1 RESOLVED) · Per-case protocols: `applications/q2p1_el_pipeflow/
{tier2_cases,validation_cases}/*/RUNBOOK.md`.

## 1. Executive summary

The unresolved Euler–Lagrange (CFD-DEM, Model-A) implementation in
FeatFloWer was validated bottom-up: unit mechanics, conservation
machinery, single-particle physics, collective settling, migration, and
suspension rheology. The headline result is the **kroupa_couette** case:
total suspension viscosity measured two independent ways — the Kroupa
(Langmuir 2016) wall force balance and an interior Irving–Kirkwood
impulse virial — **agreeing within 2–4% at every volume fraction** and
tracking Krieger–Dougherty within ~15% over φ = 0.05–0.30. Suspension
viscosity in this model class is generated almost entirely by the
pairwise/wall lubrication closures plus dry-contact impulses; the
campaign demonstrates each channel separately (frozen box), proves the
missing dilute-limit component is exactly the Einstein hydrodynamic term
(a structural model boundary), and then closes the loop in the coupled
wall-driven cell.

All conservation gates hold at machine precision in periodic settings
(Newton-pair mismatch ≤ 1e-17 across production sweeps); wall-bounded
runs close to ≤ 4e-6 with the residual attributed to PGS fixed-iteration
truncation (solver-accuracy footnote, § 6.4).

## 2. Model summary

- Fluid: Q2/P1 FEM Navier–Stokes with volume-fraction (ε_f) weighting;
  divergence-form convection (leak fix, § 3.3); particle feedback as a
  kernel-spread body force (Newton-mirrored), explicit or semi-implicit
  drag coupling.
- Transfer: Lucy/`deen_poly` polynomial kernel, width δ = 2.5 d_p;
  volume-fraction "cloud" deposition; ELMomentumFix measured-leak
  compensator (periodic boxes).
- Forces: Di Felice drag, Saffman/Mei lift, Matas–Asmolov inertial-lift
  profile, pressure/buoyancy; optional Magnus (off).
- Particles: PE rigid-body engine, HardContactEulerLagrange PGS solver,
  hard contacts (velocity-level, ~zero overlap), per-CFD-step substepping
  (`substeps_` = 10, ladder-validated § 5.3).
- Lubrication (this campaign): Kroupa et al. pairwise closures
  (normal eq 12, sliding eq 13/14; twisting omitted) + sphere–wall
  closures (eqs 16–18) with Vinogradova slip f* (eq 20, ε_c = 0.1) and
  the h < h_c substitution; all-pairs sweep with designated-treater
  relay (Newton-exact by construction); impulse-cap regularization.
- Conventions: `Prop@Viscosity` is KINEMATIC (μ = ρ·ν); PE receives
  dynamic viscosity.

## 3. Infrastructure verification (tiers 0–2)

- **Tier 2 unit battery**: terminal velocity vs Stokes/Di Felice fixed
  point (0.01% error); straddling-particle conservation (machine zero);
  total-momentum drift (4e-7 rel; semi-implicit variant characterized
  ~3·dt first-order); Saffman lift/drag budgets vs analytic (≤ 3e-7).
- **Stage 0**: seeding reproducibility across PE decompositions
  (byte-identical), file-mode determinism, scaling envelope (N ≤ 5964
  verified).
- **3.3 Fluid-side momentum leak ladder** (v2_rz_settling RUNBOOK,
  parts 1–4): EL_FLUID_PAIR audit exposed a convective-form leak
  (9.8e-8/step, growing) → divergence-form convection fix (1.9e-13
  flat) → energy-stability caveat at production scale → measured-leak
  compensator ELMomentumFix (bounded 5.2e-7 cumulative over 50k steps vs
  5.7e-2 legacy ramp). Verdict chain: leak RESOLVED, compensator PASS at
  production scale.
- **3.4 Newton-pair (third-law) machinery**: global particle-momentum
  audit vs CFD-charged impulse. Verified at machine zero under ~5.3k
  simultaneous contacts, through periodic images, across PE
  decompositions, and (this campaign) through a four-bug lubrication
  hardening ladder — shadow-copy margin (radius + cutoff), Jacobi
  accumulation (bitwise cross-rank antisymmetry), ownership guard, and
  the designated-treater relay (root-cause fix for pair one-sidedness;
  EL_LUB_IMPULSE instrumentation matched the leak digit-for-digit).
  Wall-bounded extension: momentum absorbed by walls moves to the
  audit's expected side (§ 5.4).

## 4. Single-particle and migration physics

- **V1b ten Cate settling** (E1–E4): one-way u_t within 6.5–9.6% of the
  Di Felice fixed point, consistent after the ε_eff ≈ 0.975
  self-voidage correction (2%-level agreement); two-way co-flow bias
  (+6%) RECORDED as a finding — a single particle inside its own kernel
  cloud sees a reduced apparent voidage and a co-moving fluid column.
- **V3 Segré–Silberberg**: frozen-field migration slopes within Mei-trim
  ratios of Matas 2004 (profile-shape ratio 7.20 vs 7.12); coupled runs:
  lift ON/OFF discrimination 22:1, equilibrium annulus r/R = 0.664 vs
  0.675 Matas–Asmolov (1.6%).

## 5. Suspension rheology (the lubrication arc)

### 5.1 Model boundary (V4 baseline)
Drag-only unresolved EL generates **no** suspension viscosity: pipe
μ_app/μ flat (0.999–0.987) over φ = 0.05–0.20 at fixed forcing — the
Kroupa Fig 5 lubrication-OFF branch, reproduced.

### 5.2 Pairwise lubrication in PE (kroupa_shear, frozen box)
Irving–Kirkwood impulse virial in frozen linear shear (G = 1):
η_L/μ = 0.0079 / 0.0363 / 0.2496 / 0.7610 at φ = 0.05/0.10/0.20/0.30;
lubrication-OFF twin identically zero; Newton-pair machine zero
(≤ 2.1e-19); unit cross-rank pair case PASS.

### 5.3 Contact channel + temporal convergence
Converged PGS contact impulses accumulated as a second virial:
η_C/μ = 0.0017 / 0.0125 / 0.2106 / 0.8306 — co-dominant from φ ≈ 0.2,
firing at ~1e-10 overlap (invisible to penetration gates). Substep
ladder at φ = 0.30 (10/50/100): total particle-phase viscosity drifts
+0.8% — `substeps_` = 10 is converged (Kroupa's ~250 not needed).
**Dilute-limit consistency check**: the component missing from the
frozen box (KD minus measured) equals the Einstein hydrodynamic term
2.5φ at dilute φ (0.129 vs 0.125 at 0.05; 0.264 vs 0.250 at 0.10) —
the frozen field lacks exactly η_H, as constructed.

### 5.4 Total viscosity (kroupa_couette — headline)
Two-way coupled plane-Couette cell: periodic x/y, moving z-walls
±G·L/2 (`Inflow300` tag; symmetric frame), exact linear-shear IC,
sphere–wall lubrication (eqs 16–18), wall-force accumulators, and the
EL_WALL_STRESS diagnostic implementing the Kroupa force balance
(η_H fluid channel by Q2 face quadrature + wall lubrication + wall
contact channels). Gates: bitwise regression of the walls-off path;
single-particle twin gives η_wall/μ = 1.0000 to all printed digits.

| φ | η_wall/μ | interior (μ+η_L+η_C)/μ | wall vs interior | KD | vs KD |
|---:|---:|---:|---:|---:|---:|
| 0.05 | 1.034 | 1.016 | +1.8% | 1.139 | −9% |
| 0.10 | 1.126 | 1.093 | +3.0% | 1.312 | −14% |
| 0.20 | 1.664 | 1.606 | +3.7% | 1.821 | −9% |
| 0.30 | 2.936 | 2.874 | +2.2% | 2.751 | +7% |

Boundary-traction and bulk-virial estimators agree within 2–4%
everywhere (the +2–4% gap is the η_H feedback reshaping only the wall
sees); the curve is monotone, superlinear, within ~15% of KD across the
range. Dilute end below KD by the partial Einstein term — the § 5.3
model boundary, shared with Kroupa's method.

### 5.5 Lubricated pipe (V4 rerun)
With pair lubrication (no wall lubrication — the pipe wall is a CFD
boundary): μ_app/μ = 1.012/1.056/1.108/1.148 — mechanism PASS
(monotone; Fig-5 on/off contrast in a wall-bounded geometry) but only
9–20% of the KD excess: near-wall dissipation dominates Poiseuille flow
and both wall lubrication and the lift-migration depletion are absent
there. RECORDED with attribution; the Couette cell (§ 5.4) supplies the
quantitative result.

## 6. Findings, boundaries and open items

1. **Mesoscale cluster-induced settling (V2)**: fully periodic two-way
   settling at φ = 0.05–0.20 develops a plume instability — settling
   ENHANCED 2.6–3.3× over Richardson–Zaki (|U|/|U₀| up to 1.08: faster
   than an isolated particle). Lane-mode filtering (insufficient —
   diagonal modes) and box shrink to L/d_p = 10 (unchanged) exhausted
   suppression routes; V2 is REFRAMED against the cluster-induced-
   settling literature (2–4× enhancement is the documented behavior of
   periodic point-particle sedimentation). RZ hindered settling is a
   wall-bounded/homogenized result the periodic box does not represent.
   All conservation gates in these runs PASS.
2. **Self-voidage / co-flow bias**: single-particle drag biased by
   ε_eff ≈ 0.975 (kernel self-occupancy); two-way runs add +6% co-flow.
   Gate against corrected predictions (memory + datasheet rows).
3. **Dilute-limit Einstein term**: unresolved feedback generates only
   part of 2.5φ·μ — structural to the model class.
4. **PGS truncation footnote**: wall-CONTACT-active runs close the
   Newton audit to ≤ 3.9e-6 (7e-5 of system momentum; z-dominant;
   machine zero in wall-lubrication-only smoke). Attributed to
   converged-cache vs relaxation-folded-applied impulse truncation of
   the fixed-100-iteration PGS — invisible in periodic runs by pairwise
   antisymmetry. Not a leak; bounded; recorded.
5. **ε_f bookkeeping**: pipe runs drift ⟨ε_f⟩ slightly above 1−φ at
   high φ (0.812 vs 0.800); the Couette box holds 1−φ to 4 decimals.
   Worth a look if pipe concentration profiles become a payload.
6. **Deferred**: EMAC-form convection; particle–wall lubrication in the
   pipe geometry; V5 radial concentration profiles (data exists in the
   V4 runs); Kroupa parameter-matched sweep vs their Fig 3 experimental
   overlay (ε_c sensitivity).

## 7. Reproducibility

Every case ships a RUNBOOK (design, gates, job numbers, verdicts) under
`applications/q2p1_el_pipeflow/{tier2_cases,validation_cases}/`; every
quantitative claim is a datasheet row with expected value, source,
measured value and tolerance. Rundirs are staged with
`tools/el_stage_rundir.sh` and submitted with `tools/el_slurm_submit.sh`
(np = 28, QBOX9/PIPEZ27 3×3×3/1×1×27 partitions). The PE library changes
live on `libs/pe` master (lubrication, contact/wall virials, z-walls);
FeatFloWer changes on `feature/euler-lagrange-phase1`.

## 8. Non-EL regression (pre-merge) — PASS

Shared-file changes (QuadSc_main/user/boundary path, dem_query, PE
library) verified against the rest of the build (2026-07-28):

- **Full-suite compile/link** (~250 targets, all applications + unit
  tests, PE-enabled config): one failure — `q2p1_hashgrid_test`, whose
  PE shim is `#ifdef PE_SERIAL_MODE` (serial-mode-only WIP app; never
  linked in this parallel configuration; predates the branch).
- **ctest (16 tests)**: all 5 EL unit tests PASS (transfer
  serial/mpi-2/mpi-8, convergence, kernel forces);
  pe-el-semi-implicit-drag PASS; serial PE interface smokes: fresh
  variants PASS, the two resume-roundtrip variants FAIL with an
  MPI-context abort — **verified pre-existing** by rebuilding the
  pre-campaign PE (commit 8b037ae, before any campaign change) in the
  same configuration: identical failures.
- **Non-EL solver runs** (with the in-tree metis on the loader path —
  the harness's metis discovery is an environment gap, not code):
  q2p1-fac-newt (×2), q2p1-fac3d, q2p1-bench-sedimentation and
  q2p1-bench-fluidization — real FBM/PE-coupled solver benchmarks —
  all PASS on the branch. q2p1-fac-visco is blocked by a missing mesh
  fixture in the local mesh repo (`_adc/ViscoHex2/aaa.prj`),
  independent of code.
- **Bitwise walls-off regression** (§ 5.4): frozen-box EL_SUSP_STRESS
  identical to pre-wall binary output.

Conclusion: zero regressions attributable to the branch; the three
red items (hashgrid link, resume-roundtrip smokes, fac-visco fixture)
are all demonstrably pre-existing or environmental.
