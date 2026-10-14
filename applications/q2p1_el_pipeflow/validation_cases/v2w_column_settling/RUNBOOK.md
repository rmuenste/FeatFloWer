# Runbook: v2w_column_settling — wall-bounded settling column (EL twin of DNS D3.2)

Staged 2026-10-02 on branch `feature/el-rz-postfix` (worktree `FeatFloWer-el`).
Status: **staged, G0/G1 gates pending, nothing submitted for the ladder.**

## 1. Purpose and background

The periodic Richardson–Zaki case (`v2_rz_settling`) develops a box-scale
concentration instability and settles enhanced. The DNS campaign settled the
same question in a WALLED column (D3.2, rows `d32_ut_ref`, `d32_n40_smoke`,
`d32_phi_ladder`, `d32_wide_attribution`): finite clouds of N = 20…120
resolved spheres in a closed 6d × 6d × 24d column settle stably hindered,
U/u_t collapses on (1−φ_cloud)^n with a confined exponent n = 4.57 (free
intercept 4.31), velocity spread saturates, no clustering. This case runs the
SAME column, SAME fluid, SAME particles and the byte-identical particle clouds
through the unresolved EL solver. It is therefore both (a) the wall-bounded
hindered-settling measurement RZ actually describes, and (b) the first
DNS ↔ EL same-physics twin of the campaign plan (stage D5.1), with the DNS
ladder as the in-house reference. Reference numbers: DNS u_t = 0.4061
(Re_t = 86.4), U/u_t = 0.895/0.90/0.88 (N=20), 0.80/0.84/0.83 (N=40),
0.68/0.71/0.68 (N=80), 0.60/0.56/0.61 (N=120); plateau window t = 15–25.

## 2. Configuration (mirrors the DNS deck; EL-specific keys stated)

- Mesh: `dkt_mesh_v1/coarseB` (12 × 12 × 48 coarse el, cell d/2, d = 1,
  x,y ∈ [−3,3], z ∈ [0,24], six `Wall` faces), partition `DKT_B_108`
  (x3-y3-z12, axis-uniform). Links (untracked fixture, see MESH.md):
  ```
  M=/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer
  for f in dkt_B.tri file.prj xmin.par xmax.par ymin.par ymax.par zmin.par zmax.par; do ln -sfn $M/dkt_mesh_v1/coarseB/$f mesh/coarseB/$f; done
  ln -sfn $M/dkt_mesh_v1/coarseB/_mesh/DKT_B_108 mesh/partitions/DKT_B_108
  ```
- Resolution: `MaxMeshLevel = 2` → h = d/4, kernel δ = 2.5 d = 10 h
  (campaign standard ≥ 4.5; level 3 = DNS D/h = 8 is the resolution rung).
  55 296 Q2 cells on 108 workers, np = 109 (2 nodes, like the DNS runs).
- Fluid: ρ_f = 1, μ = 0.0047 (`Prop@DynVisc`), `Prop@Gravity = (0,0,−1)` with
  `ELFluidGravity = Yes` (closed box, hydrostatic part absorbed by the
  pressure; `NoOutflow = Yes`), `ELMomentumFix = No` (walls invalidate the
  compensator). dt = 0.005 = PE `stepsize_`, 5000 steps to t = 25 (= 10 d/u_t).
- Particles: r = 0.5, ρ_p = 1.14, Di Felice drag, `ELLiftModel = none`
  (drag-only hindrance, as in v2; a lift-on rung is a later sensitivity),
  explicit coupling (τ_p = 13.5 ≫ dt), PE `substeps_ = 2`, dry contacts
  (restitution 0, friction 0 as DNS), lubrication off, PE `gravity_ = 0`
  (buoyancy is analytic in EL: (ρ_p−ρ_f)V g from `Prop@Gravity`).
- Walls on the PE side: `boxWallsEnabled_ = true` — six fixed global planes
  on the CFD box (new key, pe branch `feature/el-box-walls` on top of
  ce3bb5f; zero velocity; refused together with any periodic axis or the
  Couette z-walls).
- Seeding: `seedMode_ = file`; `particles.xyz` of the case = the DNS u_t
  reference sphere at (0,0,20); `clouds/n{20,40,80,120}_s{1,2,3}.xyz` are the
  DNS ladder clouds, copied byte-identical from the D3.2 rundirs
  (`q2p1_dns_rundir_d32_<n>_<s>/particles.xyz`; n40_s1 = `..._d32_n40`).
- Trajectory output: `ELParticleLogFreq = 10` writes
  `el_particles_rank<NNNN>.log` per worker (new, collective-free; DNS
  `particle_force.log` column layout). Merge with
  `python3 tools/el_particle_log_merge.py <rundir>` → `particle_force.log`,
  then analyse with the DNS tool unchanged:
  `python3 tools/d32_ladder_analysis.py <rundirs> --window 15 25 --ut <EL u_t> --area 36.0 --radius 0.5`.

Staging: `stage_ladder.sh <target-root>` → `q2p1_el_rundir_v2w_ut` + twelve
`q2p1_el_rundir_v2w_n<N>_s<S>`.

## 3. Gates, in order

- **G0 — regression twin.** Tier-2 `momentum_conservation` on the rebuilt
  binary (new keys default off: `ELParticleLogFreq = 0`,
  `boxWallsEnabled_ = false`) must stay BITWISE against job 137542 /
  147442 on every EL stream.
- **G1 — u_t smoke** (`v2w_ut`, 1 sphere, 5000 steps): banner
  `box walls enabled`, PE decomposition 3×3×12, `EL_SEED_POS` at (0,0,20),
  no NaN/OOB/clipping, Newton pair ≤ 1e-10 (walls: ≤ 4e-6 only once wall
  contacts fire), EL_FLUID_PAIR not a gate (walls). Settled u_t from the
  merged log; record vs DNS 0.4061 and vs the Di Felice fixed point with the
  ε_eff ≈ 0.975 self-voidage correction (V1b precedent). Wall time/step.
- **G2 — n40 seed 1** (the DNS smoke twin): hindered plateau U/u_t < 1 from
  t ≈ 14, cloud coherent, velocity spread saturates. Compare U/u_t with
  DNS 0.80.
- **G3 — ladder** (remaining eleven): fit n through the origin over
  ln(U/u_t) vs ln(1−φ_cloud); compare with DNS 4.57 (free intercept 4.31,
  rms 0.023) and with RZ/Rowe 2.7–3.0 at swarm Re 50–80.

Verdict bookkeeping: rows `v2w_*` in the EL datasheet; expected values are
the DNS rows above (in-house ground truth) with the RZ band as the
literature reference; the EL-vs-DNS difference is the D5.1 deliverable
(closure deficit), not a pass/fail on its own.

## 4. Submission (owner decision)

```bash
EL=/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-el
R=$EL/build-el-rz-postfix/applications
export EL_BIN=$EL/build-el-rz-postfix/applications/q2p1_el_pipeflow/q2p1_el_pipeflow
# 109 ranks on 2 nodes: use the job script written by stage (job_v2w.sbatch) or
# tools/el_slurm_submit.sh with EL_NP=109 EL_NODES=2 (see section 5).
```
DNS cost reference: 6 h per run at 109 ranks, level 3. EL at level 2 is
expected well under 1 h per run.

## 5. Known differences from the DNS column (declared, not hidden)

- Unresolved vs resolved: the EL particle is a point with a 2.5 d kernel;
  in a 6 d column the kernel of a centred particle reaches 0.5 d from the
  side walls (Shepard renormalisation applies near walls). Wall hindrance
  on an isolated particle is therefore represented only through the fluid
  solve, not through a resolved wall-sphere interaction.
- No particle–wall lubrication in EL (the planes are hard contacts only).
- `ncontacts` semantics on the current pe pin: hard contacts only.

## 6. Results so far (2026-10-02)

- **G0 PASS** (job 147524, 14 s): rebuilt binary (FF + pe feature/el-box-walls
  bf20411) reproduces job 137542 BITWISE on EL_MOMENTUM_ELEMINT,
  EL_FEEDBACK_CONSERVATION, EL_NEWTON_PAIR, EL_VOLUME_CONSERVATION,
  EL_MEAN_SLIP and EL_CONTACT_STATS; 702 PERIODIC_COMM; no rank logs written
  with ELParticleLogFreq = 0.
- **G1 PASS** (job 147525, nx-02/03, 109 ranks, 1 h 27 min for 5000 steps at
  level 2): banner `box walls enabled on x=[-3,3] y=[-3,3] z=[0,24]`,
  3x3x12 decomposition, seed at (0,0,20), no NaN/OOB/clipping, Newton-pair
  worst 8.5e-17 (no wall contact in 25 t.u.), five rank logs merged by
  tools/el_particle_log_merge.py into particle_force.log. Settled
  **u_t(EL) = 0.35734** (drift < 2e-6 per t.u. over the last 1.5 t.u.).
  Context: Dallavalle/Di Felice fixed point at eps = 1 is 0.368 (EL −3%,
  the V1b self-voidage band), Schiller–Naumann 0.398, DNS in the same
  column 0.4061 (EL −12%). The EL/DNS gap at the single-particle level is
  itself a D5.1 datum: the DNS sphere is wall-hindered through resolved
  flow, the EL sphere through the kernel-smoothed no-slip field plus the
  Dallavalle correlation's own offset from the standard drag curve.
  Cost: ~1.0 s/step incl. I/O at np = 109.
- **G3 ladder SUBMITTED** (owner go-ahead 2026-10-02): med partition,
  06:00:00, 109 ranks / 2 nodes / 40 GB each: n20_s1=147526 n20_s2=147527 n20_s3=147528 n40_s1=147529 n40_s2=147530 n40_s3=147531 n80_s1=147532 n80_s2=147533 n80_s3=147534 n120_s1=147535 n120_s2=147536 n120_s3=147537
  Analysis once complete:
  `for d in /data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-el/build-el-rz-postfix/applications/q2p1_el_rundir_v2w_n*_s*; do python3 tools/el_particle_log_merge.py $d; done`
  then `python3 tools/d32_ladder_analysis.py /data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-el/build-el-rz-postfix/applications/q2p1_el_rundir_v2w_n*_s* --window 15 25 --ut 0.35734 --area 36.0 --radius 0.5`.

## 7. Ladder RESULT (2026-10-03, jobs 147526/147573–147583, med, 109 ranks, 1 h 30–1 h 49 each)

Eleven first submissions died at 7 s on mixed nx+worldgames allocations (OFI
modex failure, known cluster defect); resubmitted with `--exclude=worldgames`
(helper fixed, commit 2a46e23b). All twelve runs then completed: no NaN,
Newton-pair worst 7e-15, hard contacts 0–13 per sample (N = 120), velocity
spread 0.014–0.027 saturating in time, clouds coherent — stable hindered
settling at every rung, no clustering (the walled column suppresses the
periodic-box instability for EL exactly as it did for DNS).

| N | φ_cloud | U/u_t EL (3 seeds) | U/u_t DNS (3 seeds) | EL/DNS |
|---:|---:|---|---|---:|
| 20 | 0.019–0.021 | 0.963 / 0.966 / 0.966 | 0.895 / 0.901 / 0.879 | 1.08 |
| 40 | 0.039–0.048 | 0.928 / 0.938 / 0.938 | 0.802 / 0.842 / 0.834 | 1.13 |
| 80 | 0.077–0.080 | 0.870 / 0.877 / 0.864 | 0.684 / 0.706 / 0.678 | 1.26 |
| 120 | 0.113–0.119 | 0.812 / 0.800 / 0.803 | 0.597 / 0.561 / 0.605 | 1.37 |

Fit (tools/d32_ladder_analysis.py, window t = 15–25, u_t = 0.35734):
through-origin n = 1.72; free-intercept n = 1.78, U/ut(phi->0) = 1.005 (rms ln-residual 0.005). DNS in the same column: 4.57 (free
intercept 4.31); unbounded Richardson–Zaki/Rowe 2.7–4.65.

Reading: the EL closure set (Di Felice local-voidage hindrance + kernel-smoothed
two-way coupling, no lift) produces a clean power law but with an exponent
~2.7× smaller than the resolved spheres. The EL/DNS ratio grows monotonically
with φ (1.08 → 1.37): the deficit is in the slope, not an offset, so it is a
collective-interaction deficit, not the single-particle u_t gap (which the
ratio already divides out). Candidate mechanisms, in the order the data
suggests: (i) the confinement-driven return flow that carries the DNS
hindrance (row d32_wide_attribution) is weakened in EL because the cloud's
reaction force is smeared over kernels of 2.5 d in a 6 d column, so the
resolved backflow through the cloud is smaller; (ii) Di Felice ε^−χ at local
ε ≈ 0.9–0.98 adds only a few percent; (iii) no resolved near-field
interactions. This is the first D5.1 same-physics twin and its deliverable is
exactly this closure deficit; it is NOT a pass/fail against RZ. Verdict rows:
mechanism PASS (hindered, stable, power-law collapse), quantitative RECORDED
(n = 1.72 vs DNS 4.57).

Next rungs that would discriminate (i) from (ii): a kernel-width rung
(δ = 1.5 d) and a level-3 resolution rung, both from the same clouds.

## 8. Discriminator rungs SUBMITTED (2026-10-03, owner go-ahead)

Staged by `stage_rungs.sh` from the same case and clouds:
- **k15** — kernel width δ = 1.5 d (δ/h = 6 at level 2), u_t + all twelve
  clouds, med 06:00:00, 109 ranks / 2 nodes / 40 GB.
- **l3** — MaxMeshLevel 3 (h = d/8 = the DNS cell size, δ/h = 20), u_t +
  N = 40 ×3 + N = 120 ×3 (the DNS two-rung discriminator shape), long
  1-12:00:00, 109 ranks / 2 nodes / 200 GB.
Each rung is analysed against ITS OWN u_t. Jobs: k15_ut=147639 k15_n20_s1=147640 k15_n20_s2=147641 k15_n20_s3=147642 k15_n40_s1=147643 k15_n40_s2=147644 k15_n40_s3=147645 k15_n80_s1=147646 k15_n80_s2=147647 k15_n80_s3=147648 k15_n120_s1=147649 k15_n120_s2=147650 k15_n120_s3=147651 l3_ut=147652 l3_n40_s1=147653 l3_n40_s2=147654 l3_n40_s3=147655 l3_n120_s1=147656 l3_n120_s2=147657 l3_n120_s3=147658

## 9. Kernel-width rung k15 RESULT (2026-10-04, jobs 147639–147651, all COMPLETED, 1 h 25–1 h 32 each)

δ = 1.5 d at level 2 (δ/h = 6), same twelve clouds, own u_t. Clean: no NaN,
Newton-pair worst 5.5e-15, stable hindered settling at every rung (spread
0.026–0.045, saturating).

- u_t(k15) = **0.31890** vs 0.35734 at δ = 2.5 d (−10.8%): the narrower
  kernel samples a fluid more disturbed by the particle itself (stronger
  self-voidage / self-induced velocity), as the V1b self-voidage picture
  predicts in the wrong direction for "narrower is more accurate".
- Seed means U/u_t: 0.984 / 0.964 / 0.917 / 0.860 at N = 20/40/80/120
  (δ = 2.5 d: 0.965 / 0.935 / 0.870 / 0.805; DNS 0.892 / 0.826 / 0.689 / 0.588).
- **n = 1.15** through origin (rms 0.012; seed span at N = 40 widened to
  0.040) vs 1.72 at δ = 2.5 d vs DNS 4.57.

Reading: narrowing the kernel moves the EL hindrance law AWAY from DNS. Kernel
smearing is therefore not the mechanism behind the slope deficit; the wider
kernel helps, by coupling each particle to more of the cloud-scale return
flow and by letting the cloud's reaction force act on a larger fluid volume.
The deficit sits in what the point-particle closure can represent (local
voidage factor ε^−χ, no resolved near-field / wake interactions), or in the
resolution of the return flow — which the l3 rung (section 8) tests. Row
`v2w_k15_exponent` RECORDED.

## 10. Level-3 rung l3 RESULT (completed 2026-10-05, recorded 2026-10-14; jobs 147652–147658, 3 h 27 – 5 h 43 each)

MaxMeshLevel 3 (h = d/8, the DNS cell size; δ = 2.5 d = 20 h), own u_t, the
N = 40 and N = 120 triplets. Clean: no NaN, Newton pair ≤ 7e-15.

- u_t(l3) = **0.35714** vs 0.35734 at level 2 (−0.06%).
- U/u_t: N = 40 → 0.928 / 0.938 / 0.939; N = 120 → 0.812 / 0.800 / 0.803.
  Level-2 twins of the same six clouds: 0.928 / 0.938 / 0.938 and
  0.812 / 0.800 / 0.803. Particle-sampled voidage and velocity spread also
  identical to three digits.
- Two-rung exponent: **n = 1.73 at level 3, 1.73 at level 2** (rms 0.006 both).

Reading: NULL RESULT, which is the informative outcome. Doubling the
resolution of the return flow changes neither the single-particle velocity
nor the collective hindrance. Together with section 9 (narrower kernel →
weaker hindrance) the discriminators eliminate kernel smearing and
resolution; what remains is the input the drag closure receives: the
kernel-sampled voidage (diluted across the cloud edge: ε seen 0.927 vs
1 − φ_cloud 0.888 at N = 120) and the self-disturbed slip velocity. Next:
an undisturbed-field correction at the sampling → force boundary
(el_transfer: EL_INTEGRATE_PARTICLE → EL_FORCES_FROM_PARTICLE), gated by
the tier-2 bitwise twin with the key off, V1b (self-voidage bias must
shrink) and this ladder. Row `v2w_l3_exponent` RECORDED.

## 11. Model-A pressure-force rung pf SUBMITTED (2026-10-14, owner go-ahead; deck-only, binary = build-el-rz-postfix)

Motivation (liggghts-public-34 survey, verified on our logs 2026-10-14): the EL
fluid is plain incompressible NS and the particle gets ρ_f buoyancy only, so
the suspension-weight share of the pressure gradient never reaches the
particle (one power of ε lost from Di Felice's exponent) and no
volume-averaged return flow exists (sampled u_f at particles ≈ 0 vs
−φu_p/ε = 0.008 / 0.023 at N = 40 / 120; a second power lost). Ceiling under
the current equations ≈ χ−1 = 2.1 at Re ≈ 80; measured 1.72.

This rung tests the first of the two with existing keys: `ELFluidGravity =
No` (pressure becomes dynamic; in a closed box the dropped uniform body force
is a pure gradient) + `ELPressureForce = Yes` (particle gets −V∇p_dyn on top
of the analytic (ρ_p−ρ_f)Vg). u_t + N = 40 ×3 + N = 120 ×3, level 2, δ = 2.5 d,
own u_t. Expected: n rises toward ≈ χ ≈ 3; the volume-averaged continuity
(worktree FeatFloWer-va, branch feature/el-volume-averaged) is the second
step. Jobs: pf_ut=149697 pf_n40_s1=149698 pf_n40_s2=149699 pf_n40_s3=149700 pf_n120_s1=149701 pf_n120_s2=149702 pf_n120_s3=149703
