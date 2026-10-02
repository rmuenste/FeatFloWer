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
