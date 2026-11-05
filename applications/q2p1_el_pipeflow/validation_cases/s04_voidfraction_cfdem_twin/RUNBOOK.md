# s04_voidfraction_cfdem_twin — FeatFloWer twin of CFDEM `voidfractionTest`

FeatFloWer side of case 04 of the CFDEM cross-code comparison: void-fraction deposition and
d(ε)/dt with prescribed particle motion. Spec and tools:
`/data/warehouse17/rmuenste/code/liggghts-dem/LIGGGHTS-PUBLIC/md_docs/case_specs/04_voidfraction/`
(`SPEC.md` §3 identities, §6 FeatFloWer recipe, §8 protocol, §12 CFDEM results;
`featflower/{make_box_mesh.py,make_particles.py,stage_rundir.sh,postproc_featflower.py}`,
`reference/predict_voidfraction.py`). Binary `build-va` `q2p1_el_pipeflow` (branch
`feature/el-volume-averaged`, built 2026-10-15 15:14). Run directories
`/data/warehouse17/rmuenste/runs/vft04/{VA_L2,VB_L2,VA_L3,VB_L2_ptight}`, local Slurm 2026-10-16.

Setup: closed column 0.1 × 0.1 × 1.0 m, coarse mesh 10 × 10 × 100 hex (h = 1 cm) from
`make_box_mesh.py`, partition `VFT10` = ten z-slabs (pe `processesZ_ 10`), 11 ranks. 1000 spheres
d = 3 mm, ρ_p 3000, fluid ρ 10 / μ 1e-4, prescribed w = −0.3 m/s, no gravity, no forces
(`ELApplyParticleForces No`, `ELApplyFluidFeedback No`), kernel `deen_poly` 2.5 d, dt 1 ms,
200 steps. **Positions are the CFDEM reference run's own** (`run/vft_VA-divided_r2`, LIGGGHTS dump
frame at step 5000 shifted back by w·0.05 s; identical to machine precision to the
`vft_diag_trilinear_dump` frame at step 100 shifted by w·1 ms). VA: `ELPrescribedField
linear_shear` with rate 0 (u ≡ 0, no fluid solve). VB: `ELVolumeAveraged continuity`, fluid solved.
Decks here: `q2p1_param_VA.dat` / `_VB.dat` (spec copies; level set per run), `example.json`
(with `seedMinGap_ 1e-6`, see §3), `particles.xyz` (the CFDEM positions).

## 1. Results

Identities (spec §3): G1 Σ(1−ε)V = N V_p = 1.413717e-5 m³; G2 Σ dε/dt V = 0;
G3 Σ z dε/dt V = N V_p |w| = 4.24115e-6 m⁴/s; G5 (VB) ∫ε u_z dV / V_box = 4.24115e-4 m/s.
Predictor on these positions: G3 ratio 1.00001 (L2) / 1.00000 (L3), L1_dedt 1.129e-3 (L2) /
1.332e-3 (L3) m³/s, drag diagnostic 4.41e-3 / 4.53e-3 N.

| run (job) | wall | G1 max rel. error | w | G2 max | G3 per output mean / sd | G3 slope ratio | drag diag |
|---|---|---|---|---|---|---|---|
| VA_L2 (150194) | 1 min | 4.0e-15 | −0.3 exact | 3.1e-12 | 4.24107e-6 / 7.1e-10 | 1.00000 | 4.412e-3 |
| VA_L3 (150196) | 7.6 min | 8.5e-15 | −0.3 exact | 1.2e-12 | 4.24109e-6 / 2.6e-10 | 1.00000 | 4.533e-3 |
| VB_L2 (150195) | 7.5 min | 4.0e-15 | −0.3 exact | 3.1e-12 | 4.24107e-6 / 7.1e-10 | 1.00000 | 4.473e-3 |

G3 per-output scatter: 1.7e-4 (L2), 6e-5 (L3) of the identity; CFDEM trilinear 5 % per step,
divided 0.99992 mean. All three PASS G1 (gate 1e-12), G2, G3 (gates 2e-3 / 1e-3).

**VB return flow (G5, gate 2 %):** `EL_VA_SOURCE` L1_dedt time mean 1.1318e-3 m³/s (last step 1.1291e-3; predictor 1.1291e-3 at step 100, +0.24 %),
|∫g|/L1_dedt ≤ 1.06e-3 (compatibility removal, below α_max 0.07 as expected). Superficial mean
u_z from `EL_MEAN_SLIP uf_super` (a domain mean: here exactly the quantity wanted):

| window | ⟨u_z,super⟩ | / 4.24115e-4 |
|---|---|---|
| t ≥ 0.002 (spec window) | 4.1567e-4 | 0.980 (−2.0 %) |
| t ≥ 0.1 (mean) | 4.1611e-4 | 0.981 (−1.9 %) |
| step 2 (start-up) | 4.584e-4 | 1.081 |
| min / max | 4.111e-4 / 4.584e-4 | 0.969 / 1.081 |

PASS at the gate edge (−2.0 % on the spec window, −1.9 % over the second half; single samples at t = 0.1 / 0.2 are −1.4 / −1.2 %). CFDEM VB-b-divided: +0.7 %
with write scatter 0.85–1.16; VB-a −11 %. Discriminator for the deficit submitted as
`VB_L2_ptight` (`Pres@MGCriterion1 1d-7`, `MGMaxIterCyc 60`; the deck's 1d-3 relative pressure
defect reduction is the first suspect) — result in §4.

## 2. Verdicts (spec §8.3 rows)

| quantity | expected | measured | verdict |
|---|---|---|---|
| vft_VA_L2 volume_rel_error_max | 0 (≤ 1e-12) | 4.0e-15 | PASS |
| vft_VA_L2 particle_velocity_z | −0.3 | −0.300000000000 | PASS |
| vft_VA_L2 sum_z_depsdt_V | 4.24115e-6 (2e-3) | 4.24107e-6 (−1.9e-5) | PASS |
| vft_VA_L3 sum_z_depsdt_V | 4.24115e-6 (1e-3) | 4.24109e-6 (−1.4e-5) | PASS |
| vft_VB_L2 uf_super_z | 4.24115e-4 (2e-2) | 4.1567e-4 (−2.0 %) | PASS (edge) |
| vft_VB_L2 L1_dedt_mean | 1.136e-3 (1e-2) | 1.1318e-3 (−3.7e-3) | PASS |
| twin uf_super_z EL / CFDEM VB-b-divided | 1 (2e-2) | 0.980 / 1.0137 = 0.967 (CFDEM re-run writing every step, run 90: 1.0137 ± 0.0053; the 20-write 1.007 had SE 0.019) | FAIL by 1.3 % |

## 3. Staging corrections (to feed back into the spec)

1. `featflower-partition 10 -3 1 VFT10 …` aborts ("requires more than 1 subgrid"); with
   `10 -3 10` only `sub0001` receives `GRID0001.tri` (the partition count is distributed over the
   subgrids: 1 ⇒ first subgrid only) and the solver dies on `sub0002/GRID0001.tri`. Correct call:
   **`featflower-partition 10 -3 10 VFT10 vft_box/file.prj`** (ten partitions over ten z-slabs;
   slab z-ranges exact, 1000 coarse elements each). `stage_rundir.sh` therefore never wrote
   `job.sbatch`; jobs were submitted with `tools/el_slurm_submit_np.sh`.
2. The CFDEM dumps have no populated frame at step 0 (`run 0 upto` does not dump the inserted
   particles); `make_particles.py --from-dump` crashes on the empty frame. Positions were rebuilt
   from the first populated frame by the exact back-shift w·t (velocities exactly −0.3, no images).
3. pe refuses the CFDEM packing: minimum surface gap 3.8e-5 m = 0.013 d is below `seedMinGap_`'s
   default 0.1 d. Set `"seedMinGap_": 1.0e-6` in `example.json`.
4. `postproc_featflower.py` needs `--output-dt 0.001`: our VTK file indices are step numbers
   (00000, 00009, 00019 …), not output counts, so the default 0.01 makes the G3 slope 10× low.
5. Level 3 VA takes 7.6 min on 11 ranks (the spec's 30-min worry is moot).
6. `grep NaN` on VA logs matches the multigrid timer statistics (no fluid solve); harmless.

## 4. Discriminators for the return-flow deficit (both NULL)

| run (job) | change vs VB_L2 | ⟨u_z,super⟩ spec window | t ≥ 0.1 |
|---|---|---|---|
| VB_L2 (150195) | — | 4.15670e-4 (0.9801) | 4.16109e-4 (0.9811) |
| VB_L2_ptight (150198) | `Pres@MGCriterion1 1d-7`, `MGMaxIterCyc 60` | 4.15670e-4 (0.9801) | 4.16109e-4 (0.9811) |
| VB_L2_nl3 (150199) | + `Velo@NLmin 3`, `NLmax 6` | 4.15678e-4 (0.9801) | 4.16117e-4 (0.9811) |

Neither the pressure tolerance nor the nonlinear iteration count moves the result (2e-5): the
−2.0 % is a property of the discrete continuity source, not of solver convergence. What the
discretisation implies (`el_volume_averaged.f90:86-106`): per element g = −(∂_tε + u·∇ε)/ε with
ε = max(ε_el, 1e-3); the source enters only the constant P1 pressure dof of each element as
(g − ḡ)·V/dt. Testing the discrete constraint with the P1 function z (its linear modes are
constrained homogeneously) gives ∫u_z dV = −Σ_el z_c g_el V_el exactly at convergence, and the
1/ε factor would make that an *excess* over N V_p|w| of order the |∂_tε|-weighted α, not a
deficit; ∫α u_z dV is only ≈ 0.1 % of the integral. The sign therefore points elsewhere
(candidates: the kernel-deposited ε vs the P1 constant-mode source near the top wall, where
the cloud starts within one kernel radius; the ḡ removal's interaction with the cloud's
asymmetric position; the diagnostic's own ε-weighting). Next discriminator is a code
diagnostic printing Σ z_c g V, ∫u_z dV and the per-element constraint residual (binary change,
own twin gate); not done here. Recorded as a FeatFloWer finding, gate PASS at the edge.

## 5. Peer follow-up (liggghts-public-1e, 2026-10-16)

All four staging corrections are folded into their spec and scripts (partitioner call, dump
back-shift with vz/image checks reproducing `particles_cfdem.xyz` to 5e-8 m, `seedMinGap_`,
VTK index time from the deck). Their re-derivation from our run dirs reproduces every gate.
Correction: the L1_dedt time mean is 1.1318e-3 (the 1.1291e-3 first quoted was the last step).
CFDEM VB-b-divided re-run writing every step (run 90, `run/vft_VB-b-divided_w1`): G5 = 1.0137 ±
0.0053 (t ≥ 0.002, 199 writes), 1.020 ± 0.009 (t ≥ 0.1); the earlier 1.007 from 20 writes had a
standard error of 0.019. Twin ratio 0.980 / 1.0137 = 0.967: FAIL by 1.3 % beyond the gate. The two
codes bracket the identity (CFDEM +1.4 %, FeatFloWer −2.0 %); their excess is suspected to be the
face-flux vs cell-centre reconstruction mismatch [I]. Agreed next step: the FeatFloWer constraint
diagnostic (§4), gated on the identity itself rather than on EL/CFDEM.

## 6. Return-flow cause tests T1, T2, T4 (SPEC 04 §12.3-12.4; jobs 155110-155115, 2026-10-27)

Peer request (liggghts-public-1e, at the owner's request): the −2.0 % G5 deficit is attributed to
the discrete continuity source, `el_volume_averaged.f90:77-107`: dividing the backward difference
of ε by ε^{n+1} leaves a positive O(dt) term Q2 = (dt/2) Σ|e| (∂tε)²/ε, so int_g is persistently
positive, and the uniform removal of ḡ places the compensating sink at z = 0.5 while the spurious
source sits in the cloud (z ≈ 0.74). CFDEM's +1.4 % was its cell-centre post-processing; its
conservative face flux gives 1.00002, so the twin ratio is 0.980/1.000. Tests, all clones of
`runs/vft04/VB_L2` (CFDEM positions) changing only what is listed, `versions.txt` in each run dir:

| Run (dir under `runs/vft04/`) | change | binary | job, wall |
|---|---|---|---|
| T1 `T1_VB_L3` | `MaxMeshLevel 3` | frozen copy of build-va (2026-10-15, HEAD d89f2e64 code) | 155110, 54 min, 11 ranks |
| T2a `T2a_VB_L2_dt5e-4` | dt 5e-4, 400 steps (deck + pe) | frozen | 155111, 16 min |
| T2b `T2b_VB_L2_dt2.5e-4` | dt 2.5e-4, 800 steps | frozen | 155112, 32 min |
| T2c `T2c_VB_L2_dt2e-3` | dt 2e-3, 100 steps | frozen | 155113, 4 min |
| T4 `T4_VB_L2_diag` | none (diagnostics commit 7db4650e) | build-va rebuilt 2026-10-27 20:05 | 155114, 7.5 min |

G5 = uf_super_z / (N V_p |w| / V_box) = uf_super_z / 4.24115e-4 from `EL_MEAN_SLIP` (log, not VTK),
mean ± sd over t ≥ 0.002; int_g and L1_dedt are means of `EL_VA_SOURCE` (`vft_g5.py` here).

| Run | h, dt | G5 mean ± sd | G5 t ≥ 0.1 | int_g | L1_dedt | prediction (SPEC 12.4) |
|---|---|---|---|---|---|---|
| base VB_L2 | L2, 1e-3 | 0.9801 ± 0.0088 | 0.9812 ± 0.0045 | 4.48e-7 | 1.132e-3 | |
| T2c | L2, 2e-3 | 0.9590 ± 0.0192 | 0.960 | 8.19e-7 | 1.130e-3 | 0.94-0.955 |
| T2a | L2, 5e-4 | 0.9917 ± 0.0050 | 0.992 | 2.43e-7 | 1.132e-3 | 0.99-1.00, int_g ≤ 2.3e-7 |
| T2b | L2, 2.5e-4 | 0.9979 ± 0.0049 | 0.998 | 1.34e-7 | 1.132e-3 | 1.00-1.01 |
| **T1** | **L3, 1e-3** | **0.9926 ± 0.0279** | **0.9912 ± 0.0004** | **1.68e-7** | 1.332e-3 | **0.95-0.96, int_g 0.8-1.1e-6** |

- **T2 confirms the O(dt) mechanism**: the deficit and int_g halve with each halving of dt
  (1 − G5 = 4.1 / 2.0 / 0.8 / 0.2 % and int_g 8.2 / 4.5 / 2.4 / 1.3e-7 at dt 2e-3 / 1e-3 / 5e-4 /
  2.5e-4). T2a and T2b fall inside their predicted bands; T2c is 0.4 points above its band.
- **T1 falsifies the refinement prediction**: at level 3 and the same dt the deficit does not
  double, it shrinks to −0.9 % and int_g drops by 2.7× (4.5e-7 → 1.7e-7); the L1_dedt prediction
  (1.33e-3) is met. The level-3 series is much smoother (sd 0.0004 for t ≥ 0.1 against 0.0045
  at level 2) after a large start-up transient (G5 1.38 in the first steps, int_g −6.9e-7 for
  t < 0.01). The SPEC 12.3 statement "growing with refinement at a fixed dt" is wrong in sign
  (corrected there). Peer's reading of the h-dependence (2026-10-27, from the VTK budget of
  `VA_L3`, same ε field): the O(dt) term itself RISES at level 3 (Q2 1.32e-6, I0 1.355e-6 against
  7.28e-7 / 7.42e-7 at level 2), so adv_net ≈ 1.19e-6 at level 3 cancels about 88 % of it
  against 39 % at level 2: the residual net source is the O(dt) term times an imperfect discrete
  cancellation by the advective term, and the cancellation improves with h. A T4-diagnostics run
  at level 3 (`T4_VB_L3_diag`, staged, not submitted) would measure adv_net directly.
- **T4 confirms the budget identities** (new `EL_VA_MOMENT` line and `uz_unw=` field, print-only,
  gated bitwise against `tier2_cases/momentum_conservation/reference/run_g0_twin_va_2026-10-14.log`
  by job 155115, `tools/el_stream_compare.py`: all 13 EL streams identical, EL_MEAN_SLIP with one
  extra trailing token): ∫u_z dV = −M1raw + ½ int_g to rel rms 2.0e-5 (max residual 1.1e-9) over
  199 steps (predicted ≈ 1e-3); ∫u_z / MOM 0.9867 (predicted ≈ 0.987); Q2 7.28e-7 (7.3e-7),
  adv_net +2.94e-7 (+3.0e-7), int_g + adv_net 7.42e-7 (7.5e-7), S2 2.92e-7; G5 and int_g equal
  to the base run to four digits.
- T5 (fix A: α-weighted removal of the net source; fix B: divide by (ε^{n+1}+ε^n)/2) is a physics
  change and has NOT been applied; it waits for the owner's decision.

## 7. Fix B test (T6) and the pe-rebased binary on VB (jobs 155127-155133, 2026-10-27)

Peer request (SPEC 04 §12.4 T5/T6, owner's approval relayed): test fix B as a deck switch, not
as a change to the main EL branch. Implementation on the TEST branch `test/el-va-source-eps-mid`
(worktree `FeatFloWer-va-fixb`, build `build-fixb`, commit b97b60ce on top of 4c52af80):
`SimPar@ELVASourceEpsLevel = new` (default, divide by ε^{n+1}) `| mid` (divide by
(ε^{n+1}+ε^n)/2 with ε^n = ε^{n+1} − dt·∂tε reconstructed from the backward difference); the
same ε enters the source, the L1 norms and the `EL_VA_MOMENT` terms. Fix A is NOT implemented;
fix B is NOT on `feature/el-volume-averaged`.

Gates (switch at default): tier-2 twin job 155128 bitwise vs `reference/run_g0_twin_va_2026-10-14.log`
(13 streams); VB_L2 base job 155129 bitwise vs the T4 log (14 streams incl. `EL_VA_MOMENT`).
The pe-rebased main binary (4c52af80, pin 4d04c1e) also reproduces the T4 log bitwise on VB_L2
(job 155127, `G_VB_L2_pe4d04c1e`), so the pe bump is certified on the volume-averaged case too.

| Run | dt | G5 mean ± sd (t ≥ 0.002) | G5 t ≥ 0.1 | int_g | adv_net | Q2 | M1raw | S2 | identity rel rms | ∫u_z/MOM |
|---|---|---|---|---|---|---|---|---|---|---|
| base (new) | 1e-3 | 0.9801 ± 0.0088 | 0.9812 | +4.48e-7 | +2.94e-7 | 7.28e-7 | −3.961e-6 | 2.92e-7 | 2.0e-5 | 0.9867 |
| T6a (mid) `T6a_VB_L2_epsmid` 155130 | 1e-3 | **1.0213 ± 0.0090** | 1.0199 | **−2.88e-7** | +2.88e-7 | 7.28e-7 | −4.504e-6 | 2.92e-7 | 1.9e-5 | 1.0281 |
| T6b (mid) `T6b_VB_L2_epsmid_dt5e-4` 155133 | 5e-4 | **1.0123 ± 0.0050** | 1.0119 | **−1.26e-7** | +1.26e-7 | 3.64e-7 | −4.385e-6 | 2.92e-7 | 2.8e-6 | 1.0190 |

L1_dedt 1.132e-3 in all three. Reading:
- Fix B removes the temporal O(dt) term exactly: with the centred divisor Σ|e| ∂tε/ε̄ integrates
  to int_g + adv_net = −2.5e-10 (T6a) / +2.7e-10 (T6b), i.e. the requested "exact temporal
  source integral" is −(int_g + adv_net) and is zero to 1e-10 (log-difference telescoping). So
  under fix B int_g = −adv_net.
- The advective term keeps its own O(dt) part (Codex's (dt/2)∫d² reading): adv_net 2.88e-7 →
  1.26e-7 when dt halves, and G5 − 1 = +2.1 % → +1.2 %. Bare fix B OVERSHOOTS with a first-order
  dt dependence of the opposite sign: Codex's forecast (G5 1.017-1.023, int_g −(2.5-3.2)e-7)
  is met; reading (i) (G5 ≈ 1.004, int_g ≈ 0) is not.
- Richardson limits agree: original series 2·0.9917 − 0.9801 = 1.0033; fix-B series
  2·1.0123 − 1.0213 = 1.0033. The dt → 0 value of G5 is ≈ 1.003 either way; fix B alone moves
  the first-order error from −2.0 % to +2.1 % at dt 1e-3 and does not remove it.

### 7.1 T7: T4 diagnostics at level 3 (`T4_VB_L3_diag`, job 155126, 54 min, 11 ranks)

Main binary 4c52af80 (pe 4d04c1e, diagnostics 7db4650e), deck = T1 (`MaxMeshLevel 3`, dt 1e-3).
All streams common with T1 (frozen 2026-10-15 binary) are bitwise; only `EL_VA_MOMENT` is new.

| quantity (mean, t ≥ 0.002) | T7 level 3 | prediction (peer) | level 2 (T4) |
|---|---|---|---|
| G5 | 0.9926 ± 0.0279 (0.9912 ± 0.0004 for t ≥ 0.1) | as T1 | 0.9801 |
| int_g | 1.678e-7 | 1.68e-7 | 4.48e-7 |
| Q2 | 1.321e-6 | ≈ 1.32e-6 | 7.28e-7 |
| adv_net | +1.186e-6 (1.155e-6 for t ≥ 0.1) | ≈ +1.19e-6 | +2.94e-7 |
| I0 = int_g + adv_net | 1.354e-6 | ≈ 1.355e-6 | 7.42e-7 |
| cancellation adv_net/I0 | 87.6 % | ≈ 88 % (Codex 85.2 %) | 39.6 % |
| M1raw, S2 | −4.170e-6, 4.05e-7 | | −3.961e-6, 2.92e-7 |
| identity rel rms | 1.8e-5 (max 9.1e-10) | ≈ 2e-5 | 2.0e-5 |
| ∫u_z / MOM | 1.0031 | | 0.9867 |

The peer's reading of T1 is confirmed and the earlier guess in §6 (Q2 falls at level 3) is
wrong: the O(dt) term nearly doubles with refinement at fixed dt, and the advective term cancels
88 % of it at level 3 against 40 % at level 2. The residual net source is the imperfectly
cancelled O(dt) pair, whose cancellation improves with h.

## 8. T8: endpoint BDF2 continuity source with source/projection correctors (jobs 155139-155145, 2026-10-28)

Peer request (SPEC 04 §12.7, Codex rank-1 fix, owner's choice relayed). Implemented on the TEST
branch `test/el-va-source-eps-mid` (worktree `FeatFloWer-va-fixb`, commits 4c0dd4f3 + 9fb6157c):
`SimPar@ELVASourceTimeScheme = BE | BDF2` (d_cont = 1.5 ∂tε^{n+1} − 0.5 ∂tε^n, ε and ∇ε at n+1, BE
on the first step; previous rate saved once per step at the history advance, restart v3) and
`SimPar@ELVASourceCorrectors = n` (after the projection: reassemble (1/k)Bᵀu with the corrected
velocity, recompute g from it, log the constraint residual, solve with the same pressure solver,
correct u and accumulate p; P^n kept for the P_new extrapolation; no second particle advance or
momentum predictor). New log lines `EL_VA_TSRC` (signed temporal integral Σ|e| d_cont/ε, its L1,
scheme) and `EL_VA_RESIDUAL` (dt·defp on the constant dof: L1, max, int_g; corr 1..n before each
corrector solve, corr −1 after the last correction). Defaults byte-identical: tier-2 twin job
155139 and VB_L2 job 155140 bitwise vs the reference / the T4 log. Four-step local smoke with two
correctors: residual L1 5.2e-7 → 3.0e-8 → 2.0e-10 (×20 per pass).

Means over t ≥ 0.002 (source terms from the last solve pass); G5 for t ≥ 0.1 in brackets:

| Run | h, dt, corr | G5 ± sd | int_g | adv_net | dedt_net | Q2 | residual L1 before corr / final | ∫u_z/MOM | job, wall |
|---|---|---|---|---|---|---|---|---|---|
| base BE | L2, 1e-3, 0 | 0.9801 ± 0.0088 (0.9812) | +4.48e-7 | +2.94e-7 | −7.42e-7 | 7.28e-7 | – | 0.9867 | 155114 |
| T8a BDF2 | L2, 1e-3, 0 | **1.0046 ± 0.0100 (1.0042)** | +1.6e-8 | −1.0e-8 | −6.2e-9 | 7.37e-7 | – | 1.0113 | 155141, 8 min |
| T8b BDF2 | L2, 1e-3, 1 | **1.0039 ± 0.0050 (1.0039)** | +2.8e-8 | −2.2e-8 | −6.2e-9 | 7.37e-7 | 9.1e-8 / 1.8e-9 | 1.0106 | 155142, 10 min |
| T8c BDF2 | L2, 5e-4, 1 | 1.0041 ± 0.0050 (1.0041) | +2.4e-8 | −2.3e-8 | −1.1e-9 | 3.66e-7 | 1.3e-8 / 1.2e-10 | 1.0108 | 155143, 21 min |
| T8e BDF2 | L2, 2.5e-4, 1 | 1.0041 ± 0.0049 (1.0041) | +2.4e-8 | −2.3e-8 | −6.1e-10 | 1.82e-7 | 3.2e-9 / 2.8e-11 | 1.0108 | 155145, 39 min |
| T8d BDF2 | L3, 1e-3, 1 | 0.9926 ± 0.0019 (0.9930) | +1.70e-7 | −1.57e-7 | −1.3e-8 | 1.34e-6 | 3.5e-7 / 1.2e-8 | 1.0032 | 155144, 72 min |

L1_dedt unchanged (1.132e-3 at L2, 1.332e-3 at L3); L1_dcont 1.137e-3 / 1.339e-3. Moment identity
(final assembly) rel rms 6e-5 (T8b), 3e-7 (T8c), 9e-8 (T8e), 7e-4 (T8d).

Reading:
- **BDF2 alone removes the first-order error at level 2**: T8a gives 1.0046 and the net source
  drops from 4.5e-7 to 1.6e-8 (int_g) with the temporal integral −6e-9. The time-level match is
  the whole story at L2; the corrector (T8b) only tightens the constraint (residual ×50) and
  moves G5 by −0.0007.
- **The dt dependence is gone**: T8b/T8c/T8e = 1.0039 / 1.0041 / 1.0041 (prediction: within
  0.002, O(dt²)); the limit 1.0041 agrees with the Richardson limit 1.0033 of the BE and fix-B
  series to the scatter of those extrapolations. Codex's T8b prediction (1.004-1.005, int_g a few
  1e-8) is met.
- **Level 3 is not improved**: T8d gives 0.9926, identical to BE at L3 (T1/T7 0.9926), with
  int_g 1.7e-7 again. Under BDF2 the temporal integral is ~0 and the whole net source is the
  advective term (adv_net −1.57e-7, which was +1.19e-6 under BE). So the L3 deficit of −0.7 % is an
  h-dependent bias of the discrete advective term u·∇ε, not an O(dt) effect; it sits "below 1.01"
  as predicted but it is not a dt → 0 limit of the L2 kind. Open: an L3 dt series or an L4 run to
  tell whether it converges with h.
- Corrector residual history (L2, dt 1e-3): 9.1e-8 before the corrector, 1.8e-9 after it, i.e.
  the predictor already satisfies the constraint to 2 % of int_g; one corrector is enough.

### 8.1 G4 slab profiles (base, T8b, T8c; `vft_slab_profile.py`, full tables in `slab_profiles_T8.txt`)

From the VTK pieces: `ELEpsilonF`/`ELAlphaP` are the element-wise constant CELL fields the source
uses, `Velocity` is Q1 vertex data; per cell ε_e · ⟨u_z⟩_8vertices, volume weighted, ten z-slabs
of 0.1 m. The box integral of this estimate is 0.92-0.95 of MOM (Q1 vertices miss the Q2 edge/face/
centre dofs; the log `uf_super` is the reference), so only the slab-to-slab and run-to-run
differences are meaningful. VTK index 00099 = t 0.100, the last index (00190 / 00380) = t 0.191.

Ratio ⟨ε u_z⟩ / (|w| ⟨α⟩) per slab (slabs 1-4 hold no particles; their ⟨ε u_z⟩ is given in m/s):

| slab (z) | base t 0.10 | T8b t 0.10 | T8c t 0.10 | base t 0.19 | T8b t 0.19 | T8c t 0.19 |
|---|---|---|---|---|---|---|
| 1-4 (0-0.4), ⟨ε u_z⟩ | −1.5e-6 … −1.0e-5 | +0.6e-6 … +4.0e-6 | +0.6e-6 … +4.1e-6 | −2.8e-6 … −1.9e-5 | −0.8e-6 … −5.7e-6 | −0.7e-6 … −5.2e-6 |
| 5 (0.4-0.5) | −0.883 (α 1.8e-5) | 2.58 (α 1.8e-5) | 2.61 | 0.854 | 0.922 | 0.924 |
| 6 (0.5-0.6) | 0.918 | 0.940 | 0.940 | 0.914 | 0.933 | 0.933 |
| 7 (0.6-0.7) | 0.927 | 0.942 | 0.942 | 0.932 | 0.945 | 0.946 |
| 8 (0.7-0.8) | 0.926 | 0.937 | 0.937 | 0.931 | 0.939 | 0.938 |
| 9 (0.8-0.9) | 0.939 | 0.945 | 0.945 | 0.942 | 0.944 | 0.944 |
| 10 (0.9-1.0) | 0.944 | 0.943 | 0.942 | 0.959 | 0.952 | 0.952 |
| box integral / MOM | 0.9225 | 0.9456 | 0.9456 | 0.9174 | 0.9371 | 0.9376 |

- T8b lifts the ratio in the cloud slabs 6-9 by +0.01 to +0.02 almost uniformly, which is the
  +2.4 points of G5; the top slab 10 is unchanged or slightly lower. T8c equals T8b to 1e-4 in
  every slab (dt-independence slab by slab).
- The particle-free slabs 1-4 carry a linear-in-z return flow: in the base it is negative
  (−1e-5 m/s at z 0.35, the sink of the uniform ḡ removal sitting at z = 0.5 below the cloud);
  under BDF2 it is ten times smaller and of the opposite sign at t 0.1.

## 9. BDF2 source + projection correctors ported to the VA branch (commits 14d478a5 + cafa91bc, jobs 155235-155246, 2026-10-28)

Owner decision: move `ELVASourceTimeScheme = BDF2` and `ELVASourceCorrectors = n` from the test
branch to `feature/el-volume-averaged`, deck-switchable, modular (`docs/md_docs/el_volume_averaged_stage.md`,
last section, lists the keys and the design). `ELVASourceEpsLevel` (fix B) not ported. Corrector loop
isolated in `source/src_quadLS/QuadSc_el_va.f90` (one guarded call in the fluid core, guard
`EL_VA_SOURCE_CORRECTORS_ENABLED()` on the deck keys only); BDF2 rate in `EL_VA_CONTINUITY_RATE`
with the variable-step weight a = dt/(dt+dt_prev) (= 0.5 exactly for equal steps).

Reviews (read-only, both on 14d478a5): Codex gpt-6.1-sol high and Claude Fable 5.1 xhigh
(`docs/md_docs/el_va_bdf2_port_reviews_2026-10-28.md`). Agreed findings, all fixed in cafa91bc:
(1) with `MatrixRenewal` C ≥ 2 and a direct coarse pressure solver (types 2/4/9) `MG_Solver` frees the
coarse factorization after every solve, so a corrector solve would use a freed handle — now rejected
with a message (every EL deck is C1 + type 4); (2) the source Allreduce had grown from 8 to 10 entries
at the defaults (identity only empirical) — the 8-entry campaign reduction is kept, the two `EL_VA_TSRC`
sums are reduced separately on audit steps; (3) BDF2 assumed equal steps and `rate_valid` was never
cleared — variable-step weight, flags cleared on non-genuine rates, both rates required; (4) rank-local
return ahead of the residual Allreduces removed; restart v3 also stores `history_valid` and the two
intervals; `el_stream_compare.py` now compares `EL_VA_TSRC`/`EL_VA_RESIDUAL` and reports streams absent
from the reference as "new, not compared". Both reviewers verified by reading: step ordering of the
rate history, MPI participation, restored `valP_old`/`P_new`, no reader of `defP`/`rhsP` after the fluid
core, deck-only switch with order-independent validation.

Gates, three rounds (14d478a5 binary; review-fix binary; final cafa91bc binary), all with
`tools/el_stream_compare.py`:

| Round | tier-2 twin vs `reference/run_g0_twin_va_2026-10-14.log` (13 streams) | VB_L2 defaults vs T4 log (15 streams) | T8b keys on vs test-branch T8b log |
|---|---|---|---|
| 1 (14d478a5) | 155235 bitwise | 155236 bitwise | 155237 bitwise (15 streams, residual lines identical) |
| 2 (review fixes) | 155241 bitwise | 155242 bitwise | 155243 bitwise (17 streams incl. TSRC/RESIDUAL) |
| 3 (cafa91bc) | 155244 bitwise | 155245 bitwise | 155246 bitwise (17 streams), G5 1.0039 ± 0.0050 |

Wall: VB_L2 defaults 7 min 46, T8b 9 min 51 (11 ranks). The ported scheme therefore reproduces
§8 exactly (T8b: G5 1.0039, int_g 2.83e-8, residual L1 9.06e-8 → 1.77e-9) and the defaults reproduce
the campaign binary bitwise. Unit test `test_el_kernel_forces` extended (restart v3 round trip, history
capture of `EL_BEGIN_FIELD_UPDATE`), passes. Sample deck `q2p1_param_VB_bdf2.dat` (= VB deck + the two keys).

## 10. T9: fix A (alpha-weighted compatibility removal) on top of BDF2 (jobs 155259-155264, 2026-10-28)

Peer request (SPEC 04 §12.8; owner's choice relayed). TEST branch `test/el-va-compat-alpha` (worktree
`FeatFloWer-va-fixb`, commit 8aeb26c8 on top of the reviewed VA head 19e3349c; not for merge):
`SimPar@ELVACompatibility = uniform | alpha`. In mode alpha the raw source integral is removed as
g_e ← g_e − int_g α_p(e)/W, W = Σ α_p|e| (uniform fallback for W < 1e-30), in every pass; new audit line
`EL_VA_COMPAT` (int_g_raw, int_g_corr, W, z_α = Σ z_c α|e|/W, mshift = (z_α − z_centre) int_g, mode). The
reduction for W/z_α and the line exist only in mode alpha. Gates at the default: tier-2 twin job 155259
and VB_L2 job 155260 bitwise vs the certified logs (on the epyc nodes summergames/wintergames, same
result as on nx). Runs on epyc nodes (wall: T9a 16 min, T9c 16 min, T9d 33 min, T9b 1 h 54).

Means over t ≥ 0.002 (last solve pass); G5 for t ≥ 0.1 in brackets; `vft_g5.py` now reads EL_VA_COMPAT
and closes the moment identity with the sink at z_α (uniform: box centre):

| Run | Settings | G5 ± sd | predicted | int_g raw | int_g corr (max) | adv_net | dedt_net | z_α | mshift/MOM | res. L1 corr 1 / −1 | ∫u_z/MOM | identity rel rms | job |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T8b (ref) | BDF2+c1 uniform, L2 | 1.0039 ± 0.0050 (1.0039) | – | 2.83e-8 | – | −2.2e-8 | −6.2e-9 | – | – | 9.1e-8 / 1.8e-9 | 1.0106 | 6.1e-5 | 155142 |
| T9a | BDF2+c1 alpha, L2, dt 1e-3 | **1.0055 ± 0.0018 (1.0054)** | 1.0055 | 2.85e-8 | 1.9e-21 (1.4e-19) | −2.2e-8 | −6.2e-9 | 0.7388 | 0.0016 | 7.9e-8 / 5.3e-10 | 1.0122 | 4.0e-6 | 155261 |
| T9c | BE alpha, L2, dt 1e-3 | **1.0053 ± 0.0018 (1.0052)** | 1.005 (1.0053) | 4.54e-7 | −1.0e-21 (1.3e-19) | +2.9e-7 | −7.4e-7 | 0.7388 | 0.0255 | – | 1.0120 | 2.8e-6 | 155262 |
| T9d | BDF2+c1 alpha, L2, dt 5e-4 | 1.0054 ± 0.0018 (1.0054) | T9a ± 1e-3 | 2.45e-8 | −1.4e-21 | −2.3e-8 | −1.1e-9 | 0.7388 | 0.0014 | 1.2e-8 / 9.9e-11 | 1.0121 | 9.3e-8 | 155263 |
| T8d (ref) | BDF2+c1 uniform, L3 | 0.9926 ± 0.0019 (0.9930) | – | 1.70e-7 | – | −1.57e-7 | −1.3e-8 | – | – | 3.5e-7 / 1.2e-8 | 1.0032 | 7e-4 | 155144 |
| T9b | BDF2+c1 alpha, L3, dt 1e-3 | **1.0021 ± 0.0001 (1.0021)** | ≈1.002 (1.0018) | 1.73e-7 | −1.6e-21 (2.6e-19) | −1.59e-7 | −1.3e-8 | 0.7388 | 0.0097 | 2.5e-7 / 1.7e-9 | 1.0127 | 9.1e-6 | 155264 |

Reading:
- Every prediction is met to the printed digits: T9a +0.0016 over T8b (= mshift/MOM 0.0016), T9b +0.0095
  over T8d (mshift/MOM 0.0097, predicted 0.0096), T9c 1.0053, T9d within 1e-4 of T9a. The raw int_g is
  unchanged by the removal mode (2.85e-8 / 1.73e-7), as predicted; the corrected integral is at 1e-21.
- The L3 deficit is closed: 0.9926 → 1.0021 (+0.2 %, the sd of 1e-4 is the tightest of the campaign).
  What remains at L3 (+0.21 %) and L2 (+0.55 %) is of the opposite sign and of the size of the two-way
  co-flow bias seen in 01_settling.
- The moment identity with the sink at z_α closes to 4e-6 (T9a), 9e-6 (T9b), 9e-8 (T9d) relative, the
  cleanest confirmation that the uniform-removal moment (0.5 − z_cloud) int_g was the L3 lever.
- fix A alone (T9c, BE) gives the same 1.005 as BDF2 alone (T8a 1.0046) at L2: the two defects are
  additive and of similar size at L2 (O(dt) temporal, uniform-removal moment); at L3 the second dominates.

Slab profiles (`slab_profiles_T9.txt`, same script and frames as §8.1; ratio ε u_z/(|w| α) per slab):

| slab (z) | T8d t 0.10 | T9b t 0.10 | T8b t 0.10 | T9a t 0.10 | T9c t 0.10 | T8d t 0.19 | T9b t 0.19 | T8b t 0.19 | T9a t 0.19 |
|---|---|---|---|---|---|---|---|---|---|
| 1-4 (0-0.4), ⟨ε u_z⟩ m/s | −0.8e-6 … −5.9e-6 (linear in z) | ≤ 1e-10 | ~+1e-6 | ≤ 4e-11 | ≤ 4e-11 | −0.8e-6 … −5.9e-6 | ≤ 6e-10 | – | ≤ 3e-10 |
| 5 (0.4-0.5) | −0.324 (α 1.8e-5) | 1.070 | 2.58 | 1.63 | 1.63 | 0.945 | 0.973 | 0.922 | 0.949 |
| 6 (0.5-0.6) | 0.962 | 0.971 | 0.940 | 0.934 | 0.934 | 0.966 | 0.973 | 0.933 | 0.940 |
| 7 (0.6-0.7) | 0.966 | 0.973 | 0.942 | 0.938 | 0.938 | 0.972 | 0.978 | 0.945 | 0.951 |
| 8 (0.7-0.8) | 0.967 | 0.972 | 0.937 | 0.934 | 0.934 | 0.972 | 0.975 | 0.939 | 0.942 |
| 9 (0.8-0.9) | 0.971 | 0.973 | 0.945 | 0.944 | 0.943 | 0.974 | 0.975 | 0.944 | 0.945 |
| 10 (0.9-1.0) | 0.973 | 0.973 | 0.943 | 0.943 | 0.942 | 0.983 | 0.982 | 0.952 | 0.951 |
| box integral / MOM | 0.9628 | 0.9724 | 0.9456 | 0.9391 | 0.9391 | 0.9673 | 0.9757 | 0.9371 | 0.9453 |

- Below the cloud (slabs 1-4) the spurious flow of the uniform removal (T8d: linear, −0.8e-6 to −5.9e-6
  m/s, predicted order 1e-5) vanishes under fix A at both levels (≤ 1e-9 m/s, predicted ≤ 1e-6).
- At L3 fix A lifts the cloud slabs by +0.001 to +0.009 and the box integral by +0.0096 (= mshift/MOM),
  consistent with G5. At L2 the cloud-slab ratios move by −0.001 to −0.006 and the box integral by
  −0.0065 while G5 rises by +0.0016: the slab quadrature (cell ε × vertex-mean u_z) and the code's own
  ∫ε u_z disagree on the sign of a 0.5 % change once the sink sits in the particle cells (sub-cell
  velocity structure at L2). The slab tables are therefore read for the below-cloud flow and for the
  L3 trend only; G5 is the measure. (T9c = T9a slab by slab to 1e-3.)

Not merged. Decision for the owner: fix A is a one-key change on top of the ported BDF2 scheme
(`ELVACompatibility = alpha`), default-gated bitwise; it closes the L3 deficit and leaves +0.2 % (L3) /
+0.55 % (L2).

Provenance note (peer check, 2026-10-28 14:20): the run directories of §9 (G_port*, P_port*) and §10
(T9a-d, G_t9_*) were cloned with rsync from the T8 directories and had kept the clones' `versions.txt`
(9fb6157c, build of 2026-10-27 23:00). The logs prove the right binaries ran (EL_VA_TSRC / EL_VA_COMPAT
lines, build mtimes); the files were rewritten by hand: round 1 14d478a5 (build-va 08:24), round 2
14d478a5 + uncommitted review fixes (08:38), round 3 cafa91bc (08:42:48), T9 and G_t9 8aeb26c8
(build-fixb 12:11:18), pe 4d04c1e throughout. From now on `tools/el_slurm_submit_np.sh` regenerates
`versions.txt` at every submission from the git state of the binary's source tree
(`tools/el_write_versions.sh`), so a cloned directory can no longer carry stale provenance.

## 11. fix A ported to the VA branch with a restructured source module (commits cef8eeaa + f93eb55e, jobs 155644-155653, 2026-10-30)

Owner decision: move fix A to `feature/el-volume-averaged` with the same review and gate procedure as
the BDF2 port (§9), with the explicit requirement of clean, maintainable, extensible code. Deck key
`SimPar@ELVACompatibility = uniform | alpha` (default uniform = the campaign code). The port is also a
restructuring of `source/src_el/el_volume_averaged.f90`: named slots for the three reduction blocks
(`C_*` core = the campaign 8-entry Allreduce, `T_*` temporal on audit steps, `W_*` weight block for
weighted sink modes), a sink-mode table (`el_va_compat_names` in `EL_CONFIG`, deck value = audit label,
validated against the same table), `EL_VA_SINK_WEIGHT` (one CASE per weighted mode), `EL_VA_SINK`
(uniform: gbar; weighted: int_g w_e/W, expressions verbatim from the validated runs),
`EL_VA_EFFECTIVE_SINK_MODE` (requested vs effective mode, uniform fallback without particles decided on
the reduced W), `EL_VA_WRITE_LINE` for every audit line, four labelled sections, no collective in the
audit section. A further weighted sink mode is one table entry and one CASE.

Reviews (read-only, both on cef8eeaa; `docs/md_docs/el_va_fixa_port_reviews_2026-10-30.md`): Codex
gpt-6.1-sol high and Claude Fable 5.1 xhigh. Both: arithmetic and collective pattern preserved, deck-only
switch, materially more maintainable. Fixed in f93eb55e: requested/effective mode as two variables (both),
table-driven modes so a third mode is local (Fable S3), the corrected-integral reduce hoisted out of the
audit section (Fable S1), zero moment shift under the no-particle fallback (Codex 2), header paragraphs
and extension-point docs made accurate (both), `EL_VA_COMPAT` added to the compare tool (Fable S4),
helper renamed `EL_VA_WRITE_LINE`, PURE helpers. Not changed: the port writes `EL_VA_COMPAT` after
`EL_VA_TSRC` (the test branch wrote it first; Codex called this a blocker against a byte-level log diff,
the gate criterion is per-stream identity, recorded here so the reorder is not read as a numeric change).

Gates, two rounds (cef8eeaa binary; f93eb55e binary), `tools/el_stream_compare.py` plus a byte comparison
of all `EL_VA_*` lines against the §9 round-3 logs where the keys match:

| Round | tier-2 twin, defaults (13 streams) | VB_L2 defaults vs T4 log (15) | T8b: BDF2+c1 uniform vs T8b log (17) | T9a: BDF2+c1 alpha vs T9a log (18 incl. COMPAT) | T9b: L3 alpha vs T9b log |
|---|---|---|---|---|---|
| 1 (cef8eeaa) | 155644 bitwise | 155645 bitwise, EL_VA_* bytes = G_port3 | 155646 bitwise | 155647 bitwise, G5 1.0055 | 155648 bitwise, G5 1.0021 (summergames, 1 h 55) |
| 2 (f93eb55e) | 155649 bitwise | 155650 bitwise, EL_VA_* bytes = G_port3 | 155651 bitwise, EL_VA_* bytes = P_port3_T8b | 155652 bitwise, G5 1.0055 | 155653 bitwise, G5 1.0021 (nx-05, 1 h 10) |

Wall on nx-04 (round 2): VB_L2 8 min, T8b 10 min, T9a 10 min; T9b 1 h 10 on nx-05 (1 h 55 on summergames in round 1). The level-3 alpha path (18 streams incl. EL_VA_COMPAT) is bitwise in both rounds and on both node types. The `EL_VA_COMPAT` stream (600 lines,
3 passes × 200 steps) is bitwise down to the 1e-21 `int_g_corr` field. Provenance: `versions.txt`
regenerated at submission (§10 note).

### 11.1 pe pin 4d04c1e → 11c6142: box walls rebased onto pe master f17116f (PR #54 interface layout), jobs 155657-155659, 2026-10-30

pe master moved by one merge since the pin (PR #54, `interface: one case per file, ABI headers apart,
helpers under support/`, a mechanical reorganisation; FeatFloWer binds to the C symbols only). The
box-walls commit cherry-picks onto it without conflict (git follows the file move to
`pe/interface/setup/el_terminal_velocity.h`): branch `feature/el-box-walls-on-master-f17116f` =
c720761 (the hunk) + 11c6142 (json contract test `pe-box-walls-config`, header comment). pe preset
`serial`: 45/45 tests pass. FeatFloWer at this pin (build-va, cmake regenerated for the new pe source
list - a plain ninja after the pin bump fails on the removed files): tier-2 twin job 155657 bitwise,
VB_L2 defaults job 155658 bitwise vs the T4 log with `EL_VA_*` bytes identical to §9 round 3 (box walls
on, "six fixed planes" line present), T9a alpha keys job 155659 bitwise vs the T9a log (G5 1.0055).
The pe branch is to be pushed by the owner and merged as a PR; the pin then moves to the merge commit.

### 11.2 pe pin 11c6142 → df8c8a6 (pe master, PR #55 merged 2026-11-01), jobs 156017-156019

PR #55 (CI: cgal, clang-serial, gcc-serial, mpi all pass) merged as df8c8a6; its tree is identical to
11c6142, so the regenerated build left the binary untouched (same inode, mtime 2026-10-30 21:46:25).
Gates for the record under the new pin: tier-2 twin 156017 bitwise, VB_L2 defaults 156018 bitwise vs the
T4 log with `EL_VA_*` bytes identical to §9 round 3, T9a alpha keys 156019 bitwise vs the T9a log
(G5 1.0055). The VA worktree now follows pe master; the EL box walls are upstream.

## 12. T10: conservative continuity form (form C) on the VA head (branch feature/el-va-conservative-form, commits dc2221e0 + 9634485e, jobs 157063-157113, 2026-11-05)

Design, derivation and code: `docs/md_docs/el_va_conservative_form_2026-11-05.md` (GPT-6 Astra design
consultation before coding, Codex + Fable reviews after; `el_va_conservative_form_*_2026-11-05.md`).
One deck key more than T9, `SimPar@ELVAContinuityForm = conservative` (sample deck
`q2p1_param_VB_bdf2_alpha_cons.dat`): the constraint div(εu) = −∂ε/∂t is enforced in the weak form of
all four P1 test functions with shared face fluxes of the Q1-reconstructed transport porosity ε_h, iterated
through the existing corrector passes. The audit lines `EL_VA_FLUX` (storage rate S, boundary flux Φ,
ε-weighted sink integral, unsinked residuals) and `EL_VA_FLUXMOM` (coordinate identities, ε_h-weighted and
P0-weighted superficial velocity, storage first-moment rate m1rate, dvol gap) are new; `vft_g5.py` reads
the campaign G5 as before (`EL_MEAN_SLIP`, P0 ε), the form-C numbers below come from the new lines
(`vft_formc.py` in this folder, means over t ≥ 0.002).

Smoke (VB_L2, 20 steps, 11 ranks, S10_* run dirs): the face flux equals the divergence of B to 1e-17
(`divchk`), Φ and the three identities are 1e-20 on the closed box, the unsinked row-1 residual falls
1.1e-3 → 8.7e-6 → 1.5e-7 over the predictor and one corrector (3.8e-9 with two, 3.8e-12 with four) and the
ε-weighted sink integral vanishes (7e-10 → 2e-11 → 3e-14). Two bugs found by the smoke and fixed before
the ladder: FEAT accepts left-handed hexes (only |det J| is used), so the outward normal needs sgn(det J);
and a SIGN-intrinsic misuse. Both are covered by the unit test (`test_el_va_flux`, 21 checks).

Ladder (clones of the T9 decks with the key added; `versions.txt` in each run dir; dc2221e0 binary for the
runs started 12:20, 9634485e for T10a2 and round 2):

| Run | change vs T9a | G5 ε_h-weighted (`uf_eps_z`/MOM), t ≥ 0.002 | same, settled t ≥ 0.1 | G5 campaign (`EL_MEAN_SLIP`, P0 ε) | m1rate/(N V_p\|w\|), t ≥ 0.1 | per-step identity `uf_eps_z` V/m1rate | res1_L1 max | \|sink_eps\| max | wall |
|---|---|---|---|---|---|---|---|---|---|
| T10a `T10a_VB_L2_cons_c1` (157089) | form C | 0.99994 ± 0.00033 | 1.00000 ± 0.00016 | 0.9987 ± 0.0003 | 1.000000 | 1 + 2.3e-6 ± 7.9e-6 | 1.5e-7 | 7e-10 | 18 min |
| T10c `T10c_VB_L2_cons_c1_dt5e-4` (157091) | dt 5e-4, 400 steps | 0.99995 ± 0.00032 | 1.00000 ± 0.00015 | 0.9987 ± 0.0003 | 1.000000 | 1 + 4.4e-7 ± 6.8e-7 | 1.0e-9 | 2e-12 | 30 min |
| T10d `T10d_VB_L2_cons_c2` (157092) | 2 correctors | 0.99994 ± 0.00034 | 1.00000 ± 0.00016 | 0.9987 | 1.000000 | 1 − 5e-8 ± 1.7e-7 | 3.8e-9 | 2e-11 | 20 min |
| T10d `T10d_VB_L2_cons_c4` (157093) | 4 correctors | 0.99994 ± 0.00034 | 1.00000 ± 0.00016 | 0.9987 | 1.000000 | 1 ± 8e-17 | 3.8e-12 | 3e-14 | 29 min |
| T10e `T10e_VB_L2_cons_c1_uniform` (157094) | uniform sink | 0.99988 ± 0.00065 | 0.99996 ± 0.00016 | 0.9986 ± 0.0006 | 1.000000 | 1 − 6.8e-5 ± 4.1e-4 | 4.3e-7 | 2.8e-9 | 21 min |
| T10b `T10b_VB_L3_cons_c1` (157090) | level 3 | 0.99996 ± 0.00024 | 1.00000 ± 0.00001 | 0.9986 ± 0.0002 | 1.000000 | 1 + 2.7e-6 ± 1.7e-5 | 9.9e-8 | 8e-10 | 146 min |
| T10a2 `T10a2_VB_L2_cons_c1_fixbin` (157113) | T10a with 9634485e | 0.99994 ± 0.00033 | 1.00000 ± 0.00016 | 0.9987 ± 0.0003 | 1.000000 | 1 + 2.3e-6 ± 7.9e-6 | 1.5e-7 | 7e-10 | 16 min |

- **The identity closes.** With form C the ε_h-weighted return flow equals the deposited first-moment
  rate step by step (`uf_eps_z` V/m1rate = 1 within 2e-6 with one corrector, 1e-16 with four), and
  once the kernel deposition has settled (t ≥ 0.1; it ramps from 0.9963 over the first ~50 steps, which
  is all the 6e-5 of the t ≥ 0.002 mean) m1rate equals N V_p|w| to six digits, so the settled
  ε_h-weighted G5 is 1.00000 ± 0.00016 at level 2, independent of dt (T10c), of the corrector count
  beyond one (T10d: the residual contracts by 40-1000 per pass, G5 unchanged to 1e-5) and of the sink
  mode (T10e: the sink vanishes at the fixed point, as derived; its spread is larger only because the
  uniform sink converges more slowly within the single pass).
- **Level 3 closes the same way** (T10b, 2 h 26 min on 28 ranks, 200 steps, no NaN): ε_h-weighted G5
  0.99996 ± 0.00024, settled 1.00000 ± 0.00001 (the level-3 scatter is 16 times smaller than at level 2),
  identity within 3e-6, residual 9.9e-8, sink 8e-10. T9a had +0.21 % here with the advective form; the
  h-dependence of that bias is gone, as the derivation requires.
- **Like quantities.** The campaign G5 (P0 ε times the Q2 velocity) reads 0.9987 (−0.13 %), against
  +0.55 % for T9a: the remaining gap between the two readings is the representation difference between
  the P0 storage porosity and its Q1 reconstruction at level 2, not a continuity defect (the P0-weighted
  Q2-rule integral of the same loop gives 0.99867). Compare CFDEM's face-flux G5 (1.00002) with the
  ε_h-weighted number. At level 3 the campaign G5 reads 0.9986 (−0.14 %): the representation gap
  between the P0 and the Q1 weighting of the same field does not shrink from level 2 to level 3 on
  this ladder; why it is h-independent is not explained here (follow-up, together with the
  campaign-slot dedup).
- Phi, S and the coordinate identities stay at 1e-19 to 1e-17 throughout (closed box, no clipping,
  relaxation 0); `dvol_gap` −2.5e-9 (single-precision element volumes).
- T10a2 (review-fix binary 9634485e) reproduces T10a to rounding, not bitwise: the residual rows moved
  into a pure array-expression routine and the compiler contracts it differently, so the conservative
  path differs in the last bit (worst relative difference 5e-7 in a 1e-13 slip component, 1e-21 in
  the identity residuals, m1rate now carries the Q_z term); every G5 digit, residual norm and sink
  value is identical. The default path is unaffected (round 2 below).

Default gates (`tools/el_stream_compare.py`, every certified configuration of the merged-VA set):
round 1 with dc2221e0 = VB_L2 defaults vs the T4 log (157086), T9a keys vs the T9a log (157087,
`EL_VA_*` bytes identical), s01 G2_twoway_lev2 (157088) bitwise; the six 28-rank cases of round 1 were
cancelled unstarted when the review-fix rebuild replaced the binary and resubmitted as round 2 with
9634485e (157104-157112, `runs/coupling_p4p5/formc_round2`): tier-2 twin vs the certified reference,
semi-implicit twin vs the e0164ff5 run, VB_L2 defaults vs the T4 log, T9a keys vs the T9a log (`EL_VA_*`
bytes identical, 600/600/600/600/400 lines), s01 G2_twoway_lev2 vs its log, P4 explicit and P5 No vs the
pre-fix runs, and all nine vs the merged-VA runs of e7ae87fc including EL_FLUID_SOURCE: all bitwise.
Lesson (repeated): never rebuild a binary that queued jobs point to; the T1/T2 frozen-copy rule (§6)
applies to every rebuild.
