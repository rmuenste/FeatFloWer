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
  t < 0.01). So the O(dt) term is not "fixed by dt alone": its magnitude depends on h as well
  (Q2 ∝ Σ|e| (∂tε)² — the cell-averaged kernel field is smoother at level 3), and the SPEC
  12.3 statement "growing with refinement at a fixed dt" is wrong in sign.
- **T4 confirms the budget identities** (new `EL_VA_MOMENT` line and `uz_unw=` field, print-only,
  gated bitwise against `tier2_cases/momentum_conservation/reference/run_g0_twin_va_2026-10-14.log`
  by job 155115, `tools/el_stream_compare.py`: all 13 EL streams identical, EL_MEAN_SLIP with one
  extra trailing token): ∫u_z dV = −M1raw + ½ int_g to rel rms 2.0e-5 (max residual 1.1e-9) over
  199 steps (predicted ≈ 1e-3); ∫u_z / MOM 0.9867 (predicted ≈ 0.987); Q2 7.28e-7 (7.3e-7),
  adv_net +2.94e-7 (+3.0e-7), int_g + adv_net 7.42e-7 (7.5e-7), S2 2.92e-7; G5 and int_g equal
  to the base run to four digits.
- T5 (fix A: α-weighted removal of the net source; fix B: divide by (ε^{n+1}+ε^n)/2) is a physics
  change and has NOT been applied; it waits for the owner's decision.
