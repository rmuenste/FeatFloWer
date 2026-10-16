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

**VB return flow (G5, gate 2 %):** `EL_VA_SOURCE` L1_dedt 1.1291e-3 m³/s (predictor 1.1291e-3),
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
| twin uf_super_z EL / CFDEM VB-b-divided | 1 (2e-2) | 0.980 / 1.007 = 0.973 | FAIL by 0.7 % (sum of the two codes' deviations) |

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
