# s01_settling_cfdem_twin — single-sphere settling twin of CFDEM `settlingTestMPI`

FeatFloWer side of the unresolved Euler-Lagrange cross-code comparison with CFDEM
(`cfdemSolverPiso`), case 01 of the LIGGGHTS-PUBLIC specs. Spec, case table, ODE reference
and post script: `/data/warehouse17/rmuenste/code/liggghts-dem/LIGGGHTS-PUBLIC/md_docs/case_specs/01_settling/`
(`SPEC.md` §4 cases, §6 FeatFloWer recipe, §8 gates; `settling_cases.py`, `settling_post.py`,
`featflower/stage_case.sh`). Binary: `build-va` `q2p1_el_pipeflow` (branch
`feature/el-volume-averaged`, built 2026-10-15 15:14, c39c7334 tree). Mesh: tier-2
`quiescent_box_9`, `MaxMeshLevel 2` (18³, h = 1/18), partition `QBOX9Z3`, 4 ranks. Run
directories: `/data/warehouse17/rmuenste/runs/settling01/` (local Slurm jobs 150153–150170,
short partition, 2026-10-16; level-3 checks 150171–150173). `q2p1_param.dat`, `example.json`,
`particles.xyz` here are the spec's primary deck (M-L2, one-way, difelice_std); the other 17
differ only in the keys `make_case.py` sets (dt, steps, kernel width, drag, feedback, viscosity).

Variants: **M** d = 0.005, `ELKernelWidthFactor 20` (kernel radius 0.1 = 20 d, self-voidage
5.4e-5, the like-for-like twin of CFDEM's 25 d cell); **N** d = 0.04, factor 2.5 (campaign
standard, kernel radius 0.1 = 2.5 d, d/h = 0.72). All runs `ELPressureForce No`,
`ELFluidGravity Yes`, `ELVolumeAveraged none`, `ELDragCoupling explicit`.

## 1. One-way, variant M (gate G1: |mean u/u_ode − 1| ≤ 1 %)

| run | drag | plateau u | prediction at ε_self | plateau/pred − 1 | u/u_ode − 1 | verdict |
|---|---|---:|---:|---:|---:|---|
| M-G2 (ρ_r 300, Ar 29.3) | difelice_std | 0.917827 | 0.918171 | −3.7e-4 | −6e-5 | PASS |
| M-L1 (Re 0.1) | difelice_std | 0.0326617 | 0.0326649 | −1.0e-4 | −9e-5 | PASS |
| M-L2 (Re 1.29) | difelice_std | 0.106190 | 0.106199 | −0.8e-4 | −8e-5 | PASS |
| M-L3 (Re 80) | difelice_std | 0.442612 | 0.442631 | −0.4e-4 | −4e-5 | PASS |
| M-L2 | difelice (campaign ε^−χ) | 0.106183 | 0.106199 | −1.5e-4 | −1.4e-4 | PASS |
| M-L2 | schiller_naumann | 0.113405 | 0.113405 | +0.0 | +0.0 | PASS |

G3 law discrimination: schiller_naumann / difelice_std = 1.0680 at M-L2 (spec: 6.8 % above). The
two Di Felice forms coincide at ε ≈ 1 (−0.7e-4 apart). Force budget at the plateau closes to
1e-4 of the weight (`EL_FORCE_BUDGET total`).

## 2. Cross-code (gate G2: Re_t FF M vs CFDEM at equal Ar within 1 %) and two-way (R1)

CFDEM one-way plateaus from their RUNBOOK runs 33–41 (SI, d = 1e-4 m); Re_t = ε d u/ν on each side.

| case | FF u_t | Re_t FF | Re_t CFDEM | ratio − 1 | FF two-way/one-way − 1 | CFDEM two-way − 1 |
|---|---:|---:|---:|---:|---:|---:|
| M-G2 | 0.917827 | 1.28561 | 1.28603 | -3.2e-04 | +0.85% | +0.31% |
| M-L1 | 0.032662 | 0.09999 | 0.10000 | -5.6e-05 | +0.19% | +0.12% |
| M-L2 | 0.106190 | 1.28600 | 1.28610 | -7.8e-05 | +0.17% | +0.06% |
| M-L3 | 0.442612 | 79.99673 | 79.99781 | -1.3e-05 | +0.05% | +0.02% |

G2 PASS on all four (|ratio − 1| ≤ 3.2e-4; spec expected ≤ 2e-4 from the self-voidage difference,
G2 is 1.6× that). R1: the FeatFloWer two-way bias is 1.5–3× CFDEM's at the same ordering
(largest for the ρ_r = 300 case), 0.05–0.85 % for the 20 d kernel. Level-3 reruns (M-L2 two-way,
N-L2 one- and two-way; jobs 150171–150173) are the resolution check — see §5.

## 3. Variant N (gate G4: plateau within 2 % of the prediction at the logged ε)

**The gate cannot be applied as written.** `settling_post.py` reads `eps=` from
`EL_MEAN_SLIP_TAVG` as the void fraction at the particle; in `el_diagnostics.f90:157-163` that
line's `uf_intr`, `uf_super` and `eps` are *domain* means (∫ρ ε u dV / ρ∫ε dV and ∫ε dV / V),
which is why `eps=` prints 0.999966 (= 1 − V_p/V_box) for variant N and `uf_intr` ≈ 1e-10 in
every run, two-way included. The kernel-sampled ε and u_f at the particle, which the drag uses,
are not logged (per-particle log has position, velocity, force only). Spec §6.6 and the post
script need correcting on the peer's side; on ours a per-particle `eps_p`/`uf_p` column in the
particle log is the fix (not implemented here; it changes the binary and needs its own twin gate).

What can be done instead: invert the drag law at the plateau for the ε the drag saw.

| run | drag | plateau u | ε inverted | u/u_ode − 1 at ε 0.9724 (continuum) | at ε 0.980 |
|---|---|---:|---:|---:|---:|
| N-L1 | difelice_std | 0.0878584 | 0.98035 | +2.09 % | +0.09 % |
| N-L2 | difelice_std | 0.288091 | 0.98002 | +1.64 % | +0.01 % |
| N-L2 | difelice | 0.277807 | 0.97992 | +3.03 % | −0.03 % |
| N-L3 | difelice_std | 1.21815 (box too short, u/u_t 0.996) | — | +0.76 % | −0.00 % |

The sampled ε at level 2 for d = 0.04, factor 2.5 is **0.980**, consistent across Re 0.1–80 and
across both Di Felice forms to 1e-4 (the ε^−2 ratio of the forms is confirmed by the df/std pair).
It is above the continuum estimate 0.9724 (spec §4.3) because the deposited volume is
cell-averaged before it is sampled back; v1b measured 0.974–0.978 on its meshes. Against the
a-priori ε_self = 0.9724 the 2 % tolerance holds for L2/L3 and is missed by 0.1 % at L1;
against the inverted ε every run is at 1e-3. Verdicts recorded as such.

## 4. Two-way, variant N (R1, informative)

| run | one-way u | two-way u | two-way/one-way − 1 |
|---|---:|---:|---:|
| N-L1 | 0.0878584 | 0.107863 | +22.8 % |
| N-L2 | 0.288091 | 0.35757 | +24.1 % |
| N-L3 | 1.21815 | 1.26009 | +3.4 % |

At d/h = 0.72 the particle rides its own deposited momentum: +23 % at Re ≤ 1.3 (v1b E4 at
δ/h 4.5: +6 %), +3.4 % at Re 80 where the particle outruns its wake (crossing time of the kernel
0.08 s vs viscous time 16 s). The drag at the plateau is the weight in every case (budget), so the
speed-up is entirely a fluid velocity at the particle of u_two − u_one (7.8e-3 for M-G2,
0.069 for N-L2), which the current diagnostics do not print (§3).

## 5. Level-3 resolution check (jobs 150171-150173, `MaxMeshLevel 3`, 36³, h = 1/36)

| run | level 2 | level 3 | change |
|---|---:|---:|---:|
| M-L2 two-way plateau | 0.106373 | 0.106372 | −1e-5 (two-way bias +0.17 % at both levels) |
| N-L2 one-way plateau | 0.288091 | 0.285013 | −1.1 % |
| N-L2 one-way ε inverted | 0.980 | 0.975 | → continuum 0.9724 (δ/h 1.8 → 3.6) |
| N-L2 two-way plateau | 0.35757 | 0.354616 | −0.8 % |
| N-L2 two-way/one-way − 1 | +24.1 % | +24.4 % | resolution-independent |

Variant M is mesh-converged at level 2 (one-way at 1e-4, two-way bias unchanged to 1e-5). For
variant N the sampled ε moves toward the continuum kernel estimate as the cell shrinks
(0.980 → 0.975 → 0.9724 at δ/h 1.8 / 3.6 / ∞): the level-2 value is a cell-averaging artefact
of d/h = 0.72, which is why the a-priori G4 prediction is 1.6–2.1 % low at level 2. The N two-way
bias is set by the kernel width (2.5 d), not the mesh: +24 % at both levels.

## 6. Pitfalls

- `settling_post.py`'s "eps at particle (logged)" and "uf_intr co-flow" are domain means (§3).
- `grep -i error` on the logs matches `EL_VOLUME_CONSERVATION ... rel_error=`; the runs are clean.
- Spec deck keeps `OutputLevel = MAX+1` (drops EL VTK fields); harmless here, no VTK written.

## 7. Peer confirmation (liggghts-public-1e, 2026-10-16)

The domain-mean reading of `EL_MEAN_SLIP` was confirmed in our source and fixed in their spec
(SPEC 4.3/6.6/8.1/8.2, pitfall 9.11) and in `settling_post.py` (prediction at the kernel
self-voidage or `--eps`, inverted ε printed with a non-plateau warning, N gated at 2 %, brentq
fallback without scipy). Their re-run over `runs/settling01` reproduces every number above:
inverted ε 0.98035 / 0.98003 / 0.97993 (level 2), 0.97499 (level 3); G4 +2.09 / +1.64 / +0.76 %;
two-way +0.85 / +0.19 / +0.17 / +0.05 % (M), +22.8 / +24.1 / +3.4 % (N). Cross-code Re_t with each
side's ε inside Re: −3.4e-4 / −7.7e-5 / −9.8e-5 / −3.4e-5 (ours −3.2e-4 / −5.6e-5 / −7.8e-5 /
−1.3e-5; the difference is the ε factor). All SPEC 8.3 rows are filled on their side (SPEC 11.2).

Open on R1 (M two-way bias 1.5–3× CFDEM at matched scale), candidate causes not yet separated:
(a) the kernel's effective volume V_p/α_self is 1.6× smaller than CFDEM's cell; (b) our closed
box against CFDEM's open bottom. Requested from us when the per-particle u_f column exists: the
M-G2 two-way u_f at the particle. CFDEM model-A (ε^(2−χ)) variants for 04/05/06 are not run yet.

## 8. Same-mesh twin (SPEC 01 §12; jobs 155064-155101, 2026-10-27)

Peer request `LIGGGHTS-PUBLIC/md_docs/case_specs/01_settling/twin/REQUEST_featflower.md`: the CFDEM
column with y and z swapped (0.01 × 0.01 × 0.05 m, 4 × 4 × 20 coarse hex, h = 2.5 mm = 25 d), SI
units of the CFDEM case table, the same dt and end time, particle at the cell centre
(0.00625, 0.00625, 0.02625), `difelice_std`, kernel width 21.4 d (the inviscid response of
CFDEM's centred cell), `featflower-partition 4 -3 4` (four z-slabs), pe `processesZ_ 4`,
`mpirun -np 5`. Staged with the peer's `twin/stage_ff_twin.py`, post-processed with
`twin/twin_post.py` (u_f inverted from `EL_FORCE_BUDGET` with the one-way twin's ε(t)).
Run dirs `/data/warehouse17/rmuenste/runs/settling01_twin/<case>_<oneway|twoway>_lev<1|2|3>` and
`..._lev2_wf10p7`; `submit_all.sh`, `post_all.sh`, `summary.txt` there and in `twin/` here.
Binary `FeatFloWer-va/build-va` q2p1_el_pipeflow (built 2026-10-15; `git diff fc63d601 HEAD --
source libs/pe applications/*.f90` is empty, so the binary is HEAD 3c2e7055 code).

### 8.1 Checks requested by the peer

1. **`Velo@MinDef` is absolute.** `QuadSc_main.f90:169`: the nonlinear stop is
   `max(RhsUVW·defCrit, MinDef)`. Set to 1d-16 in all 32 decks (the MG criteria scale with the
   velocity norm and were left alone). Levels 2/3 are unchanged by it to the digits printed.
2. **pe box = mesh.** `app_init.f90:454` passes the global mesh extents; the log prints
   `Global CFD domain extents 0 0.01 0 0.01 0 0.05`, six box walls, `PE decomposition 1x1x4`
   (slabs of 0.0125 = the four 80-element partition slabs).
3. **Level 1 as staged diverges.** With `Pres@MGMinLev = MGMedLev = MaxMeshLevel = 1` the driver
   clamps the pressure cycle count to 1 (`QuadSc_mg.f90:52-61`) and coarse-solver type 4 is a
   single approximate sweep (defect ×0.48 per step): the hydrostatic imbalance grows ×1.3 per
   step, the fluid reaches 30 m/s, NaN at step 45 (kept as `L1_oneway_lev1.as_staged_NaN`).
   Fix in the eight level-1 decks only: `Pres@MGCrsSolverType = 1` (BiCGStab on the P1 system)
   and `Pres@MGDefImprCoarse = 1d-6`. Patched level 1: pressure defect 1e-17 in 45 iterations,
   residual fluid speed 1.4e-9 m/s (level 2: 4e-10), plateau equal to level 2 to 7 digits.
   The `.as_staged` deck is kept next to every patched deck.

### 8.2 Results (two-way bias = mean u/u_ode − 1 over the last 20 %; u_f window mean / last)

| Case | lev | one-way u/u_ode−1 | ε (window) | two-way plateau | two-way bias | u_f [m/s] | CFDEM 4×20×4 |
|---|---|---|---|---|---|---|---|
| G2 | 1 | +2.2e-8 | 0.999993 | 0.1294284 | +0.639 % | 7.57e-4 / 7.67e-4 | +0.659 % |
| G2 | 2 | −1.7e-6 | 0.999995 | 0.1296380 | +0.802 % | 9.48e-4 / 9.46e-4 | |
| G2 | 3 | −5.9e-7 | 0.999987 | 0.1296472 | +0.811 % | 9.61e-4 / 9.61e-4 | |
| L1 | 1 | −1.3e-7 | 0.999974 | 4.582223e-3 | +0.155 % | 7.58e-6 / 8.46e-6 | +0.210 % |
| L1 | 2 | −7.3e-8 | 0.999974 | 4.584047e-3 | +0.195 % | 9.53e-6 / 1.06e-5 | |
| L1 | 3 | −8.4e-8 | 0.999961 | 4.583885e-3 | +0.195 % | 9.52e-6 / 1.06e-5 | |
| L2 | 1 | −6.1e-7 | 0.999976 | 1.489542e-2 | +0.138 % | 2.19e-5 / 2.43e-5 | +0.193 % |
| L2 | 2 | +2.3e-7 | 0.999972 | 1.490123e-2 | +0.178 % | 2.82e-5 / 3.15e-5 | |
| L2 | 3 | −6.2e-8 | 0.999961 | 1.490063e-2 | +0.176 % | 2.80e-5 / 3.11e-5 | |
| L3 | 1 | +5.1e-7 | 0.999983 | 6.202164e-2 | +0.037 % | 2.25e-5 / 2.24e-5 | +0.035 % |
| L3 | 2 | +1.4e-7 | 0.999973 | 6.202669e-2 | +0.046 % | 2.84e-5 / 2.88e-5 | |
| L3 | 3 | +1.6e-8 | 0.999961 | 6.202622e-2 | +0.047 % | 2.87e-5 / 2.90e-5 | |

One-way plateaus: G2 0.1286061 / 0.1286063 / 0.1286044, L1 4.575137e-3 / 4.575137e-3 /
4.574987e-3, L2 1.487491e-2 / 1.487481e-2 / 1.487445e-2, L3 6.199869e-2 / 6.199804e-2 /
6.199728e-2 (lev 1/2/3). All one-way runs match the ODE to ≤ 3.3e-6.

Level 2 with kernel width 10.7 d (counterpart of CFDEM's 8×40×8 mesh, δ/h = 0.86 again):

| Case | one-way u/u_ode−1 | ε | two-way plateau | two-way bias | u_f [m/s] | CFDEM 8×40×8 |
|---|---|---|---|---|---|---|
| G2 | −3.3e-6 | 0.999975 | 0.1326053 | +3.11 % | 3.33e-3 / 3.34e-3 | +2.578 % |
| L1 | +7.9e-8 | 0.999966 | 4.645468e-3 | +1.54 % | 7.53e-5 / 8.40e-5 | +1.599 % |
| L2 | +3.5e-6 | 0.999947 | 1.507506e-2 | +1.35 % | 2.15e-4 / 2.39e-4 | +1.490 % |
| L3 | +1.8e-6 | 0.999952 | 6.211752e-2 | +0.19 % | 1.15e-4 / 1.22e-4 | +0.172 % |

Reading: at matched smoothing volume (δ/h 0.86, level 1 at 21.4 d and level 2 at 10.7 d) the
FeatFloWer two-way bias is within 0.06 points of CFDEM for G2/L1/L2/L3 at mesh 1 and within
0.5 / 0.06 / 0.14 / 0.02 points at the mesh-2 counterpart. Refining the mesh under a fixed kernel
(lev 2 → 3) changes the bias by ≤ 0.01 points: the FeatFloWer bias is set by the kernel, CFDEM's
by the cell (its bias grows 3.9-7.6× from mesh 1 to mesh 2, ours 4.9-9.9× from 21.4 d to 10.7 d
at level 2). The co-flow u_f is still ramping at the end of the L1/L2 runs (last > window mean).
ε ≥ 0.99995 everywhere: self-voidage is irrelevant at d/h = 0.04.

### 8.3 Follow-ups on the unit-box decks (jobs 155064-155068; u_f by `twin/ff_unitbox_uf.py`)

| Run | result | compared with |
|---|---|---|
| a. M-L2 two-way, 540 steps (`fu_a_M-L2_twoway_std_540`, one-way partner `fu_a_M-L2_oneway_std_540`) | bias +0.49 %, u_f 5.27e-4 window mean / 5.67e-4 last; one-way plateau 0.1061901 (+3e-5 vs ODE) | 180 steps: +0.17 %, u_f 1.97e-4 / 2.18e-4; Oseen prediction u_f ≈ 5.6e-4, +0.53 % |
| b. N-L1 at z0 = 0.5167 (`fu_b_N-L1_{oneway,twoway}_std_z0.5167`) | one-way plateau 8.78668e-2, inverted ε 0.98042 (z0 0.85: 8.78584e-2, 0.98041); two-way plateau 0.1095177, bias +24.6 %, u_f 2.215e-2 / 2.282e-2 | z0 0.85: +22.8 %, u_f 2.040e-2 / 2.090e-2; prediction u_f → ≈ 0.0224 |
| c. N-L2 one-way, level 4, 30 steps (`fu_c_N-L2_oneway_std_L4_30steps`, 14.5 min, 22.6 GB RSS) | inverted ε 0.97309 at step 1, 0.97315 at step 30 (window mean 0.97315) | prediction 0.9731; ladder 0.980 / 0.975 / 0.9731 at levels 2/3/4, continuum 0.9724 |

Wall times (short partition, 1 node): level 1 10 s, level 2 30 s, level 3 2.4-4.8 min (5 ranks).
