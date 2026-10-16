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

## 5. Level-3 resolution check

Pending: jobs 150171 (M-L2 two-way), 150172 (N-L2 one-way), 150173 (N-L2 two-way) at
`MaxMeshLevel 3` (36³, h = 1/36, δ/h = 3.6, d/h = 1.44 for N).

## 6. Pitfalls

- `settling_post.py`'s "eps at particle (logged)" and "uf_intr co-flow" are domain means (§3).
- `grep -i error` on the logs matches `EL_VOLUME_CONSERVATION ... rel_error=`; the runs are clean.
- Spec deck keeps `OutputLevel = MAX+1` (drops EL VTK fields); harmless here, no VTK written.
