# M2 — Periodic sphere arrays with the Chimera component (Phase 5)

Design: `chimera-integration-design.md` §10 row 5 ("Arrays + periodicity").
Sibling of the FBM ladder `applications/q2p1_dns_drag/validation_cases/
d11_hasimoto/RUNBOOK.md` (same case, same conventions). First executed
2026-09-04 on the single-core sandbox (9 oversubscribed ranks).

## Case design

- Unit cell [0,1]³, one **fixed** sphere, d/L = 1/3, φ = π/6·(1/3)³ =
  0.0193925 (R = 0.1666665 from `seed_array.py --phi`). Stokes regime:
  ρ = 1, μ = 1 (`Prop@Viscosity = 1d0`), body-force driving
  `SimPar@UseConstantForcing = YES`, `ConstantForcing = 0,0,1d-2`
  (Re ≈ 3e-3), backward Euler dt = 0.01 to t = 4.01 (≈ 4 diffusive
  times; the force plateau is flat to all printed digits from t ≈ 3).
- Periodicity: `SimPar@PeriodicLength = 1,1,1` (new shared parser key →
  `dPeriodicity`), background `_data/CHIMERA_BOX6` (6³ unit cube from
  `channel_tri.py --periodic`, all faces tagged `Periodic`), axis-aligned
  2×2×2 partition (`partition_periodic_box.py` = `PyPartitioner 1 -123 8`
  + `flatten_axis_partition.py`), 9 ranks, `SubMeshNumber = 1`.
- Chimera: cubed-sphere shell `sphere_shell_coarse.tri` (192 hexes,
  level 2 → 1536 hexes, 13842 Q2 dofs, 47670 unknowns), fitted to
  [R, 2R] (H = R, the `H_k = min(H_max, ½ gap)` rule with H_max = R;
  gap to the periodic image = 2d), Robin outer condition (α = 1),
  `ChimeraSubStokes = Yes` (linear submesh, factorized once; the Picard
  variant gives the same forces to 8 digits and costs 2 min/step).
  Strong (Chimera-S) deck `_data/q2p1_param_hasimoto.dat`; weak
  (Chimera-W) deck `_data/q2p1_param_hasimoto_weak.dat` with the Phase-4
  settings (nodal lumped penalty, plain projection, ramp 0.25/0.5) and
  γ = 1e5 (viscous rate μ/h² = 144 at L2).
- D/h ladder via `SimPar@MaxMeshLevel`: L2/L3 → D/h = 4/8 (h = 1/(6·2^(L−1))).
  The atmosphere stays at submesh level 2 (radial spacing H/4 = 1/24 =
  the L3 background spacing).
- Periodic-donor check: the same sphere placed on the box CORNER
  (`_data/hasimoto_corner.dat`, `--offset 0,0,0`): its atmosphere
  straddles all six faces, every Robin sample point is wrapped, and the
  hole/fringe markers come from minimum-image distances. The flow is a
  pure translation of the centred case, so F and U must agree to the
  discretisation-equivalence of the two node sets (the Q2 lattice is
  invariant under the half-cell shift, so: to round-off).

## Staging checklist

1. Build `q2p1_chimera` (any Chimera-capable build: the component is
   always compiled). The CMake target stages the shell fixture into
   `_data/CHIMERA_BOX6/`; copy the vendored `_data/CHIMERA_BOX6/*`,
   `_data/hasimoto_*.dat`, `_data/random8_phi005.dat` and the decks into
   the run directory (the FAC staging also applies: `start/`, `_data/MG.dat`).
2. Partition: `python3 tools/chimera_meshgen/partition_periodic_box.py
   CHIBOX6 _data/CHIMERA_BOX6/box6.prj` (the flat `sub0001/GRID0001..8`
   layout; the solver maps ranks to subgrids by host name, so the
   subgrid layout of the partitioner cannot be used on one host).
3. Deck comments must not contain `=` (the key/value splitter takes the
   LAST `=` of the line).
4. Run: `mpirun -np 9 [--oversubscribe] ./q2p1_chimera`; post-process:
   `python3 validation_cases/m2_hasimoto/hasimoto_k.py _data/prot.txt
   --phi 0.0193925 --every 100`.

## Gates and conventions

- K_ref = 1/(1 − 1.7601 φ^{1/3} + φ − 1.5593 φ²) = **1.8317** at
  φ = 0.0193925 (the D1.1 runbook quotes 1.8322 for the same formula).
- U = the k = 0 Fourier mode = volume average over the WHOLE cell
  (superficial velocity). Chimera measures it as the composite average
  (`ChimeraBulk:` line): 0 inside the body, the submesh solution in the
  atmosphere, the background elsewhere (27-point Gauss per element).
- F: Hasimoto's sphere force balances the mean pressure gradient over the
  whole cell. A body force acting on the FLUID only gives the surface
  traction F_s = f V_fluid at steady state; the mean-gradient force is
  F_H = F_s + f V_solid (the FBM constraint force contains this
  "buoyancy" automatically because the hole fluid is forced). Reported:
  **K_meas** = F_H/(6πμRU), **K_bal** = f V_cell/(6πμRU) (exact momentum
  balance), and the balance residual F_s/(f V_fluid) − 1 (a conservative
  coupling gives 0; the Chimera value measures the momentum leak of the
  overlap).
- Momentum balance |F_s/(f V_fluid) − 1| and the fluid fraction of the
  composite quadrature (0.980442 vs exact 0.980608) are the sanity gates.

## Results ladder (2026-09-04)

| Case | Variant | D/h | F_s | U_sup | K_meas | vs 1.8317 | K_bal | vs 1.8317 | balance | wall |
|---|---|---|---|---|---|---|---|---|---|---|
| H_S_L2 | S | 4 | 9.7261000e-3 | 1.7523269e-3 | 1.8020 | −1.63 % | 1.8165 | −0.83 % | −0.82 % | 4.4 min |
| H_W_L2 | W (γ 1e5) | 4 | 9.7242570e-3 | 1.7531250e-3 | 1.8008 | −1.69 % | 1.8157 | −0.88 % | −0.83 % | 4.3 min |
| H_Scorner_L2 | S, sphere on the box corner | 4 | 9.7260995e-3 | 1.7523268e-3 | 1.8020 | −1.63 % | 1.8165 | −0.83 % | −0.82 % | 4.4 min |
| H_S_L3 | S | 8 | TBD | | | | | | | |
| H_W_L3 | W | 8 | TBD | | | | | | | |

FBM ladder for comparison (D1.1, post-fix, raw K): −6.05 % at D/h 6,
−5.01 % at 12, −2.29 % at 24, −1.35 % at 48.

## Random array (Beetstra / Tenneti band)

`_data/random8_phi005.dat`: 8 spheres, φ = 0.05, random sequential
addition with a minimum surface gap of 0.74 d (`seed_array.py --seed 1`),
R = 0.11427, H_k ∈ [0.0865, 0.1143] (≥ 2 background cells at L3).
Reference closures for the superficial-velocity normalised drag at
Re → 0: Beetstra et al. (2007) F = 10φ/(1−φ)² + (1−φ)²(1 + 1.5√φ) =
1.759; Tenneti et al. (2011) F = 1/(1−φ)³ + 5.7φ/(1−φ)³ + 10φ²/(1−φ)⁴ =
1.530 (both per sphere, U = superficial). The mean over the 8 spheres
of K_meas (F_H per sphere; f V_solid split equally) is the number to
compare; one 8-sphere realisation scatters by several per cent.

| Case | Variant | D/h | mean K_meas | K_bal | balance | wall |
|---|---|---|---|---|---|---|
| R8_S_L3 | S | 5.5 | TBD | | | |
| R8_W_L3 | W | 5.5 | TBD | | | |

## Reading of the ladder

- **Periodic images work.** The corner-placed sphere (atmosphere across
  all six faces, all Robin points wrapped, hole/fringe from minimum-image
  distances) reproduces the centred case to 5e-8 in F and 6e-8 in U —
  the residual is the different position of the partition faces relative
  to the body, not the coupling.
- **Both variants agree** on the periodic array: strong and weak differ
  by 0.02 % in F and 0.05 % in U at D/h = 4 (the FAC lift disagreement of
  Phase 4 was the cut-cell staircase of a 2-D body at h = D/4; a sphere in
  creeping flow is benign).
- **Accuracy at a coarse background.** At D/h = 4 the Chimera K is within
  −1.6 % (measured) / −0.8 % (balance) of Hasimoto; the FBM ladder reaches
  −2.3 % only at D/h = 24 and −1.35 % at D/h = 48. The gap between K_meas
  and K_bal is the 0.8 % momentum leak of the overlap coupling (surface
  traction below f V_fluid); it is the number to watch on the L3 rung.
- **Cost.** With the linear submesh operator a coupling update is one
  back-substitution (400 steps ≈ 4.4 min on one oversubscribed core for
  9 ranks); the Picard path refactorises twice per step (≈ 2 min/step)
  and gives the same forces to 8 digits at Re 3e-3.

## Open rungs

- L3 atmosphere (`ChimeraSubmeshLev 3`, 12288 hexes, ~350k unknowns per
  body) at the L3/L4 background — memory is available (377 GB) but the
  factorisation is minutes per body on this core.
- Halo (neighbourhood) exchange upgrade of `chi_exchange` for many-body
  arrays on many ranks (the replicated collective is O(N_points) per
  rank; irrelevant for 9 ranks).
- Variationally consistent forces (the balance residual is the target).
