# D1.1 — Periodic SC-array drag vs Hasimoto (q2p1_dns_drag)

Campaign: `dns-validation-campaign-plan.md` §D1.1. First executed 2026-08-02.

## Case design

- Unit cell [0,1]³, single **fixed** sphere at (0.5, 0.5, 0.5), d/L = 1/3,
  φ = π/6·(1/3)³ = 0.0193925. Stokes regime: ρ = 1, μ = 1 (deck
  `Prop@Viscosity = 1d0` + `example.json` `fluidViscosity_` matching),
  body-force driving `SimPar@ConstantForcing = 0,0,1d-2` → Re ≈ 3e-3.
- Sphere placement trick (no code change): `setupDNSDragSerial` places a
  lattice with `wallCenterOffset = 3r`; `volumeFraction_ = 0` forces a
  single sphere, and **r = 1/6 makes 3r = 0.5** → exactly centered.
- Periodicity: EL-proven path (plan D1.1 note) — QBOX9 fixture from
  `q2p1_el_pipeflow/tier2_cases/momentum_conservation/mesh/` (9³ unit
  cube, all faces `Periodic`-tagged, 3×3×3/27-subdomain partition,
  **28 ranks**). `q2p1_dns_drag` hard-codes `dPeriodicity = 1.0`, which
  equals the unit-box length — correct here by construction; the
  deck-key port remains a D1.1 code deliverable.
- D/h ladder via `SimPar@MaxMeshLevel`: L2/L3/L4 → D/h = 6/12/24
  (h = 1/(9·2^(L−1)), d = 1/3). 400 steps at dt = 0.01 (BE) ≈ 4 diffusive
  times — force plateau flat to 5 digits.

## Staging checklist (learned via failures 137360–137362)

1. Partition goes under **`_mesh/<MeshFolder>/`** (`_mesh/QBOX9/sub00NN/
   GRID0001.tri`; rank 0 reads `_mesh/QBOX9/GRID.tri` — QBOX9 ships it).
   Project dir (`quiescent_box_9/`) stays at the rundir root.
2. **`_data/MG.dat` is required** by `CreateDumpStructures` (level-
   dependent, mesh-independent reference table; copy from any bench
   rundir — 6 levels = 41739 lines).
3. **`start/` dir required** by `init_fc_rigid_body` (FullC0ntact init:
   `data.TXT` + `sampleRigidBody.xml`; copy from the app dir).
4. Plus the usual: `_data` deck, `example.json` (pinned per case),
   empty `_dump/_gmv/_vtk`.

## Gates

- Momentum balance (estimator 2): steady F_drag vs f·V_cell.
- Hasimoto (1959): F/(6πμ r U_superficial) vs 1/(1 − 1.7601φ^{1/3} + φ
  − 1.5593φ²) = **1.8322** at φ = 0.0193925. U_superficial = bulk_flow
  column 1 (column 2 = interstitial; col1 = col2 × fluid_frac).
- Ladder: error vs D/h → measured convergence order of the α-integral
  force on a stationary interface.

## Results

| Job | Level | D/h | F_z steady | K = F/(6πμrU_sup) | vs Hasimoto 1.8322 | Momentum balance |
|---|---|---|---|---|---|---|
| 137363 | L2 | 6 | 1.00137e-2 | 1.623 | −11.4% | +0.14% PASS |
| 137364 | L3 | 12 | — | — | — | — |
| 137365 | L4 | 24 | — | — | — | — |

Symmetry check (L2): transverse forces ~1e-12 vs axial 1e-2.
Indicator fluid fraction 0.98005 vs exact 0.98061 (L2 volume error
0.06% of cell).
