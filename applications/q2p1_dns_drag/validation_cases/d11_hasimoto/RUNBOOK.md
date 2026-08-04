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

> **All pre-fix numbers below were measured on a box whose periodic faces
> were solved traction-free** (D1.1b, datasheet `d11_periodicity_absent`
> / `d11_uxtest`): the Q2/E013 communicator contained zero periodic DOF
> pairs. The fix landed 2026-08-03; **use the post-fix table.** Anything
> derived from the pre-fix ladder (the r_h calibration curve, the −8.9%
> asymptote, the φ-scaling anomaly) has to be re-measured.

### Post-fix (periodic Q2 coupling active)

| Job | Level | D/h | F_z steady | K = F/(6πμrU_sup) | vs Hasimoto 1.8322 | Momentum balance |
|---|---|---|---|---|---|---|
| 137539 | L2 | 6 | 1.00087e-2 | 1.7213 | −6.05% | +0.09% PASS |
| 137540 | L3 | 12 | 1.00123e-2 | 1.7404 | −5.01% | +0.12% PASS |
| 137541 | L4 | 24 | 1.00052e-2 | 1.7902 | −2.29% | +0.05% PASS |
| 137611 | L5 | 48 | 9.98535e-3 | 1.8075 | −1.35% | −0.15% (still creeping +0.01%/0.1s at t=4) |

Both tables quote the **raw** K = F/(6πμr·U_sup), so pre- and post-fix are
directly comparable. Ladder increments no longer flatten at ≈1.70 — the
pre-fix convergence to a wrong asymptote is gone. Transverse forces stay
at 5e-12 (L2) / 2.8e-10 (L4).

**Closure analysis (datasheet `d11_rh_collapse`).** The residual ladder is
explained by two measurable discretization terms, no unexplained physics
left:

1. **Indicator-volume error** per rung (measured fluid_frac vs analytic
   0.9806075): the L2 discrete sphere is +0.95% in radius, L3 −0.43%,
   L4/L5 ≲0.05%. Dividing each K by (1+Δr/r)^1.87 removes it.
2. **O(h) interface widening**: the volume-corrected superficial ladder
   is near-first-order (increments +0.063/+0.037/+0.017, orders
   0.76/1.18) with a level-constant grid coefficient c = (r_h−1)·a/h ≈
   0.12–0.17, i.e. **a_eff ≈ a + 0.14h**. Richardson: K_∞ = 1.821–1.825
   (−0.4…−0.6% vs 1.8322) — closes within resolution uncertainty.

**Convention PINNED from the paper** (`literature/hasimoto_1959.pdf`):
eq. (2.14) P₀ = F/τ₀ — the sphere force balances the mean pressure
gradient over the WHOLE cell, so with body-force driving F = f·V_cell,
which is exactly the measured force (no f·V_fluid "exterior" correction);
and eqs. (2.5)+(4.7–4.8) — the U in F = 6πμaU/(1−κa) is v₀, the k=0
Fourier mode = volume average over the entire unit cell = the superficial
(Darcy) velocity (bulk_flow col 2). The raw K_sup table above IS
Hasimoto's convention on both counts. The earlier `d11_convention`
"exterior drag = f·V_fluid" reading (decided on the traction-free box)
is superseded; the interstitial-velocity variant (−2.2% residual) is the
wrong comparison, not a discretization defect.

Superseded intermediate reading: the 3-point (L2–L4) Richardson quoted
below commit 63482ea1 (order ≈0.33, K_∞ ≈ 1.815) mixed conventions and
lacked the volume correction; the collapse analysis above replaces it.
Re-measured hydrodynamic-radius calibration: r_h/r =
**1.034 / 1.028 / 1.013** at D/h = 6/12/24 — roughly half the pre-fix
values and clearly shrinking with h, i.e. most of the apparent effective-
radius offset was the traction-free box, not the FBM interface.

Verification that the coupling is real: `PERIODIC_COMM` lines in the log
report 702 = 27×26 neighbour pairs on the QBOX9 3×3×3 torus (every rank
neighbours every other rank), 386 of them periodic-only. The u_x
symmetry-plane test dropped from 0.102 to 1.75e-9 (job 137535).

### Pre-fix (traction-free faces — historical)

| Job | Level | D/h | F_z steady | K = F/(6πμrU_sup) | vs Hasimoto 1.8322 | Momentum balance |
|---|---|---|---|---|---|---|
| 137363 | L2 | 6 | 1.00137e-2 | 1.623 | −11.4% | +0.14% PASS |
| 137364 | L3 | 12 | 1.00173e-2 | 1.640 | −10.5% | +0.17% PASS |
| 137365 | L4 | 24 | 1.00116e-2 | 1.684 | −8.1% | +0.12% PASS |

Pre-fix implied hydrodynamic-radius calibration (ten Cate §IV-C analogue,
sensitivity dlnK/dlnr ≈ 1.87): r_h/r = 1.067 / 1.061 / 1.046 at
D/h = 6/12/24 — **superseded**, see above. Indicator VOLUME error is
only −0.4% (fluid fraction), so the residual deficit sits in the
effective no-slip surface, not the represented volume.

## Partitioning prerequisite (periodic decks)

Periodic runs require an **axis-uniform Cartesian partition with ≥2 ranks
per periodic axis**, so that subdomain faces coincide with their periodic
images. METIS graph partitions are invalid. This is now enforced: the
coarse-face neighbour decode in `parentcomm.f90` aborts with an
explanatory message (`CheckFaceClaimDecode`) if a coarse face is not
claimed by exactly two distinct subdomains. QBOX9 (3×3×3/27) satisfies
this by construction.

MUMPS (`*@MGCrsSolverType = 5`) is refused on periodic decks —
`Create_GlobalNumbering` swaps periodic DOF labels instead of merging
them, so the global numbering handed to MUMPS is wrong.

Symmetry check (L2): transverse forces ~1e-12 vs axial 1e-2.
Indicator fluid fraction 0.98005 vs exact 0.98061 (L2 volume error
0.06% of cell).
