# D5.1 numerical viscometer — case specification (for sign-off)

Status: DRAFT for owner agreement, 2026-08-20. Nothing built yet.
Reference: Prignitz & Bänsch, Kybernetika 46 (2010) 281–293
(`literature/baensch_num_visc.pdf`), §4.2 rotational viscometer; their
2D torque-ratio protocol, done here in 3D with resolved FBM spheres.

## 1. The instrument

Searle configuration: stationary outer vessel, rotating inner bob.
The bob is NOT meshed and NOT an FBM body — it is a blind coaxial bore
in the mesh whose surfaces carry a rotating Dirichlet condition
u = Ω×r = (−Ωy, Ωx, 0), steady in the lab frame.

Geometry (units of particle diameter d = 1; axis = z through origin):

| item | value |
|---|---|
| vessel (mesh outer boundary) | cylinder r_a = 10, z ∈ [0, 12], closed |
| bore (not meshed) | coaxial cylinder r_i = 4, from top face z=12 down to z=4 |
| bob immersion depth | 8 (bottom clearance 4, "partially hollow") |
| annular gap | 6 = 6 particle diameters |
| sheared annulus volume (bore height) | π(r_a²−r_i²)·8 ≈ 2111 d³ |

## 2. Boundary conditions (.par plan)

| component | .par type | geometry line | motion |
|---|---|---|---|
| outer wall | `Wall` | type 7: `0 0 0 10 1 1 0` | 0 |
| bottom z=0 | `Wall` | type 4 plane | 0 |
| top annulus z=12 | `Wall` | type 4 plane | 0 (closed lid; a real device's free surface — absorbed by the ratio protocol) |
| bore disc face z=4 | `Inflow14` | type 4 plane `0 0 1 -4` | u = Ω×r |
| bore side wall | `Inflow14` | type 7: `0 0 0 4 1 1 0` | u = Ω×r |
| bore side+disc (overlay) | `WallF` | same vertex sets | torque mask |

Authoring rules (from the code scout, all verified file:line):
- `.prj` order: the disc-plane .par MUST precede the bore-cylinder .par
  (components overwrite shared ring vertices in order; plane sets z,
  cylinder then sets r and leaves z — reversed order breaks the ring).
- The overlay trick (second .par on the same vertex set) is in-tree
  precedent (`q2p1_sse/_ianus/ZIM/ANNULAR/meshDir/in.par` + `inT.par`).
- `WallF` must tag ONLY the bore (BndrForce is a single global mask:
  anything else tagged WallF pollutes the torque sum).
- Shared bore/vessel edge nodes: `Inflow` wins over `Wall` in
  `Boundary_QuadScalar_Val` — acceptable (the contact ring rotates);
  if unwanted, exclude ring vertices from the Inflow14 files.

## 3. Physics and operating point

| quantity | value | rationale |
|---|---|---|
| fluid ν | 0.2 | see Ta/Re budget |
| ρ_f | 1.0 | |
| particles | spheres d = 1, ρ_p/ρ_f = 1.0 (neutrally buoyant) | isolate rheology |
| gravity | 0 | kill residual buoyancy numerics |
| Ω | 0.1 rad/t.u. (deck: SimPar@RPM = 60·Ω/2π ≈ 0.9549) | see budget |
| nominal gap shear rate γ̇ ≈ Ω·r_i/gap | 0.067 | |
| particle Re_p = γ̇ d²/ν | 0.33 | Stokes-ish, Einstein-comparable |
| gap Reynolds Ω r_i (r_a−r_i)/ν | 12 | laminar |
| Taylor number Ω² r_i gap³ / ν² | 216 | ≥8× below the narrow-gap critical ~1712 (itself a lower bound at radius ratio 0.4) — rotating-inner is the unstable direction, this is the safety margin |
| viscous transient τ_ν = gap²/ν | 180 t.u. | sets run length |
| bob period 2π/Ω | 63 t.u. | average torque over ≥1 period |
| dt | 0.05 (u_max = Ω r_i = 0.4, h_min = d/8 → CFL ≈ 0.16) | |
| steps | 5000 → t = 250 ≈ 1.4 τ_ν + averaging window | d32-like budget |
| steadiness gate | endpoint torque drift < 5×10⁻³ over final t.u. (D3.1 gate, applied to T not U) | |

## 4. Resolution

Coarse cell ≈ d/2, production at MG level with h = d/8 (D/h = 8) —
the DKT/D3.2-proven operating point. Bob radius 4d → 32 cells/radius.
Known biases at D/h=8 that the RATIO protocol cancels or bounds:
end effects of the blind bore, closed-lid drag, the uniform +2%
force-functional bias (fbm_functional_closure); what it does NOT
cancel: per-particle grid noise (D1.2: ~4% rms at D/h=8 per particle —
mitigated by time-averaging the *boundary* torque over ≥1 period) and
near-contact under-resolution at high φ (out of scope for the gate).

## 5. Measurement and gates

Observable: η◆(φ) = T_z(φ)/T_z(0) at fixed Ω, T_z from GetDNATorque
on the WallF bore surface, time-averaged over the final period after
the steadiness gate passes.

| run | φ | N (annulus basis) | gate |
|---|---|---|---|
| V0 baseline | 0 | 0 | T_z(0) plateau; sanity vs infinite-annulus analytic per-height torque 4πμΩ/(r_i⁻²−r_a⁻²)·H_imm (expect O(10%) end corrections, recorded not gated); bob-vs-vessel torque balance |
| V1 dilute gate | 0.05 | ≈ 200 (RSA in annulus, margins: 0.5d walls, 1d bob-bottom/lid) | η◆ vs Einstein 1+2.5φ = 1.125; PASS band ±0.03 (grid noise + dilute finite-φ corrections O(φ²)≈0.006) |

Decision point after V1: extend ladder (φ = 0.10/0.20 ×2–3 seeds) vs
Krieger–Dougherty + the EL campaign's plane-Couette values at matched φ.

Secondary diagnostics: φ(r) profile drift (shear-induced migration
toward the outer wall — physics, must be quantified over the averaging
window), particle angular velocity vs local fluid rotation.

## 6. Code changes (all verified against source; torque routine read line-by-line)

Torque strategy — TWO estimators, campaign two-estimator discipline:

- PRIMARY: new residual-based torque, the variationally consistent
  analog of the FAC BenchForce method. Extend the machinery of
  `EvaluateDragLift_old` (`source/src_quadLS/QuadSc_def.f90:3856-3904`,
  benchmark-certified) with the rigid-rotation test function
  v = e_z×r: T_z = Σ over marked bob DOFs of (x_i·R_y,i − y_i·R_x,i),
  lever arms from the Q2 DOF coordinates. ~30-50 lines in proven code.
- SECONDARY (cross-check): `GetDNATorque`
  (`source/src_quadLS/QuadSc_torque.f90:637-896`). Assessed B−:
  integrand correct (full deformation tensor + proper element-local P1
  pressure; T_z exact about z; parallel reduction fine) but it is the
  ∇α-band volume-form integral — the technique family measured at
  +1.9-2.8% systematic (`fbm_functional_closure`) — and it has never
  been validated in-tree (all call sites commented out). Suspect at
  the ~2% level in absolute value; fine as an independent check.
- V0 GATE: the two estimators must agree on T_z(0) within 3%.
  Disagreement quantifies the α-band bias on a boundary — a finding
  either way.

| change | size | where |
|---|---|---|
| `iT=14` rotating-BC block (u=Ω×r via RPM) | ~6 lines | `source/src_quadLS/QuadSc_user.f90:586` (imitate iT=771 at :407-413; NOT iT=13 — that reads uninitialized extruder deck fields) |
| residual-based torque (primary) | ~30-50 lines | new routine beside `EvaluateDragLift_old`, `QuadSc_def.f90:3856` pattern, BndrForce-masked DOFs |
| enable GetDNATorque call (secondary) | uncomment | `source/src_quadLS/QuadSc_main.f90:822` (or app postprocessing) |
| raw-torque output (drop hardcoded 1-RPM normalization) | 1 line | `source/src_quadLS/QuadSc_force_torque_calc.f90:463-464` — print raw T_z per estimator, one labeled line each, normalize offline |
| deck key for Ω | none | `SimPar@RPM` exists (default 12.0, parser param_parser.f90:1095) |
| torque axis | none | z-axis through origin = our axis; Center offset not needed |

Application: DEDICATED CLONE `applications/q2p1_viscometer` (owner
decision 2026-08-20; some older applications will be retired soon, the
slot is acceptable). Clone from q2p1_dkt (serial PE,
HardContactAndFluid, external particles.xyz packing,
gravity_=[0,0,0] supported); strip DKT-specific pieces, add the torque
postprocessing calls app-side. Build in the post-merge tree. The PE
z=0 ground plane in the inherited setup is harmless (particles live at
z ≥ 5).

## 7. Mesh order (for the meshing agent, after sign-off)

Hex mesh of cylinder-minus-blind-bore: butterfly/O-grid cross-section
for the annulus, solid O-grid core below the bore; target coarse cell
≈ d/2 (expect ~25-30k coarse elements); .par set per §2 with type-7
cylinders (flags 1 1 0) and type-4 planes; .prj order plane-before-
cylinder for the bore pair; partition ~108 subdomains (walled case —
METIS acceptable, Cartesian not required); verify projection by
refining once and checking r-histograms of boundary nodes at both
radii. Deliverables mirror dkt_mesh_v1: mesh dir + README + REPRODUCE.

## 8. Sign-off state

1. Setup/geometry: owner approved 2026-08-20 ("Setup seems ok")
   — r_i=4, r_a=10, H=12, immersion 8.
2. Operating point ν=0.2, Ω=0.1: proposed, not yet vetoed.
3. Scope: gate-first (V0+V1) — agreed.
4. App: dedicated clone q2p1_viscometer — owner decision 2026-08-20.
5. Torque: two-estimator strategy (residual primary, GetDNATorque
   secondary, 3% agreement gate) — added after owner questioned
   GetDNATorque quality; assessment in §6.
6. Case id: D5.1 (bridges DNS → EL closure validation).
