# Phase 4 V&V Case: Single-Particle Terminal Velocity (Stokes)

Verification case for the Euler-Lagrange scheme (plan: Tier-2 test #8). Verifies
the full Phase-2 chain — drag closure × driver-side `grav_buoy` × semi-implicit PE
integrator — against the **exact analytic** terminal velocity and relaxation time.
References are analytic, not prior output.

## Design decision: run ONE-WAY, quiescent fluid

A single settling particle in a fully periodic box with **two-way** coupling never
reaches a lab-frame terminal velocity: gravity injects momentum into the particle,
the particle injects it into the fluid, and with no wall to absorb it the whole
mixture accelerates indefinitely. So this case uses **one-way** coupling:

- `ELApplyParticleForces = Yes`, `ELApplyFluidFeedback = No`.
- Fluid: zero initial condition, no body force, no feedback source ⇒ `u_f ≡ 0`
  for the whole run (trivial steady state of the momentum solve).
- The particle settles against the quiescent fluid; `u_slip = −u_p` and the
  particle reaches the analytic terminal velocity.

Consequence: with `u_f ≡ 0` the **mesh is not the sensitive ingredient** here — it
only samples a trivial zero field. A modest uniform hex cube is sufficient; no fine
or path-aligned resolution is needed. (Transfer/sampling accuracy is covered by the
Tier-1 convergence tests; two-way momentum behavior by conservation test #7.)

## Configuration

| Item | Value | Notes |
|---|---|---|
| Domain | unit cube `[0,1]³` | periodic in all directions (app default, set in `applications/q2p1_dns_drag/app_init.f90`) |
| Mesh | structured hex, **18³** (h = 1/18 ≈ 0.0556) | uniform; resolution is non-critical here (see note). d_p/h ≈ 0.72, kernel width / h ≈ 1.8 — particle stays sub-grid |
| Fluid density | `ρ_f = 1.0` | `example.json: fluidDensity_` |
| Fluid viscosity | `ν = 0.1` (kinematic) ⇒ `μ = ρ_f ν = 0.1` | `Properties%Viscosity(1)` is kinematic; `example.json: fluidViscosity_ = 0.1` |
| Particle | one **free** sphere, `r = 0.02` (d = 0.04) | `example.json: benchRadius_ = 0.02` |
| Particle density | `ρ_p = 5.0` | `example.json: particleDensity_ = 5.0` (density ratio 5) |
| Initial particle velocity | `0` | starts from rest, at a known interior position |
| Gravity | `g = [0, 0, −10]` | `example.json: gravity_ = [0.0, 0.0, -10.0]` |
| Drag model | `ELDragModel = stokes` | exact `F = 6πμr·slip`; matches analytic `u_t`. (DiFelice also reduces to this at Re→0, ε=1) |
| Lift | `ELLiftModel = none` | no shear; pure axial settling |
| Pressure force | irrelevant (`∇p ≈ 0`) | leave default |
| Coupling | `ELApplyParticleForces = Yes`, `ELApplyFluidFeedback = No` | one-way |
| BCs | periodic; fluid `u = 0` IC, no forcing | keeps `u_f ≡ 0` |
| Time step | `dt = 5e-4` | resolves the transient (~11 steps per τ_p) |
| Run length | ~100 steps to `t ≈ 0.05` (> 10 τ_p) | particle reaches terminal well before this |

**Mesh note.** The mesh resolution is **not critical** for this one-way case: with
`u_f ≡ 0` the grid only samples a trivial zero field, so the result is set by the
force closure + grav_buoy + semi-implicit integrator, not by `h`. The analytic
targets below depend only on `r, ρ_p, ρ_f, μ, g` — **they do not change with the
mesh**. `18³` (h = 1/18) is used because it is convenient to generate for the fully
periodic box; `16³` or other uniform resolutions give the same `u_t`/`τ_p`. Keep
`r = 0.02` (do not re-tune it to the mesh — `u_t ∝ r²` and `τ_p ∝ r²`, so changing
`r` would shift the baselines). `dt` is set by `τ_p`, not the mesh (no CFL
constraint, since the fluid is at rest). **Caveat:** this mesh-insensitivity is
specific to one-way coupling; the optional two-way variant (below) genuinely
depends on `h` and should be checked for mesh convergence.

## Analytic targets (the baseline)

Stokes regime (`ε_f = 1`, dilute single particle). Net submerged weight balances
Stokes drag at terminal velocity:

```
(ρ_p − ρ_f) · V · g = 6 π μ r · u_t,     V = (4/3) π r³
  ⇒  u_t = (2/9) (ρ_p − ρ_f) g r² / μ
  ⇒  τ_p = m / (6 π μ r) = (2/9) ρ_p r² / μ
```

With the configuration above:

| Quantity | Formula | Value |
|---|---|---|
| Terminal velocity `u_t` (z) | `(2/9)(ρ_p−ρ_f) g r²/μ` | **0.03556** |
| Relaxation time `τ_p` | `(2/9) ρ_p r²/μ` | **4.44e-3** |
| Drag coefficient `B = 6πμr` | — | 0.03770 |
| Particle mass `m = ρ_p V` | — | 1.676e-4 |
| Particle Reynolds `Re_p = ρ_f d u_t/μ` | — | **0.0142** (≪ 1 ⇒ Stokes exact) |

Transient (from rest): `u_p,z(t) = u_t · (1 − e^{−t/τ_p})`.

Sanity: the particle reaches terminal within ≈ 5 τ_p ≈ 0.022, having fallen
`< 1e-3` (≪ one cell) — so the run is short and the path needs no resolution.

## Acceptance criteria

- **Terminal velocity:** steady `u_p,z` within **1%** of `u_t = 0.03556`.
- **Relaxation:** fit of `u_p,z(t)` to `u_t(1 − e^{−t/τ_p})` recovers `τ_p` within
  a few % (or check `u_p,z(τ_p) ≈ 0.632 u_t` and `u_p,z(3τ_p) ≈ 0.950 u_t`).
- **Regime guard:** assert `Re_p < 0.1` at terminal (confirms Stokes-exact).

## Diagnostic / harness wiring

Mirror the existing `SED_BENCH_VEL` pattern that `tools/featflower_test` already
consumes. Print one keyword line per step on the representative rank (gate like the
existing `el_write_diagnostics` / `showid` block in `el_transfer.f90`):

```
EL_TERMINAL_VEL time= <t> ip= <idx> <vx> <vy> <vz>
```

`featflower_test` YAML (model on
`tools/featflower_test/testcases/definitions/q2p1_bench_sedimentation_pe_serial.yaml`):
`keyword_columns` parser on `EL_TERMINAL_VEL`, `numeric_format: fortran_d_or_e`,
extract `vz`, `compare: tolerance` against the analytic baseline `u_t` (and
checkpoints at `t = τ_p, 3τ_p, 5τ_p` for the transient). Baseline values are the
analytic numbers above, not a captured run.

## Responsibilities

- **User provides:** the unit-cube structured hex mesh (`.tri` / `GRID` for the
  partitioner), 16³ or similar, periodic-compatible. This is the only mesh-side
  action.
- **Implementation (remote agent):**
  - PE setup creating **one mobile sphere** at a known interior position with zero
    initial velocity. NB: the frozen-trace `example.json` places a *fixed sphere
    array* — this case needs a single *free* sphere, so the PE setup needs a small
    tweak (or a dedicated setup variant).
  - `example.json` + `q2p1_param.dat` keys as tabulated above.
  - The `EL_TERMINAL_VEL` per-step printout.
  - The `featflower_test` YAML definition + analytic baseline.

## Optional two-way variant (follow-up)

To exercise the feedback while still reaching a steady state, add a constant fluid
body force equal and opposite to the particle's weight per unit domain volume (via
the existing `UseConstantForcing` / `ConstantForcing`), so mixture momentum stays
bounded and the particle settles at `u_t` measured as the steady **slip**
`u_p − ⟨u_f⟩`. Same mesh and fluid/particle parameters; `ELApplyFluidFeedback = Yes`
and the balancing force added. Keep as a follow-up; the one-way case above is the
primary check.
