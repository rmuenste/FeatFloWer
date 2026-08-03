# ATC slug training case — RPM 40 (D6-ATC rehearsal)

> Portable staging instructions (for fresh clusters / other agents):
> `docs/md_docs/guide_05_q2p1_atc_case_staging.md`.

Staged 2026-08-02 from the old FF-ATC-NEW setup, to be run with the new
`build-atc-ninja-release-eigen` binary (branch `feature/dns-validation-phase1`,
head e4bd41b3; guide 04 configuration: PE serial + CGAL + JSON + Eigen +
FBM/KVEL). Purpose: training for the DNS campaign D6-ATC capstone —
suspension regime vs RPM in the Archimedes tube crystallizer slug.

## Provenance

| Piece | Source | Notes |
|---|---|---|
| `_data/q2p1_param.dat` | `FF-ATC-NEW/build/applications/q2p1_ATC/_data/` | verbatim: RPM=40, dt=0.25 ms, levels 1–4, CN, rotating BC `Inflow771` |
| `_data/MG.dat` | same | READ by OutputProfiles (multigrid ke/jv map) — input, not output |
| `atc_35_mesh/` | same | ProjectFile: coarse mesh + 3 region .par (Slip / Inflow771 / Slip) |
| `CASE_090_771/` | same | `PrectFile` target; key is NOT parsed by current param_parser — kept for fidelity |
| `_mesh/NEWFAC/` | same | partitioned tree, sub0001 = 63 subdomains → np = 64 |
| `sorted_vertices_by_x_world.txt` | same (rundir root) | REQUIRED by new `setupATCSerial` (centerline for particle seeding + stuck-particle diagnostics) |
| `atc_boundary_param_zero.obj` | same | PE domain-boundary triangle mesh (slug wall for particle containment) |
| `example.json` | MERGED — see deviations | old-case physics + new-schema blocks |
| `start/` (sampleRigidBody.xml, data.TXT) | same | REQUIRED: FullC0ntact loads the XML unconditionally (crash 137403 without it). nBodies=0 → its RigidBody entry (missing `CASES/geo090/FD_090.off`) is skipped; only the BoundaryDescription is read |
| `OFF/` (atc.off, xp_sphere.off, xm_sphere.off, +3 unused) | same | REQUIRED: the three BoundaryShape meshes referenced by sampleRigidBody.xml (slug wall + sphere caps) |
| `q2p1_ATC` | symlink → `../build-atc-ninja-release-eigen/...` | new binary, git stamp acebb68f |

Production reference: `FF-ATC-NEW/d273_simulations_and_vids/CASE_20260131090504_atc_r60_config_3/`
(same physics; differs: RPM=25, dt=0.5 ms, tilted gravity (-848.7,-490,0) = 60° tilt,
`packingMethod_="grid"`).

## Deviations from the old staged rundir (all deliberate)

1. **`stepsize_` 0.0005 → 0.00025** — synced to CFD `SimPar@TimeStep`. The old
   pair was a PE/CFD timestep mismatch, exactly the `pe_stepsize_mismatch`
   defect class found in the DNS campaign; the new binary's fatal first-step
   guard (`fbm_updateDefaultFC2`) would abort on the old value.
2. **`packingMethod_` "none" → "grid"** — the old staged rundir seeded NO
   particles (its `spheres_*.txt` outputs are empty). "grid" activates the
   centerline seeding in the new `setupATCSerial` (~7000 spheres, r=0.0137,
   as in the Jan/Feb 2026 production batches). Set back to "none" for a
   fluid/mixer-only twin.
3. **`domainBoundary_` block added** (new schema, from FF-ATC-NEW root
   `example.json` template): `atc_boundary_param_zero.obj` at
   (0.000006, -2.25096, 0.08719), distance map res=128 tol=6 inverted,
   `writeVti` off. The old json predates this schema; without the block the
   new setup runs with NO particle containment mesh.
4. **`vtk_` false** (was true) — matches production; avoids per-rank PE
   paraview output in serial mode. CFD VTK output is unaffected.
5. Old `xyzFilePath_` "sphere_centers_50.xyz" kept but inert (only read when
   `packingMethod_="external"`; the file does not exist anywhere in FF-ATC-NEW).
6. Slurm script: explicit module loads (gcc 13 + OpenMPI 4.1.6) instead of
   `~/.bashrc` aliases; otherwise mirrors old `atc.sh` (1 node × 64 tasks).

7. **`fluidViscosity_` 0.058 → 0.01** (2026-08-02, user-confirmed): the case
   is CGS and the fluid is ultra-pure water — `Prop@Viscosity = 0.01` cm²/s
   is water's kinematic viscosity at 20 °C; PE takes DYNAMIC viscosity
   (μ = ρ·ν = 1·0.01 = 0.01 P). The old 0.058 was stale (≈5.8× water) and
   would have rescaled PE lubrication once particles are enabled.
8. **`StartingProc` 3 → 0** — mode 3 = restart from a repartitioned solution
   and reads it UNCONDITIONALLY (`init_sol_repart` → EOF crash, job 137404);
   there is no fresh-start fallback for modes 1–3
   (`app_initialization.f90:55-70`: 0 fresh, 1 same-level, 2 lower-level,
   3 repartition). The old rundir was resuming mid-campaign; a clean rundir
   must use 0.

## DNS particle resolution (measured from atc_35_mesh, MaxMeshLevel=4)

Coarse hex edges 107–976 µm (mean 688 µm), strongly graded — finest cells in
the near-wall layer. `benchRadius_=0.0137` cm = the 273 µm production
particle (d273 batches); 365 µm ⇒ r=0.01825 cm.

| Level | h_mean | h_coarsest | elements | D/h (365 µm) | D/h (273 µm) |
|---|---|---|---|---|---|
| **L4 (case)** | 86 µm | 122 µm | 1.23 M | **4.2** (worst 3.0, near-wall up to 27) | **3.2** (worst 2.2, near-wall up to 20) |
| L5 | 43 µm | 61 µm | 9.83 M | 8.5 | 6.4 |

Q2 velocity (27 DOFs/hex) samples at h/2 ⇒ ~8 (365 µm) / ~6 (273 µm) DOF
spacings per diameter at L4 — well above the FictKNPR transparency floor
(D/h ≲ 1.5, D1.1) and enough for regime/mechanism physics, but NOT force
metrology (campaign guidance: D/h ≳ 8 usable, ~16–24 quantitative — and
those guidelines were measured WITH Q2, so no extra credit). Verdict:
L4 = mechanism-level (matches the D6-ATC gate as written); L5 buys 365 µm
into the usable band at 8× cost; D/h ≈ 3–4 is prime DNS↔EL crossover
territory (D5).

## Gates (probe-first — do NOT go straight to the 48 h run)

Smoke (short partition, ~25 steps: set `SimPar@MaxNumStep=25`, submit with
`--partition=short --time=02:00:00`):
- [ ] PE/CFD dt guard passes (no FATAL at step 1); `ATC SETUP` banner shows
      dt = 0.00025, fresh start, ~7000 particles, volume fraction printed.
- [ ] `DOMAIN BOUNDARY SETUP` banner: distance map 128³ built, escape-detection
      sanity check = PASS on all ranks.
- [ ] Centerline loaded (vertex count > 0, total length > 0).
- [ ] No rank asymmetry hangs; `_data/prot.txt` advances; CFL sane.
- [ ] Rotating BC active: velocity field responds to RPM=40 (Mixer field in VTK).

Full run: restore `MaxNumStep=15000`, submit `atc_slurm.sh` (long, 48 h).
Observable: vertical concentration profile / suspension regime at RPM=40 vs
the prior in-house study (settling ↔ even suspension boundary).

## Job log

| Date | Job ID | What | Verdict |
|---|---|---|---|
| 2026-08-02 | 137403 | smoke v1 (`q2p1_atc_rundir_smoke`, 25 steps, no particles) | FAIL: missing `start/sampleRigidBody.xml` (+ `OFF/`); fileset amended |
| 2026-08-02 | 137404 | smoke v2 | FAIL: `StartingProc=3` restart read on empty `_sol` (`solution_io.f90:1395` EOF); set 0 |
| 2026-08-02 | 137405 | smoke v3 | **PASS**: dt guard OK (PE dt=0.00025 synced), fresh start, 0 particles, domain-boundary distance map 128³ + escape check PASS, 25/25 steps to t=6.25 ms, velocity solver converging (defect ~1e-9), dumps + VTK written, clean exit in 4:22 |
| 2026-08-02 | 137406 | fluid-only spin-up v1 (`q2p1_atc_rundir_fluid_rpm40`, 15000×0.25 ms = 3.75 s) | CANCELLED before start: 3.75 s < spin-up τ≈6 s; also OutputFreq 0.025 s → 150 VTK frames too heavy (space) |
| 2026-08-02 | 137407 | fluid-only spin-up v2: dt=1 ms (PE synced), 18000 steps = 18 s ≈ 3τ, OutputFreq 0.5 s (36 frames) | submitted — target: fully developed Dean/secondary vortices as the particle-injection base state; restart with particles at dt=0.25 ms from its final dump |
