# FeatFloWer Guide 05: Stage and Run the ATC Slug Case (`q2p1_ATC`)

Smoke-validated 2026-08-02 (job 137405, 25/25 steps, clean exit) on the
`feature/dns-validation-phase1` branch with the guide-04 build
(PE serial + CGAL + JSON + Eigen + FBM/KVEL). This guide exists so any agent
can stage the case from scratch on ANY cluster — e.g. when the home queues
are full. It is the training case for the DNS-campaign D6-ATC capstone
(suspension regime vs RPM in the Archimedes tube crystallizer).

Reference rundirs (this repo, warehouse17):
- `q2p1_atc_rundir_smoke/` — validated 25-step fluid-only smoke (PASS).
- `q2p1_atc_rundir_rpm40/` — full case incl. particles; its `RUNBOOK.md`
  holds provenance, the deviation list, and the job log. READ IT FIRST.
- `q2p1_atc_rundir_fluid_rpm40/` — fluid-only spin-up variant + `compute_Re.py`.
Original (historical) source: `FF-ATC-NEW/build/applications/q2p1_ATC/`.

## 0) The case in one paragraph

A liquid slug (ultra-pure water, CGS units) trapped between two gas menisci
inside one turn of the rotating ATC helix. The mesh is fixed; the drum
rotation enters as a velocity BC on the tube wall (`Inflow771` →
`GetVeloBCVal`, `source/src_quadLS/QuadSc_user.f90`: swirl about the global
z-axis at `SimPar@RPM` + axial advance 0.772 cm/rev). Slip caps = menisci;
`NoOutflow=Yes` → closed volume. Characteristic flow: wall-driven conveyor
cell (forward annulus, core return) + Dean secondary vortices;
Re_wall ≈ 340 (RPM 25) … ≈ 550 (RPM 40), De ≈ 100 — steady laminar
(Ito Re_crit ≈ 9800 at this curvature). Particles: PE serial mode,
centerline-seeded spheres r=0.0137 cm, ρ_p=1.42 g/cm³, contained by a
triangle-mesh domain boundary with CGAL distance map.

## 1) Prerequisites

1. Binary built per `guide_04_q2p1_atc_pe_serial_fbm_kvel_from_scratch.md`
   (NOTE: the old `--target cgal` pre-step is obsolete — build `q2p1_ATC`
   directly). Any GCC13/14 + OpenMPI stack; no other runtime deps
   (binary links METIS statically — no `libmetis.so` needed in the rundir).
2. MPI able to run 64 ranks (1 master + 63 workers; see §4 to change).

## 2) Input manifest (complete — nothing else is read)

```
rundir/
├── _data/
│   ├── q2p1_param.dat        # solver config (see §3 for required edits)
│   └── MG.dat                # multigrid ke/jv map — an INPUT (OutputProfiles reads it)
├── atc_35_mesh/              # ProjectFile: coarse mesh, 2400 el / 2907 vt
│   ├── file.prj  atc_35_mesh.tri
│   └── region_0001.par (Slip)  region_0002.par (Inflow771)  region_0003.par (Slip)
├── CASE_090_771/             # older alt mesh (1600 el); NOT read — see PrectFile trick, §5
├── _mesh/NEWFAC/             # partition tree; sub0001/GRID0001..0063.tri = 63 subdomains
│                             # (NEWFAC/GRID.tri is byte-identical to atc_35_mesh.tri)
├── start/
│   ├── sampleRigidBody.xml   # REQUIRED — FullC0ntact loads it unconditionally.
│   │                         # nBodies=0 => its RigidBody entry (missing FD_090.off) is SKIPPED;
│   │                         # only the 3 BoundaryShape OFF meshes are read
│   └── data.TXT
├── OFF/                      # atc.off, xp_sphere.off, xm_sphere.off (+3 unused)
├── atc_boundary_param_zero.obj   # PE domain-boundary mesh (particle containment)
├── sorted_vertices_by_x_world.txt # centerline (193 pts, L=4.401 cm) — REQUIRED
│                             # by setupATCSerial even in fluid-only runs
├── example.json              # PE serial config — see §3
├── q2p1_ATC                  # binary (symlink or copy)
└── empty dirs: _dump _gmv _ns _vtk _sol/1 spheres checkpoints paraview
```

Geometry hard numbers (measured): bore a=0.2375 cm (d=0.475), centerline
L=4.401 cm, coil radius R_c=2.521 cm about the drum z-axis, V=0.78 cm³;
slug sits at y≈−2.2 (gravity −y). Everything CGS.

## 3) Settings that MUST be right (each one cost us a crashed/wrong job)

| Setting | Value | Why |
|---|---|---|
| `SimPar@StartingProc` | `0` for a fresh rundir | modes 1–3 are restarts and read the solution UNCONDITIONALLY (mode 3 crashed job 137404: EOF in `solution_io.f90`). 0=fresh, 1=same-level, 2=lower-level, 3=repartitioned |
| `example.json stepsize_` | == `SimPar@TimeStep` | serial PE integrates with the JSON dt, NOT the CFD dt; the branch's fatal first-step guard aborts on mismatch. The old setup shipped a 2× mismatch |
| `example.json fluidViscosity_` | `0.01` | DYNAMIC viscosity: fluid is water, ν=0.01 cm²/s (CFD `Prop@Viscosity`, kinematic), ρ=1 → μ=0.01 P. Old value 0.058 was stale and rescales PE lubrication |
| `example.json timesteps_` | == `SimPar@MaxNumStep` | keeps PE checkpointer/vis spacing consistent |
| `packingMethod_` | `"none"` fluid-only / `"grid"` particles | "grid" = centerline seeding (~7000 spheres from `sorted_vertices...txt` + `benchRadius_`); `"external"`+`xyzFilePath_` needs an .xyz that does NOT exist — leave inert |
| `domainBoundary_` block | enabled, `atc_boundary_param_zero.obj`, position `[0.000006,-2.25096,0.08719]`, distanceMap res 128 tol 6 inverted, writeVti false | without it particles have no containment mesh (old JSON schema predates it) |
| `SimPar@OutputFreq` | budget it! | VTK per frame = 63 piece files, ~10 M points at level 4. 0.5 s spacing → 36 frames for an 18 s run. 0.025 s → 720 frames = disk killer |
| dt for FLUID-ONLY spin-up | 1 ms is safe (CFL≈0.7 at level 4, CN) | 0.25 ms is the particle-run dt; use it only when particles are on |

Gravity `[0,-980,0]` = straight-down variant (RPM40 training case). The
production tilted variant used `[-848.7,-490.0,0]` with RPM=25, dt=0.5 ms
(see `FF-ATC-NEW/d273_simulations_and_vids/CASE_*_atc_r60_config_3/`).

## 3.1) DNS particle resolution at this mesh (know before you claim)

Measured from `atc_35_mesh.tri` (coarse edges 107–976 µm, mean 688 µm,
graded fine near the wall). Production particle: d = 273 µm
(`benchRadius_=0.0137` cm — the "d273" batches); alternative size 365 µm.

| Level | h_mean | h_coarsest | elements | D/h (365 µm) | D/h (273 µm) |
|---|---|---|---|---|---|
| **L4 (case default)** | 86 µm | 122 µm | 1.23 M | **4.2** (near-wall up to 27) | **3.2** (near-wall up to 20) |
| L5 | 43 µm | 61 µm | 9.83 M | 8.5 | 6.4 |

Q2 velocity (27 DOFs/hex, spacing h/2) gives ~8 / ~6 DOF spacings per
diameter at L4 — above the FictKNPR transparency floor (D/h ≲ 1.5) and fine
for regime/mechanism observables, but not for per-particle force metrology
(campaign guidance D/h ≳ 8 usable, ~16–24 quantitative, measured with Q2
already). L4 supports the D6-ATC gate as written (regime boundary +
concentration profiles); use L5 only for a single confirmation run
(8× elements, ~10× wall-clock); at these D/h the case is also a natural
DNS↔EL crossover candidate (campaign stage D5).

## 4) Rank count and other clusters

np = (number of `GRID00NN.tri` in `_mesh/NEWFAC/sub0001/`) + 1 = **64**.
The partition tree is data, not config — to run a different np you must
repartition the coarse mesh (see
`docs/md_docs/featflower_partitioner_usage_guide.md`; serial-PE mode has no
CFD/PE decomposition-alignment constraint, unlike parallel PE).
Cluster-agnostic job essentials: load a GCC+OpenMPI stack explicitly in the
job script (never rely on `~/.bashrc` aliases), `cd $SLURM_SUBMIT_DIR`,
`mpirun -np 64 ./q2p1_ATC | tee run_slurm.log`, ~3 GB/rank. Wall-time data
point: ~5.5–7 s/step on one 64-core node (25 steps + setup = 4:22).
Setup includes a 128³ CGAL distance-map build on EVERY rank (~1–2 min).

## 5) Known traps (verified on this codebase)

1. **PrectFile "trick"**: the parser ignores unknown keys silently. Lines
   like `SimPar@PrectFile = 'CASE_090_771/file.prj'` are PARKED meshes —
   you switch by renaming which line says `ProjectFile`. Retired practice;
   do not imitate. Only `ProjectFile` + `MeshFolder`/`SubMeshNumber` count.
2. **Double `iT==771` block** in `GetVeloBCVal` (QuadSc_user.f90:388 and
   :407): the first hardcodes 4 RPM, the second (a re-labeled 773 branch)
   overwrites with `SimPar@RPM` — LAST ONE WINS, so RPM works today, but do
   not reorder these blocks, and audit before trusting a new BC id
   (D6.2 rotation-BC recertification covers this).
3. `spheres_*.txt` in the old rundir are OUTPUTS (empty headers = the old
   staged run had zero particles). Not inputs.
4. `_adc` symlink, `parametrization/spheres.txt`, `libmetis.so`, `main.vtu`
   in the old rundir are NOT needed at runtime.
5. `.gitignore` excludes `*.tri`/`*.par` globally — committing mesh
   fixtures needs `git add -f`.
6. Restart-for-particle-injection plan: run fluid-only to a developed state
   (steadiness gate: ||u(t)−u(t−Δt)||/||u(t)|| < ~1e-3 between frames, see
   `compute_Re.py`), then restart `StartingProc=1` (same level & np) with
   `packingMethod_="grid"` and dt back to 0.25 ms, PE stepsize synced.

## 6) Smoke gates (run these before ANY long submission)

25-step smoke (`MaxNumStep=25`, short queue): expect in `run_slurm.log`:
- `ATC SETUP` banner: correct dt, `fresh start`, expected particle count;
- `DOMAIN BOUNDARY SETUP` + `ESCAPE DETECTION ... Result: PASS`;
- 25/25 steps in `_data/prot.txt` (`time: ... itns: 25/25`), no FATAL;
- `PP3D_LES has successfully finished.` + dumps in `_dump/processor_*/`.
Then check physics with `compute_Re.py` on a late frame: max|u| must equal
Ω·r_wall (2π·RPM/60 × ~2.76 cm) — this catches wrong-RPM/BC regressions and
even identifies which RPM a dataset was run at.
