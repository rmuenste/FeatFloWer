# FeatFloWer Guide 04: Configure and Build `q2p1_bench_sedimentation` with PE Parallel Mode

This guide is the parallel-mode counterpart to Guide 02.
It runs the **same** sphere sedimentation benchmark, but with the PE rigid-body library in its
standard **MPI domain-decomposition mode** instead of serial mode.

Scope of this guide:

- load the required RHEL 9.7 module environment
- configure FeatFloWer with PE enabled in **parallel mode** (no `USE_PE_SERIAL_MODE`)
- build `q2p1_bench_sedimentation`
- install the **Cartesian 1×1×12** partitioned mesh (this is the critical difference)
- stage the PE runtime JSON config with a matching processor grid
- run on 13 MPI ranks and validate against the serial reference
- keep logs in files (avoid terminal/context flooding)

> **Read Guide 02 first.** This guide only spells out what *differs* from serial mode in
> full detail. Everything else (JSON/Eigen FetchContent, `MG.dat`, physical parameters,
> benchmark background) is identical and is covered there.

---

## 1) The One Thing That Makes Parallel Mode Different

In serial mode every CFD rank runs its own independent, complete PE instance, so PE never
needs to know how the CFD mesh was split. **In parallel mode it does.** PE builds an MPI
Cartesian decomposition and each PE process owns a rectangular box; that box must coincide
with the CFD subdomain of the same rank. Two consequences follow, and both bite immediately.

### 1a) The mesh directory layout is different

`app_init.f90` includes a *different* partition reader per mode:

| Mode | Reader | Layout expected |
|---|---|---|
| Serial (`PE_SERIAL_MODE`) | `source/include/PartitionReader2.f90` | `_mesh/<Mesh>/sub0001/GRID<myid>.tri` — **one** subdir, N grids |
| **Parallel** | `source/include/PartitionReader.f90` | `_mesh/<Mesh>/sub<myid>/GRID0001.tri` — **N** subdirs, one grid each |

The parallel reader hardcodes `iPart = 1; iSubpart = myid`. So Guide 02's METIS output is not
merely the wrong partition *count* here — it is the wrong directory *shape*, and reusing it
fails at mesh load.

### 1b) METIS cannot produce the required partitioning

`commf2c` dispatches to `setupParticleBench()` (`libs/pe/pe/interface/setup_part_bench.h`),
which builds the decomposition as:

```cpp
const real dz( 0.16 / config.getProcessesZ() );
decomposeDomain(center, 0.0, 0.0, 0.0, dx, dy, dz,
                config.getPx(), config.getPy(), config.getPz());
```

`decomposeDomain` gives the process at Cartesian coordinate `center` the box
`[b + c·d, b + (c+1)·d]` on each axis. With `MPI_Cart_create(dims={1,1,12}, reorder=false)`,
rank *r* owns z-slab *r*, ascending.

METIS is a **graph** partitioner: it produces arbitrarily shaped, arbitrarily ordered
subdomains. It cannot express ordered, axis-aligned slabs, so neither `tools/PyPartitioner.py`
nor `featflower-partition` can generate a valid mesh for this mode. You need a **Cartesian**
partitioning (see section 5).

---

## 2) Why the Quarter Domain and `bx = by = 0.0` Are Correct

This trips people up when they compare the mesh bounds against the PE domain, so it is worth
stating explicitly.

The CFD mesh is a **quarter domain**: `x ∈ [-0.05, 0]`, `y ∈ [-0.05, 0]`, `z ∈ [0, 0.16]`,
with symmetry planes at `x = 0` and `y = 0`. The sphere sits exactly on the symmetry axis at
`(0, 0, 0.1275)`.

The CFD therefore integrates the hydrodynamic force over only a quarter of the sphere surface.
`SED_BENCH` reconstructs the full-sphere force in `source/src_quadLS/QuadSc_force_extension.f90`:

```fortran
DResForceX = 0.0
DResForceY = 0.0
DResForceZ = 4.0 * DResForceZ
theParticles(ip)%torque(:) = (/0.0, 0.0, 0.0/)
```

The `4.0` accounts for the four quadrants; the transverse components and the torque are zeroed
because they cancel identically for a sphere on the symmetry axis.

So **PE receives the full physical force and carries a full sphere** (full `benchRadius_`, full
mass). PE is not modelling the quarter — it models the real sphere. That is why
`decomposeDomain(center, 0.0, 0.0, 0.0, ...)` is correct even though PE's x/y quadrant is the
mirror image of the mesh's:

- with `px = py = 1` nothing is subdivided in x or y, so the x/y extent partitions nothing;
- the sphere lies on the corner line `x = y = 0`, shared by either quadrant box;
- only the **z** slabs correspond to a real decomposition, and those match exactly.

**Do not "fix" the origin to `-0.05`.** It is not a bug.

A useful corollary: `v_x` and `v_y` are zeroed by construction and never computed, so they
print as exactly `0.000000E+00` for the whole run. That is expected, not a symmetry artefact
of the solver.

---

## 3) Environment Setup on RHEL 9.7 (module-based)

```bash
source /etc/profile.d/modules.sh
module purge
module load gcc/14.3.0
module load openmpi/4.1.6
module list
```

> **Load the modules as separate `module load` commands.** Combining them
> (`module load gcc/14.3.0 openmpi/4.1.6`) returns exit code 0 but silently leaves the MPI
> wrappers unset, because the OpenMPI modulefile resolves its compiler-specific prefix from
> the GCC module at evaluation time. Verify before configuring:

```bash
which mpicc mpicxx mpifort
```

**Agent/automation form** — each command runs in a fresh shell, so include the loads every time:

```bash
bash -c '
  source /etc/profile.d/modules.sh
  module purge
  module load gcc/14.3.0
  module load openmpi/4.1.6
  # actual command here
'
```

---

## 4) Configure for PE Parallel Mode

From the repository root:

```bash
cmake -S . -B build-pe-parallel -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_APPLICATIONS=ON \
  -DUSE_PE=ON \
  -DPE_USE_JSON=ON \
  -DSED_BENCH=ON \
  -DPE_VERIFY_HASHGRID=OFF \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpifort
```

**The critical difference from Guide 02: there is no `-DUSE_PE_SERIAL_MODE=ON`.** Everything
else is the same. `-DPE_USE_JSON=ON` remains **mandatory** — the processor grid is read from
`example.json`, and without JSON support it silently stays at its default.

Expected configure banner:

- `Adding pe library to build`
- `SED_BENCH: sedimentation benchmark output enabled`
- `nlohmann/json will be downloaded to: ...`
- **no** `PE SERIAL MODE ENABLED` line

Verify the cache:

```bash
grep -E "^USE_PE|^SED_BENCH|^PE_USE_JSON" build-pe-parallel/CMakeCache.txt
```

Expected:

```
PE_USE_JSON:BOOL=ON
SED_BENCH:BOOL=ON
USE_PE:BOOL=ON
USE_PE_SERIAL_MODE:BOOL=OFF
```

### `HAVE_MPI` is the real switch

Serial mode builds the PE library **without** MPI. Confirm the parallel build did the opposite:

```bash
grep "HAVE_MPI" build-pe-parallel/libs/pe/config.h
# Expected: #define HAVE_MPI 1     (serial mode gives 0)
```

This matters more than it looks: large parts of `libs/pe/src/interface/sim_setup.cpp` sit
inside `#if HAVE_MPI`, so **parallel mode compiles code paths that serial mode never touches.**
A tree that builds cleanly in serial mode can still fail to build in parallel mode.

---

## 5) Build the Application

```bash
cmake --build build-pe-parallel --target q2p1_bench_sedimentation -- -j8 \
  > build_q2p1_bench_sed_parallel.log 2>&1
```

Failure scan:

```bash
rg -n "error:|undefined reference|FAILED:|ninja: build stopped" build_q2p1_bench_sed_parallel.log
```

Verify:

```bash
ls -lh build-pe-parallel/applications/q2p1_bench_sedimentation/q2p1_bench_sedimentation
# Expected: ~50 MB executable
```

> **Note on `-Wfatal-errors`.** The PE build uses `-Wfatal-errors`, so each translation unit
> reports only its **first** error. If a build fails on a missing collision-system or body
> method, fixing it may simply reveal the next one. Do not assume a single error is the whole
> story — rebuild and re-scan until clean.

---

## 6) Install the Cartesian 1×1×12 Mesh

This replaces Guide 02's METIS partitioning step entirely.

### 6a) Choosing the subdivision

`1 × 1 × 12` is used because it divides the coarse mesh cleanly. The coarse mesh
(`benchSym/mesh.tri`, 876 elements / 1239 vertices) has **36 uniform element layers** in
`z ∈ [0, 0.16]`. Since `12 | 36`, each of the 12 slabs is exactly 3 element layers, and every
slab boundary falls on an existing mesh plane:

```
slab width = 0.16 / 12 = 0.0133333   (identical to PE's dz = 0.16 / getProcessesZ())
cuts at z = 0, 0.01333, 0.02667, ..., 0.16
```

Pick any `pz` that divides the z-layer count. `1×1×12` also keeps `px = py = 1`, which avoids
subdividing the quarter-domain symmetry directions.

### 6b) Use the pre-partitioned mesh

> **Known gap — the mesh is not in the repository yet.** The validated run used a ready-made
> partitioning at `benchSym/mesh12/NEWFAC`, which is currently **untracked**: meshes are
> deliberately kept out of the repo. Until that is resolved you must either obtain
> `benchSym/mesh12/` separately or generate an equivalent mesh yourself (section 6d). A
> recipe for building these Cartesian partitionings with the maintained FeatFloWer
> partitioner is planned and will replace section 6d.

The expected layout is:

```
benchSym/
├── mesh.tri, bench.prj, top.par, bot.par, xwall.par, ywall.par, x.par, y.par   # source case
└── mesh12/NEWFAC/
    ├── GRID.tri, GRID.prj          # coarse mesh, read by the master rank
    ├── <boundary>.par
    └── sub0001/ ... sub0012/       # one directory per CFD worker rank
        ├── GRID0001.tri            # that rank's z-slab
        └── <boundary>_0001.par
```

Install it over whatever the staging target produced:

```bash
cmake --build build-pe-parallel --target q2p1_bench_sedimentation_stage

cd build-pe-parallel/applications/q2p1_bench_sedimentation/_mesh
mv NEWFAC NEWFAC_metis_unused          # staging emits the serial METIS layout
cp -r /path/to/FeatFloWer/benchSym/mesh12/NEWFAC ./NEWFAC
```

> The staging target still runs the METIS partitioner and writes the **serial** layout. That
> output is unusable here; run staging for the other assets it provides (`_data/MG.dat`,
> `example.json`, `_adc/`, `libmetis.so`) and then overwrite `_mesh/NEWFAC`.

### 6c) Verify the mesh

```bash
cd build-pe-parallel/applications/q2p1_bench_sedimentation

ls -d _mesh/NEWFAC/sub* | wc -l          # Expected: 12
ls _mesh/NEWFAC/sub0001/GRID0001.tri     # Expected: exists
ls _mesh/NEWFAC/sub0012/GRID0001.tri     # Expected: exists
ls _mesh/NEWFAC/GRID.tri                 # Expected: exists (master rank)
```

Stronger check — the slabs must be a complete, non-overlapping, **ascending** cover:

```bash
python3 - <<'EOF'
def hdr(f):
    L=open(f).read().split("\n"); p=L[2].split(); return int(p[0]),int(p[1])
tot=sum(hdr(f"_mesh/NEWFAC/sub{i:04d}/GRID0001.tri")[0] for i in range(1,13))
cn,_=hdr("_mesh/NEWFAC/GRID.tri")
print(f"sum of slabs NEL={tot}  coarse NEL={cn}  complete cover={tot==cn}")
EOF
# Expected: sum of slabs NEL=876  coarse NEL=876  complete cover=True
```

### 6d) If you need to regenerate it

The maintained `featflower-partition` (METIS 5) cannot do this. Use the legacy METIS-4
`tools/partpy` variant, whose method `-4` performs uniform Cartesian splitting from the
bounding box and numbers the parts in ascending order:

```bash
python3 tools/partpy/PyPartitioner.py 12 -4 x1-y1-z12 NEWFAC _adc/benchSym/bench.prj
```

`NPart` must equal the product of the axis counts (`1·1·12 = 12`). Note `tools/partpy` links
**METIS 4** (`extern/libraries/metis-4.0.3`), not the METIS 5 used elsewhere.

---

## 7) Stage the PE JSON Config

The processor grid must match the mesh **and** the rank count. Edit `example.json` in the
application run directory:

```json
{
    "processesX_": 1,
    "processesY_": 1,
    "processesZ_": 12,

    "timesteps_": 1300,
    "stepsize_": 0.001,
    "benchRadius_": 0.0075,
    "particleDensity_": 1120.0,
    "fluidViscosity_": 0.058,
    "fluidDensity_": 960.0,
    "gravity_": [0.0, 0.0, -9.81]
}
```

The physical parameters are unchanged from Guide 02 (Nylon sphere in silicone oil). Only the
`processes*_` triple is parallel-specific — the stock `example.json` ships `3, 3, 3`, which
will abort this run.

Verify:

```bash
grep -E "processesX_|processesY_|processesZ_" example.json
```

---

## 8) Rank Count: Why 13 and Not 12

`setupParticleBench` asserts that the processor grid matches the PE communicator:

```cpp
if( config.getProcessesX() * config.getProcessesY() * config.getProcessesZ()
    != mpisystem->getSize() ) { /* abort */ }
```

`mpisystem` is set to `ex0`, the communicator that **excludes CFD rank 0** (the master, which
holds the coarse problem and owns no PE subdomain). So:

```
total MPI ranks = px·py·pz + 1 = 1·1·12 + 1 = 13
```

The rank mapping chains through cleanly, and you can read it straight out of the log:

```
CFD rank n  →  ex0 rank n-1  →  Cart coords (0,0,n-1)  →  z-slab n-1  →  _mesh/NEWFAC/sub{n}
```

---

## 9) Run

### Direct run

```bash
cd build-pe-parallel/applications/q2p1_bench_sedimentation
mpirun -np 13 ./q2p1_bench_sedimentation > run_sed_parallel_np13.log 2>&1
```

### SLURM submission

Create `run_sed_parallel.sh` in the application directory:

```bash
#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=13
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --job-name=q2p1_sed_par

source /etc/profile.d/modules.sh
module purge
module load gcc/14.3.0
module load openmpi/4.1.6

cd "$SLURM_SUBMIT_DIR"

echo "=== Job started at $(date) ==="
ls -lh q2p1_bench_sedimentation example.json _data/MG.dat
ls -d _mesh/NEWFAC/sub* | wc -l

mpirun -np 13 ./q2p1_bench_sedimentation

echo "=== Job completed at $(date) ==="
```

```bash
sbatch run_sed_parallel.sh
```

---

## 10) Verifying the Run

### Startup markers (check these first)

```bash
grep -E "Configuration Sed bench|3D communicator|3D coordinates|Total number of MPI" run_sed_parallel_np13.log
```

Expected:

```
1> C) Configuration Sed bench with 13 processes.
> 3D communicator created
3D coordinates were created
 Total number of MPI processes           = 12
```

Note both numbers: **13** total ranks, **12** PE processes. If you instead see

```
 Invalid number of MPI processes: 12!=27
```

the `processes*_` triple in `example.json` does not match your rank count.

### Rank → subdomain mapping

```bash
grep "is assigned to the mesh" run_sed_parallel_np13.log | sort -k2 -n
```

Expected — rank *n* to `sub{n}`, all 12 present:

```
Rank:    1 is assigned to the mesh:_mesh/NEWFAC/sub0001/GRID0001.tri
...
Rank:   12 is assigned to the mesh:_mesh/NEWFAC/sub0012/GRID0001.tri
```

### Particle trajectory — different output channel than serial

**This is a real observability difference.** `SED_BENCH_VEL` lines are emitted only by
`source/src_quadLS/QuadSc_force_serial.f90`, the serial-mode force path. In parallel mode the
trajectory comes from PE's own printer instead:

```
Position: <systemID> <z>  <time>
Velocity: <systemID> <v_z> <time>
```

Extract it with:

```bash
grep "^Velocity:" run_sed_parallel_np13.log | tail -5
```

Do **not** conclude the run failed just because `grep SED_BENCH_VEL` comes back empty — that
is expected here. `Force with SED BENCH settings!` is still printed, confirming the ×4
quarter-domain force path is active:

```bash
grep -c "Force with SED BENCH settings" run_sed_parallel_np13.log
```

---

## 11) Expected Results and Serial Cross-Validation

A completed run advances all 1300 steps to `t = 1.30` and ends with
`PP3D_LES has successfully finished.`

| Quantity | Parallel (1×1×12, np 13) | Serial (np 32) | Publication |
|---|---|---|---|
| Terminal velocity | **−0.132877 m/s** @ t=0.924 | −0.132921 m/s @ t=0.892 | 0.128 m/s |
| Deviation vs serial | **0.033%** | — | +3.8% |
| Wall impact | t ≈ 1.09 | t ≈ 1.09 | ~1.2 s settling |

Agreement across the full trajectory (1299 common sample times):

| Window | Samples | max abs deviation | mean relative |
|---|---|---|---|
| Settling phase, t ∈ [0.05, 1.05] | 1001 | **8.4e-04 m/s** | **0.227%** |
| Impact/rebound, t ∈ [1.05, 1.31] | 250 | 3.9e-02 m/s | 10.4% |

The settling phase — the physically meaningful part, and the part the benchmark is about —
matches serial to under 1e-3 m/s. The larger spread after impact is expected and **not** a
defect: a sub-millisecond difference in wall-contact timing produces a large instantaneous
velocity difference during the bounce, while the pre-impact trajectory and the terminal
velocity are unaffected.

Regression criterion: **terminal velocity within ~0.1% of the serial reference.**

---

## 12) Troubleshooting

**Issue:** `Invalid number of MPI processes: N!=M`
- **Cause:** `processesX_·processesY_·processesZ_` ≠ (total ranks − 1)
- **Fix:** for `1×1×12`, run with `-np 13`. Remember `ex0` excludes the master rank.

**Issue:** Mesh load fails, or ranks read the wrong/missing grid file
- **Cause:** serial (METIS) mesh layout in `_mesh/NEWFAC` — one `sub0001/` with many `GRID*.tri`
- **Fix:** install the Cartesian mesh (section 6). Parallel mode needs `sub0001…sub0012`, each
  with a single `GRID0001.tri`.

**Issue:** Particles leak across subdomains, or contacts are missed at slab boundaries
- **Cause:** PE's z-cuts do not coincide with the CFD mesh partition boundaries
- **Fix:** ensure `pz` divides the z element-layer count (36 here) and that the mesh was cut
  uniformly on the bounding box in ascending z. Re-run the section 6c cover check.

**Issue:** `grep SED_BENCH_VEL` returns nothing
- **Cause:** not a failure — that output exists only in the serial force path
- **Fix:** use `grep "^Velocity:"` instead (section 10)

**Issue:** Build fails on a missing collision-system/body method, e.g.
`has no member named 'resetElLubricationVirial'`
- **Cause:** the coupling interface reaching an API that the selected `pe_CONSTRAINT_SOLVER`
  does not implement, in a path guarded by `#if HAVE_MPI` (so serial builds never see it)
- **Fix:** guard the call with the `std::void_t` detector + `if constexpr` pattern in
  `libs/pe/pe/interface/el_optional_api.h`. Because of `-Wfatal-errors`, rebuild after each
  fix — the next missing method only appears once the previous one is resolved.

**Issue:** All parameters at defaults despite `example.json`
- **Cause:** PE built without JSON support
- **Fix:** reconfigure with `-DPE_USE_JSON=ON`; check `HAVE_JSON 1` in `libs/pe/config.h`

---

## 13) Guide Series Context

- Guide 01: baseline `q2p1_fc_ext` cylinder benchmark from scratch
- Guide 02: `q2p1_bench_sedimentation` with **PE serial mode**
- Guide 03: `q2p1_sse` workflow for SSE/TSE execution
- Guide 04 (this guide): `q2p1_bench_sedimentation` with **PE parallel mode** and a Cartesian
  1×1×12 decomposition

Serial vs parallel at a glance:

| | Guide 02 (serial) | Guide 04 (parallel) |
|---|---|---|
| CMake | `+ -DUSE_PE_SERIAL_MODE=ON` | (omit it) |
| `HAVE_MPI` in PE | `0` | `1` |
| Partition reader | `PartitionReader2.f90` | `PartitionReader.f90` |
| Mesh layout | `sub0001/GRID<myid>.tri` | `sub<myid>/GRID0001.tri` |
| Partitioner | METIS (`PyPartitioner.py`) | Cartesian (`partpy -4`) |
| Decomposition | none (PE replicated) | MPI Cartesian, must match CFD |
| Ranks | `-np 32` (31 parts + master) | `-np 13` (12 parts + master) |
| Trajectory output | `SED_BENCH_VEL` | `Position:` / `Velocity:` |
