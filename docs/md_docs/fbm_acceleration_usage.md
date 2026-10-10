# FBM Acceleration Runtime Control

## Quick Reference

Add to your `_data/q2p1_param.dat`:

```
SimPar@UseHashGridAccel = Yes
SimPar@UseKVELAccel = Yes
```

## Parameters

### `SimPar@UseHashGridAccel`

Controls HashGrid-accelerated alpha field computation (DOF solid/fluid classification).

- **Values:** `Yes` | `No`
- **Default:** `Yes` in serial PE builds (`USE_PE_SERIAL_MODE=ON`), `No` in
  parallel PE builds (see [Parallel PE mode](#parallel-pe-mode) below)
- **Effect when `Yes`:** Uses spatial HashGrid for O(1) particle containment queries (after timestep 1)
- **Effect when `No`:** Uses brute-force `verifyAllParticles` (linear search over all particles)
- **Performance impact:** ~100-1000x speedup for many-particle systems

### `SimPar@UseKVELAccel`

Controls KVEL/KEEL/KAAL candidate element acceleration for force integration.

- **Values:** `Yes` (default) | `No`
- **Effect when `Yes`:** Only integrates over elements adjacent to particle DOFs
- **Effect when `No`:** Integrates over all mesh elements (brute-force)
- **Performance impact:** ~10,000-20,000x reduction in element evaluations for many-particle systems
- **Scope:** serial PE mode only; the parallel PE force path always integrates
  over all elements and no KVEL cache is built there. Details:
  `kvel_force_acceleration.md`.

## Parallel PE Mode

In parallel PE builds (`USE_PE=ON`, `USE_PE_SERIAL_MODE=OFF`) the HashGrid
query is **off by default** and must be requested explicitly with
`SimPar@UseHashGridAccel = Yes`. Reason: the accelerated query
(`pointInsideParticlesAccelerated` in PE) only sees bodies that are stored in
the HashGrids data structure. PE's `simulationStep()` runs `findContacts()`
first and `synchronize()` last; shadow copies created in that
`synchronize()` (a particle entering a rank's halo) are added to the coarse
detector's `bodiesToAdd_` list and only inserted into the grid at the next
`findContacts()`. `HashGrids::getBodiesNearPoint` does not scan
`bodiesToAdd_`, so during the following fluid step the DOFs of such a particle
are classified as fluid on that rank, while the baseline linear search
(`pointInsideParticles`, local bodies then shadow copies) finds them. The
HashGrid path has not been verified in parallel PE mode.

With the default, a parallel PE build with `ENABLE_FBM_ACCELERATION=ON`
classifies DOFs exactly like a build without acceleration.

## Known Limitations of the HashGrid Query (All Modes)

These apply to serial PE mode as well and are not fixed:

1. **Bodies added between PE steps are invisible.** Any body added to PE after
   the last `findContacts()` (e.g. particle insertion during the run) waits in
   `bodiesToAdd_` and is not returned by the accelerated query until the next
   PE step. The first timestep is safe: the baseline is used until the
   collision pipeline has run once (`collision_pipeline_initialized`).
2. **Stale hashing after integration.** Bodies are hashed by their AABB in
   `findContacts()`, before the positions are integrated. The query visits the
   query point's cell plus its 26 neighbours, which covers a body only while
   its displacement since the last hashing plus its AABB size stays below the
   cell span of its grid level. Large displacements per PE step, or bodies
   whose size is close to their grid's cell span (polydisperse systems), can
   be missed.
3. **Overlapping bodies.** If a point lies inside two bodies, the HashGrid and
   the baseline may return different owners (first hit in different
   iteration orders), so the system id stored in `FictKNPR_uint64` can differ.

Use `-DPE_VERIFY_HASHGRID=ON` (`hashgrid_verification.md`) to measure
mismatches for a given case.

## Build Requirements

`ENABLE_FBM_ACCELERATION` defaults to `ON` whenever `USE_PE=ON` and is always
`OFF` with `USE_PE=OFF`. A serial PE build therefore compiles both
accelerations without extra flags:
```bash
cmake -S . -B build \
  -DUSE_PE=ON \
  -DUSE_PE_SERIAL_MODE=ON
```

- Switch it off explicitly with `-DENABLE_FBM_ACCELERATION=OFF` (pure baseline
  build, see below).
- `ENABLE_FBM_ACCELERATION=ON` also forces PE's `PE_USE_ACCELERATED_POINT_QUERY=ON`.
- **Existing build trees keep their cached value.** Trees configured before the
  default changed have `ENABLE_FBM_ACCELERATION:BOOL=OFF` in `CMakeCache.txt`
  and stay unaccelerated; CMake prints a hint. Reconfigure with
  `-DENABLE_FBM_ACCELERATION=ON` or use a fresh build directory.
- With `USE_PE=OFF` the cache entry (if any) is hidden as
  `ENABLE_FBM_ACCELERATION:INTERNAL` and the effective value is `OFF`.

## Testing Workflow

### 1. Correctness Verification (Same Binary)

Build once with acceleration compiled in:
```bash
cmake -S . -B build \
  -DUSE_PE=ON \
  -DUSE_PE_SERIAL_MODE=ON
cmake --build build -- -j8
```

**Run A — Accelerated:**
```
SimPar@UseHashGridAccel = Yes
SimPar@UseKVELAccel = Yes
```

**Run B — Brute-force:**
```
SimPar@UseHashGridAccel = No
SimPar@UseKVELAccel = No
```

Compare forces/positions/velocities. They agree to round-off, not bitwise:
KVEL sums the element contributions in a different order than the
brute-force loop; the HashGrid may differ from the baseline in
the cases listed under *Known Limitations*.

### 2. Final Verification (Pure Baseline Build)

Rebuild without acceleration code:
```bash
cmake -S . -B build-baseline \
  -DUSE_PE=ON \
  -DUSE_PE_SERIAL_MODE=ON \
  -DENABLE_FBM_ACCELERATION=OFF
cmake --build build-baseline -- -j8
```

The runtime flags have no effect. Compare this pure baseline against the accelerated run.

## Output

When enabled (serial PE mode), rank 0 prints one line per force evaluation:

```
KVEL: 1.3932E+05 candidates vs 3.0818E+09 brute-force (22121.0x speedup)
```

The parameter echo at startup shows `UseHashGridAccel = ON|OFF` and
`UseKVELAccel = ON|OFF`.

When disabled via parameter file:
- No cache building
- No candidate building
- Full element loop for every particle

When disabled at compile time (`-DENABLE_FBM_ACCELERATION=OFF`):
- No log lines
- Acceleration code never compiled
- Pure baseline behavior guaranteed

## Performance Examples

### Sphere Sedimentation (Single Particle, 18432 elements/rank, 63 ranks)
- **Brute-force:** 18432 × 1 particle = 18,432 element evaluations/rank
- **KVEL:** ~55 elements/rank (only boundary elements)
- **Speedup:** ~335x per rank

### Many-Particle System (2508 particles, 18432 elements/rank, 63 ranks)
- **Brute-force:** 18432 × 2508 ≈ 46.2 million element evaluations/rank
- **KVEL:** ~2,211 element evaluations/rank (only particles on this subdomain)
- **Speedup:** ~21,000x per rank

The speedup scales with:
1. Number of particles (more particles = more wasted work in brute-force)
2. Mesh refinement (finer mesh = more elements to skip)
3. Particle size relative to domain (smaller particles = fewer boundary elements)
