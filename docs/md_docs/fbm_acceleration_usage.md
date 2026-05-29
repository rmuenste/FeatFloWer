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

- **Values:** `Yes` (default) | `No`
- **Effect when `Yes`:** Uses spatial HashGrid for O(1) particle containment queries (after timestep 1)
- **Effect when `No`:** Uses brute-force `verifyAllParticles` (linear search over all particles)
- **Performance impact:** ~100-1000x speedup for many-particle systems

### `SimPar@UseKVELAccel`

Controls KVEL/KEEL/KAAL candidate element acceleration for force integration.

- **Values:** `Yes` (default) | `No`
- **Effect when `Yes`:** Only integrates over elements adjacent to particle DOFs
- **Effect when `No`:** Integrates over all mesh elements (brute-force)
- **Performance impact:** ~10,000-20,000x reduction in element evaluations for many-particle systems

## Build Requirements

Both accelerations require:
```bash
cmake -DUSE_PE=ON -DUSE_PE_SERIAL_MODE=ON ..
```

The acceleration code is compiled in by default (`-DENABLE_FBM_ACCELERATION=ON`).

## Testing Workflow

### 1. Correctness Verification (Same Binary)

Build once with acceleration compiled in:
```bash
cmake -DUSE_PE=ON -DUSE_PE_SERIAL_MODE=ON ..
make -j8
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

Compare forces/positions/velocities — should match to machine precision.

### 2. Final Verification (Pure Baseline Build)

Rebuild without acceleration code:
```bash
cmake -DENABLE_FBM_ACCELERATION=OFF ..
make -j8
```

The runtime flags have no effect. Compare this pure baseline against the accelerated run.

## Output

When enabled, you'll see log lines like:

```
KVEL cache: 354699 DOFs cached (all ranks)
KVEL: 1.3932E+05 candidates vs 3.0818E+09 brute-force (22121.0x speedup)
```

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
