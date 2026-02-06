# HashGrid Verification System

## Overview

The HashGrid verification system provides runtime validation of the accelerated HashGrid-based point-inside queries against the baseline linear search implementation. This ensures correctness of the HashGrid acceleration by comparing results for every single query during simulation.

## How It Works

When enabled, the verification system:

1. **Calls both implementations** - For each point query, both the accelerated HashGrid method and the baseline linear search are executed
2. **Compares results** - Checks if both methods agree on:
   - Whether the point is inside any particle (boolean match)
   - Which particle contains the point (particle ID match)
3. **Reports discrepancies** - Prints detailed warnings when results differ
4. **Tracks statistics** - Maintains running counts of total queries and mismatches
5. **Periodic reporting** - Prints summary every 100,000 queries

## Building with Verification

### Enable HashGrid Acceleration + Verification

```bash
cd /path/to/ff-accel
mkdir build_verify
cd build_verify

# Load required modules
module load gcc/latest-v13 openmpi/options/interface/ethernet openmpi/4.1.6

# Configure with both acceleration and verification
cmake -DCMAKE_BUILD_TYPE=Release \
      -DUSE_PE=ON \
      -DUSE_PE_SERIAL_MODE=ON \
      -DUSE_ACCELERATED_POINT_QUERY=ON \
      -DVERIFY_HASHGRID=ON \
      -DBUILD_APPLICATIONS=ON \
      ..

make -j8
```

### CMake Options

- **`USE_ACCELERATED_POINT_QUERY=ON`** - Enable HashGrid acceleration (required for verification)
- **`VERIFY_HASHGRID=ON`** - Enable runtime verification
- **`VERIFY_HASHGRID=OFF`** - Disable verification (default, for production runs)

## Output Format

### Immediate Mismatch Detection

When a discrepancy is detected, you'll see output like:

```
HASHGRID MISMATCH at vertex 12345 point:   0.523000    0.412000    0.789000
  Accelerated: T Baseline: F
```

Or for particle ID mismatches:

```
HASHGRID PARTICLE MISMATCH at vertex 12345 point:   0.523000    0.412000    0.789000
  Accelerated particle: 42 Baseline particle: 87
```

### Periodic Statistics

Every 100,000 queries, you'll see a summary:

```
HashGrid verification: 100000 queries, 0 mismatches (  0.0000%)
HashGrid verification: 200000 queries, 2 mismatches (  0.0010%)
  Particle ID mismatches: 1 (both found particle but different IDs)
HashGrid verification: 300000 queries, 2 mismatches (  0.0007%)
  Particle ID mismatches: 1 (both found particle but different IDs)
```

### Interpreting Results

**Perfect match (0.0000% mismatch rate):**
```
HashGrid verification: 500000 queries, 0 mismatches (  0.0000%)
```
✅ HashGrid implementation is correct!

**Non-zero mismatch rate:**
```
HashGrid verification: 500000 queries, 150 mismatches (  0.0300%)
```
❌ There's a bug in the HashGrid implementation that needs investigation.

**Particle ID mismatches only:**
```
HashGrid verification: 500000 queries, 0 mismatches (  0.0000%)
  Particle ID mismatches: 5 (both found particle but different IDs)
```
⚠️ Rare edge case: Both methods detect the point is inside a particle, but disagree on which particle. This can happen at boundaries between overlapping particles.

## Performance Impact

**WARNING:** Verification mode significantly impacts performance because:

1. Every query calls **both** the accelerated and baseline methods
2. The baseline method uses O(N) linear search over all particles
3. Additional comparison and logging overhead

**Expected performance:**
- Normal mode: ~1000-10000 queries/sec (depending on particle count)
- Verification mode: ~50-500 queries/sec

**Use verification only for:**
- ✅ Testing and debugging
- ✅ Validating correctness with small test cases
- ✅ Regression testing after code changes
- ❌ NOT for production simulations
- ❌ NOT for large-scale runs

## Example Workflow

### 1. Test with Verification Enabled

```bash
# Build with verification
cd build_verify
cmake -DUSE_ACCELERATED_POINT_QUERY=ON -DVERIFY_HASHGRID=ON ..
make -j8 q2p1_hashgrid_test

# Run test (will be slower)
cd applications/q2p1_hashgrid_test
mpirun -np 4 ./q2p1_hashgrid_test > verification_output.txt 2>&1

# Check results
grep "MISMATCH" verification_output.txt
grep "HashGrid verification" verification_output.txt | tail -5
```

### 2. If Verification Passes, Build Production Version

```bash
# Build without verification for production
cd ../..
mkdir build_production
cd build_production

cmake -DUSE_ACCELERATED_POINT_QUERY=ON -DVERIFY_HASHGRID=OFF ..
make -j8

# Now run at full speed
```

## Implementation Details

### Code Location

**Verification code:** `source/src_fbm/fbm_main.f90` in `fbm_getFictKnprFC2()` function

**Key sections:**
```fortran
#ifdef VERIFY_HASHGRID
  ! Call both methods
  accel_result = checkAllParticles(cvidx, key, point, longFictId%bytes)
  baseline_result = pointInsideParticles(cvidx, baseline_key, point, baseline_bytes)

  ! Compare and report
  if (accel_result .neqv. baseline_result) then
    ! Print mismatch warning
  endif
#endif
```

### Statistics Tracking

Statistics are tracked using `SAVE` variables:
- `total_queries` - Total number of point queries executed
- `mismatch_count` - Number of inside/outside disagreements
- `particle_mismatch_count` - Number of particle ID disagreements

### CMake Integration

**Configuration:** `libs/pe/CMakeLists.txt`

The build system:
1. Adds `-DVERIFY_HASHGRID` preprocessor definition when enabled
2. Warns if verification is enabled without acceleration
3. Prints configuration status during CMake

## Troubleshooting

### "Verification enabled but acceleration is OFF"

**Error:**
```
CMake Warning:
  VERIFY_HASHGRID is enabled but USE_ACCELERATED_POINT_QUERY is OFF.
  Verification requires acceleration to be enabled.
```

**Solution:** Enable both options:
```bash
cmake -DUSE_ACCELERATED_POINT_QUERY=ON -DVERIFY_HASHGRID=ON ..
```

### No verification output during run

**Possible causes:**
1. Verification not enabled during build
2. No particles in simulation
3. No point queries being executed

**Check:**
```bash
# Verify preprocessor definition is set
strings ./q2p1_hashgrid_test | grep VERIFY_HASHGRID

# Check CMake configuration
cmake -L | grep VERIFY
```

### High mismatch rate

**If you see > 0.01% mismatch rate:**
1. Investigate the specific query points causing mismatches
2. Check if there are precision issues near particle boundaries
3. Review HashGrid cell size and resolution
4. Check for coordinate transformation bugs

## Advanced Usage

### Custom Reporting Interval

To change the reporting frequency, modify the periodic check in `fbm_getFictKnprFC2`:

```fortran
! Change 100000 to desired interval
if (mod(total_queries, 100000) == 0) then
```

### Per-Rank Reporting

By default, only rank 1 prints output. To see output from all ranks, remove the `if (myid == 1)` guards:

```fortran
#ifdef VERIFY_HASHGRID
  ! Remove this line to see output from all ranks:
  ! if (myid == 1) then
  write(*,'(A,I0,A)') 'Rank ', myid, ' verification stats...'
```

### Export Statistics to File

To save statistics to a file instead of stdout:

```fortran
integer :: verify_unit = 99
open(verify_unit, file='hashgrid_verification.log', status='replace')
write(verify_unit, '(A,I0)') 'Total queries: ', total_queries
close(verify_unit)
```

## Best Practices

1. **Always verify after code changes** - Run with verification after any modifications to HashGrid or collision detection code

2. **Use small test cases first** - Start with 10-100 particles to quickly identify issues

3. **Disable for production** - Never run production simulations with verification enabled

4. **Check logs carefully** - Even a small mismatch rate can indicate bugs

5. **Test with different particle counts** - Verify with various numbers of particles to test different HashGrid configurations

6. **Test edge cases** - Verify with particles near domain boundaries, overlapping particles, etc.

## Summary

The HashGrid verification system provides:
- ✅ **Per-query validation** - Every single query is verified
- ✅ **Detailed diagnostics** - Exact coordinates and particle IDs for mismatches
- ✅ **Statistical tracking** - Running counts and percentages
- ✅ **Easy to enable/disable** - Simple CMake flag
- ✅ **Minimal code changes** - Contained within one function

Use it during development and testing to ensure HashGrid acceleration is correct, then disable it for production runs.
