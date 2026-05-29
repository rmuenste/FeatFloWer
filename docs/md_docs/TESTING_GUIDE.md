# KVEL Acceleration Testing Guide

## Quick Start Testing

### 1. Build the Code

```bash
cd /data/warehouse17/rmuenste/code/FF-ATC-NEW/ff-accel/build-ninja-baseline
ninja  # Full build
```

### 2. Run a Test Case

Use any existing sedimentation or FBM benchmark that uses PE_SERIAL_MODE:

```bash
cd /data/warehouse17/rmuenste/code/FF-ATC-NEW/ff-accel
# Run your favorite test case - the acceleration is automatic!
```

### 3. Check the Output

Look for this line in the output (from rank 1):
```
KVEL cache: 12345 vertices cached
```

And later (after force computation):
```
KVEL: 5432 candidates vs 2000000 brute-force (368.3x speedup)
```

---

## Enabling/Disabling Acceleration

### Enable (Default)
The acceleration is **ON by default**. No action needed.

### Disable for Comparison
Edit `source/src_quadLS/QuadSc_var.f90` line ~321:
```fortran
LOGICAL :: bUseKVEL_Accel = .FALSE.  ! Changed from .TRUE.
```

Rebuild and run to compare performance.

---

## Verification Tests

### Test 1: Bit-Identical Forces

**Goal**: Verify KVEL produces identical forces to brute-force

**Steps**:
1. Run test with `bUseKVEL_Accel = .TRUE.` → save forces
2. Run test with `bUseKVEL_Accel = .FALSE.` → save forces
3. Compare forces (should match to machine precision)

**Expected**: Forces identical to ~1e-14 relative error

### Test 2: Performance Scaling

**Goal**: Measure speedup vs mesh refinement

**Test Cases**:
- NLMAX=4 (32k elements)
- NLMAX=5 (262k elements)
- NLMAX=6 (2M elements)

**Expected Speedup**:
- NLMAX=4: ~8x
- NLMAX=5: ~16x
- NLMAX=6: ~33x

### Test 3: Multiple Particles

**Goal**: Verify cache works correctly with many particles

**Test Cases**:
- 1 particle
- 10 particles
- 50 particles

**Expected**: Speedup should increase with more particles (more reuse of cached data)

---

## Debugging

### No Cache Output?

**Problem**: Don't see "KVEL cache: X vertices cached"

**Check**:
1. Are you using PE_SERIAL_MODE? (Check build flags)
2. Is myFBM%nParticles > 0?
3. Are any vertices marked in FictKNPR?

### Warning: No KVEL Candidates

**Problem**: See "WARNING: No KVEL candidates for particle X"

**Possible Causes**:
1. Particle completely outside domain
2. Cache not built (check alpha computation)
3. Particle has zero cached vertices

**Fix**: This should be rare. Fallback to brute-force is automatic.

### Segmentation Fault

**Problem**: Code crashes during force computation

**Check**:
1. KVEL array bounds: Is `ivt <= NVT`? ✅ (we check this)
2. Particle indices: Is `IP <= size(ParticleVertexCache)`? ✅ (we check this)
3. NLMAX defined: Is MGPAR common block accessible? ✅ (we added this)

---

## Performance Profiling

### Measure Force Computation Time

Add timing around force computation in your application:

```fortran
call ztime(t_start)
call ForcesLocalParticlesSerial(...)
call ztime(t_end)
if (myid == 1) write(*,*) 'Force time:', t_end - t_start
```

### Expected Results

**Baseline (brute-force)**:
- NLMAX=5: ~10 seconds per timestep (1 particle)

**KVEL Accelerated**:
- NLMAX=5: ~0.6 seconds per timestep (1 particle)
- NLMAX=6: ~0.5 seconds per timestep (vs ~20 seconds brute-force)

---

## Common Issues

### Cache Not Deallocated

**Symptom**: Memory usage grows over timesteps

**Check**: Verify cleanup code runs (end of ForcesLocalParticlesSerial)

**Fix**: Should not happen - we deallocate every timestep

### Statistics Show 0x Speedup

**Symptom**: `KVEL: X candidates vs Y brute-force (0.0x speedup)`

**Cause**: Division by zero protection

**Check**: Are candidate elements actually being found?

---

## Success Indicators

✅ Code compiles without errors
✅ "KVEL cache: X vertices" message appears
✅ "KVEL: X candidates vs Y brute-force (Z.Zx speedup)" shows Z > 10
✅ Forces match brute-force to machine precision
✅ Simulation runs faster than baseline

---

## Contact Info for Issues

If you encounter problems:

1. Check build log for errors: `/tmp/build_kvel3.log`
2. Verify PE_SERIAL_MODE is enabled
3. Check that KVEL mesh structure is built (mg_mesh)
4. Review implementation summary: `KVEL_IMPLEMENTATION_SUMMARY.md`

---

## Advanced: Verification Mode (Optional)

To add comprehensive verification, add to QuadSc_force_serial.f90:

```fortran
#ifdef VERIFY_KVEL_CACHE
  ! Compare candidate set with brute-force
  LOGICAL :: bBruteForce(NEL)
  bBruteForce = .FALSE.

  DO IEL = 1, NEL
    ! [existing boundary check code]
    IF (NJALFA /= 27 .AND. NIALFA /= 27) bBruteForce(IEL) = .TRUE.
  END DO

  ! Check for missed elements
  INTEGER :: nMissed
  nMissed = 0
  DO IEL = 1, NEL
    IF (bBruteForce(IEL) .AND. .NOT. bCandidateElement(IEL)) THEN
      WRITE(*,*) 'ERROR: Missed boundary element', IEL
      nMissed = nMissed + 1
    END IF
  END DO

  IF (nMissed > 0) THEN
    WRITE(*,*) 'VERIFICATION FAILED: Missed', nMissed, 'elements'
  ELSE
    WRITE(*,*) 'VERIFICATION PASSED: All boundary elements found'
  END IF
#endif
```

Build with `-DVERIFY_KVEL_CACHE` to enable.
