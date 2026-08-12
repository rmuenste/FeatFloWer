# Runbook: unit_straddling_contact (cross-rank contact Newton-pair audit)

**Purpose:** fast unit-style reproducer for the Newton's-third-law violation that
the PE hard-contact solver exhibited for contacts whose two bodies are owned by
different MPI ranks (fixed in libs/pe commit "HardContactEulerLagrange: restore
exact Newton-pair momentum balance (contacts + migration)"). Two spheres
(r = 0.005) are seeded IN CONTACT (center distance 0.008, overlap 2e-3)
straddling the x = 1/3 PE/CFD partition face of the shared QBOX9 3x3x3
periodic box (mesh symlinked from `tier2_cases/momentum_conservation/mesh`).
Gravity acts in z (Prop@Gravity), fluid feedback is OFF, difelice drag,
`ELMomentumAuditFreq = 1`, 20 steps of dt = 1e-3 with 2 PE substeps.

The intentional seed overlap requires the PE-side json flag
`"seedAllowContact_": true` (added together with the fix), which skips the
min-gap sanity guard of the EL terminal-velocity setup in file seed mode.

**np = 28** (27 workers from the 3x3x3 decomposition + 1 master).

## Run

```bash
cd /path/to/FeatFloWer
tools/el_stage_rundir.sh unit_straddling_contact /path/to/fresh-rundir
BIN=$PWD/build-el-phase2-pe-gcc14/applications/q2p1_el_pipeflow/q2p1_el_pipeflow
MPIPREFIX=$(ldd "$BIN" | grep -m1 'libmpi\.so' | awk '{print $3}' | sed 's#/lib/libmpi.*##')
export PATH="$MPIPREFIX/bin:$PATH"; export LD_LIBRARY_PATH="$MPIPREFIX/lib:$LD_LIBRARY_PATH"
export OMPI_MCA_rmaps_base_oversubscribe=1
mpirun --oversubscribe --wdir /path/to/fresh-rundir -np 28 "$BIN" | tee run.log
grep -E "EL_NEWTON_PAIR|EL_CONTACT_STATS" run.log
```

## Acceptance

- Step 1 must show `EL_CONTACT_STATS ... ncontacts= 2` (the seeded pair is in
  contact; the contact is resolved on the rank owning the left sphere while the
  right sphere is a shadow copy there) and the Baumgarte impulse separates the
  pair along x.
- `EL_NEWTON_PAIR mismatch` must be at machine zero (|mismatch| < 1e-15,
  typically <= 2e-19) for EVERY step, including step 1 with the active
  cross-rank contact.

## Reference results

- Broken PE (before fix): step-1 `mismatch_x = -7.7e-5` (elRelax asymmetry),
  or `-1.5e-12` with only the relaxation fix applied (epsilon-tolerant
  send-guard dropping the PGS convergence tail).
- Fixed PE: step-1 `mismatch = (0.0, 0.0, -5.1e-21)`; all 20 steps < 2e-19.
