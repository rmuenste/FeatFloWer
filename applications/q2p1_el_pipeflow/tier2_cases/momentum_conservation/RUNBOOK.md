# Runbook: Tier-2 np=28 momentum-conservation case

**Goal:** run the no-wall periodic-box momentum case end-to-end and confirm the
post-cleanup acceptance metric (`EL_MOMENTUM_ELEMINT drift_rel < 1e-5`). This is the
runtime validation that needs the staged mesh + a hand-launched MPI job (no
`module load` in the shell).

## Fixed facts (this environment)
- **Repo root:** `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer`
- **Binary (must be built):** `build-el-phase2-pe-gcc14/applications/q2p1_el_pipeflow/q2p1_el_pipeflow`
- **Case source (committed):** `applications/q2p1_el_pipeflow/tier2_cases/momentum_conservation/`
  — `q2p1_param.dat` (MaxNumStep=60, MaxMeshLevel=2, **MeshFolder="QBOX9"**),
  `example.json` (processesX/Y/Z=3 -> **27 workers + 1 master = np 28**),
  `particles.xyz`, and `mesh/partitions/QBOX9` (27 subdomains).
- **np = 28.** Fixed by the 3x3x3 decomposition. Not 4, not any other value.
- The np=4 terminal-velocity case uses `QBOX9Z3` (3 subs) instead -- do not use it here.

## Step 1 - resolve MPI from the binary (no modules)
The shell has no MPI on `PATH`. Don't `module load`; read it from the binary's rpath:
```bash
cd /data/warehouse17/rmuenste/code/FF-EL/FeatFloWer
BIN=$PWD/build-el-phase2-pe-gcc14/applications/q2p1_el_pipeflow/q2p1_el_pipeflow
ldd "$BIN" | grep -q 'not found' && { echo "MISSING LIBS"; ldd "$BIN" | grep 'not found'; exit 1; }
MPIPREFIX=$(ldd "$BIN" | grep -m1 'libmpi\.so' | awk '{print $3}' | sed 's#/lib/libmpi.*##')
# Expected here: /sfw/openmpi/gcc14.3.x/4.1.6/ucx-threaded-noverbs
export PATH="$MPIPREFIX/bin:$PATH"
export LD_LIBRARY_PATH="$MPIPREFIX/lib:$LD_LIBRARY_PATH"
echo "Using MPI: $MPIPREFIX"
```

## Step 2 - stage the committed case into the run dir
Re-stage from source so the run uses the current committed param (incl. the
`ELMomentumAuditFreq` key and MaxNumStep=60):
```bash
SRC=$PWD
CASE=$SRC/applications/q2p1_el_pipeflow/tier2_cases/momentum_conservation
RUN=$SRC/build-el-phase2-pe-gcc14/applications/q2p1_el_pipeflow
cp -f  "$CASE/q2p1_param.dat" "$RUN/_data/q2p1_param.dat"
cp -f  "$CASE/example.json"   "$RUN/example.json"
cp -f  "$CASE/particles.xyz"  "$RUN/particles.xyz"
mkdir -p "$RUN/_mesh"
cp -rf "$CASE/mesh/partitions/QBOX9" "$RUN/_mesh/QBOX9"   # 27 subdomains
ls -d "$RUN/_mesh/QBOX9"/sub* | wc -l   # must print 27
```

## Step 3 - launch 28 ranks
On an interactive/login node, OpenMPI must be told to oversubscribe (28 ranks >
advertised slots). `--wdir` makes each rank find `_data/`, `_mesh/`, `example.json`
by relative path:
```bash
LOG="$RUN/run_mom_$(date +%s).log"
mpirun --oversubscribe --wdir "$RUN" -np 28 "$RUN/q2p1_el_pipeflow" > "$LOG" 2>&1
echo "exit=$? log=$LOG"
```
> If you land on a **Slurm allocation** (`$SLURM_JOB_ID` set), drop
> `--oversubscribe` and use a real allocation: `srun -n 28 "$RUN/q2p1_el_pipeflow"`.
> The `--oversubscribe` flag is only a login-node workaround.

## Step 4 - acceptance checks
```bash
echo "=== crashes? (must be empty) ==="
grep -E 'signal (6|11)|Aborted|Segmentation|EL_KDFG|out of range' "$LOG"

echo "=== momentum drift (last step) ==="
grep 'EL_MOMENTUM_ELEMINT' "$LOG" | tail -1

echo "=== feedback Newton pair (last step) ==="
grep 'EL_FEEDBACK_CONSERVATION' "$LOG" | tail -1
```
**Pass criteria:**
- No crash/abort lines, and the **NDFGL bounds assertion must NOT fire** (it replaced
  the `EL_OOB` probe -- if it trips, the `/TRIAD/` pin regressed).
- `EL_MOMENTUM_ELEMINT ... drift_rel=` **< 1e-5** (observed peak with the fix was
  ~4.3e-7, so expect much less than tol). Value is whitespace-token **column 16,
  0-indexed** on that line.
- `EL_FEEDBACK_CONSERVATION ... residual=` <= ~1e-10 (post-fix ~1e-19). Column 12,
  0-indexed.

## Step 5 (optional) - run it through the harness instead
To exercise the YAML metric extraction + baseline compare rather than eyeballing,
drive the definition at
`tools/featflower_test/testcases/definitions/q2p1_el_pipeflow_tier2_momentum.yaml`
(builds the `*_stage` target, launches np=28, parses `EL_MOMENTUM_ELEMINT` vs
baseline). Verify the definition keys the metric on `EL_MOMENTUM_ELEMINT` at column
16 with tol 1e-5, and that the baseline file key/value match. This is the end-to-end
check that the harness change actually parses real output.

## Common gotchas
- **Wrong mesh:** `QBOX9Z3` (3 subs) with np=28 fails partition/rank matching. Use
  `QBOX9` (27 subs).
- **`mpirun: command not found`** -> Step 1 wasn't sourced into the current shell.
- **"There are not enough slots"** -> missing `--oversubscribe` on a login node.
- **Missing inputs at runtime** -> forgot `--wdir "$RUN"`, or `_mesh/QBOX9` /
  `_data/q2p1_param.dat` not staged.
- **Stale binary:** if source changed, rebuild first:
  `cmake --build build-el-phase2-pe-gcc14 --target q2p1_el_pipeflow -j8`.
