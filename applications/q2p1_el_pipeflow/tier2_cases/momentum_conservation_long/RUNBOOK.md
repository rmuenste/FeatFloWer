# Runbook: Tier-2 long-horizon momentum-conservation case

Goal: run the opt-in 10000-step periodic-box momentum case on 28 ranks and
confirm that the permanent momentum audit stays bounded over a longer horizon.

## Harness Run

The case is disabled by default because it is a long runtime check. Run it
explicitly through the harness when validating long-horizon behavior:

```bash
cd /data/warehouse17/rmuenste/code/FF-EL/FeatFloWer

BIN=$PWD/build-el-tier2-pe/applications/q2p1_el_pipeflow/q2p1_el_pipeflow
MPIPREFIX=$(ldd "$BIN" | grep -m1 'libmpi\.so' | awk '{print $3}' | sed 's#/lib/libmpi.*##')
export PATH="$MPIPREFIX/bin:$PATH"
export LD_LIBRARY_PATH="$MPIPREFIX/lib:$LD_LIBRARY_PATH"
export OMPI_MCA_rmaps_base_oversubscribe=1
export PYTHONPATH="$PWD/tools/featflower_test/src:$PYTHONPATH"

/usr/bin/python -m featflower_test run \
  tools/featflower_test/testcases/definitions/q2p1_el_pipeflow_tier2_momentum_long.yaml \
  --repo-root "$PWD"
```

## Acceptance Checks

The harness compares the final audited lines:

- `EL_MOMENTUM_ELEMINT` column 16 as `drift_rel`, tolerance `2.0e-5`.
- `EL_FEEDBACK_CONSERVATION` column 12 as `residual`, tolerance `1.0e-10`.

The first Slurm validation run (`20260702-211527-1fafec3a`) measured
`drift_rel=1.000013e-5` and `residual=1.229441e-27`. The tolerance is set
to roughly 2x the measured value to leave room for accumulated linear-solver
noise over this 10000-step audit while still catching a real drift regression.

Also confirm no void-fraction clipping occurred. The absence check is not a YAML
metric, so inspect the run log:

```bash
grep -c 'void-fraction clipping' results/runs/<run-id>/tests/q2p1_el_pipeflow_tier2_momentum_long/stages/run-level-2.log
```

The command must print `0`.

Archive the raw workdir log before launching another harness case, because the
shared application workdir is overwritten by the next staged run:

```bash
cp build-el-tier2-pe/applications/q2p1_el_pipeflow/simulation_output_level_2.log \
  results/runs/<run-id>/tests/q2p1_el_pipeflow_tier2_momentum_long/
```
