# Runbook: pipe_hp_check — single-phase Hagen–Poiseuille validation (W3.4)

Goal: validate the axially-periodic pipe infrastructure (z-periodicity via
`PeriodicAxis`/`PeriodicLength`, analytic-cylinder wall parametrization,
constant axial body force) against the exact Hagen–Poiseuille solution before
any particle-laden pipe run (V3/V4/V5).

## Configuration

- Mesh: PIPEZ27 (see MESH.md), R = 0.5, L = 2.0, np = 28 (27 subs + master).
- Fluid: ρ = 1, ν = 0.02 (μ = 0.02), axial body force f_z = 0.02.
- One PASSIVE probe particle (forces not applied, no feedback) — the EL PE
  setup requires ≥1 particle; it does not influence the flow.
- Time: dt = 0.02 (PE `stepsize_` matched — the startup assertion enforces it),
  750 steps to t = 15. Startup transient time constant τ = R²/(λ₁²ν) ≈ 2.16
  (λ₁ = 2.405), so t = 15 ≈ 7τ (residual transient < 0.1%).
- Time averaging: `ELTAvgWindow = 0.2` → TAVG over t ∈ [12, 15].

## Analytic targets

```
u_max  = f R² / (4 μ) = 0.02·0.25/0.08  = 0.0625
u_mean = f R² / (8 μ)                    = 0.03125
Q      = π R⁴ f / (8 μ) = π R² u_mean    = 0.024544
```

The volume-average axial velocity over the pipe equals u_mean, so the
acceptance metric is the `uf_super`/`uf_intr` z-component of the
`EL_MEAN_SLIP_TAVG` line (ε_f ≈ 1 everywhere; the probe's α_p is O(1e-4)).

## Run

```bash
BUILD=/home/user/rmuenste/nobackup/code/FF-EL/FeatFloWer/build-el-phase2-pe-gcc14
BIN=$BUILD/applications/q2p1_el_pipeflow/q2p1_el_pipeflow
MPIPREFIX=$(ldd "$BIN" | grep -m1 'libmpi\.so' | awk '{print $3}' | sed 's#/lib/libmpi.*##')
export PATH="$MPIPREFIX/bin:$PATH"; export LD_LIBRARY_PATH="$MPIPREFIX/lib:$LD_LIBRARY_PATH"

cmake --build "$BUILD" --target q2p1_el_pipeflow_val_pipe_hp_check_stage
RUN=$BUILD/applications/q2p1_el_pipeflow
mpirun --oversubscribe --wdir "$RUN" -np 28 "$BIN" > "$RUN/run_hp.log" 2>&1
grep EL_MEAN_SLIP_TAVG "$RUN/simulation_output_level_2.log" | tail -1
```

For the mesh-convergence pair, rerun with `SimPar@MaxMeshLevel = 2` in the
staged `_data/q2p1_param.dat` (halved resolution; expect a larger deficit).

## Acceptance

- Level 3: |uf_z − 0.03125| / 0.03125 ≤ 3%.
- Level 2 recorded alongside (convergence direction must be toward analytic).
- No periodic-plane anomalies (pressure/velocity fields smooth across z = 0/2,
  spot-check VTK output); startup dt assertion passes; run completes cleanly.

## Results (2026-07-04, np=28, build-el-phase2-pe-gcc14)

| Level | elements | uf_intr_z (TAVG) | error vs 0.03125 |
|---|---|---|---|
| 2 | 10,368 | 0.0313647 | +0.37% |
| 3 | 82,944 | 0.0312310 | **−0.06%** |

- PASS: level-3 error 0.06% (gate 3%); |error| drops 6.2× from level 2 —
  clean mesh convergence toward analytic.
- eps = 0.9999583 = 1 − V_probe/V_pipe exactly (probe bookkeeping correct);
  transverse velocity components O(1e-14); no fatal/OOB/clipping lines;
  startup dt assertion passed; textbook startup transient
  (u(t): 0.0072 @ t=0.5 → 0.03113 @ t=10 → converged by t=12).
- Two defects found and fixed by this case (commit da536c8e): unquoted
  type-7 hull descriptor in the generated .par files (list-directed read
  truncation → EOF abort in InitParametrization), and the unsliced
  dvol(nel+1) total-volume slot doubling the EL volume integrals.
