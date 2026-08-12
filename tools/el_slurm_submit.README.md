# Slurm cluster notes for EL campaign runs

Field knowledge collected while running the Euler–Lagrange validation
campaign (2026-07). For any agent submitting FeatFloWer jobs here.

## Partitions (shared node pool, 33 nodes)

| Partition | Time limit | Use for |
|---|---|---|
| `short`* | 2 h | smokes, benchmarks |
| `med` | 8 h | diagnostics, U₀-class runs |
| `long` | 2 days | production sweeps (V2/W4.2 class) |
| `ultralong` | 28 days | only if truly necessary |

## Nodes

- **`nx-01..06`: 64 cores, 1.4 TB RAM — the fast nodes.** Reference: EL
  coupled box, level 3 (36³ cells, ~3k particles, np=28): **~1.65 s/step**.
- 48-core boxes (`worldgames`, `wintergames`, `82ndalltheway`,
  `christmastruce`, `gottmituns`, …): same workload ~2.0–2.1 s/step,
  **only ~29 GB RAM**.
- Use **`--prefer=nx` (soft), never `--constraint=nx` (hard)**: nx nodes
  are frequently all allocated or drained, and a hard constraint then
  blocks the job forever (this bit the Stage-0 battery; the stale
  constraint was removed from the momentum_long harness yaml).
- Memory: request **`--mem=25G` total** for np=28 EL runs. Do NOT use
  `--mem-per-cpu=2G` with `--ntasks=28` (=56 G — unschedulable on the
  48-core/29 GB boxes).

## Submitting EL runs

Always: isolated run dir + wrapper script (never run in the shared staged
build dir — concurrent runs clobber configs):

```bash
tools/el_stage_rundir.sh <case> <rundir>     # stage from validation_cases/
# per-run sed edits on <rundir>/_data/q2p1_param.dat and <rundir>/example.json
tools/el_slurm_submit.sh <rundir> <jobname> <walltime e.g. 1-18:00:00> [partition]
```

`el_slurm_submit.sh` writes `<rundir>/job.sbatch` (nodes=1, ntasks=28,
mem=25G, `--prefer=nx` built in), loads `gcc/latest-v13` +
`openmpi/options/interface/ethernet` + `openmpi/4.1.6` inside the job,
resolves the MPI prefix from the binary rpath, and runs
`mpirun -np 28 <BIN>` with output to `<rundir>/run_slurm.log`
(Slurm log: `<rundir>/slurm-%j.out`). Binary override: `EL_BIN=...`.

## Sizing rules of thumb (np=28)

| Workload | s/step (nx) | s/step (48-core) |
|---|---|---|
| Coupled box L3, ~3k particles | 1.65 | 2.0–2.1 |
| Coupled box L2, ~750 particles | ~0.16 (30k steps ≈ 80 min) | — |
| Coupled pipe L3 (PIPEZ27), ~30 particles | ~0.75 (50k ≈ 11 h) | — |
| Frozen-field runs (CFD solve skipped) | ≪ coupled; often local-feasible | — |

Local dev node has 6 physical cores: np=28 needs
`OMPI_MCA_rmaps_base_oversubscribe=1` and runs ~10× slower than nx —
fine for smokes/canaries, not production. On cluster nodes the
oversubscribe flag is unnecessary (but harmless).

## Queue behavior

- The queue is usually near-empty; jobs typically start within minutes.
  Cancel + resubmit is cheap (used e.g. to add VTK output to an
  already-launched sweep).
- Monitor: `squeue -u $USER -o '%.10i %j %T %P %l %N'`;
  post-mortem: `sacct -S <date> -o JobID,JobName%14,State,Elapsed,NodeList -X`.
- No user quota issues observed up to 7 concurrent np=28 jobs.

## VTK output in batch runs

- PE particle series: json `"vtk_": true` + `"visspacing_": N` → frames in
  `<rundir>/paraview/` (writer emits N−1 frames, known off-by-one, harmless).
- Fluid fields: `SimPar@OutputFreq = <T>` is in PHYSICAL TIME units →
  dumps to `<rundir>/_vtk/` (fields per `SimPar@OutputFields`).
- Both add negligible cost at modest cadence (e.g. visspacing 250,
  OutputFreq 10 t.u. on a 100 t.u. run).
