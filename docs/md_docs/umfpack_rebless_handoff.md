# Handoff: Benchmark Re-Bless of UMFPACK 7.12.2 (FetchContent) vs Phase-1 Numbers

Goal: the Phase-1 performance numbers for pressure coarse types 2 and 4 were
measured against **site UMFPACK 5.7.8** (`USE_EXTERNAL_SUITESPARSE` +
`EXTERNAL_SUITESPARSE_DIR`). The FetchContent integration now builds
**SuiteSparse 7.12.2** as the fetch provider. Nobody has benchmarked 7.12.2.
Run the standard harness against the fetch-provider build and produce a
verdict: does UMFPACK 7.12.2 match/beat the 5.7.8 Phase-1 numbers for t2/t4,
and are the forces still correct?

## Background you must read first

- `docs/md_docs/solver_baseline_phase0.md` — the harness
  (`tools/cluster_scripts/solver_baseline_phase0.slurm`), its environment
  knobs (`FF_P0_BUILD` selects the build tree; also `FF_P0_VELO_TYPES`,
  `FF_P0_VELO_MINLEV`, backend knob `FF_P0_BACKEND`), and the summarizer.
  Read the slurm script itself for the authoritative knob list — it may have
  a pressure-type knob; if it does not, adapt minimally (see Constraints).
- `docs/md_docs/solver_phase1_direct_solvers.md` — how Phase 1 measured
  t2/t4/t9 and where its results live.
- `docs/md_docs/fetchcontent_dependencies.md` — the provider pattern and the
  `nm | grep '_gfortran_st_'` sanity check (all of `st_open`, `st_read`,
  `st_close` must be `T` in the binary; a split runtime hangs at startup).

## Reference numbers to beat/match

- Vendored baseline: `stage1-work/results/phase0-baseline-v2/summary.tsv`.
- Phase-1 external 5.7.8: `stage1-work/results/phase1-direct/`.
- Headline Phase-1 facts (medians per coarse call): external UMFPACK t2 was
  1.9–2.3x faster than vendored on cylinder, 1.5x on FAC3D; t4 unchanged
  (extraction-dominated); t9 (PARDISO) fastest everywhere. Forces were
  bitwise identical per type across builds.

## The build under test

`stage1-work/build-fetchdeps-cpu` — already configured and fully built
(fetch-provider hypre 2.33.0 + SuiteSparse 7.12.2, `USE_MKL_PARDISO=ON`,
`ENABLE_SOLVER_TIMING=ON`); `q2p1_fc_ext` and `q2p1_fac3d` binaries exist and
the cylinder case has been run to completion on the login node. Do NOT
reconfigure it. If you must rebuild a target, use the environment below.
Before submitting anything, run the `nm` libgfortran sanity check on both
binaries and confirm `ldd` shows no stray site SuiteSparse.

## Task

1. Run the Phase-0/1 harness on **both** benchmarks — cylinder
   (`q2p1_fc_ext`) and undeformed FAC3D (`q2p1_fac3d`) — with
   `FF_P0_BUILD` pointing at `build-fetchdeps-cpu`, covering at minimum
   pressure coarse types **2 and 4** plus **9** as a cross-build control
   (t9 does not use UMFPACK; its timings should be unchanged, which
   validates the comparison is apples-to-apples). The SSE benchmark is
   EXCLUDED by standing user decision — do not add it.
2. Mirror the Phase-1 measurement setup (same np, same levels, same number
   of timesteps) so medians are directly comparable. The Phase-1 doc and
   the `phase1-direct` results directory record what was used.
3. Collect results into `stage1-work/results/rebless-umfpack712/` with the
   same summarizer/summary.tsv format as previous phases.
4. Compare, per benchmark and per type: median time per coarse call vs
   phase1-direct (5.7.8) and vs phase0-baseline-v2 (vendored); drag/lift
   forces vs the Phase-1 runs (bitwise identical expected per type; report
   the max relative deviation if not bitwise).
5. Write the verdict into a new section appended to
   `docs/md_docs/fetchcontent_dependencies.md` (it already flags the
   re-bless as pending) — table of numbers + one-paragraph conclusion.
   No new doc file, no README.md index edit needed.

## Slurm — ALLOWED for this task (unlike previous handoffs)

- You MAY submit jobs with `sbatch` and poll with `squeue`/`sacct`. This is
  the point of the task.
- Cluster gotchas (tardis): the partition only accepts jobs that request its
  GPU gres even for CPU jobs (copy the gres line from
  `tools/cluster_scripts/solver_baseline_phase0.slurm` or a previous phase
  job script); the short partition caps at 2 h. Previous full harness runs
  fit comfortably within that.
- The phase0 slurm script already stages `libmetis.so` next to the app for
  PyPartitioner — keep that logic.
- Poll politely (sleep ≥60 s between squeue checks). If a job pends for a
  long time, keep waiting — do not cancel and resubmit repeatedly.

## Constraints

- Work directly in `/data/warehouse17/rmuenste/code/FF-GPU/FeatFloWer`.
  No git worktree, no branch, no commits, no checkout/restore/reset/stash —
  the tree holds uncommitted work that must survive.
- Do NOT modify `source/` or any CMake files. This is a measurement task.
  If the harness script needs a small extension (e.g. a pressure-type knob
  or a build-dir knob it lacks), prefer environment-variable-driven,
  backward-compatible edits to
  `tools/cluster_scripts/solver_baseline_phase0.slurm` only, and document
  them in your report. Alternatively copy it to a new script
  `tools/cluster_scripts/rebless_umfpack712.slurm` and leave the original
  untouched (preferred if changes are more than trivial).
- Do NOT reconfigure or rebuild `build-phase2`, `build-gpu-phase3`,
  `build-typevalid-*`, or `build-fetchdeps-gpu`.
- Do not overwrite `phase0-baseline-v2` or `phase1-direct` results.

## Environment (no modules needed)

```
export PATH=/sfw/gcc/13.2.0-static-gmp-mpfr-mpc-isl/bin:/sfw/openmpi/gcc13.2.x/4.1.6/ucx-threaded-noverbs/bin:$PATH
export LD_LIBRARY_PATH=/sfw/gcc/13.2.0-static-gmp-mpfr-mpc-isl/lib64:$LD_LIBRARY_PATH
```

## Verification (definition of done)

1. Jobs completed (sacct state COMPLETED, exit 0) for both benchmarks.
2. `stage1-work/results/rebless-umfpack712/summary.tsv` exists and contains
   t2/t4/t9 rows for both benchmarks.
3. The t9 control medians agree with Phase-1 t9 within run-to-run noise
   (~10–20%); if they do not, flag the whole comparison as suspect and say
   why rather than papering over it.
4. Forces compared and reported.
5. `docs/md_docs/fetchcontent_dependencies.md` updated with the verdict
   section.

## Deliverables

- Results directory + summary.tsv.
- Doc section with the comparison table and verdict
  (7.12.2 vs 5.7.8 vs vendored, per benchmark, per type).
- Final report: job IDs, exact sbatch lines and env, any harness-script
  changes, the numbers, forces outcome, and a clear one-line verdict
  ("re-blessed" or "regression found: ...").
