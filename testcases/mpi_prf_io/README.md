# MPI-PRF collective-I/O benchmark

This benchmark measures the MPI-PRF checkpoint path used by `q2p1_sse` and
checks each checkpoint before the next repetition can overwrite it. Run it
inside a Slurm allocation with a prepared SSE run directory and a binary built
from this branch:

```bash
PRF_BENCH_RANKS=8 PRF_BENCH_REPETITIONS=3 \
  ./benchmark_romio_matrix.sh /path/to/run ./q2p1_sse -a 0
```

The default matrix contains:

- untouched OpenMPI defaults;
- the detected ROMIO component with its built-in defaults;
- 4, 16, and 64 MiB collective buffers with data sieving enabled;
- a 16 MiB collective buffer with data sieving disabled.

The harness records the selected environment, full solver logs, raw samples in
`summary.tsv`, aggregates in `profile_summary.tsv`, and a review-ready table in
`profile_summary.md`. Checkpoint time comes from the `MPI_WTIME` record emitted
by `postprocessing_sse`; payload bytes are measured from the key/component PRF
files, so throughput excludes visualization and solver work.

The active MPI module normally puts `ompi_info` on `PATH`. If it does not, set
`PRF_BENCH_OMPI_INFO=/path/to/openmpi/bin/ompi_info`; this ensures component
detection matches the OpenMPI installation used by the solver.

## Reload and bitwise validation

Provide separate parameter decks to reload every produced checkpoint before
the next measurement:

```bash
PRF_BENCH_WRITE_DECK=_data/q2p1_param.write.dat \
PRF_BENCH_RESTART_DECK=_data/q2p1_param.restart.dat \
  ./benchmark_romio_matrix.sh /path/to/run ./q2p1_sse -a 0
```

The write deck must produce at least one MPI-PRF checkpoint. The restart deck
must use `StartingProc = 1`, select MPI-PRF (`DumpFormat = 3` in the current
SSE input), and address the checkpoint slot selected by the SSE angle. A
successful restart is recorded as `reload_ok`.

For a true bitwise write→read→write comparison, configure the restart run to
read a different staged slot and regenerate the writer slot, then set:

```bash
PRF_BENCH_STAGED_RESTART_SLOT=2
```

The harness archives the original, moves it to slot 2, performs the reload,
and compares all `*_key*.prf` and `*_comp*.prf` files against the regenerated
writer slot. Only bitwise equality passes. This is the mode used to detect the
pressure/coordinate corruption fixed on this branch. For SSE, the restart
deck and invocation must select the staged angle/slot consistently.

By default only exit code 0 is successful. Cases that intentionally finish
with another documented application status can set a comma-separated list,
for example `PRF_BENCH_SUCCESS_CODES=0,57`.

Use `PRF_BENCH_PROFILES` to run a subset after screening at 8 ranks, for
example:

```bash
PRF_BENCH_RANKS=16 \
PRF_BENCH_PROFILES='openmpi_default romio_4m_enable' \
  ./benchmark_romio_matrix.sh /path/to/run ./q2p1_sse -a 0
```

Keep the screening runs small (4 or 8 ranks for correctness, 8 ranks for the
full performance matrix). Large 72-to-288-rank experiments are deliberately
outside this benchmark workflow.
