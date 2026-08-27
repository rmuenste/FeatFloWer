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
- 4, 16, and 64 MiB `cb_buffer_size` hints with `romio_ds_write` enabled;
- a 16 MiB `cb_buffer_size` hint with `romio_ds_write` disabled.

The tuned ROMIO profiles generate one profile-specific hints file under
`romio_hints/` and pass it through the real `ROMIO_HINTS` mechanism. The two
default profiles explicitly unset `ROMIO_HINTS`; their `cb_buffer_size` and
`romio_ds_write` columns mean “no requested hint”, not a claimed ROMIO
default. `OMPI_MCA_io` remains the component selector for the detected ROMIO
component. `profile_settings.tsv` records each profile's component, hints
path, source, and requested keys; these same provenance fields are copied to
each row in `summary.tsv`.

The harness records the selected environment, full solver logs, raw samples in
`summary.tsv`, aggregates in `profile_summary.tsv`, and a review-ready table in
`profile_summary.md`. Checkpoint time comes from the `MPI_WTIME` record emitted
by `postprocessing_sse`; payload bytes are measured from the key/component PRF
files, so throughput excludes visualization and solver work.

The active MPI module normally puts `ompi_info` on `PATH`. If it does not, set
`PRF_BENCH_OMPI_INFO=/path/to/openmpi/bin/ompi_info`; this ensures component
detection matches the OpenMPI installation used by the solver.

To inspect effective ROMIO settings, run a focused, untimed verification with
`PRF_BENCH_PRINT_HINTS=1`. It uses the profile named by
`PRF_BENCH_PRINT_HINTS_PROFILE` (default `romio_4m_enable`) and writes
`<profile>.romio_print_hints.log`. The probe runs immediately before that
profile's timed samples and checks that ROMIO printed the requested values.
The selected profile must be present in `PRF_BENCH_PROFILES`. Leave the option
unset for normal measurements.

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

If `PRF_BENCH_RESTART_COMPONENT` selects a different ROMIO component, it uses
that component's built-in defaults. Set `PRF_BENCH_RESTART_CB_BUFFER` and/or
`PRF_BENCH_RESTART_DS_WRITE` to generate a separate restart hints file with
those real ROMIO keys.

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

`summary.tsv` keeps a row for a failed write, including its solver exit code.
Solver failures use `write_status=failed`; successful solver exits without a
valid timed checkpoint use `write_status=checkpoint_failed`. Missing timing or
payload fields are recorded as `not_available`; the summarizer excludes those
values from timing statistics and reports `timed_samples`.

Keep the screening runs small (4 or 8 ranks for correctness, 8 ranks for the
full performance matrix). Large 72-to-288-rank experiments are deliberately
outside this benchmark workflow.
