# MPI-PRF checkpoint tests

Use a small prepared application case whose `q2p1_param.dat` selects
`SimPar@DumpFormat = MPI_PRF` and produces at least one checkpoint. The
recommended correctness run writes with 4 MPI ranks and restarts the same slot
with 8 ranks (`StartingProc = 3`). Compare all restored fields bit-for-bit on
the same architecture and confirm that time and completed step match
`checkpoint_meta.prf`.

For performance screening, submit an interactive or batch allocation and run:

```bash
PRF_BENCH_RANKS=8 PRF_BENCH_REPETITIONS=3 \
  ./benchmark_romio_matrix.sh /path/to/prepared/run ./q2p1_fc_ext
```

To validate and reload every checkpoint before the next repetition can
overwrite it, prepare separate writer and restart parameter decks and set:

```bash
PRF_BENCH_WRITE_DECK=_data/q2p1_param.write.dat \
PRF_BENCH_RESTART_DECK=_data/q2p1_param.restart.dat \
  ./benchmark_romio_matrix.sh /path/to/prepared/run ./q2p1_fc_ext
```

The restart deck should select `StartingProc = 1`, point `StartFile` at the
written slot, and use `MaxNumStep = 0` so the application validates and loads
the checkpoint without advancing the simulation. Some applications still
produce a final dump; use staged-slot validation below for those cases. A
missing timing record, writer failure, or restart failure stops the matrix
immediately. The active deck defaults to `_data/q2p1_param.dat` and can be
overridden with `PRF_BENCH_ACTIVE_DECK`.

By default, reload uses the profile being measured. To benchmark writers while
validating every artifact through one fixed reader configuration, set
`PRF_BENCH_RESTART_COMPONENT` and optionally
`PRF_BENCH_RESTART_CB_BUFFER`/`PRF_BENCH_RESTART_DS_WRITE`. For example:

```bash
PRF_BENCH_RESTART_COMPONENT=romio321 \
PRF_BENCH_RESTART_CB_BUFFER=16777216 \
PRF_BENCH_RESTART_DS_WRITE=enable \
  ./benchmark_romio_matrix.sh /path/to/prepared/run ./q2p1_fc_ext
```

For applications that always write during finalization, stage the writer slot
under the slot selected by the restart deck and compare the regenerated dump:

```bash
PRF_BENCH_WRITE_SLOT=1 PRF_BENCH_STAGED_RESTART_SLOT=2 \
  ./benchmark_romio_matrix.sh /path/to/prepared/run ./q2p1_fc_ext
```

The harness moves slot 1 to slot 2, reloads slot 2, and compares every key and
payload file in the regenerated slot 1 bit-for-bit. A failure preserves both
slots and stops the matrix. Successful originals are archived below the result
directory before the next repetition starts.

The matrix includes the untouched OpenMPI default, ROMIO's own defaults, and
4/16/64 MiB collective buffers with data sieving enabled, plus the colleague's
16 MiB setting with data sieving disabled. It records total elapsed time and
preserves each solver log; the solver log also reports checkpoint-only elapsed
time and payload bytes. With restart validation enabled, the summary also
records reload elapsed time and status. Screen all profiles at 8 ranks, then
rerun only the fastest two at 16 ranks. The former 72-to-288-rank experiment is
intentionally deferred until queue conditions justify it.

Run each profile from an equivalent prepared state. For overwrite coverage,
reuse a slot after reducing the selected field or chunk count and verify no
orphan `*_comp*_chunk_*.prf` files remain. Also retain cases for a marker-less
legacy MPI-PRF slot, a legacy SmartDump `time.prf`, a deliberately incomplete
generation, same-level restart, lower-level prolongation, and 4-to-8-rank
repartition restart.
