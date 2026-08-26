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

The matrix includes the untouched OpenMPI default, ROMIO's own defaults, and
4/16/64 MiB collective buffers with data sieving enabled, plus the colleague's
16 MiB setting with data sieving disabled. It records total elapsed time and
preserves each solver log; the solver log also reports checkpoint-only elapsed
time and payload bytes. Screen all profiles at 8 ranks, then rerun only the
fastest two at 16 ranks. The former 72-to-288-rank experiment is intentionally
deferred until queue conditions justify it.

Run each profile from an equivalent prepared state. For overwrite coverage,
reuse a slot after reducing the selected field or chunk count and verify no
orphan `*_comp*_chunk_*.prf` files remain. Also retain cases for a marker-less
legacy MPI-PRF slot, a legacy SmartDump `time.prf`, a deliberately incomplete
generation, same-level restart, lower-level prolongation, and 4-to-8-rank
repartition restart.
