# q2p1_sse Output And Finalization

Primary files:

- `source/postprocessing/solution_io.f90`
- `applications/q2p1_sse/q2p1_sse.f90`
- `source/src_quadLS/QuadSc_utilities.f90`
- `applications/q2p1_sse/CMakeLists.txt`

## Main Output Entry

The application calls:

```text
postprocessing_sse(dout, inonln_u, inonln_t, ufile)
```

from `source/postprocessing/solution_io.f90`.

The call happens after each successful momentum solve and optional scalar solve.
There is also a pressure-convergence shortcut path that advances `timens` by
`dtgmv`, calls `postprocessing_sse`, and exits the main loop.

## Runtime Directories

The CMake staging/install logic creates and carries the application-specific
output directories:

- `_mesh`
- `_vtk`
- `_data`
- `_dump`
- `_1D`
- `_hist`
- `_prot0`
- `_prot1`
- `_RTD`
- `_ianus`

The launcher and solver also use `_data/prot.txt` and copy angle-specific
protocol files such as `_data/prot_0000.txt`.

Checkpoint backend selection, MPI-PRF metadata, restart behavior, and ROMIO
tuning are documented in [MPI-PRF checkpoints](../mpi_prf_checkpoints.md).

## Progress And Status

During the solve loop, `q2p1_sse.f90` calls:

- `print_time`
- `handle_statistics`
- `MemoryPrint`

The Python launcher watches `_data/prot.txt` and updates `e3d.log` with angle
iteration, inner iteration, heat iteration, and status information.

## Finalization

At the end of the Fortran program:

```text
DetermineIfGoalsWereReached
sim_finalize_sse
exit(0) or exit(myErrorCode%RUNOUTOFTIMESTEPS)
```

The launcher interprets selected exit codes:

- `55`: deformation failure path eligible for retry when `--retry-deformation`
  is enabled.
- `57`: run out of iterations.
- `88`: screw could not be created, usually a wrong-angle geometry failure.
