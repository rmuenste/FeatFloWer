# q2p1_sse Solver Loop

Primary files:

- `applications/q2p1_sse/q2p1_sse.f90`
- `source/src_quadLS/QuadSc_transport_extensions.f90`
- `source/src_quadLS/QuadSc_utilities.f90`

## Program-Level Loop

`Q2P1_SSE` performs a thin lifecycle:

```text
parse command-line options
init_q2p1_ext
start notification
for itns = 1, nitns:
  update time state
  Transport_q2p1_UxyzP_sse
  if diverged: exit
  if pressure-convergence mode:
    postprocess and exit
  optional Transport_LinScalar
  postprocessing_sse
  print_time
  handle_statistics
  stop if timemx reached
finish notification
abort on divergence
DetermineIfGoalsWereReached
sim_finalize_sse
exit with success or RUNOUTOFTIMESTEPS code
```

## Momentum Solver Boundary

The main application calls:

```text
Transport_q2p1_UxyzP_sse(ufile, inl_u, itns)
```

This routine lives in `source/src_quadLS/QuadSc_transport_extensions.f90` and is
the main momentum/pressure solve boundary for the SSE application.

## Scalar Transport Boundary

If `bTracer` is enabled, the loop calls:

```text
Transport_LinScalar(ufile, inonln_t)
```

Otherwise `inonln_t` is set to `2` and postprocessing still runs.

## Goal And Failure Handling

- `DivergedSolution` exits the loop and triggers `MPI_Abort` with the divergence
  error code.
- `DetermineIfGoalsWereReached` performs the final convergence/goal check.
- The current source has the MPI abort for run-out-of-timesteps commented out;
  the program instead exits with `myErrorCode%RUNOUTOFTIMESTEPS`.

