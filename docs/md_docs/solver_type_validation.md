# Early Validation of Solver Types vs Compiled Capabilities

FeatFloWer selects the coarse-grid (and, for the pressure, optionally the
fine-level) linear solver through two integers in `_data/q2p1_param.dat`:

```
Velo@MGCrsSolverType = 1
Pres@MGCrsSolverType = 4
```

Several of those types are backed by optional third-party libraries that are
only present if the binary was configured with the matching CMake option.
Before this validation existed, a mismatch was discovered *at the first
solve* — minutes into a large job, after startup, partitioning and matrix
assembly — and in the worst case not at all.

## What is validated

`ValidateSolverTypes` (in `source/src_util/param_parser.f90`) checks two
things for each component:

1. **Is the type dispatched at all?** The dispatchers are chains of
   independent `IF` blocks; a value that matches no branch used to be a
   *silent no-op*: the coarse solve did nothing and multigrid ran without a
   coarse correction, producing degraded or diverging solutions with no error
   message. This is now a hard error.

2. **Is the backing library compiled in?** Checked via the `#ifdef` symbols
   that the dispatchers themselves use.

### Supported types per component

The two components do **not** support the same set. The sets below were read
off the dispatchers, not assumed.

| Type | Solver                                   | Pressure (`mgCoarseGridSolver_P`) | Velocity (`mgCoarseGridSolver_U`) |
|-----:|------------------------------------------|:---------------------------------:|:---------------------------------:|
|  1   | SSOR / block-Jacobi (master BiCGStab)     | yes | yes |
|  2   | UMFPACK on the gathered coarse matrix     | yes | yes |
|  3   | BiCGStab on the /16 coarsened matrix      | yes | no  |
|  4   | UMFPACK on the /16 coarsened matrix       | yes | no  |
|  5   | MUMPS                                     | yes | yes |
|  7   | Hypre GMRES+BoomerAMG, full coarse matrix | yes | no  |
|  8   | Hypre with /16 geometric coarsening       | yes | no  |
|  9   | MKL PARDISO on the gathered coarse matrix | yes | no  |
| 10   | Hypre GMRES+BoomerAMG on the *finest* level (bypasses geometric MG, see `Solve_General_LinScalar`) | yes | no |

So:

- `Pres@MGCrsSolverType` accepts **1, 2, 3, 4, 5, 7, 8, 9, 10**
- `Velo@MGCrsSolverType` accepts **1, 2, 5**

### Library requirements

| Type(s) | Library     | Preprocessor symbol  | CMake option            |
|--------:|-------------|----------------------|-------------------------|
| 1–4     | vendored UMFPACK/AMD, always built | — | — |
| 5       | MUMPS       | `MUMPS_AVAIL`        | `-DUSE_MUMPS=ON`        |
| 7, 8, 10| Hypre       | `HYPRE_AVAIL`        | `-DUSE_HYPRE=ON`        |
| 9       | MKL PARDISO | `MKL_PARDISO_AVAIL`  | `-DUSE_MKL_PARDISO=ON`  |

## When it runs

`ValidateSolverTypes` is called from `Init_QuadScalar`
(`source/src_quadLS/QuadSc_handlers.f90`), immediately after
`GetVeloParameters` and `GetPresParameters` have populated
`QuadSc%prm%MGprmIn%CrsSolverType` and `LinSc%prm%MGprmIn%CrsSolverType`.

`Init_QuadScalar` is the single shared entry point that every `q2p1_*`
application uses to read the `Velo@`/`Pres@` sections, so no per-application
wiring is needed. It runs:

- on **all ranks**, on identical values (every rank parses the same file);
- **before** the type-7/8/10 multigrid hierarchy reconfiguration in
  `Init_QuadScalar_Structures*` (`NLMAX = NLMIN`, `bMasterTurnedON = .FALSE.`),
  so a doomed run never reshapes the master level hierarchy;
- long before partitioned matrix assembly and the first solve.

It is also re-run after a runtime parameter reload (`Reload_Velo` /
`Reload_Pres` in `_data/ProcCtrl.txt`, see `source/ProcCtrl.f90`), so a live
reload cannot smuggle in an unsupported type.

## What a failure looks like

The master rank writes the block to stdout, the `showID` rank mirrors it into
the protocol file, both flush, and every rank then calls
`MPI_Abort(MPI_COMM_WORLD, myErrorCode%SOLVER_TYPE_INVALID)`. The job exits
with status **58**.

Library missing:

```
==============================================================================
 FATAL: invalid solver type requested in _data/q2p1_param.dat
  Pres@MGCrsSolverType = 9 requires MKL PARDISO, which is not compiled into this binary.
  Rebuild FeatFloWer with -DUSE_MKL_PARDISO=ON or select a different solver type.
==============================================================================
```

```
==============================================================================
 FATAL: invalid solver type requested in _data/q2p1_param.dat
  Pres@MGCrsSolverType = 10 requires Hypre, which is not compiled into this binary.
  Rebuild FeatFloWer with -DUSE_HYPRE=ON or select a different solver type.
==============================================================================
```

Unknown / unsupported type (both components are reported in one go):

```
==============================================================================
 FATAL: invalid solver type requested in _data/q2p1_param.dat
  Velo@MGCrsSolverType = 4 is not a supported solver type.
  Valid Velo@MGCrsSolverType values are 1,2,5.
  Pres@MGCrsSolverType = 6 is not a supported solver type.
  Valid Pres@MGCrsSolverType values are 1,2,3,4,5,7,8,9,10.
==============================================================================
```

## Defence in depth at the dispatch sites

The pre-existing "`… is not available!`" messages at the dispatch sites are
now followed by `MPI_Abort` instead of `STOP`, and name the CMake option:

| File | Context | Was | Is |
|------|---------|-----|----|
| `source/src_quadLS/QuadSc_mg.f90` | `mgCoarseGridSolver_U`, type 5 without MUMPS | `STOP` | `MPI_Abort` |
| `source/src_quadLS/QuadSc_mg.f90` | `mgCoarseGridSolver_P`, type 5 without MUMPS | `STOP` | `MPI_Abort` |
| `source/src_quadLS/QuadSc_mg.f90` | `mgCoarseGridSolver_P`, type 9 without PARDISO (**master-only branch**) | `STOP` | `MPI_Abort` |
| `source/src_quadLS/QuadSc_mg.f90` | `mgCoarseGridSolver_P`, types 7/8 without Hypre | `STOP` | `MPI_Abort` |
| `source/src_quadLS/QuadSc_def.f90` | `Solve_General_LinScalar`, type 10 without Hypre | `STOP` | `MPI_Abort` |
| `source/src_quadLS/QuadSc_solver_coarse.f90` | `Setup_UMFPACK_CoarseSolver`, type 9 factorize without PARDISO (**master-only routine**) | `STOP` | `MPI_Abort` |

The two master-only conversions are the important ones: a bare `STOP` on rank
0 left every worker blocked inside MPI, and the job only died when the MPI
runtime eventually noticed — not a clean abort, and prone to hanging on MPI
stacks other than Open MPI.

These paths should now be unreachable, since validation rejects the same
configurations at startup; they remain as a safety net.

## Related change: the missing `Pres@` default

`GetVeloParameters` has always defaulted `Velo@MGCrsSolverType` to `1` when
the key is absent from `q2p1_param.dat`; `GetPresParameters` had no such
default, so the field kept its static-storage value `0`, which matches no
branch of `mgCoarseGridSolver_P`. Five shipped applications
(`q2p1_cc`, `q2p1_fac_bench3D`, `q2p1_hashgrid_test`, `q2p1_particle_tracer`,
`q2p1_particle_tracer_xse`) ship a `q2p1_param.dat` without
`Pres@MGCrsSolverType` and were therefore running the silent-no-op case.
`GetPresParameters` now mirrors the velocity default (`= 1`), which both
fixes that latent case and keeps those applications startable under the new
validation.

## Adding a new solver type

1. Implement the branch in `mgCoarseGridSolver_P` and/or
   `mgCoarseGridSolver_U` (`source/src_quadLS/QuadSc_mg.f90`), or in
   `Solve_General_LinScalar` for a fine-level solver.
2. Extend `CheckOneSolverType` in `source/src_util/param_parser.f90`: add the
   value to the `bKnown` expression of the right component, and — if it needs
   an optional library — add an availability block guarded by that library's
   `#ifdef` symbol together with its CMake option string.
3. Update the tables above.
