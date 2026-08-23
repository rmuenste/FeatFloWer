# Coupled CC solver (`q2p1_cc`)

## Purpose and relation to Turek's CC approach

`applications/q2p1_cc` implements a coupled velocity-pressure defect-correction
method for the three-dimensional Q2/P1 discretization. In the terminology of
Turek's book, this is the CC (coupled) alternative to the PP
(pressure-correction/projection) solver used by `q2p1_fc_ext`:

- the nonlinear outer iteration assembles an Oseen problem;
- velocity and discontinuous P1 pressure remain in one saddle-point system;
- each local MPSC/Vanka patch contains the 27 Q2 velocity nodes of one
  hexahedron (81 velocity unknowns) and its four local P1 pressure unknowns;
- an 85-by-85 coupled patch is factorized with UMFPACK and used for defect
  correction;
- velocity and pressure are updated together, rather than by a separate
  pressure-Poisson projection.

The book discusses local/global MPSC and nonlinear defect correction as the
core ingredients of CC2D/CC3D. This application applies those ideas to the
codebase's conforming Q2/P1 element; it is not a literal copy of the book's
rotated-Q1/Q0 implementation.

## Code map

| File | Responsibility |
| --- | --- |
| `q2p1_cc.f90` | time loop and application lifecycle |
| `app_init.f90` | MPI, mesh hierarchy, parameters, boundary data |
| `q2p1_transport_cc.f90` | nonlinear Oseen/defect-correction loop |
| `q2p1_def_cc.f90` | Newtonian/generalized-stress matrices, coupled defect, 85-by-85 patch extraction, coarse matrix |
| `q2p1_mg_cc.f90` | coupled Vanka smoother and optional multilevel cycles |
| `q2p1_cc_umfpacksolver.f90` | local and coarse direct-solver adapters |
| `assemblies_cc.f`, `iso_assemblies.f` | Q2/P1 coupling and cylinder-force integration |
| `postprocessing.f90` | VTK, dump, timing, and process control |

The algebraic system used in a Newtonian time step is, up to the sign convention
of the stored `B` blocks,

```text
[ M + dt (nu D + K(u^k))   -dt B ] [delta u] = [r_u]
[             B^T                0 ] [delta p]   [r_p]
```

`AAij` holds the matrix used by the coupled correction. `Aij` is used to
evaluate the nonlinear defect. Dirichlet rows are filtered before patch
extraction.

### Newtonian Oseen defect that was repaired

This equation is directly connected to the central defect in the original
implementation. The complete Newtonian branch that should have assembled its
velocity operator was commented out. Consequently, choosing
`FlowType = Newtonian` did not form the diagonal velocity blocks at all: the
coupled solver received empty/zero blocks instead of the Oseen operator. The
only apparent workaround was to enter the non-Newtonian path with a power-law
exponent of one.

The restored Newtonian Navier--Stokes branch in `q2p1_def_cc.f90` always forms
the pure Oseen defect operator, for each velocity component,

```text
Aii  = M + dt (D + K(u^k))
```

Here `M` is the transient mass matrix, `D` is the Newtonian viscous operator,
and `K(u^k)` is the convection operator linearized at the current nonlinear
iterate. `Aii` evaluates the nonlinear defect, so the root of the defect
correction is always the plain Navier--Stokes solution, independent of the
correction operator chosen below. The Stokes branch similarly forms `M + dt D`
when convection is disabled.

### Newton treatment of the correction operator

The correction matrix `AA` (extracted into the coupled patch system and the
coarse matrix) is selected by `CCuvwp@NewtonTreatment`:

| Value | Correction operator | Turek reference |
| --- | --- | --- |
| `Off` | `AAii = Aii`, no reactive blocks | favored fixed-point/Oseen preconditioner `S^F`, eq. (3.174) |
| `Diagonal` (default) | `AAii = Aii + dt alpha barMii` | hybrid: diagonal part of the Newton reaction term only |
| `Full` | additionally `AAij = dt alpha barMij` for all six cross blocks | full Newton derivative, eqs. (3.167)--(3.169) |

`barMij` is the assembled reaction tensor `int rho (du_i/dx_j) phi phi` — the
discretization of the Newton term `(delta u . grad) u`. These blocks are *not*
zero for a Newtonian fluid; `Off` and `Diagonal` deliberately omit some or all
of them from the correction operator, which is legitimate in the
defect-correction framework (the converged root is defined by `A`, not `AA`)
but affects the nonlinear contraction rate. Turek (Sec. 3.3.1) warns that the
reactive blocks resist robust multigrid smoothing; in this solver the Vanka
patches are solved by direct LU, so `Full` is worth measuring. All nine `AA`
blocks are always allocated because patch extraction, coarse assembly, and
Vanka smoothing share a uniform nine-block velocity layout.

The blending factor `alpha` (`CCuvwp@Alpha`, adapted each nonlinear iteration)
scales the retained reactive blocks inside `AA`. Note that this is an operator
blending between Oseen (`alpha=0`) and Newton (`alpha=1`); it is not Turek's
adaptive step-length damping `omega` of eqs. (3.175)--(3.177) — the computed
correction is applied undamped. With `NewtonTreatment = Off`, `alpha` only
influences the multigrid stopping criterion.

## Rehabilitated implementation

The original application was a WIP and could not run the cylinder case. The
following faults were corrected:

1. The application is now built independently of MUMPS. UMFPACK provides the
   local patch solves and a pressure-pinned coarse-grid fallback when MUMPS is
   unavailable.
2. Core FeatFloWer sources are no longer compiled a second time into the
   executable. They come from `FF_APPLICATION_LIBS`, avoiding duplicate module
   state and COMMON blocks.
3. The missing Newtonian Oseen operator assembly described above was restored.
   Previously selecting `FlowType = Newtonian` left the velocity blocks empty
   and required the accidental workaround `nonNewtonian` plus power-law
   exponent one.
4. All nine velocity blocks are allocated. The extraction and Vanka code use a
   uniform 4-by-4 block layout even when the six Newtonian off-diagonal blocks
   are zero.
5. Newtonian cylinder forces are evaluated; the former code only assigned the
   force vector in the non-Newtonian branch.
6. Mesh-volume arrays now include the sentinel entry required by `SETARE`,
   removing an initialization heap overwrite.
7. Operator cleanup checks the parent matrix allocation before dereferencing
   component arrays, so `S=0` is valid.
8. Level-dependent boundary masks are used during multilevel operations.
9. Non-Newtonian viscosity projection is forced onto the solver level instead
   of the worker-only output level.
10. A clean-build Fortran module dependency hidden by `defs_include.h` is now
    explicit.
11. The build stages the parameter file, `MG.dat`, METIS, and partitioner. The
    obsolete FAC3Ds Python test driver was replaced by a current launcher.

## Cycle choice

`CCuvwp@MGCycType` accepts the historical `F`, `V`, and `W` multilevel modes.
It also accepts `S`, the supported local MPSC/Vanka-only cycle. Use `S` for the
cylinder benchmark.

On the two-level cylinder mesh, the legacy coarse transfer does not provide a
useful contraction: an isolated exact coarse correction increased the initial
defect from approximately `1.074e-3` to `1.510e-3`. The UMFPACK solution of the
coarse matrix was checked against the assembled coordinate matrix and had a
residual below `1e-8`; the remaining issue is the legacy intergrid correction,
not the direct solve. `S` therefore preserves the working coupled algorithm
without allowing the experimental coarse transfer to undo local Vanka work.

### Coarse pressure handling in the UMFPACK fallback

`SolveCoarseWithUMFPACK` pins the first coarse pressure DOF *unconditionally*.
This was investigated against the shared solver's practice
(`Setup_UMFPACK_CoarseSolver` applies singular-pressure handling only for
`bNoOutflow`): the two operators differ. The CC coarse saddle-point system
couples the pressure only through the element-wise discontinuous-P1
gradient/divergence blocks, so a globally constant pressure is a nullspace
vector regardless of outflow boundaries — an unpinned factorization fails on
the very first coarse solve (verification residual `~1e15` at rhs scale
`~1e-3`), while the pinned solution satisfies every unmodified equation of
the original system to relative accuracy below `1e-8`. The PP
pressure-Poisson matrix is a different operator in which the outflow boundary
condition does enter, which is why the shared solver's conditional handling
does not transfer to the CC system.

The fallback's verification residual is measured *relative to the incoming
right-hand side*: in a diverging `F`/`V`/`W` cycle the rhs grows without
bound, and an absolute check would misattribute that outer divergence to the
direct solve. The routine also counts structurally zero matrix rows and
reports them, so a defective coarse assembly is diagnosed explicitly instead
of surfacing as an unexplained solver failure.

The shipped cylinder deck solves the *stationary* Re-20 problem directly:

```text
SimPar@MatrixRenewal = M1D1K3S0C0
SimPar@FlowType = Newtonian
SimPar@SteadyState = Yes
SimPar@MaxNumStep = 1
CCuvwp@NewtonTreatment = Diagonal
CCuvwp@StrictConvergence = Yes
CCuvwp@NLmax = 25
CCuvwp@Stopping = 1d-8
CCuvwp@MGMaxIterCyc = 100
CCuvwp@MGSmoothSteps = 1
CCuvwp@MGCycType = S
CCuvwp@MGRelaxPrm = 0.2
```

`SteadyState = Yes` does two things: the shared mass-matrix assembly zeroes
`MMat`, so the assembled operator is the stationary Oseen system `dt (D + K)`
(with `dt` a pure scale factor), and the nonlinear stopping test uses the
*absolute* defect `CCuvwp@Stopping = 1d-8`. On the level-1/2 `2Dbench` mesh
with three worker partitions this converges from a cold start in about 13
nonlinear iterations (defect `1.07e-3` down to below `1e-8`, contraction
roughly `0.4` per iteration) in under a minute, exits with status 0, and
prints no warnings. `NewtonTreatment = Full` saves about one iteration at a
higher cost per iteration; `Off` and `Diagonal` are equivalent within one
iteration on this case.

Each nonlinear loop that reaches `CCuvwp@NLmax` without meeting
`CCuvwp@Stopping` prints a per-step warning with the achieved and required
criteria, and `sim_finalize` summarizes how many loops were exhausted together
with the worst achieved criterion. In strict mode
(`CCuvwp@StrictConvergence = Yes`, the default) such a run additionally
replaces the success banner with a failure message and exits with status 1;
set `CCuvwp@StrictConvergence = No` to keep the warnings but exit with
status 0.

A pseudo-transient run remains possible (`SteadyState = No`, e.g. `dt=1` with
ten implicit steps, or `dt=0.01` when the transient history itself is the
comparison target). Note that its stopping test is a *relative per-step*
reduction that is re-baselined every time step; it becomes increasingly hard
to meet as the flow approaches steady state, so NLmax warnings are expected
in the tail of such runs, the absolute defect is the meaningful convergence
measure there, and strict mode should be disabled or the tolerance chosen
accordingly.

## Build and run the cylinder benchmark

Configure exactly as for `q2p1_fc_ext`; `Q2P1_MESH_DIR` must contain
`2D_FAC/2Dbench.prj`:

```bash
cmake -S . -B build-cc-check \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_APPLICATIONS=ON \
  -DQ2P1_MESH_DIR=/path/to/mesh_repo
cmake --build build-cc-check --target q2p1_cc -- -j8
cd build-cc-check/applications/q2p1_cc
./q2p1_cc.py --num-processors 4
```

The launcher creates three worker partitions and runs four MPI ranks (one
master plus three workers). If `_mesh/NEWFAC` is already current:

```bash
python3 ./q2p1_cc.py -n 4 --skip-partition
```

Manual equivalents are:

```bash
LD_LIBRARY_PATH="$PWD:$LD_LIBRARY_PATH" \
python3 ./PyPartitioner.py 3 1 1 NEWFAC _adc/2D_FAC/2Dbench.prj
mpirun -np 4 ./q2p1_cc > run_q2p1_cc_np4.log 2>&1
```

A fully successful run ends with exit status 0 and
`CC3D_iso_adaptive has successfully finished.` — with the shipped stationary
deck this is the out-of-the-box outcome. If any nonlinear loop was exhausted,
strict mode (see above) prints a failure message instead and the run exits
with status 1. The normalized cylinder values are printed as
`BenchForce: time drag lift` in either case.

## Validation status and remaining work

Verified on the level-1/2 `2Dbench` mesh with three worker partitions:

- clean build with GCC/OpenMPI and no MUMPS;
- partitioning and MPI initialization;
- Newtonian `D`-matrix assembly and coupled block extraction;
- repeated 85-by-85 Vanka factorization/solve;
- coupled defect reduction;
- Newtonian drag/lift integration;
- VTK and dump output;
- successful finalization without heap corruption or segmentation faults.

The first successful runs recorded during rehabilitation were:

| Run | Time | Drag | Lift | Interpretation |
| --- | ---: | ---: | ---: | --- |
| Physical-step smoke test, `dt=0.01`, four nonlinear corrections | `0.01` | `139.42` | `0.476` | Intentionally under-converged impulsive first step; not a benchmark endpoint |
| Pseudo-transient smoke test, `dt=1`, four nonlinear corrections | `1` | `6.9194472` | `0.013909134` | First successful complete CC step |
| Ten-step pseudo-transient validation, `dt=1` | `10` | `5.5880793` | `0.0099092682` | Under-converged endpoint (NLmax=4 exhausted each step); approximately 163 seconds |
| Stationary solve, shipped deck, absolute defect below `1d-8` | steady | `5.5981453` | `0.0099736961` | Converged benchmark result; 13 nonlinear iterations, about 50 seconds, exit status 0 |

The PP reference from the from-scratch guide is drag `5.601296` and lift
`0.00994712` near `t=10`. The converged stationary CC result differs by
approximately `0.06%` in drag and `0.26%` in lift; at absolute defect `1d-8`
the forces are settled to about four digits, so these residual differences
reflect the discretization/solver differences between the CC and PP paths,
not under-convergence. Tightening from `1d-6` to `1d-8` still moved the drag
by about `0.08%`, which is why the shipped tolerance is `1d-8`.

Known follow-up items are:

- repair or replace the legacy Q2/P1 intergrid transfer before recommending
  `F`, `V`, or `W` cycles;
- reduce the cost of repeated patch solves (factor reuse/coloring or a modern
  block preconditioner);
- convert the external-mesh cylinder run into a bounded CTest fixture;
- reconcile legacy FEAT COMMON-block size warnings emitted by `feat2d`. They
  also occur in old application code and should not be ignored in a broader
  modernization effort.
