# q2p1_sse Initialization Flow

Primary files:

- `applications/q2p1_sse/q2p1_sse.f90`
- `applications/q2p1_sse/app_init.f90`
- `source/src_util/ReadExtrud3DParameters.f90`
- `source/src_quadLS/QuadSc_Sigma_User.f90`

## Entry Points

`q2p1_sse.f90` calls:

```text
init_q2p1_ext(ufile)
```

`init_q2p1_ext` calls:

```text
General_init_ext(79, log_unit)
Init_QuadScalar_Structures_sse
Init_GenLinSc_MULTIMATALPHA_Q1 or Init_LinScalar/InitCond_LinScalar
restart or initial-condition branch
InitOperators
```

## General_init_ext Responsibilities

`General_init_ext` is the largest application-local setup routine. It currently
combines several concerns:

- MPI initialization through `INIT_MPI` and `FindNodes`.
- Parameter loading through `GDATNEW("SimPar", 0)`.
- Partition reader include selection depending on PE/serial mode.
- Q2/P1 and scalar parameter initialization through `Init_QuadScalar`.
- E3D setup loading from `_data/Extrud3D.dat` via `ReadS3Dfile`.
- STL/OFF segment registration through `Setup_STL_Segments`.
- Automatic timestep estimation from geometry, speed, and rheology.
- Multigrid mesh allocation, `readTriCoarse`, and `refineMesh`.
- Parallel communication setup for coarse, refined, Q2, and scalar structures.
- CGAL/FBM rigid-body initialization through `init_fc_rigid_body`,
  `FBM_GetParticles`, and `FBM_ScatterParticles`.
- Boundary parametrization, coordinate prolongation, and initial smoothing.
- Dump-structure creation and coordinate transfer to the master rank.
- Optional PE communicator setup excluding rank 0.

## Restart Branches

`init_q2p1_ext` branches on `istart`:

- `istart == 0`: fresh start, create dump structures, deform mesh, initialize
  operators, initialize Q2 scalar fields.
- `istart == 1`: restart from same level and same partition count, optionally
  using list-file or MPI dump formats.
- `istart == 2`: restart from lower level, either through MPI dump prolongation
  or scalar solution read/prolongation.
- `istart == 3`: restart from same level with repartitioning through
  `SolFromFileRepart`.

`RestartFromMPILowerLevelSSE` handles the MPI dump lower-level path. It loads
dumped fields, creates prolongation/restriction structures, prolongates scalar
fields, rebuilds operators, updates FBM geometry, recomputes non-Newtonian
viscosity, and communicates screw/shell scalar fields back to the master rank.

## Important Couplings

- `General_init_ext` assumes `_data/q2p1_param.dat`, `_data/Extrud3D.dat`, and
  `_data/meshDir/file.prj` already exist. The Python launcher creates these.
- Rank 0 is treated specially. Several geometry and dump operations only run on
  compute ranks, while master-rank structures are reconstructed through
  communication helpers.
- Mesh and parametrization setup depends on values parsed into `mySigma`,
  `myProcess`, `mySetup`, and `myMultiMat`.

