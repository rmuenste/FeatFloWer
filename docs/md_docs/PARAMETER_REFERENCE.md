# FeatFloWer Parameter Reference

**Version**: 1.2
**Last Updated**: 2026-01-21

## Introduction

This document provides a comprehensive reference for all parameters in the `q2p1_param.dat` configuration file used by FeatFloWer applications. Understanding these parameters is essential for correctly configuring simulations of fluid dynamics problems using the Q2/P1 finite element discretization.

## Parameter File Format

Parameters in `q2p1_param.dat` follow the format:
```
Category@ParameterName = value
```

### Parameter Categories

- **SimPar@**: Simulation control, I/O, and time stepping
- **Velo@**: Velocity solver configuration (multigrid)
- **Pres@**: Pressure solver configuration (multigrid)
- **Prop@**: Physical properties (fluid properties, materials)

### Common Value Types

- **Integer**: `1`, `100`, `1000000`
- **Real (Fortran double precision)**: `1d0`, `0.01d0`, `1d-6`
- **String**: `"filename"`, `Newtonian`, `BE`
- **Yes/No**: `Yes`, `No` (case-sensitive)
- **Array**: `0d0,0d0,-9.81d0` (comma-separated, no spaces)

---

## Quick Reference Tables

### SimPar@ - Simulation Control Parameters

**Source**: `source/src_util/param_parser.f90` (GDATNEW subroutine, lines 542-1076)
**Note**: GDATNEW was migrated from `source/Init.f90` to the centralized `param_parser` module for better maintainability.

| Parameter | Type | Units | Description | Example |
|-----------|------|-------|-------------|---------|
| **Mesh & Domain** |
| MeshFolder | string | - | Name of the mesh directory | `"NEWFAC"` |
| SubMeshNumber | integer | - | Number of subpartitions for multilevel partitioning | `1` |
| ProjectFile | string | - | Path to project file | `"_adc/2D_FAC/2Dbench.prj"` |
| MinMeshLevel | integer | - | Minimum mesh refinement level | `1` |
| MaxMeshLevel | integer | - | Maximum mesh refinement level | `2` |
| **Time Discretization** |
| TimeScheme | string | - | Time scheme: `BE` (Backward Euler), `CN` (Crank-Nicolson), `FE` (Forward Euler) | `BE` |
| TimeStep | real | s | Time step size | `0.01d0` |
| TimeAdaptivity | Yes/No | - | Enable adaptive time stepping | `No` |
| MinTimeAdapt | real | s | Minimum adaptive time step | `0.00001d0` |
| MaxTimeAdapt | real | s | Maximum adaptive time step | `0.050d0` |
| StartSimTime | real | s | Starting simulation time | `0d0` |
| MaxSimTime | real | s | Maximum simulation time | `10.0001d0` |
| MaxNumStep | integer | - | Maximum number of time steps | `1000000` |
| **File I/O** |
| ParticleFile | string | - | Path to particle configuration file | `"_data/particles.dat"` |
| ProtocolFile | string | - | Path to output protocol/log file | `"_data/prot.txt"` |
| StartFile | string | - | Restart file path for continuing simulation | `"_dump/01"` |
| SolFile | string | - | Solution output file path | `"_sol/1/mySol"` |
| BackUpFreq | integer | - | Backup/checkpoint frequency (timesteps) | `1` |
| BackUpNum | integer | - | Number of backup files to keep | `10` |
| OutputFreq | real | s | Output file write frequency | `0.50d0` |
| OutputFields | string | - | Comma-separated list of fields to export | `"Pressure_V,Velocity,Viscosity"` |
| OutputFormat | string | - | Output file format: `VTK` or `GMV` | `"VTK"` |
| OutputLevel | string/int | - | Mesh level for output: `1`-`N`, `MAX`, `MAX+1`, `MAX-1` | `MAX+1` |
| **Solver Control** |
| MatrixRenewal | string | - | Matrix renewal scheme (`M#D#K#S#C#` format) | `M1D1K3S0C1` |
| FlowType | string | - | Flow model: `Newtonian` or non-Newtonian | `Newtonian` |
| NoOutflow | Yes/No | - | Apply no-outflow boundary condition | `No` |
| NS_Stabilization | Yes/No | - | Enable Navier-Stokes stabilization | `No` |
| **Physics Models** |
| Tracer | Yes/No | - | Enable tracer transport equation | `No` |
| bViscoElastic | Yes/No | - | Enable viscoelastic model | `No` |
| **Fictitious Boundary Method (FBM)** |
| skipFBMForce | Yes/No | - | Skip FBM force computation (for testing/debugging) | `No` |
| skipFBMDynamics | Yes/No | - | Skip FBM dynamics update (for testing/debugging) | `No` |
| **Mesh Processing** |
| Umbrella | integer | - | Umbrella smoothing steps in `InitCond_QuadScalar()` (set to 0, see §Umbrella Mesh Smoothing) | `0` |
| InitUmbrella | integer | - | Umbrella smoothing steps in `General_init_ext()` (set to 0 on restart, see §Umbrella Mesh Smoothing) | `20` |
| UmbrellaStepM | integer | - | Main umbrella smoothing steps | - |
| UmbrellaStepL | integer | - | Umbrella steps per level | - |
| LoadAdaptedMesh | string | - | Path to pre-adapted mesh file | `"_dump/msh"` |
| **Advanced** |
| StartingProc | integer | - | Restart mode: 0=fresh start, 1=same level/partitions, 2=lower level, 3=repartition (see §Umbrella Mesh Smoothing for restart workflow) | `0` |
| ElemTrans | integer | - | Element transformation type (1 or 2) | `2` |
| ProlongationDirection | integer | - | Mesh prolongation: 0=standard, 1=X-cyl, 2=Y-cyl, 3=Z-cyl | `0` |
| CGALtoRealFactor | real | - | Conversion factor for CGAL geometry | - |
| BoundaryCheck | Yes/No | - | Enable boundary checking | `No` |
| **MPI Communication** |
| aSynchComm | Yes/No | - | Enable asynchronous MPI communication | - |
| CommSwitch | integer | - | Communication switching parameter | - |
| **Benchmarks** |
| GammaDot | real | 1/s | Characteristic shear rate | - |
| AlphaRelax | real | - | Relaxation parameter | - |
| RadParticle | real | m | Particle radius | - |
| Bench_U_mean | real | m/s | Benchmark: mean velocity | - |
| Bench_H | real | m | Benchmark: height parameter | - |
| Bench_D | real | m | Benchmark: diameter parameter | - |
| Bench_Sc_U | real | - | Benchmark: velocity scaling | - |
| Bench_Sc_Mu | real | - | Benchmark: viscosity scaling | - |
| Bench_Sc_a | real | - | Benchmark: scaling parameter | - |
| ViscoElasticBench | integer | - | Viscoelastic benchmark type (2=2D, 3=3D) | - |
| SteadyState | Yes/No | - | Steady-state computation mode | - |
| ReferenceFrame | Yes/No | - | Enable reference frame | - |

---

### Velo@ - Velocity Solver Parameters

**Source**: `source/src_quadLS/QuadSc_def.f90` (GetVeloParameters, lines 4235-4409)

| Parameter | Type | Default | Description | Example |
|-----------|------|---------|-------------|---------|
| **Solver Selection** |
| iMass | integer | - | Mass matrix type | - |
| SolvType | string | - | Solver type: `Jacobi`, `BiCGSTAB` | - |
| **Convergence Criteria** |
| defCrit | real | - | Defect criterion for convergence | `1d-6` |
| epsCrit | real | - | Epsilon criterion | `1d-5` |
| MinDef | real | - | Minimum defect threshold | `1d-10` |
| **Nonlinear Iteration** |
| NLmin | integer | - | Minimum nonlinear iterations | `1` |
| NLmax | integer | - | Maximum nonlinear iterations | `5` |
| Alpha | real | - | Nonlinear relaxation parameter | `1d0` |
| **Multigrid Hierarchy** |
| MGMinLev | integer | - | Multigrid minimum level | `2` |
| MGMedLev | integer | - | Multigrid medium level | `2` |
| **Multigrid Cycling** |
| MGMinIterCyc | integer | - | Minimum MG iterations per cycle | `1` |
| MGMaxIterCyc | integer | - | Maximum MG iterations per cycle | `1` |
| MGCycType | char | - | MG cycle type: `F` (full), `V`, `W` | `F` |
| **Multigrid Smoothing** |
| MGSmoothSteps | integer | - | MG smoothing steps (pre/post smoother) | `2` |
| MG_VANKA | integer | - | VANKA smoother parameter | - |
| MGRelaxPrm | real | `0.66` | MG relaxation parameter (ω) | `0.7d0` |
| **Coarse Grid Solver** |
| MGIterCoarse | integer | - | MG iterations on coarsest grid | `50` |
| MGCrsSolverType | integer | `1` | Coarse solver: 1=SSOR+blockJacobi, 5=MUMPS | `1` |
| MGCrsRelaxPrm | real | - | Coarse grid relaxation parameter | - |
| MGCrsRelaxParPrm | real | - | Coarse grid parallel relaxation parameter | - |
| **Multigrid Convergence** |
| MGDefImprCoarse | real | - | Defect improvement on coarse grid | `0.1d0` |
| MGCriterion1 | real | - | MG convergence criterion 1 | `0.9d-1` |
| MGCriterion2 | real | - | MG convergence criterion 2 | `0.9d-1` |

---

### Pres@ - Pressure Solver Parameters

**Source**: `source/src_quadLS/QuadSc_def.f90` (GetPresParameters, lines 4447-4569)

| Parameter | Type | Default | Description | Example |
|-----------|------|---------|-------------|---------|
| **Multigrid Hierarchy** |
| MGMinLev | integer | - | Multigrid minimum level | `1` |
| MGMedLev | integer | - | Multigrid medium level | `1` |
| **Multigrid Cycling** |
| MGMinIterCyc | integer | - | Minimum MG iterations per cycle | `2` |
| MGMaxIterCyc | integer | - | Maximum MG iterations per cycle | `20` |
| MGCycType | char | - | MG cycle type: `F` (full), `V`, `W` | `F` |
| **Multigrid Smoothing** |
| MGSmoothSteps | integer | - | MG smoothing steps | `16` |
| MGRelaxPrm | real | `0.66` | MG relaxation parameter (ω) | - |
| **Coarse Grid Solver** |
| MGIterCoarse | integer | - | MG iterations on coarsest grid | `799` |
| MGCrsSolverType | integer | - | Coarse grid solver type | `4` |
| MGCrsRelaxPrm | real | - | Coarse grid relaxation parameter | - |
| **Multigrid Convergence** |
| MGDefImprCoarse | real | - | Defect improvement on coarse grid | `1d-3` |
| MGCriterion1 | real | - | MG convergence criterion 1 | `1d-3` |
| MGCriterion2 | real | - | MG convergence criterion 2 | `1d-10` |

---

### Prop@ - Physical Properties

**Source**: `source/src_quadLS/QuadSc_def.f90` (GetPhysiclaParameters, lines 4573-4689)

> ⚠️ **CRITICAL NOTE**: `Prop@Viscosity` specifies **KINEMATIC viscosity ν [m²/s]**, NOT dynamic viscosity μ [Pa·s]!

| Parameter | Type | Units | Description | Example |
|-----------|------|-------|-------------|---------|
| **Basic Fluid Properties** |
| Gravity | real(3) | m/s² | Gravity vector [x, y, z] | `0d0,0d0,-9.81d0` |
| Density | real(2) | kg/m³ | Fluid densities (phase 1, phase 2) | `1000d0,1d0` |
| **Viscosity** |
| Viscosity | real(2) | **m²/s** | **KINEMATIC viscosity ν** (phase 1, phase 2) | `1d-6,1.5d-5` |
| DynVisc | real(2) | **Pa·s** | **DYNAMIC viscosity μ** (auto-converted to kinematic) | `1d-3,1d0` |
| DiffCoeff | real(2) | m²/s | Diffusion coefficients (phase 1, phase 2) | `1d-4,1d-4` |
| **Interfacial Properties** |
| Sigma | real | N/m | Surface tension coefficient | `0.073d0` |
| DiracEps | real | - | Dirac epsilon for interface smoothing | - |
| **Non-Newtonian Flow** |
| PowerLawExp | real | - | Power law exponent for non-Newtonian fluids | `1.0` |
| **Viscoelastic Properties** |
| ViscoLambda | real | s | Viscoelastic relaxation time λ | - |
| ViscoAlphaImp | real | - | Viscoelastic alpha (implicit) | - |
| ViscoAlphaExp | real | - | Viscoelastic alpha (explicit) | - |
| ViscoModel | integer | - | Viscoelastic model type | - |
| **Navier-Stokes Stabilization** |
| NS_StabAlpha_Imp | real | - | NS stabilization alpha (implicit) | - |
| NS_StabAlpha_Exp | real | - | NS stabilization alpha (explicit) | - |
| **Force Scaling** |
| ForceScale | real(6) | - | Force scaling factors | - |
| **Transport Equation** |
| nTPSubSteps | integer | - | Transport equation substeps | - |
| nTPFSubSteps | integer | - | Transport fluid substeps | - |
| nTPIterations | integer | - | Transport iterations | - |
| **Material** |
| Material | string | - | Material name | - |

---

## Critical Parameter Details

### ⚠️ Prop@Viscosity: KINEMATIC VISCOSITY!

**This is the most commonly misunderstood parameter in FeatFloWer.**

`Prop@Viscosity` specifies **kinematic viscosity ν [m²/s]**, not dynamic viscosity μ [Pa·s].

**Recommendation**: Use **`Prop@DynVisc`** to provide the dynamic viscosity μ [Pa·s] directly. The solver will automatically convert it to kinematic viscosity using the provided `Prop@Density`.

**Relationship**:
```
ν = μ / ρ
```
where:
- ν = kinematic viscosity [m²/s]
- μ = dynamic viscosity [Pa·s = kg/(m·s)]
- ρ = density [kg/m³]

**Example for water at 20°C**:
- Dynamic viscosity: μ = 1.002×10⁻³ Pa·s
- Density: ρ = 998.2 kg/m³
- Kinematic viscosity: ν = μ/ρ = 1.004×10⁻⁶ m²/s

**Preferred parameter file (using DynVisc)**:
```
Prop@Density = 998.2d0,1d0
Prop@DynVisc = 1.002d-3,1d0  ! μ in [Pa·s], auto-converted to ν
```

**Alternative (using manual conversion)**:
```
Prop@Density = 998.2d0,1d0
Prop@Viscosity = 1.004d-6,1d0  ! This is ν, NOT μ!
```

**Evidence from source code**:
- `source/src_pp3d/init1.f:1723` uses parameter **NU** (standard notation for kinematic viscosity)
- `source/src_quadLS/QuadSc_force_torque_calc.f90:464` computes dynamic viscosity as:
  ```fortran
  ! daux = 1/(mu*rho*OMEGA*R_i^2*H)
  daux = 1d0/(Properties%Viscosity(1)*Properties%Density(1)*...)
  ```
  This confirms: `Viscosity × Density = μ`, therefore `Viscosity = ν`

---

## Common Fluid Property Values

### Water (at 20°C, 1 atm)
```
Prop@Density = 998.2d0,1d0
Prop@Viscosity = 1.004d-6,1d0  ! Kinematic viscosity
```

### Air (at 20°C, 1 atm)
```
Prop@Density = 1.204d0,1d0
Prop@Viscosity = 1.516d-5,1d0  ! Kinematic viscosity
```

### Glycerin (at 20°C)
```
Prop@Density = 1260d0,1d0
Prop@Viscosity = 1.19d-3,1d0  ! Kinematic viscosity (highly viscous!)
```

### Oil (typical mineral oil at 20°C)
```
Prop@Density = 900d0,1d0
Prop@Viscosity = 1.0d-4,1d0  ! Kinematic viscosity
```

---

## Multigrid Parameters Explained

### Understanding Multigrid Levels

- **MinLev**: Coarsest grid level (smallest number of DOF)
- **MedLev**: Working level for the solver
- **MaxLev**: Finest grid level (largest number of DOF)

Typical configuration:
```
SimPar@MinMeshLevel = 1
SimPar@MaxMeshLevel = 3
Velo@MGMinLev = 1
Velo@MGMedLev = 3
```

### Cycle Types

- **V-cycle** (`V`): Cheaper per iteration, more iterations needed
- **W-cycle** (`W`): More expensive per iteration, fewer iterations needed
- **F-cycle** (`F`): Full multigrid, good balance (recommended for most cases)

### Coarse Grid Solvers

**Velo@MGCrsSolverType** and **Pres@MGCrsSolverType**:
- `1`: SSOR + block Jacobi (default, fast, works on all systems)
- `5`: MUMPS direct solver (requires Intel compiler, more robust for difficult problems)

---

## Time Discretization Schemes

### SimPar@TimeScheme

- **`BE`** (Backward Euler): θ = 1.0
  - **Unconditionally stable**
  - First-order accurate in time
  - **Recommended for most simulations**
  - Good for stiff problems

- **`CN`** (Crank-Nicolson): θ = 0.5
  - Conditionally stable
  - Second-order accurate in time
  - Can produce oscillations for high Reynolds numbers
  - Use with caution for convection-dominated flows

- **`FE`** (Forward Euler): θ = 0.0
  - **Conditionally stable** (very restrictive timestep)
  - First-order accurate in time
  - Rarely used (only for special cases)

---

## Matrix Renewal Scheme

### SimPar@MatrixRenewal

Format: `M#D#K#S#C#` where each number controls renewal frequency:

- **M**: Mass matrix
- **D**: Diffusion matrix (viscous term)
- **K**: Convection matrix (advection term)
- **S**: Stokes matrix
- **C**: Gradient/Divergence matrices

**Values**:
- `0`: Never renew (fastest, only for linear problems)
- `1`: Renew every timestep (safest for nonlinear problems)
- `N`: Renew every N timesteps

**Example**: `M1D1K3S0C1`
- Mass matrix: every timestep
- Diffusion: every timestep
- Convection: every 3 timesteps
- Stokes: never
- Gradient/Divergence: every timestep

**Recommendations**:
- **Nonlinear problems**: `M1D1K1S0C1` (safest)
- **Linear/Stokes flow**: `M0D0K0S0C0` (fastest)
- **Moderate nonlinearity**: `M1D1K3S0C1` (good balance)

---

## Umbrella Mesh Smoothing

### SimPar@InitUmbrella and SimPar@Umbrella

**⚠️ CRITICAL**: These two parameters control mesh smoothing at different stages of execution. Understanding when each is applied is essential for correct restart behavior.

#### Execution Stages

**Fresh Start (StartingProc=0)**:
1. **`General_init_ext()`** (mesh initialization):
   - Applies mesh parametrization (projects nodes to boundary surfaces)
   - Applies umbrella smoothing: **`SimPar@InitUmbrella`** steps
   - Generates Q2 mesh nodes (edges, faces, element centers) from P1 vertices
   - Writes `initial_mesh.pvtu` (if debug output enabled)

2. **`InitCond_QuadScalar()`** (solution initialization):
   - Sets up initial velocity/pressure fields
   - Applies umbrella smoothing: **`SimPar@Umbrella`** steps
   - First output written with this final mesh

**Restart (StartingProc=1,2,3)**:
1. **`General_init_ext()`** (mesh re-initialization):
   - Re-applies parametrization and `InitUmbrella` smoothing ← **PROBLEM IF >0**
   - Re-generates Q2 mesh nodes from P1 vertices

2. **`SolFromFile()`** (load saved state):
   - Loads **P1 vertex coordinates only** from dump files
   - Loads velocity and pressure fields
   - **Overwrites** P1 vertices in mesh
   - **Q2 nodes (edges/faces/elements) remain as generated in step 1**

3. **`InitCond_QuadScalar()`**:
   - **Skipped during restart** (solution already loaded)
   - `Umbrella` smoothing is NOT applied

#### Understanding Q2 Mesh Coordinate Storage

**Key Insight**: Dump files store **only P1 vertex coordinates**, not the full Q2 nodal structure.

- **P1 vertices**: Corner nodes of hexahedral elements (stored in dumps)
- **Q2 nodes**: Edge midpoints, face centers, element centers (regenerated on each run)

On restart, Q2 nodes are **regenerated** from P1 vertices via `SetUp_myQ2Coor()`:
```fortran
! Edge nodes: midpoint of two vertices
PX = 0.5 * (vertex1 + vertex2)

! Face nodes: center of four vertices
PX = 0.25 * (vertex1 + vertex2 + vertex3 + vertex4)

! Element nodes: center of eight vertices
PX = 0.125 * (sum of 8 vertices)
```

**Source**: `source/src_quadLS/QuadSc_def.f90:4015`

#### Problem: Mesh Smoothing Mismatch on Restart

If smoothing parameters are not adjusted for restart, the mesh undergoes different processing:

**Fresh Start:**
```
P1 vertices → parametrization → InitUmbrella → Q2 generation → Umbrella → FINAL MESH
```

**Restart (with InitUmbrella>0):**
```
Loaded P1 → parametrization → InitUmbrella → Q2 generation → NO Umbrella
                                   ↑
                           EXTRA SMOOTHING APPLIED!
```

The restarted mesh has undergone **different smoothing** than the dumped mesh, causing the loaded velocity field to be inconsistent with mesh geometry → solver divergence.

#### ✅ Recommended Workflow

**Stage 1: Fresh Start**
```fortran
SimPar@StartingProc = 0
SimPar@InitUmbrella = 20    ! Apply smoothing during initialization
SimPar@Umbrella = 0         ! No additional smoothing (recommended)
```

**Stage 2: Restart from Dump**
```fortran
SimPar@StartingProc = 1
SimPar@StartFile = "_dump/10"
SimPar@InitUmbrella = 0     ! ⚠️ CRITICAL: Disable, coordinates in dump already smoothed
SimPar@Umbrella = 0         ! Already disabled
```

#### Alternative: Two-Stage Smoothing (Advanced)

If you need smoothing at both stages (rare):

**Fresh Start**:
```fortran
SimPar@InitUmbrella = 10    ! First stage smoothing
SimPar@Umbrella = 10        ! Second stage smoothing
```

**Restart**:
```fortran
SimPar@InitUmbrella = 0     ! Disable both stages
SimPar@Umbrella = 0
```

The dump files contain P1 vertex coordinates **after both smoothing stages**. Since Q2 nodes are regenerated from these P1 vertices, any smoothing during restart must be disabled to match the original mesh.

#### Why Two Smoothing Stages?

Historically, umbrella smoothing was applied:
1. **`InitUmbrella`**: After parametrization, before multilevel hierarchy generation
2. **`Umbrella`**: After full solver setup, with complete Q2 coordinate structures

Modern practice: **Use only `InitUmbrella`** and set `Umbrella = 0` to avoid confusion and ensure restart consistency.

#### Source Code Locations

- **InitUmbrella**: Applied in `applications/*/app_init.f90` (`General_init_ext`), line ~334
- **Umbrella**: Applied in `source/src_quadLS/QuadSc_initialization.f90` (`InitCond_QuadScalar`), line ~1221
- **Q2 regeneration**: `source/src_quadLS/QuadSc_def.f90:4015` (`SetUp_myQ2Coor`)
- **Coordinate restoration**: `source/initialization/app_initialization.f90:151,189,233,305` (restores P1 vertices only)

#### Common Pitfalls

❌ **Wrong**: Using same parameters for fresh start and restart
```fortran
! Fresh start
SimPar@InitUmbrella = 20
SimPar@Umbrella = 0

! Restart (WRONG - will re-smooth mesh!)
SimPar@InitUmbrella = 20    ← This will apply smoothing before loading dump!
SimPar@Umbrella = 0
```

✅ **Correct**: Adjust parameters for restart
```fortran
! Fresh start
SimPar@InitUmbrella = 20
SimPar@Umbrella = 0

! Restart (CORRECT)
SimPar@InitUmbrella = 0    ← Disabled, dump already contains smoothed mesh
SimPar@Umbrella = 0
```

---

## Output Configuration

### SimPar@OutputFields

Comma-separated list of fields to export. Common options:
- `Pressure_V`: Vertex-based pressure
- `Pressure_E`: Element-based pressure
- `Velocity`: Velocity vector field
- `Viscosity`: Viscosity field (for non-Newtonian flows)
- `Temperature`: Temperature field (if heat transfer enabled)

**Example**:
```
SimPar@OutputFields = "Pressure_V,Velocity,Viscosity"
```

### SimPar@OutputLevel

Controls which mesh level is exported:
- `1`, `2`, `3`, ...: Specific level number
- `MAX`: Maximum (finest) level
- `MAX+1`: One level above maximum (for visualization of finer grid)
- `MAX-1`: One level below maximum

**Recommendation**: Use `MAX` for production runs, `MAX+1` for detailed visualization.

---

## Adaptive Time Stepping

When **SimPar@TimeAdaptivity = Yes**:

```
SimPar@TimeAdaptivity = Yes
SimPar@MinTimeAdapt = 1d-5    ! Smallest allowed timestep
SimPar@MaxTimeAdapt = 0.1d0   ! Largest allowed timestep
SimPar@TimeStep = 0.01d0      ! Initial timestep
```

The solver will automatically adjust the timestep based on:
- Convergence behavior
- Truncation error estimates
- Stability criteria

**When to use**:
- Transient problems with varying time scales
- Problems with sharp transients
- Uncertainty about optimal timestep

**When NOT to use**:
- Periodic flows (use fixed timestep matching period)
- Benchmarking (need reproducible timesteps)
- Already know optimal timestep

---

## Troubleshooting

### Problem: Solver Diverges

**Check**:
1. **Timestep too large**: Reduce `SimPar@TimeStep`
2. **Insufficient iterations**: Increase `Velo@MGMaxIterCyc` or `Velo@NLmax`
3. **Wrong time scheme**: Switch to `BE` (Backward Euler)
4. **Matrix renewal**: Use `M1D1K1S0C1` for safety

### Problem: Solver Too Slow

**Try**:
1. **F-cycle**: Set `Velo@MGCycType = F`
2. **Fewer smoothing steps**: Reduce `Velo@MGSmoothSteps` (but not below 2)
3. **Coarser working level**: Increase `Velo@MGMedLev`
4. **Matrix renewal**: Use `M1D1K3S0C1` or less frequent renewal

### Problem: Wrong Viscosity Behavior

**Most likely cause**: Using dynamic viscosity μ instead of kinematic viscosity ν!

**Check**: `Prop@Viscosity` should be ν = μ/ρ in units of [m²/s]

**For water at 1 mPa·s and 1000 kg/m³**:
- ❌ **WRONG**: `Prop@Viscosity = 0.001d0` (this is μ in Pa·s!)
- ✅ **CORRECT**: `Prop@Viscosity = 1d-6` (this is ν = μ/ρ in m²/s)

### Problem: Pressure Oscillations

**Solutions**:
1. **Use Backward Euler**: `SimPar@TimeScheme = BE`
2. **Enable NS stabilization**: `SimPar@NS_Stabilization = Yes`
3. **Refine mesh**: Increase `SimPar@MaxMeshLevel`
4. **More pressure iterations**: Increase `Pres@MGMaxIterCyc`

### Problem: Non-Convergence of Pressure Solver

**Try**:
1. **More iterations**: Increase `Pres@MGIterCoarse` (try 1000+)
2. **Tighter tolerance**: Reduce `Pres@MGCriterion1` and `Pres@MGCriterion2`
3. **More smoothing**: Increase `Pres@MGSmoothSteps` (try 20-30)
4. **MUMPS solver**: Set `Pres@MGCrsSolverType = 5` (if available)

---

## Fictitious Boundary Method (FBM) Control

### SimPar@skipFBMForce and SimPar@skipFBMDynamics

These boolean parameters provide fine-grained control over FBM (Fictitious Boundary Method) computations for particle-laden flows and fluid-structure interaction simulations.

**SimPar@skipFBMForce** (default: `No`):
- When set to `Yes`: Skips the computation of hydrodynamic forces acting on FBM particles
- **Use case**: Testing pure fluid flow without particle feedback
- **Warning**: Particles will not affect the fluid when this is enabled

**SimPar@skipFBMDynamics** (default: `No`):
- When set to `Yes`: Skips the rigid body dynamics update for FBM particles
- **Use case**: Fixed particles/bodies (stationary obstacles)
- **Note**: Forces are still computed (unless `skipFBMForce = Yes`), but particle positions/velocities are not updated

**Example configurations**:

```fortran
! Standard particle simulation (default)
SimPar@skipFBMForce = No
SimPar@skipFBMDynamics = No

! Stationary particles (forces computed but particles don't move)
SimPar@skipFBMForce = No
SimPar@skipFBMDynamics = Yes

! Pure fluid flow (particles present but non-interacting)
SimPar@skipFBMForce = Yes
SimPar@skipFBMDynamics = Yes
```

**Debugging workflow**:
1. First run: Both `No` (full coupling)
2. If unstable: Set `skipFBMDynamics = Yes` to check if particle motion is the issue
3. For one-way coupling tests: Set `skipFBMForce = Yes` (fluid affects particles, but not vice versa is not implemented - use both `Yes` for decoupling)

---

## Source Code Reference

All parameters are parsed from `q2p1_param.dat` in these subroutines:

| Category | File | Subroutine | Lines |
|----------|------|------------|-------|
| **SimPar@** | `source/src_util/param_parser.f90` | GDATNEW | 542-1076 |
| **Velo@** | `source/src_util/param_parser.f90` | GetVeloParameters | 104-249 |
| **Pres@** | `source/src_util/param_parser.f90` | GetPresParameters | 254-353 |
| **Prop@** | `source/src_util/param_parser.f90` | GetPhysiclaParameters | 359-500 |

**Note**:
- `GetPhysiclaParameters` has a spelling error in the source code (missing 'a' in "Physical").
- **IMPORTANT**: All parameter parsing functions were migrated from `source/Init.f90` and `source/src_quadLS/QuadSc_def.f90` to the centralized `source/src_util/param_parser.f90` module for better code organization and maintainability.

---

## Example Configuration Files

### Example 1: Laminar Channel Flow

```fortran
----------- Simulation Parameters --------------
SimPar@TimeScheme = BE
SimPar@TimeStep = 0.01d0
SimPar@MaxSimTime = 10.0d0
SimPar@OutputFreq = 0.1d0
SimPar@OutputFormat = "VTK"
SimPar@FlowType = Newtonian

----------------- Velocity --------------------
Velo@defCrit = 1d-6
Velo@NLmax = 5
Velo@MGCycType = F
Velo@MGSmoothSteps = 2
Velo@MGMaxIterCyc = 5

----------------- Pressure ---------------------
Pres@MGMaxIterCyc = 20
Pres@MGSmoothSteps = 16
Pres@MGCriterion1 = 1d-3

------------- Physical parameters --------------
Prop@Gravity = 0d0,0d0,0d0
Prop@Density = 1000d0,1d0
Prop@Viscosity = 1d-6,1d0  ! Water at 20°C (kinematic!)
```

### Example 2: Particle-Laden Flow

```fortran
----------- Simulation Parameters --------------
SimPar@TimeScheme = BE
SimPar@TimeStep = 0.001d0
SimPar@ParticleFile = "_data/particles.dat"
SimPar@BackUpFreq = 100
SimPar@OutputFields = "Pressure_V,Velocity"

----------------- Velocity --------------------
Velo@defCrit = 1d-7
Velo@NLmax = 10
Velo@MGCrsSolverType = 5  ! MUMPS for robustness

----------------- Pressure ---------------------
Pres@MGMaxIterCyc = 30
Pres@MGIterCoarse = 1000

------------- Physical parameters --------------
Prop@Gravity = 0d0,0d0,-9.81d0
Prop@Density = 1000d0,2500d0  ! Fluid and particle
Prop@Viscosity = 1d-6,1d0
```

---

## Notes and Conventions

1. **Fortran double precision**: Always use `d0` notation (e.g., `1d-6`, not `1e-6`)
2. **Case sensitivity**: Parameter names are case-insensitive, but values like `Yes/No` are case-sensitive
3. **Comments**: Use `!` for comments in parameter files
4. **Arrays**: No spaces in comma-separated arrays (e.g., `1d0,2d0,3d0`)
5. **String values**: Quotes are optional for strings without spaces

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.2 | 2026-01-21 | Added `Prop@DynVisc` parameter for intuitive specification of dynamic viscosity μ with automatic conversion to kinematic viscosity ν. |
| 1.1 | 2025-01-12 | Added comprehensive §Umbrella Mesh Smoothing section; clarified that dump files store P1 vertices only and Q2 nodes are regenerated via `SetUp_myQ2Coor()` |
| 1.0 | 2025-01-27 | Initial comprehensive documentation of all ~85 parameters |

---

## References

For more information about FeatFloWer and the Q2/P1 discretization, see:
- FeatFloWer documentation: `docs/`
- Example parameter files: `applications/*/\_data/q2p1_param.dat`
- CLAUDE.md: Repository root for build and usage instructions
