# q2p1_el_frozen_trace Parameter Reference

This note documents the user-facing runtime parameters that currently matter for the frozen-field Euler-Lagrange application `q2p1_el_frozen_trace`.

The application loads a CFD solution from dump/checkpoint output, keeps the carrier field fixed, and runs particle-side Euler-Lagrange kernels against that frozen field.

## Required Startup Parameters

These parameters must be set correctly or the application aborts during initialization.

- `SimPar@StartingProc`
  - must be nonzero
  - current supported frozen-field use is checkpoint/dump based startup
  - typical value: `1`

- `SimPar@StartFile`
  - path prefix of the dump/checkpoint to load
  - must be nonempty
  - typical value: `"_dump/01"`

## Frozen-Field Euler-Lagrange Parameters

The following `SimPar@...` entries control the current frozen-field particle pass.

### `SimPar@ELForceKernel`

Selects the hydrodynamic closure kernel.

Supported values:

- `none`
  - read-only sampling mode
  - sampled carrier quantities are reconstructed across MPI
  - no hydrodynamic force is generated

- `tracer`
  - placeholder kernel for tracer-style operation
  - currently behaves like a no-force kernel in the frozen-field pass
  - intended as the bridge toward pure advection mode

- `stokes_drag`
  - translational Stokes drag
  - uses carrier velocity minus particle velocity as slip
  - force is written back to PE

- `schiller_naumann`
  - drag law with Schiller-Naumann Reynolds correction
  - extends `stokes_drag`
  - force is written back to PE

Recommended current value:

- `none` for validating particle visibility, ownership, and `G2P` sampling

### `SimPar@ELWriteDiagnostics`

Controls whether the frozen-field particle pass writes a CSV diagnostics file.

Supported values:

- `Yes`
- `No`

Current output:

- `el_frozen_particles_stepXXXXXX.csv`

The file is written by the representative rank (`myid == showid`) and contains, per particle:

- particle id
- position
- radius
- particle velocity
- sampled carrier velocity
- slip velocity
- hydrodynamic force
- particle Reynolds number
- `found_count`
- ownership/debug sums

Recommended current value:

- `Yes`

### `SimPar@ELApplyForces`

Controls whether the reconstructed hydrodynamic force/torque is written back into PE.

Supported values:

- `Yes`
- `No`

Notes:

- in `PE_SERIAL_MODE`, every CFD rank holds the full PE particle replica
- after MPI reconstruction, every rank has the same particle-wise force values
- with `ELApplyForces = Yes`, each worker rank writes the same force data into its local PE replica

Recommended current value:

- `Yes` for closure-path validation
- `No` if only sampling diagnostics are being tested

### `SimPar@ELEnableBuoyancy`

Controls whether the FeatFloWer-side closure layer adds the net gravity/buoyancy body force to
force-producing kernels.

Supported values:

- `Yes`
- `No`

Implemented force contribution:

- `F_body = (rho_p - rho_f) * V_p * ELGravity`
- `rho_p` comes from `SimPar@ELParticleDensity`
- `rho_f` comes from `SimPar@ELFluidDensity`
- `V_p` is the spherical particle volume computed from the PE particle radius

Notes:

- this contribution is added on the FeatFloWer side, not inside PE
- PE should therefore keep its own gravity/body-force contribution neutral for this E-L mode
- the current implementation adds this term to drag-producing kernels such as `stokes_drag` and `schiller_naumann`

Recommended current value:

- `No`

### `SimPar@ELGravity`

Gravity/acceleration vector used by the FeatFloWer-side E-L closure layer when
`SimPar@ELEnableBuoyancy = Yes`.

Format:

- three real values: `gx gy gz`

Examples:

- `0.0d0 0.0d0 0.0d0` keeps buoyancy/gravity force inactive even if the switch is enabled
- `0.0d0 0.0d0 -981.0d0` is the usual CGS gravity vector for a `z`-vertical setup

Default sample value:

- `0.0d0 0.0d0 0.0d0`

### `SimPar@ELFluidDensity`

Fluid density used by the frozen-field closure models.

Typical CGS value for the current benchmark family:

- `1.0d0`

Current use:

- used by drag-based kernels when computing Reynolds number and force scaling

### `SimPar@ELKinematicViscosity`

Fluid kinematic viscosity used by the frozen-field closure models.

Typical CGS value for the current benchmark family:

- `1.0d-3`

Current use:

- dynamic viscosity is formed internally as
  - `mu = rho_f * nu`

### `SimPar@ELParticleDensity`

Reference particle density for future closure extensions.

Current status:

- parsed and reported
- not yet used by the present closure implementation

Typical current value:

- `1.0d0`

## Current MPI / PE_SERIAL_MODE Semantics

The current frozen-field implementation is designed for `PE_SERIAL_MODE`.

That means:

- every CFD worker rank has a full PE copy of all particles
- CFD sampling is local to each mesh partition
- sampled particle quantities are reconstructed with `MPI_Allreduce(..., SUM, ...)`
- after reconstruction, every worker rank holds the same hydrodynamic particle data
- the same force contribution is then written into every PE replica

The representative diagnostic rank is:

- `myid == showid`

But all worker ranks must still participate in:

- local particle sampling
- MPI reduction of sampled quantities
- PE force assignment

## Diagnostics Interpretation

The most important validation field in the CSV output is:

- `found_count`

Meaning:

- `0`
  - particle was not found in any local CFD partition
  - indicates a point-location or domain-coverage problem

- `1`
  - expected case for a well-owned particle in the current frozen-field path

- `>1`
  - more than one rank contributed a local sample
  - indicates duplicate ownership or overlap in the particle sampling path

For `ELForceKernel = none`, expect:

- `fx = fy = fz = 0`

For drag kernels, expect:

- force direction opposite to particle slip
- nonzero Reynolds number when slip is nonzero

## Example Block

Example `SimPar` block for the current validation stage:

```text
SimPar@StartingProc = 1
SimPar@StartFile = "_dump/01"

SimPar@ELForceKernel = none
SimPar@ELWriteDiagnostics = Yes
SimPar@ELApplyForces = Yes
SimPar@ELEnableBuoyancy = No
SimPar@ELFluidDensity = 1.0d0
SimPar@ELKinematicViscosity = 1.0d-3
SimPar@ELParticleDensity = 1.0d0
SimPar@ELGravity = 0.0d0 0.0d0 0.0d0
```

## Planned Extensions

The following parameters are likely to appear later, but are not part of the current implementation yet:

- pressure-gradient usage toggle
- velocity-gradient / strain-rate sampling toggle
- torque-capable closure selection
- closure-specific relaxation or damping parameters
- tracer-advection switch
- particle output frequency controls dedicated to frozen-field diagnostics

For the broader development roadmap, see:

- [euler_lagrange_development.md](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/docs/md_docs/euler_lagrange_development.md)
- [README.md](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-euler-lagrange/el-frozen/README.md)
