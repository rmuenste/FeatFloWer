# Euler-Lagrange Development Notes

This note summarizes the current repository state with respect to a future two-way coupled Euler-Lagrange (E-L) solver in FeatFlower.

Short answer: the current code base has strong pieces for resolved FBM/PE coupling, but it does not yet contain a full unresolved two-way Euler-Lagrange coupling layer.

## Terminology

- `G2P` means `Grid-to-Particle`.
  In practice, this is the operation that samples Eulerian fluid data from the mesh at particle positions. Typical examples are:
  - fluid velocity at a particle center
  - pressure or pressure gradient at a particle center
  - velocity gradients or strain rate at a particle center

- `P2G` means `Particle-to-Grid`.
  In practice, this is the reverse operation: a particle quantity is distributed back to the Eulerian mesh. For two-way E-L, this usually means:
  - spreading the reaction force from each particle to nearby fluid degrees of freedom or elements
  - accumulating particle volume into a porosity or void-fraction field

These are core coupling operators in an Euler-Lagrange method.

## Main Conclusion

FeatFlower currently has resolved FBM/PE coupling, not yet a true unresolved two-way Euler-Lagrange implementation.

The existing workflow is essentially:

1. solve the fluid system
2. compute resolved hydrodynamic forces from the fluid field
3. send those forces to PE
4. advance particle motion in PE

That is useful infrastructure, but it is not yet the same as unresolved two-way E-L with drag-law-based force closures and particle reaction-force feedback into the Navier-Stokes equations.

## Stage 0 Decisions

Stage 0 is the interface-cleanup and baseline-freeze stage. The implementation choices for this stage are:

- primary regression case: `q2p1_dns_drag`
- canonical CFD-to-PE force write path: struct-based particle setter
- legacy index-based force setter behavior: normalized to the same convention as the struct-based path
- force convention: FeatFlower writes hydrodynamic force and torque acting on the particle
- hidden interface scaling: not allowed
- future unresolved Euler-Lagrange subsystem location: `source/src_el/`

The intended coupling boundary for Stage 0 is:

`fluid solve -> force evaluation -> PE force write -> PE step -> particle state readback`

Ownership is:

- FeatFlower: fluid solution and hydrodynamic force evaluation
- `pe`: particle integration, contact handling, and authoritative particle state evolution

Stage 0 does not yet implement particle reaction-force feedback into the fluid equations.

## Missing Components For Two-Way E-L

### 1. Unresolved particle-to-fluid feedback

I did not find a `P2G` force-spreading kernel that takes particle hydrodynamic reaction forces and assembles them into the Navier-Stokes momentum right-hand side.

The existing generic forcing path still looks like a standard body-force interface. For example, the user RHS function still returns zero in:

- [source/src_pp3d/coeff.f](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-dns-drag/source/src_pp3d/coeff.f):193

and existing correction helpers cover gravity or constant forcing, e.g.:

- [source/src_quadLS/QuadSc_corrections.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-dns-drag/source/src_quadLS/QuadSc_corrections.f90):142

For a real two-way E-L solver, the following are still missing:

- particle force spreading to Q2 velocity DOFs or element/cell RHS vectors
- conservative sign convention: force on fluid equals minus hydrodynamic force on particle
- MPI-consistent ownership and reduction for particles near subdomain boundaries
- a stable hook into momentum assembly or defect formation

### 2. Unresolved Euler-Lagrange force closures

The present force evaluation is FBM/DNS-style: compute force from resolved stresses around the particle interface and pass that force to PE.

The transport path calls:

- [source/src_quadLS/QuadSc_main.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-dns-drag/source/src_quadLS/QuadSc_main.f90):568
- [source/src_quadLS/QuadSc_main.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-dns-drag/source/src_quadLS/QuadSc_main.f90):597

This is not yet the same as unresolved E-L force closure evaluation.

I did not find a general drag-law layer for point particles or unresolved finite-size particles, such as:

- Schiller-Naumann drag
- pressure-gradient force
- added mass
- lift
- hydrodynamic torque
- wall corrections
- hindered-settling or finite-volume-fraction corrections

These closures need to exist as a dedicated E-L module, separate from the resolved FBM force integration path.

### 3. A clean Grid-to-Particle sampling API

There is useful partial infrastructure for `G2P`.

Velocity sampling and point location already exist in reusable form, for example:

- [source/src_fbm/fbm_loc.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-dns-drag/source/src_fbm/fbm_loc.f90):61
- [source/src_fbm/fbm_particle_reynolds.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-dns-drag/source/src_fbm/fbm_particle_reynolds.f90):27

There is also particle tracing logic using element lookup and velocity evaluation:

- [source/src_fbm/fbm_ptrace.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-dns-drag/source/src_fbm/fbm_ptrace.f90):16

However, for E-L this still needs to be turned into a clean public coupling interface that can reliably provide:

- velocity at particle positions
- pressure or pressure gradient at particle positions
- velocity gradient or strain-rate tensor at particle positions
- element id and reference coordinates for reuse across substeps
- robust handling of remote particles and partition-boundary points

### 4. Void fraction or porosity accumulation

The code already has resolved solid-marking via `FictKNPR`, which is part of the FBM machinery. That is not the same as unresolved E-L porosity coupling.

Current handler setup is still FBM-oriented:

- [source/src_quadLS/QuadSc_handlers.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-dns-drag/source/src_quadLS/QuadSc_handlers.f90):69

For unresolved two-way E-L, the following are still missing:

- particle volume spreading to cells or DOFs
- a fluid void-fraction field such as `epsilon_f`
- optional solids fraction field such as `alpha_s`
- integration of porosity into momentum or constitutive closures

### 5. Two-way coupling time integration strategy

The current order is effectively:

1. solve fluid
2. compute hydrodynamic forces
3. advance particles

For example:

- [applications/q2p1_dns_drag/q2p1_dns_drag.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-dns-drag/applications/q2p1_dns_drag/q2p1_dns_drag.f90):41

and then in the transport code:

- [source/src_quadLS/QuadSc_main.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-dns-drag/source/src_quadLS/QuadSc_main.f90):570
- [source/src_quadLS/QuadSc_main.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-dns-drag/source/src_quadLS/QuadSc_main.f90):599

For two-way E-L, the reaction force must feed back into the fluid equations. Even with an explicit first version, the infrastructure must define:

- when particle reaction forces are accumulated
- whether the feedback enters the current or next fluid step
- whether drag is treated explicitly or semi-implicitly
- whether particles subcycle relative to the fluid time step

### 6. PE force interface cleanup

The PE bridge is already useful and exposes a good amount of particle data:

- [source/src_particles/dem_query.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-dns-drag/source/src_particles/dem_query.f90):11

The struct-based setter reaches PE here:

- [libs/pe/src/interface/object_queries.cpp](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-dns-drag/libs/pe/src/interface/object_queries.cpp):1313

and writes force and torque here:

- [libs/pe/src/interface/object_queries.cpp](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-dns-drag/libs/pe/src/interface/object_queries.cpp):1402

This is a good base for E-L. However, the older index-based setter still contains special behavior:

- [libs/pe/src/interface/object_queries.cpp](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-dns-drag/libs/pe/src/interface/object_queries.cpp):1126

That path applies a hardcoded `1.025` multiplier and zeroes torque instead of using the incoming torque. For a clean E-L solver, all force-setting paths should use a single, well-defined convention for:

- units
- scaling
- sign
- torque handling

### 7. Application-level separation of FBM and E-L

The main applications such as `q2p1_ATC` and `q2p1_dns_drag` still run through the FBM/FC extension path:

- [applications/q2p1_ATC/q2p1_ATC.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-dns-drag/applications/q2p1_ATC/q2p1_ATC.f90):67
- [applications/q2p1_dns_drag/q2p1_dns_drag.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-dns-drag/applications/q2p1_dns_drag/q2p1_dns_drag.f90):42

There is also ongoing work in:

- [applications/q2p1_bench_fluidization/app_init.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-dns-drag/applications/q2p1_bench_fluidization/app_init.f90)

but this still looks like PE/FBM-oriented infrastructure and benchmark staging, not yet a separate unresolved E-L source assembly path.

## What Already Exists And Can Be Reused

The following pieces are already useful for an E-L implementation:

- PE particle state and mechanics backend
- PE stepping and particle force/torque interfaces
- particle state exchange via the Fortran/C++ bridge
- point location and field interpolation infrastructure
- resolved hydrodynamic force calculation for FBM
- particle Reynolds-number diagnostics

In particular:

- particle bridge definitions: [source/src_particles/dem_query.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-dns-drag/source/src_particles/dem_query.f90)
- point location: [source/src_fbm/fbm_loc.f90](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-dns-drag/source/src_fbm/fbm_loc.f90)
- PE simulation stepping: [libs/pe/src/interface/sim_setup.cpp](/data/warehouse17/rmuenste/code/FF-ATC-NEW/feature-dns-drag/libs/pe/src/interface/sim_setup.cpp):151

## Bottom Line

The central missing layer for two-way E-L is:

`Grid-to-Particle sampling -> E-L force closure -> set force in PE -> Particle-to-Grid reaction spreading -> Navier-Stokes RHS/source assembly`

At the moment, the repository has parts of:

- `Grid-to-Particle sampling`
- `set force in PE`
- an alternative resolved FBM force path

but it does not yet have:

- an unresolved E-L force-closure layer
- a `Particle-to-Grid` reaction-force spreading operator
- a momentum-source assembly path for particle feedback into the fluid solve
- unresolved void-fraction or porosity coupling

That is the main gap between the current code base and a true two-way coupled Euler-Lagrange solver.

## Stage 0 Status

The following Stage 0 interface cleanup has been applied:

- the canonical CFD-to-PE write path remains the struct-based particle setter used by the force-extension code
- the legacy index-based PE setter has been normalized to the same convention:
  - no hidden force scaling
  - incoming torque is preserved instead of being discarded
- the Stage 0 baseline target is `q2p1_dns_drag`

This keeps the interface behavior internally consistent before Stage 1 begins.

## Proposed Implementation Roadmap

The roadmap below is intentionally conservative. Each stage ends with a test gate. The next stage should only begin after the current stage passes.

The main principle is:

1. build the coupling infrastructure first
2. validate it in one-way mode
3. add two-way feedback only after interpolation and particle forcing are trustworthy
4. add denser-flow corrections only after the dilute solver is stable

## Stage 0: Interface Cleanup And Baseline Freeze

### Goal

Create a stable baseline before adding new physics.

### Work Items

- identify and document the current FBM/PE and transport call chain
- normalize PE force-setting paths so unit handling and sign conventions are explicit
- remove or isolate special-case scaling in legacy PE force setter paths
- define a single E-L particle data access interface on top of `dem_query`
- decide where the future E-L module will live, for example:
  - `source/src_el/`
  - `applications/q2p1_el_*`

### Deliverables

- documented force/sign convention
- documented particle ownership convention
- documented insertion points for `G2P` and `P2G`
- stable baseline application for future regression testing

### Tests

- build test for current PE-enabled configurations
- baseline regression run for one existing resolved case:
  - `q2p1_dns_drag` or `q2p1_bench_sedimentation`
- verify that particle positions, velocities, and forces remain unchanged relative to the pre-cleanup baseline within tolerance

### Exit Criteria

- no behavior change in existing benchmark results beyond agreed numerical tolerance
- no ambiguous force scaling path remains in the PE interface

## Stage 1: Grid-to-Particle Sampling API

### Goal

Build a reusable `G2P` coupling layer for unresolved particles.

### Work Items

- create a public sampling API for:
  - velocity at particle positions
  - pressure at particle positions
  - pressure gradient at particle positions
  - velocity gradient or strain-rate tensor at particle positions
- expose element id plus reference coordinates so repeated evaluations can reuse local information
- support remote particles and points near partition boundaries
- make the API independent of FBM-specific logic where possible

### Suggested Implementation Shape

- `el_sample_velocity(particles, U, V, W, ...)`
- `el_sample_pressure(particles, P, ...)`
- `el_sample_gradients(particles, U, V, W, ...)`
- particle sample record containing:
  - particle id
  - owner rank
  - element id
  - reference coordinates
  - sampled values
  - status flag

### Tests

- manufactured interpolation test:
  - prescribe an analytic velocity field
  - sample at particle positions
  - compare against exact value
- manufactured gradient test:
  - prescribe an analytic field with known gradient
  - sample gradient at particle positions
  - compare against exact gradient
- MPI boundary test:
  - place sample points near partition boundaries
  - ensure local and distributed search return the same result within tolerance
- regression test for existing particle Reynolds diagnostic path using the new sampling API

### Exit Criteria

- velocity interpolation error is within agreed tolerance on a static mesh
- gradient reconstruction error is within agreed tolerance
- MPI point-location and sampling are stable and deterministic

## Stage 2: One-Way Euler-Lagrange Force Closure Layer

### Goal

Build the first unresolved particle force module without fluid feedback.

### Work Items

- implement a modular force-closure framework
- first closures should be:
  - drag
  - gravity
  - buoyancy
- use sampled fluid values from Stage 1
- write forces and torques into PE through the cleaned particle interface
- keep this stage strictly one-way:
  - fluid affects particles
  - particles do not affect fluid yet

### Recommended First Closure Set

- drag:
  - start with Schiller-Naumann for spheres
- gravity and buoyancy:
  - apply consistently with particle density and fluid density

### Optional But Useful Extensions Within This Stage

- pressure-gradient force
- rotational drag torque for spheres

### Tests

- single particle settling in quiescent fluid:
  - compare transient and terminal velocity against analytic or semi-analytic reference
- particle advection in prescribed uniform flow:
  - verify particle asymptotically approaches carrier velocity
- particle response-time test:
  - compare velocity relaxation against exponential reference
- one-way dilute benchmark:
  - compare against old tracer/LPT-style behavior where applicable

### Exit Criteria

- terminal settling velocity matches reference within agreed tolerance
- particle relaxation test matches expected response-time behavior
- particle motion is stable for the target time-step range

## Stage 3: Particle-To-Grid Reaction Force Spreading

### Goal

Build the `P2G` operator for conservative particle reaction-force feedback.

### Work Items

- design a force-spreading kernel from particle force to fluid DOFs or element RHS
- ensure force conservation:
  - sum of distributed force equals total reaction force
- decide spatial support:
  - containing element only
  - neighboring support region
  - shape-function-weighted projection
- add accumulation buffers for particle source terms
- define ownership and reduction rules for distributed particles

### Design Requirements

- conservation of total force
- stable distribution on distorted meshes
- deterministic MPI reductions
- explicit sign convention:
  - if drag on particle is `F_p`, source to fluid is `-F_p`

### Tests

- algebraic conservation test:
  - spread a known particle force
  - verify sum over fluid RHS equals expected reaction force
- symmetry test:
  - symmetric particle arrangement should produce symmetric source distribution
- MPI consistency test:
  - same setup on 1 rank and multiple ranks gives matching integrated source terms
- null test:
  - zero particle force produces zero fluid source everywhere

### Exit Criteria

- total spread force is conservative within machine or assembly tolerance
- distributed and serial runs agree on integrated source terms
- no spurious source appears in the null test

## Stage 4: Two-Way Coupled Momentum Source Assembly

### Goal

Insert the `P2G` source term into the Navier-Stokes momentum equations.

### Work Items

- identify the correct momentum RHS insertion point
- add accumulated particle feedback to the momentum defect or forcing vectors
- decide first coupling strategy:
  - explicit lagged coupling is the safest first target
- add source reset and accumulation lifecycle management per time step
- verify interaction with boundary conditions and existing body-force terms

### Scope Recommendation

Start with:

- explicit two-way coupling
- no porosity correction yet
- dilute suspension regime only

Do not combine this stage yet with:

- hindered settling
- lubrication corrections
- porosity-modified momentum equations

### Tests

- momentum conservation test:
  - particle receives `F`
  - fluid receives `-F`
  - integrated coupled momentum budget is correct
- single settling particle with two-way feedback:
  - compare against expected reduction in local fluid velocity
- low-volume-fraction suspension test:
  - verify stable two-way feedback without solver divergence
- time-step sensitivity test:
  - confirm expected explicit-coupling time-step limits

### Exit Criteria

- fluid RHS contains the expected particle reaction force
- coupled runs remain stable in dilute regimes
- integrated momentum budget is consistent

## Stage 5: Void Fraction And Porosity Coupling

### Goal

Add unresolved particle volume displacement effects.

### Work Items

- introduce particle volume accumulation onto mesh support
- build fields such as:
  - fluid void fraction `epsilon_f`
  - solids fraction `alpha_s`
- define how porosity enters:
  - drag correction
  - fluid coefficients
  - optional porosity-modified equations
- ensure consistency between force closure and occupancy model

### Tests

- volume conservation test:
  - integrated solids fraction matches total particle volume
- uniform cloud test:
  - homogeneous particle distribution yields nearly homogeneous porosity field
- mesh refinement test:
  - porosity field converges under refinement
- low-concentration regression:
  - porosity-enabled solver reduces to Stage 4 behavior as particle volume fraction approaches zero

### Exit Criteria

- integrated particle volume is conserved in the porosity field
- porosity behaves smoothly under mesh refinement
- low-volume-fraction limit matches the dilute two-way solver

## Stage 6: Dense-Regime Hydrodynamic Corrections

### Goal

Add the first set of dense-suspension corrections required for more realistic ATC behavior.

### Work Items

- add hindered-settling or finite-volume-fraction drag correction
- add near-wall corrections
- add lubrication correction for close particle-particle and particle-wall approach
- optionally add lift and hydrodynamic torque closures if required by validation targets

### Scope Recommendation

Do this in the following order:

1. hindered-settling or volume-fraction drag correction
2. wall correction
3. lubrication correction
4. lift and torque

### Tests

- wall-settling test:
  - compare settling near a wall against reference trends
- pair interaction test:
  - confirm close-approach correction behaves qualitatively correctly
- dilute-to-dense transition test:
  - verify drag increases with local solids fraction as expected
- stability test in clustered initial conditions

### Exit Criteria

- added closures improve agreement in the intended regime
- no regression in dilute cases from earlier stages
- clustered and near-wall cases remain numerically stable

## Stage 7: Parallel Ownership, Migration, And Production Hardening

### Goal

Make the solver operational for larger distributed runs.

### Work Items

- finalize particle ownership rules
- implement particle migration across CFD subdomains if needed by the E-L layer
- verify compatibility with PE ghost and shadow particles
- optimize `G2P` and `P2G` communication patterns
- add diagnostics for:
  - missed samples
  - lost particles
  - force conservation
  - porosity conservation

### Tests

- rank-scaling consistency test:
  - same case on different rank counts gives matching integrated results
- particle migration test:
  - particles crossing many subdomains remain valid and coupled
- long-run stability test:
  - no drift in particle count, total volume, or total spread force
- restart test:
  - restart preserves particle state and coupling fields

### Exit Criteria

- serial and MPI runs agree within tolerance
- migration and restart work reliably
- conservation diagnostics remain bounded over long runs

## Stage 8: Validation Against DNS And Experiments

### Goal

Decide whether the E-L model is credible for the intended ATC regimes.

### Suggested Validation Ladder

1. single particle settling
2. particle in uniform flow
3. particle near wall
4. two-particle interaction
5. dilute suspension
6. two-way coupled suspension
7. ATC-relevant benchmark against DNS or experiment

### Tests

- compare bulk slip velocity
- compare pressure drop
- compare particle concentration profiles
- compare near-wall accumulation trends
- compare transient startup behavior and clustering trends

### Exit Criteria

- documented agreement envelope for each validation case
- clear statement of where the model is valid and where it is not

## Recommended Stage Boundaries

The safest stage boundaries are:

1. `Stage 0 -> Stage 1`
   only after baseline cleanup is behavior-neutral
2. `Stage 1 -> Stage 2`
   only after interpolation and gradients are verified
3. `Stage 2 -> Stage 3`
   only after one-way particle dynamics match reference behavior
4. `Stage 3 -> Stage 4`
   only after force spreading is conservative
5. `Stage 4 -> Stage 5`
   only after dilute two-way momentum coupling is stable
6. `Stage 5 -> Stage 6`
   only after porosity transport and accumulation are reliable
7. `Stage 6 -> Stage 7`
   only after dense closures pass small-scale validation
8. `Stage 7 -> Stage 8`
   only after MPI consistency and long-run robustness are established

## Suggested First Practical Milestone

The most defensible first milestone is:

- Stage 0 completed
- Stage 1 completed
- Stage 2 completed

That yields a one-way unresolved Euler-Lagrange solver on top of PE, with a reusable `G2P` sampling layer and validated drag-driven particle motion.

The next milestone should then be:

- Stage 3 completed
- Stage 4 completed

That yields the first true two-way coupled dilute E-L solver.
