# Euler-Lagrange Frozen-Field Project

## Purpose

This folder tracks the frozen-field Euler-Lagrange development path.

The goal is not yet a full two-way coupled Euler-Lagrange solver. The goal is to build a near-real application path that:

- loads a CFD solution from dump/checkpoint output
- keeps the carrier field fixed
- initializes the usual mesh, FE, MPI, and PE structures
- runs particle-side kernels against that frozen carrier field

This gives a controlled environment for validating important Euler-Lagrange kernels before introducing live fluid stepping and two-way feedback.

## Current Idea

The frozen-field application should behave like a real Euler-Lagrange driver, except that the carrier field is not advanced.

Conceptually:

1. initialize the application normally
2. load `U,V,W,P` from a dump/checkpoint
3. do not call the Navier-Stokes transport solve in the main loop
4. use the loaded field as the static carrier field for particle-side operations

This makes it possible to validate:

- point location
- `G2P` interpolation
- sampled pressure / gradients
- particle time stepping
- PE interaction
- later force closures such as drag and buoyancy

without mixing in:

- fluid time integration
- nonlinear flow solver behavior
- `P2G` feedback
- porosity coupling

## Current State

The first scaffold is working.

Implemented so far:

- a dedicated application: `applications/q2p1_el_frozen_trace`
- dedicated PE bridge entry:
  - `commf2c_el_frozen_trace`
- dedicated serial PE setup:
  - `setupELFrozenTraceSerial`
- mandatory checkpoint-load guard:
  - the app aborts if `SimPar@StartingProc == 0`
  - the app aborts if `SimPar@StartFile` is empty
- main program scaffold with no live transport solve
- successful configure/build path with:
  - `USE_PE=ON`
  - `USE_PE_SERIAL_MODE=ON`
  - `USE_JSON=ON`
  - `USE_CGAL=OFF`
  - `EIGEN=ON`
  - `ENABLE_FBM_ACCELERATION=ON`

Important implementation correction already identified:

- `q2p1_dns_drag` was not a safe one-to-one initialization template because it assumes a periodic configuration
- the hardcoded `dPeriodicity(:)=1.0d0` setup was removed
- partition setup was aligned with `q2p1_ATC`

At the moment, the application is an initialization harness that reaches the point just before a live application would enter its transport loop.

## Architectural Decision

This application should not become a one-off test driver.

The long-term structure should be:

- application driver:
  - `applications/q2p1_el_frozen_trace`
- reusable Euler-Lagrange kernel layer:
  - future location: `source/src_el/`

The frozen-field app should call the same public sampling and particle-side kernels that the later live Euler-Lagrange solver will call.

## What Is Not Implemented Yet

Not implemented yet:

- explicit frozen-field time loop with particle-side work
- post-load validation of the dumped carrier field
- particle position queries and diagnostics in the new app loop
- `G2P` sampling API for the frozen-field app
- tracer advection mode
- one-way force closure mode
- `P2G` feedback
- porosity / void-fraction coupling

## Immediate Next Steps

### Step 1: Validate the Loaded Carrier Field

After the checkpoint load path completes, add a small validation summary:

- simulation time `timens`
- norms or min/max of `U,V,W,P`
- one or more sampled probe values

Purpose:

- confirm that the dump was really loaded
- make frozen-field regressions easier to spot

### Step 2: Introduce a Frozen-Field Main Loop

Add a loop structure that advances application time bookkeeping but does not call:

- `Transport_q2p1_UxyzP_fc_ext`

This loop should become the execution shell for all later frozen-field particle work.

### Step 3: First Read-Only Particle Kernel Pass

Add a first particle-side diagnostic pass with no particle motion yet:

- query particle positions from PE
- perform element lookup
- sample carrier velocity at particle positions
- write diagnostic output for a small set of particles

Purpose:

- validate the first real `G2P` usage in the new app

### Step 4: Tracer Advection Mode

Before adding drag laws, implement the simplest dynamic mode:

- particles move with sampled carrier velocity
- no force closure yet

This isolates:

- point location
- interpolation
- particle stepping
- MPI ownership issues

### Step 5: One-Way Force Closures on Frozen Field

After tracer mode works, add:

- slip velocity evaluation
- drag
- gravity
- buoyancy

Forces should be written to PE through the cleaned interface, while the carrier field remains frozen.

## Recommended Development Order

The recommended whole-project order is:

1. make frozen-field initialization and field loading observable
2. add frozen-field loop without fluid solve
3. add read-only particle sampling diagnostics
4. add tracer advection
5. add one-way closures on the frozen field
6. extract or stabilize reusable `src_el` kernel interfaces
7. transition to a live one-way Euler-Lagrange application
8. add `P2G` reaction-force spreading
9. add true two-way coupled momentum feedback

## Practical Short-Term Milestone

The next practical milestone should be:

- frozen-field app runs through initialization
- frozen-field app loads a dump solution
- frozen-field app enters a loop without solving Navier-Stokes
- frozen-field app samples carrier velocity at particle positions
- frozen-field app writes deterministic diagnostics

Once that works, the project has moved from "scaffold" to "usable kernel validation harness".
