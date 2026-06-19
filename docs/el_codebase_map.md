# Euler-Lagrange Codebase Map

This map records the integration points for the unresolved Euler-Lagrange
runtime. The implementation uses a dedicated `q2p1_el_pipeflow` executable
and `Transport_q2p1_UxyzP_el`; no global simulation-mode switch is required.

## Fluid time integration

- `source/src_quadLS/QuadSc_main.f90`
  - `Transport_Q2P1_UxyzP` and `Transport_q2p1_UxyzP_fc_ext` drive the
    fractional-step-theta momentum and pressure solves.
  - `Matdef_General_QuadScalar(QuadSc, 1)` assembles the explicit momentum
    right-hand side.
  - `Matdef_General_QuadScalar(QuadSc, -1)` assembles the implicit defect and
    fine-grid operator.
  - `AddPressureGradient`, `AddGravForce`, and `AddConstantForce` add the
    existing volume-force terms.
- `source/src_quadLS/QuadSc_boundary.f90`
  - `Boundary_QuadScalar_Def`, `Boundary_QuadScalar_Val`,
    `Boundary_QuadScalar_Mat`, and `Boundary_QuadScalar_Mat_9` impose physical
    and fictitious-boundary constraints.
  - The E-L transport must keep `FictKNPR` clear and must not call
    `updateFBMGeometry`, so these routines apply only the physical boundary
    data in the dedicated executable.

## Resolved-particle force path

- `source/src_quadLS/QuadSc_main.f90`
  - `fbm_updateForces` evaluates the resolved FBM hydrodynamic force after the
    fluid correction.
  - `fbm_updateFBM` advances the resolved rigid bodies.
- `source/src_quadLS/QuadSc_force_torque_calc.f90` contains the stress/volume
  force integration used by the FBM path.

The E-L transport bypasses both calls. It supplies closure forces through
`dem_query` and advances PE through the E-L step entry point.

## CFD to PE exchange

- `source/src_particles/dem_query.f90`
  - `getAllParticles` returns locally owned particles.
  - `getAllRemoteParticles` returns PE shadow particles.
  - `setForcesMapped` and `setRemoteForcesMapped` write complete particle
    records back to PE.
- `libs/pe/src/interface/c2f_interface.cpp`
  - `commf2c_el_frozen_trace_` initializes the existing E-L-capable PE setup.
  - `step_el_frozen_trace_` advances PE and preserves its contact,
    lubrication, and synchronization behavior.

The first transfer implementation broadcasts owned particles over the CFD
communicator and normalizes their kernel support collectively. This is
conservative across partition boundaries without depending on the lifetime
or overlap width of PE shadow particles.

## Point location and field evaluation

- `applications/q2p1_el_frozen_trace/el_frozen_driver.f90` contains a
  prototype containing-element search and direct Q2 interpolation.
- `source/src_fbm/fbm_aux.f90` provides `fbmaux_PointInHex`.
- Q2 velocity unknowns are ordered as vertices, edges, faces, and element
  centers. The first transfer implementation samples the element-center Q2
  value and applies positive compact kernel weights. This avoids using
  sign-changing Q2 basis values as accumulation weights.

An optimized support search and quadrature-based Q2 projection remain later
Phase 1 work; the collective element-center implementation is the reference
conservative path.

## Time-step and PE stepping

- Application loops, for example
  `applications/q2p1_fc_ext/q2p1_fc_ext.f90`, update `timens`, call one
  transport routine, then postprocess.
- PE subcycling is implemented inside the PE E-L setup selected by
  `commf2c_el_frozen_trace_`; `step_el_frozen_trace_` is called once per fluid
  step from the dedicated E-L transport.

## Pipe boundary configuration

- `SimPar@UseConstantForcing` and `SimPar@ConstantForcing` configure the
  existing body-force driver.
- `applications/q2p1_dns_drag/app_init.f90` enables periodicity in all three
  directions and provides a suitable initializer to share with the new
  executable.

The initial smoke case reuses the periodic unit-cube data. A circular,
axially periodic pipe case and PE wall description are separate validation
assets.

## Plan adaptations

- There is no shared runtime dispatcher. A separate executable and transport
  routine are less invasive than adding `simulation_mode` throughout common
  initialization.
- The current boundary callbacks combine physical and FBM constraints.
  Keeping `FictKNPR` empty in the dedicated E-L runtime reuses the physical
  boundary implementation while bypassing particle constraints.
- The existing frozen-field prototype is application-local. Reusable
  transfer and force code therefore belongs in `source/src_el`.
- The first milestone uses P0 element fields for `alpha_p`, `epsilon_f`,
  `epsilon_f_old`, and `d epsilon_f/dt`. Q2 force feedback is diagnostic
  until two-way coupling is enabled.
- The transfer cubature (`EL_Q1_MAP` in `source/src_el/el_quadrature.f90`)
  evaluates the geometric Jacobian with the trilinear (Q1) hex map even
  though the velocity fields are triquadratic (Q2/E013). This is exact for
  the straight-edged hexes the structured mesh hierarchy produces, so the
  cubature weights and kernel-volume conservation are correct. Curved/Q2
  geometry would require upgrading this to the full Q2 isoparametric
  Jacobian.

## Phase 1 transfer tests

- `source/src_el/tests/test_el_kernel_forces.f90` — serial unit test for
  kernel positivity, Stokes drag, and field restart round-trip
  (`el-kernel-forces-restart`).
- `source/src_el/tests/test_el_transfer.f90` — particle-mesh transfer
  conservation/interpolation suite on a synthetic structured Q2/P1 mesh,
  registered as `el-transfer-serial` and `el-transfer-mpi-{2,8}`. Covers
  constant/linear interpolation, central + wall-adjacent volume
  conservation, componentwise force-spread conservation, zero-force
  behaviour, halo overlap predicates, and rank-independence of the transfer
  math. TODO (documented in the file header): the genuinely distributed
  straddling-particle path (`EL_REFRESH_COUPLING_HALO` ->
  `EL_REDUCE_TO_OWNERS` -> `E013Sum3`) needs the partitioned production mesh
  and PE ownership, so it must be driven from a `q2p1_el_pipeflow`
  regression harness on the staged `unit_cube_27_case` (face/edge/corner
  straddling on 2 and 8 ranks).
