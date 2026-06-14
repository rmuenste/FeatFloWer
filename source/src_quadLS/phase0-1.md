• Proposed Plan


  # Euler–Lagrange Phase 0-1 Implementation

  ## Summary

  Create a separate live E-L application, q2p1_el_pipeflow, driven by Transport_q2p1_UxyzP_el. Existing FBM transport routines remain unchanged. Phase 1
  provides conservative particle-mesh transfer and provisional one-way Stokes/Schiller–Naumann particle dynamics, but no particle feedback in the fluid
  equations yet.

  ## Key Changes

  - Add docs/el_codebase_map.md documenting momentum assembly, FBM constraints and forces, PE exchange, point location, time stepping, and periodic
    forcing.

  - Add common configuration with simulation_mode = DNS_FBM | EULER_LAGRANGE; existing applications default to DNS_FBM, while the new application
    requires EULER_LAGRANGE.

  - Build reusable modules under source/src_el:
      - Configuration and field ownership.
      - Positive compact kernels: deen_poly and truncated gaussian.
      - Particle/element support search and wall renormalization.
      - P0 volume-fraction accumulation.
      - Kernel-weighted fluid sampling.
      - Conservative Q2 force projection.
      - Provisional Stokes and Schiller–Naumann closures migrated from the frozen-field prototype.

  - Store authoritative fields per element:
      - alpha_p
      - epsilon_f
      - previous epsilon_f
      - d epsilon_f/dt

  - Store particle feedback as three Q2 RHS vectors. Phase 1 computes these for diagnostics/tests but does not inject them into momentum assembly.
  - Add checkpoint/restart and VTK cell output for the P0 fields.
  - Refactor q2p1_el_frozen_trace to consume the shared sampling and provisional closure APIs rather than retaining an application-local implementation.

  ## Runtime Flow

  Transport_q2p1_UxyzP_el will:

  1. Retrieve owned and shadow PE particles.
  2. Determine kernel-support elements and normalized positive weights.
  3. Assemble alpha_p, clipped epsilon_f, and d epsilon_f/dt.
  4. Kernel-sample the current carrier velocity at particles.
  5. Evaluate provisional one-way drag and write forces to PE.
  6. Advance PE using its configured substeps.
  7. Rebuild transfer fields at updated particle positions.
  8. Compute conservative reaction-force projection for diagnostics only.
  9. Execute a self-contained Q2/P1 momentum and pressure solve without:
      - updateFBMGeometry
      - FictKNPR particle constraints
      - FBM stress-integral forces
      - FBM velocity boundary conditions

  10. Apply normal wall, inflow, slip, periodic, gravity, and constant-force behavior unchanged.

  The E-L routine will reproduce the existing low-level fluid assembly sequence directly. Shared FBM/E-L solver extraction is deferred to avoid changing
  validated FBM behavior.

  ## Transfer Rules

  - Kernel width is delta = el.kernel_width_factor * particle_diameter, default 2.5.
  - Kernel support is independent of mesh size.
  - Raw nonnegative quadrature weights are normalized over all fluid elements intersecting a particle’s support.
  - Wall and partition truncation use the same normalization path.
  - Contributions from owned and shadow particles are deduplicated by PE system ID.
  - Per-particle normalization is reduced between owner and shadow ranks before local assembly.
  - alpha_p conservation uses element quadrature:
    sum_K volume(K) * alpha_p(K) = sum_i volume(i).

  - Force spreading uses:
    b_a = -sum_i F_i integral(N_a W_i dx).

  - Shared Q2 contributions use the existing velocity-DOF communication/reduction machinery.
  - Clipping uses epsilon_f = max(1-alpha_p, el.eps_f_min) with default 0.4; every clipping event increments a global warning counter.
  - Q2 basis values are never used as particle-volume weights.

  ## Public Interfaces

  Introduce stable APIs equivalent to:

  call el_initialize(config, mesh, velocity_ndof)
  call el_update_volume_fields(particles, dt, fields)
  call el_sample_velocity(particles, velocity, samples)
  call el_evaluate_provisional_drag(particles, samples, forces)
  call el_spread_forces(particles, forces, force_rhs)
  call el_finalize()

  Configuration keys for this increment:

  simulation_mode
  el.kernel
  el.kernel_width_factor
  el.eps_f_min
  el.eps_f_relax
  el.drag_model = stokes | schiller_naumann
  el.apply_particle_forces

  Legacy frozen-trace force-kernel parameters remain supported as aliases with a deprecation message.

  ## Test Plan

  Add a dedicated MPI CTest executable, el_transfer_test, using the production mesh, quadrature, PE ownership, and shadow-particle path.

  - Kernel values are finite and nonnegative.
  - Constant-field interpolation is exact within assembly tolerance.
  - Linear-field interpolation has the expected second-order kernel error away from walls.
  - A central particle conserves volume.
  - A wall-adjacent particle conserves volume after renormalization.
  - A partition-straddling particle gives rank-independent volume.
  - Spread force sums componentwise to the negative particle force.
  - Symmetric particle arrangements produce symmetric fields.
  - Zero particles and zero force produce zero fields.
  - epsilon_f remains positive and clipping warnings are deterministic.
  - Restart reproduces epsilon_f, its history, and d epsilon_f/dt.

  Runtime regressions:

  - Existing q2p1_dns_drag FBM results remain unchanged.
  - q2p1_el_pipeflow completes a periodic constant-force smoke run without any FictKNPR particle constraints.
  - Frozen-field Stokes and Schiller–Naumann trajectories remain within existing tolerances after migration to shared APIs.

  ## Assumptions

  - Phase 1 is one-way coupled at runtime; fluid reaction-force insertion and semi-implicit drag belong to Phase 3.
  - Provisional drag is explicit and will be replaced or extended by Di Felice drag in Phase 2.
  - PE remains authoritative for particle positions, velocities, contacts, and substeps.
  - Fluid viscosity remains the carrier-fluid viscosity; no effective-viscosity model is introduced.

  • Proposed Plan


  # Complete Euler–Lagrange Phase 1

  ## Summary

  Replace the current global particle broadcast and element-center approximation with a scalable, quadrature-based particle–
  mesh transfer implementation. PE remains authoritative for particle ownership and dynamics; FeatFloWer maintains a
  transient coupling halo containing only particles whose kernels overlap another CFD partition.

  Phase 1 remains one-way coupled at runtime. The conservative fluid reaction vector is assembled for diagnostics and
  testing but is not inserted into the momentum equation.

  ## Implementation Changes

  ### Distributed Coupling Halo

  - Build a static MPI distributed-graph communicator from CFD partition bounding boxes and the global maximum kernel radius
    δmax.

  - Obtain δmax from the globally reduced maximum PE particle radius and el.kernel_width_factor.
  - Exchange only PE-owned particle records with graph neighbors whose CFD bounding boxes intersect the particle’s kernel-
    support sphere.

  - Use compact records containing system ID, owner rank, position, velocity, radius, and density.
  - Never insert coupling-halo records into PE or treat them as PE collision shadows.
  - Remove all per-particle MPI_Bcast and global particle iteration from EL_TRANSFER.
  - Deduplicate records by PE system ID and assert that every received record has exactly one owner.
  - Support periodic partition images when constructing overlap tests.

  At initialization, compare δmax against each partition’s minimum bounding-box width. If it exceeds any width:

  - Emit a prominent warning to terminal and protocol output.
  - Report δmax, the minimum partition width, affected rank count, and maximum number of partition layers crossed.
  - Explain that transfer remains correct because the graph directly connects all intersected ranks, but the decomposition
    is inefficient for the selected kernel.

  - Emit the warning once, not every timestep.

  ### Quadrature Transfer

  - Replace element-center sampling and nearest-element fallback with production hexahedral cubature using the existing
    CB3H, E013, Jacobian, and Q2 local-to-global machinery.

  - Reject a particle with no positive integrated kernel support as a fatal transfer error; do not silently deposit it into
    the nearest cell.

  - Compute per-particle raw kernel integrals on every participating rank.
  - Return partial normalization and field numerators to the PE owner.
  - Normalize over the complete distributed fluid support, providing wall and partition renormalization through the same
    path.

  - Kernel-sample Q2 velocity and its gradient, P1 pressure and gradient, and P0 void fraction.
  - Assemble:
      - P0 alpha_p, epsilon_f, epsilon_f_old, and deps_f_dt.
      - Q2 diagnostic reaction force using -F_i ∫N_a W_i dx.

  - Use E013Sum3 for shared Q2 force-vector contributions.
  - Count and globally report void-fraction clipping events.

  ### Runtime Data Flow

  For each fluid step:

  1. Read PE-owned particles.
  2. Refresh the sparse coupling halo.
  3. Integrate local kernel support and field numerators.
  4. Reduce partial samples to each PE owner.
  5. Evaluate provisional Stokes or Schiller–Naumann drag on the owner.
  6. Apply force only to the PE-owned particle and advance PE.
  7. Refresh the halo after particle motion.
  8. Reassemble volume fields and diagnostic reaction force.
  9. Execute Transport_q2p1_UxyzP_el without FBM geometry, constraints, stress-force evaluation, or particle stepping.

  Refactor the E-L transport path into a self-contained low-level Q2/P1 solve rather than controlling
  Transport_q2p1_UxyzP_fc_ext through el_transport_active. Existing FBM transport behavior must remain unchanged.

  ### Fields and I/O

  - Add explicit initialization, update, restart, output, and finalization APIs to EL_FIELDS.
  - Save and restore epsilon_f_old so deps_f_dt is restart-consistent.
  - Add VTK cell fields:
      - ELAlphaP
      - ELEpsilonF
      - ELDepsFDt

  - Add optional Q2 vector output ELReactionForce.
  - Update the E-L application’s documented output-field list and sample configuration.
  - Refactor q2p1_el_frozen_trace to use the shared sampling and provisional drag APIs.

  ## Public Interfaces

  Provide stable calls equivalent to:

  CALL EL_INITIALIZE(mesh, ilev, communicator, particles)
  CALL EL_REFRESH_COUPLING_HALO(particles)
  CALL EL_UPDATE_VOLUME_FIELDS(mesh, particles, dt)
  CALL EL_SAMPLE_FLUID_FIELDS(mesh, particles, samples)
  CALL EL_EVALUATE_PROVISIONAL_DRAG(particles, samples, forces)
  CALL EL_SPREAD_FORCES(mesh, particles, forces, force_rhs)
  CALL EL_WRITE_RESTART(output_index)
  CALL EL_READ_RESTART(start_path)
  CALL EL_FINALIZE()

  Retain these Phase 1 configuration keys:

  ELKernel                 = deen_poly | gaussian
  ELKernelWidthFactor      = 2.5
  ELEpsFMin                = 0.4
  ELEpsFRelax              = 0.0
  ELDragModel              = stokes | schiller_naumann
  ELApplyParticleForces    = Yes | No
  ELWriteDiagnostics       = Yes | No

  No global simulation_mode option is introduced.

  ## Test Plan

  - Unit tests for positive compact support, kernel cutoff, configuration validation, Stokes drag, and Schiller–Naumann
    drag.

  - Production-mesh transfer tests for:
      - Exact constant-field interpolation.
      - Expected second-order linear-field error away from walls.
      - Central and wall-adjacent particle-volume conservation.
      - Componentwise conservative Q2 force spreading.
      - Zero-particle and zero-force behavior.
      - Deterministic clipping counts.
      - Restart preservation of void-fraction history.

  - MPI CTests on 1, 2, and 8 ranks with particles crossing face, edge, and corner partition boundaries.
  - Add a case where δmax crosses more than one partition layer; verify conservation and the one-time decomposition warning.
  - Assert rank-independent volume, samples, and force sums within floating-point reduction tolerance.
  - Assert that no transfer path uses global particle broadcasts.
  - Run an E-L periodic constant-force smoke simulation and verify that no FBM particle constraints or stress-integral
    forces execute.

  - Build and run the existing DNS/FBM application to confirm unchanged behavior.

  ## Assumptions

  - Particle radii are fixed after PE initialization; changing particle size requires rebuilding the coupling graph.
  - CFD and PE rank ownership correspond on the PE worker communicator, excluding FeatFloWer rank zero.
  - Partition bounding boxes conservatively contain every local fluid element.
  - A large δmax relative to partition width is permitted with a clear warning because direct graph connections preserve
    correctness.

  - Fluid feedback and semi-implicit drag coupling remain Phase 3 work.