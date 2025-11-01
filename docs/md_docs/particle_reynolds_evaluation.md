# Particle Reynolds Number Evaluation

## Motivation

During coupled fluid–particle simulations it is often useful to track the instantaneous Reynolds number associated with each rigid body.  The value
$$\mathrm{Re}_p = \frac{\rho_f \, \lVert \mathbf{u}_f(\mathbf{x}_c) - \mathbf{u}_p \rVert \, 2R}{\mu_f}$$
provides a quick indicator of the drag regime experienced by a particle of radius $R$ translating with velocity $\mathbf{u}_p$ through a carrier fluid of density $\rho_f$ and dynamic viscosity $\mu_f$.  Having these values available inside the time-integration loop helps with online diagnostics (e.g. slip-flow monitoring, post-processing thresholds, model tuning).

## Implementation Summary

* **Module:** `source/src_fbm/fbm_particle_reynolds.f90`
* **Entry point:** `fbm_compute_particle_reynolds`
* **Call site:** `Transport_q2p1_UxyzP` in `source/src_quadLS/QuadSc_main.f90` immediately after `fbm_updateForces`.

### Data Flow

1. **Inputs**
   - Velocity DoF arrays `QuadSc%valU,V,W` (Q2 field) and elementwise fluid viscosity `Viscosity(:)`.
   - Particle state `myFBM%ParticleNew(:)%Position`, `Velocity`, `sizes(1)` (radius), and current fluid density `Properties%Density(1)`.
   - Mesh description (`mg_mesh%level`) already in scope for FBM operations.

2. **Per-Particle Sampling**
   - For each particle on non-master MPI ranks, the routine loops over local elements, quickly eliminating obvious misses via an AABB test.
   - A candidate element is confirmed with `fbmaux_PointInHex`, which returns the reference coordinates of the particle centroid inside that hexahedron.
   - Using `ELE`/`NDFGL`, the Q2 basis functions at the reference point are evaluated and the surrounding velocity is interpolated.

3. **Reynolds Calculation**
   - Slip velocity is formed from the interpolated fluid velocity and the particle translational velocity.
   - The diameter is taken as `2 * sizes(1)`; the local dynamic viscosity comes directly from the element-indexed `Viscosity` array.
   - If viscosity or radius is zero the routine assigns `Re_p = 0` to keep the value well-defined.

4. **Parallel Reduction**
   - Each rank accumulates partial results (`re_local`, `re_weight`) which are reduced with `COMM_SUMMN`.
   - The global per-particle result is stored in the new array `myFBM%ParticleRe(:)` (added to `tFBM`).
   - The element index used for sampling is optionally captured through `myFBM%iel_ug` for re-use by other diagnostics.

5. **Logging**
   - Rank 1 prints the minimum and maximum Reynolds numbers to the protocol file to provide quick feedback during runs.

## Usage Notes

* The routine expects the viscosity array to be sized at least `NEL` for the finest level; the call in `Transport_q2p1_UxyzP` uses the `Viscosity` workspace already populated before the FBM force evaluation.
* Reynolds numbers are zeroed automatically when no particles exist or when the reduction detects missing samples.
* The values are stored alongside the FBM force buffer, so they can be exported with custom routines just like hydrodynamic forces (e.g. via the C++ PE interface or Fortran post-processing).

## Related Documentation

* [Velocity Evaluation at element midpoints](velocity_midpoint_evaluation.md) – illustrates the same Q2 interpolation workflow.
* [Hydrodynamic force computation](hydrodynamic_force_computation.md) – describes the volume integration that precedes Reynolds evaluation.
* [Strain-rate dissipation calculation](strain_rate_dissipation_calculation.md) – another example of post-processing that uses the Q2 field and MPI reductions.

### 1. Sample on (or just outside) the particle surface

  Instead of sampling at the centre, first project to the particle surface along a direction n (for example the local outward normal) and evaluate the Q2 field at
  x_surf = x_p + (R + ε) · n.
  You can reuse fbmaux_PointInHex to locate the element containing each sample point; only the point location changes. Using a small ε>0 keeps you in the fluid
  region, so the interpolated velocity reflects the true boundary-layer flow. Slip is then ‖u_f(x_surf) − u_p‖. You can pick a few quadrature directions (e.g. 6
  face normals or a spherical design) and average the resulting slips to reduce directional noise.

  ———

  ### 2. Nearest-fluid DOF averaging

  Leverage the FBM distance data (Distance, FictKNPR) to identify surrounding Q2 nodes that belong to the fluid (distance>0). Gather a small ball of those DOFs,
  interpolate/average their velocities (e.g. inverse-distance weighting), and call that u_f(x_p). Because the DOFs themselves lie outside the particle, they still
  carry the fluid velocity. This avoids an additional point-location step, at the cost of a small unstructured stencil.

  ———

  ### 3. Extrapolation from the interface element

  Given the drag force F on a particle and its Stokes drag coefficient 3πμd (or an empirical correction for finite Re), you can recover an effective slip speed from
  ‖Δu‖ ≈ ‖F‖ / (3πμd). This ties the Reynolds number to the actually transferred momentum rather than the velocity field. It is attractive when you trust the force
  evaluation more than local velocity interrogation, but it depends on a model for the drag coefficient.