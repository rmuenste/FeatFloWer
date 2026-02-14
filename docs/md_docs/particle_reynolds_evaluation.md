# Particle Reynolds Number Evaluation

## Motivation

During coupled fluid–particle simulations it is often useful to track the instantaneous Reynolds number associated with each rigid body.  The value
$$\mathrm{Re}_p = \frac{\rho_f \, \lVert \mathbf{u}_f(\mathbf{x}_c) - \mathbf{u}_p \rVert \, 2R}{\mu_f}$$
provides a quick indicator of the drag regime experienced by a particle of radius $R$ translating with velocity $\mathbf{u}_p$ through a carrier fluid of density $\rho_f$ and dynamic viscosity $\mu_f$.  Having these values available inside the time-integration loop helps with online diagnostics (e.g. slip-flow monitoring, post-processing thresholds, model tuning).

## Available Methods: Center, Interface, and Extended Interface Sampling

FeatFloWer provides **three distinct approaches** for computing particle Reynolds numbers, each suited to different simulation regimes:

| Method | Subroutine | Best For | Key Requirement |
|--------|-----------|----------|-----------------|
| **Center Sampling** | `fbm_compute_particle_reynolds` | Under-resolved, point-particle simulations | Particle position & velocity |
| **Interface Sampling** | `fbm_compute_particle_reynolds_interface` | DNS, fully-resolved simulations | Particle indicator field `ALPHA` |
| **Extended Interface** | `fbm_compute_particle_reynolds_interface_extended` | DNS with enhanced accuracy | Particle indicator field `ALPHA` + `kadj` neighbors |

**Why three methods?**

- In **under-resolved simulations**, particles are smaller than the mesh spacing and the fluid velocity field is meaningful at the particle center. Center sampling directly evaluates the Q2 velocity field at the particle centroid.

- In **DNS/fully-resolved simulations**, particles are explicitly resolved on the mesh using the fictitious boundary method. The velocity at the particle center is constrained to match the rigid body motion, making center-based sampling invalid. Interface sampling addresses this by averaging velocities from elements that cross the particle boundary, capturing the actual fluid motion near the surface.

- The **extended interface method** improves upon standard interface sampling by including a second layer of elements beyond the immediate interface. This provides a larger sampling region, smoother velocity averages, and reduced sensitivity to local mesh artifacts. Distance-based weighting ensures interface elements have more influence than second-layer elements.

## Implementation Summary

* **Module:** `source/src_fbm/fbm_particle_reynolds.f90`
* **Available methods:**
  - `fbm_compute_particle_reynolds` – Center sampling (suitable for under-resolved simulations)
  - `fbm_compute_particle_reynolds_interface` – Interface sampling (suitable for DNS/fully-resolved simulations)
  - `fbm_compute_particle_reynolds_interface_extended` – Extended interface with second-layer sampling (enhanced DNS accuracy)
* **Call site:** `Transport_q2p1_UxyzP` in `source/src_quadLS/QuadSc_main.f90` immediately after `fbm_updateForces`.

### Data Flow (Center-Based Method)

This section describes the original center-based sampling approach. For the interface-based method suitable for DNS, see the dedicated section below.

1. **Inputs**
   - Velocity DoF arrays `QuadSc%valU,V,W` (Q2 field) and elementwise fluid viscosity `Viscosity(:)`.
   - Particle state `myFBM%ParticleNew(:)%Position`, `Velocity`, `sizes(1)` (radius), and current fluid density `Properties%Density(1)`.
   - Element evaluation function `ELE` (typically `E013` for Q2 elements) passed as function pointer.
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

## Interface-Based Reynolds Calculation (DNS Method)

**Status: FULLY IMPLEMENTED** via `fbm_compute_particle_reynolds_interface`

The interface method addresses a fundamental limitation of center-based sampling in DNS/fully-resolved simulations: when particles are resolved on the computational mesh, the velocity at the particle center is constrained to match the rigid body motion, making center-based slip velocity calculations meaningless.

### Algorithm Overview

The interface method identifies and samples velocities at **interface elements** – elements that cross the particle boundary, having some degrees of freedom (DOFs) inside the particle and others in the surrounding fluid. By averaging velocities from these boundary-crossing elements, the method captures the actual fluid motion near the particle surface.

### Detailed Implementation

#### 1. Interface Element Detection

For each particle `ip`, the algorithm loops through all mesh elements and uses the particle indicator field `ALPHA` to classify DOFs:

- `ALPHA(ig) = 0` → DOF is in the fluid region (outside all particles)
- `ALPHA(ig) = ip` → DOF belongs to particle `ip` (inside the particle)

For each element with `IDFL` DOFs (27 for Q2 elements), the routine counts:
- `NJALFA` = number of DOFs with `ALPHA = 0` (outside particle)
- `NIALFA` = number of DOFs with `ALPHA = ip` (inside particle `ip`)

**Interface element criterion:**
```fortran
IF (NJALFA .EQ. 27 .OR. NIALFA .EQ. 27) CYCLE  ! Skip fully inside/outside
! Otherwise: 0 < NJALFA < 27 AND 0 < NIALFA < 27 → Interface element!
```

#### 2. Velocity Sampling at Element Centers

For each identified interface element:

1. **Setup element geometry**: Extract vertex coordinates and establish the mapping from reference space to physical space
2. **Evaluate basis functions**: Call the element evaluation function at the reference element center `(ξ₁, ξ₂, ξ₃) = (0, 0, 0)`
3. **Interpolate velocity**: Compute the Q2 velocity field at the element center:
   ```fortran
   vel_sample = Σᵢ uᵢ · φᵢ(0, 0, 0)
   ```
   where `uᵢ` are the velocity DOF values and `φᵢ` are the Q2 basis functions

4. **Accumulate**: Add `vel_sample` to `vel_sum` and increment the counter `weight_sum`

#### 3. Reynolds Number Computation

After scanning all elements for a given particle:

1. **Average interface velocities**:
   ```
   vel_avg = vel_sum / weight_sum
   ```

2. **Compute slip velocity**: The difference between the averaged fluid velocity and the particle's rigid body velocity:
   ```
   slip = vel_avg - particle_velocity
   speed = ||slip||
   ```

3. **Calculate Reynolds number**:
   ```
   Re_p = (ρ_fluid × speed × diameter) / μ_avg
   ```
   where `μ_avg` is the domain-averaged dynamic viscosity

#### 4. Parallel Reduction

Since particles may have interface elements on multiple MPI ranks:

- Each rank computes partial `re_local(ip)` and `re_weight(ip)` values
- `COMM_SUMMN` performs a global sum across all ranks
- Final Reynolds number is normalized: `ParticleRe(ip) = re_local(ip) / re_weight(ip)`

### Key Advantages

1. **DNS-appropriate**: Samples velocity in the fluid region near the particle surface, avoiding the rigid-body-constrained interior
2. **Automatic boundary detection**: Uses the existing `ALPHA` field without requiring explicit surface reconstruction
3. **Robust averaging**: Averages over multiple interface elements reduces sensitivity to local mesh artifacts
4. **MPI-compatible**: Natural domain decomposition with straightforward reduction

### Input Requirements

The interface method requires one additional input compared to center sampling:

- `ALPHA(:)` – Integer array mapping each DOF to its particle ID (0 for fluid, `ip` for particle `ip`)
  - Typically named `FictKNPR` in FeatFloWer applications
  - Must be sized to cover all Q2 DOFs on the finest mesh level

## Extended Interface Method (Enhanced DNS Approach)

**Status: FULLY IMPLEMENTED** via `fbm_compute_particle_reynolds_interface_extended`

The extended interface method builds upon the standard interface approach by incorporating a **second layer** of elements beyond the immediate particle interface. This "nice-to-have" accuracy improvement provides more robust Reynolds number estimates in DNS simulations.

### Algorithm Overview

The extended method uses a **two-pass classification algorithm** with intelligent neighbor expansion:

1. **Pass 1:** Classify all elements as inside (1) or interface (2) using the ALPHA field
2. **Pass 2:** Expand from interface elements via face neighbors using `kadj` connectivity

### Detailed Implementation

#### Element Classification Array

Instead of maintaining separate lists, a single `elem_class(nel_local)` array provides O(1) lookups:

```fortran
! Classification values:
! 0 = unclassified (fluid, far from particle)
! 1 = fully inside particle
! 2 = interface element (crosses particle boundary)
! 3 = second layer (fluid neighbor of interface element)
```

#### Pass 1: Identify Interface Elements

Same as standard interface method - count DOFs inside/outside particle:

```fortran
IF (NJALFA == 0 .AND. NIALFA > 0) THEN
  elem_class(iel) = 1  ! Fully inside
ELSE IF (NJALFA > 0 .AND. NIALFA > 0) THEN
  elem_class(iel) = 2  ! Interface
END IF
```

#### Pass 2: Expand to Second Layer

For each interface element, check all 6 face neighbors via the `kadj` connectivity array:

```fortran
DO iface = 1, 6  ! 6 faces of hexahedron
  neighbor = mg_mesh%level(ilev)%kadj(iface, iel)

  IF (neighbor > 0 .AND. elem_class(neighbor) == 0) THEN
    ! Check if neighbor is fluid-dominated (≥20 fluid DOFs out of 27)
    IF (NJALFA_neighbor >= 20) THEN
      elem_class(neighbor) = 3  ! Add to second layer
    END IF
  END IF
END DO
```

**Fluid-dominated filter:** Only elements with ≥20 fluid DOFs are added to the second layer. This prevents sampling elements that are mostly inside other particles or still heavily influenced by the solid boundary.

#### Distance-Based Weighting

Elements closer to the particle surface have more influence on the final Reynolds number:

```fortran
! Compute distance from element center to particle center
dist_to_particle = ||elem_center - particle_position||

! Weight inversely proportional to distance
weight = 1.0 / (1.0 + dist_to_particle)

! Accumulate weighted velocity
vel_sum += weight * vel_sample
weight_sum += weight
```

This ensures that:
- Interface elements (closer to particle) contribute more
- Second-layer elements (further away) contribute less
- The transition is smooth and continuous

### Key Advantages Over Standard Interface Method

1. **Larger sampling region:** More elements provide better statistical averaging
2. **Reduced mesh sensitivity:** Less dependent on exact interface element location
3. **Smoother results:** Distance weighting provides gradual influence decay
4. **Robust filtering:** Fluid-dominated criterion prevents contamination from solid regions
5. **Efficient implementation:** O(1) classification lookups via integer array

### Mesh Connectivity Requirements

The extended method leverages FeatFloWer's `kadj` array:
- `mg_mesh%level(ilev)%kadj(iface, iel)` stores the neighbor element through face `iface`
- Value = 0 for boundary faces (properly handled in the algorithm)
- Available on all mesh levels after mesh initialization

### Configurable Parameters

```fortran
SUBROUTINE fbm_compute_particle_reynolds_interface_extended(
    U1, U2, U3, ALPHA, DVISC, RHOFLUID, mfile, ELE, n_layers)
    INTEGER, INTENT(IN), OPTIONAL :: n_layers  ! Default = 2
```

- `n_layers = 1`: Equivalent to standard interface method (no expansion)
- `n_layers = 2`: Interface + second layer (default, recommended)
- `n_layers ≥ 3`: Future extension for additional layers

### Performance Considerations

- **Memory overhead:** One integer per element (`elem_class` array)
- **Computational cost:** ~2× element loops compared to standard interface method
- **MPI compatibility:** Natural domain decomposition, same reduction pattern
- **Recommended use:** DNS simulations where accuracy is critical and computational cost is acceptable

## Stokes Number Calculation (Particle Inertia Characterization)

**Status: FULLY IMPLEMENTED** via `fbm_compute_particle_stokes`

The Stokes number characterizes the ratio of particle inertia to fluid viscous forces, providing critical insight into particle response to flow variations. It is particularly useful for understanding settling behavior, flow following, and turbulent dispersion.

### Physical Meaning

The Stokes number is defined as:

$$\mathrm{St} = \frac{\rho_p}{\rho_f} \cdot \frac{\mathrm{Re}_p}{18}$$

where:
- $\rho_p$ = particle material density (from PE library)
- $\rho_f$ = fluid density
- $\mathrm{Re}_p$ = particle Reynolds number (from any of the three methods above)

**Physical interpretation:**
- **St ≪ 1**: Particle closely follows fluid motion (tracer behavior)
- **St ~ O(1)**: Particle response time comparable to flow time scales (complex dynamics)
- **St ≫ 1**: Particle dominated by inertia, weakly coupled to fluid (ballistic motion)

### Algorithm

The Stokes number calculation is straightforward once Reynolds numbers are available:

1. **Prerequisite check**: Ensures `ParticleRe` array has been populated by one of the Reynolds calculation methods
2. **Particle data retrieval**: Fetches particle density from the PE library via the extended `particleData_t` interface
3. **Formula application**: Computes `St = (ρ_p/ρ_f) * (Re_p/18)` for each particle
4. **Storage**: Results stored in `myFBM%ParticleSt(:)` array

### Implementation Details

The subroutine follows the same MPI communication pattern as Reynolds calculation:

```fortran
SUBROUTINE fbm_compute_particle_stokes(RHOFLUID, mfile)
  ! Input: RHOFLUID - fluid density
  ! Input: mfile - protocol file handle for logging

  ! Loop over all particles
  DO ip = 1, numParticles
    IF (myid .NE. 0) THEN
      ! Fetch particle data including density
      CALL getPartStructByIdx(ip-1, particle)
      rho_particle = particle%density
    END IF

    ! Compute Stokes number from Reynolds number
    IF (RHOFLUID .GT. 0D0) THEN
      myFBM%ParticleSt(ip) = (rho_particle/RHOFLUID) * (myFBM%ParticleRe(ip)/18D0)
    END IF
  END DO
END SUBROUTINE
```

### Data Flow

The Stokes number calculation relies on data from two sources:

1. **From CFD solver**: Fluid density `RHOFLUID` (typically from `Properties%Density(1)`)
2. **From PE library**: Particle density `ρ_p` via `particleData_t` struct
3. **From Reynolds calculation**: Particle Reynolds number `ParticleRe(ip)` (must be called first)

### Material Density Integration

Particle densities are automatically fetched from the PE library through the extended interface:

- **C++ side** (`libs/pe/src/interface/object_queries.cpp`): Casts rigid bodies to specific geometry types and extracts material density via `Material::getDensity(mat)`
- **Fortran side** (`source/src_particles/dem_query.f90`): Extended `tParticleData` type includes `density` field
- **Automatic transfer**: No manual density specification required - values come directly from PE's material database

### Prerequisites

Before calling `fbm_compute_particle_stokes`:

1. **Initialize particle system**: PE library must be initialized with particles having defined materials
2. **Compute Reynolds numbers**: Call one of the Reynolds calculation methods:
   - `fbm_compute_particle_reynolds` (center)
   - `fbm_compute_particle_reynolds_interface` (interface)
   - `fbm_compute_particle_reynolds_interface_extended` (extended interface)
3. **Set fluid density**: Ensure `RHOFLUID` is properly set in the properties

### Output and Logging

The routine outputs summary statistics to the protocol file:

```
Particle Stokes numbers computed:
  Min St = 0.152
  Max St = 2.847
```

Results are stored in `myFBM%ParticleSt(:)` and can be exported alongside Reynolds numbers and hydrodynamic forces.

## Usage Notes

* Both routines expect the viscosity array to be sized at least `NEL` for the finest level.
* Reynolds numbers are zeroed automatically when no particles exist or when the reduction detects missing samples.
* The Reynolds values are stored in `myFBM%ParticleRe(:)` and Stokes values in `myFBM%ParticleSt(:)`. Both can be exported with custom routines just like hydrodynamic forces.
* **Element function:** All Reynolds methods require the element evaluation function `E013` to be passed as the last parameter.
* **Stokes number dependency:** `fbm_compute_particle_stokes` must be called AFTER one of the Reynolds calculation methods.
* **Choosing the method:**
  - **Center sampling** (`fbm_compute_particle_reynolds`): Use for under-resolved or point-particle simulations where the particle interior is not explicitly resolved on the mesh
  - **Interface sampling** (`fbm_compute_particle_reynolds_interface`): Use for DNS or fully-resolved simulations where particles are represented on the mesh with the fictitious boundary method and the `ALPHA` field is available
  - **Extended interface sampling** (`fbm_compute_particle_reynolds_interface_extended`): Use for DNS simulations where enhanced accuracy and smoother results are desired, particularly with fine meshes or when reducing sensitivity to mesh artifacts
  - **Stokes number** (`fbm_compute_particle_stokes`): Always called after Reynolds calculation to characterize particle inertia

### Example Usage

**Center Sampling (under-resolved):**
```fortran
CALL fbm_compute_particle_reynolds( &
     QuadSc%valU, QuadSc%valV, QuadSc%valW, &
     Viscosity, Properties%Density(1), mfile, &
     E013)  ! Element evaluation function
```

**Interface Sampling (DNS):**
```fortran
CALL fbm_compute_particle_reynolds_interface( &
     QuadSc%valU, QuadSc%valV, QuadSc%valW, &
     FictKNPR,  & ! ALPHA array (particle indicator)
     Viscosity, Properties%Density(1), mfile, &
     E013)  ! Element evaluation function
```

**Extended Interface Sampling (DNS with enhanced accuracy):**
```fortran
CALL fbm_compute_particle_reynolds_interface_extended( &
     QuadSc%valU, QuadSc%valV, QuadSc%valW, &
     FictKNPR,  & ! ALPHA array (particle indicator)
     Viscosity, Properties%Density(1), mfile, &
     E013, &      ! Element evaluation function
     2)           ! n_layers (optional, default=2)
```

**Stokes Number (must follow Reynolds calculation):**
```fortran
! First compute Reynolds numbers using any method above
CALL fbm_compute_particle_reynolds_interface( &
     QuadSc%valU, QuadSc%valV, QuadSc%valW, &
     FictKNPR, Viscosity, Properties%Density(1), mfile, E013)

! Then compute Stokes numbers
CALL fbm_compute_particle_stokes( &
     Properties%Density(1), & ! Fluid density
     mfile)                   ! Protocol file handle

! Results available in:
!   myFBM%ParticleRe(:) - Reynolds numbers
!   myFBM%ParticleSt(:) - Stokes numbers
```

## Related Documentation

* [Velocity Evaluation at element midpoints](velocity_midpoint_evaluation.md) – illustrates the same Q2 interpolation workflow.
* [Hydrodynamic force computation](hydrodynamic_force_computation.md) – describes the volume integration that precedes Reynolds evaluation.
* [Strain-rate dissipation calculation](strain_rate_dissipation_calculation.md) – another example of post-processing that uses the Q2 field and MPI reductions.

## Alternative Approaches (For Future Development)

While the center-based and interface-based methods cover most practical scenarios, alternative approaches may be useful for specialized applications:

### 1. Surface Sampling

Instead of sampling at the particle center or element centers, project to the particle surface along multiple directions and evaluate the Q2 field at `x_surf = x_p + (R + ε) · n`. Using a small `ε > 0` keeps sampling points in the fluid region. This approach provides directionally-resolved slip velocities but requires additional point-location calls to `fbmaux_PointInHex`.

### 2. Nearest-Fluid DOF Averaging

Leverage the FBM distance field to identify Q2 nodes in the fluid region near the particle surface. Use inverse-distance weighting to interpolate velocities from these DOFs. This avoids point-location overhead but requires careful stencil selection to maintain accuracy.

### 3. Force-Based Extrapolation

Infer slip velocity from the computed drag force `F` using a drag correlation: `‖Δu‖ ≈ ‖F‖ / (C_d · πμd)`. This approach ties Reynolds number directly to momentum transfer but requires an appropriate drag coefficient model for the current flow regime.