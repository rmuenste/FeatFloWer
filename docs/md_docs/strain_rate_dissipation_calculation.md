# Strain Rate Dissipation Integral Calculation

## Mathematical Background

This implementation calculates the strain rate dissipation integral:

$$\int_{\Pi_o} \mathbf{D(u)} : \mathbf{D(u)} \, dx$$

Where:
- $\mathbf{D(u)} = \frac{1}{2} (\nabla \mathbf{u} + (\nabla \mathbf{u})^T)$ is the rate of strain tensor
- `:` denotes the Frobenius inner product: $\mathbf{A} : \mathbf{B} = \sum_i \sum_j A_{ij} B_{ij}$
- $\Pi_o$ is the fluid domain (excluding particle interiors)

In 3D, this expands to:
$$\mathbf{D(u)} : \mathbf{D(u)} = D_{11}^2 + D_{22}^2 + D_{33}^2 + 2(D_{12}^2 + D_{13}^2 + D_{23}^2)$$

Where the rate of strain tensor components are:
- $D_{11} = \frac{\partial u_1}{\partial x_1}$, $D_{22} = \frac{\partial u_2}{\partial x_2}$, $D_{33} = \frac{\partial u_3}{\partial x_3}$
- $D_{12} = \frac{1}{2}(\frac{\partial u_1}{\partial x_2} + \frac{\partial u_2}{\partial x_1})$
- $D_{13} = \frac{1}{2}(\frac{\partial u_1}{\partial x_3} + \frac{\partial u_3}{\partial x_1})$
- $D_{23} = \frac{1}{2}(\frac{\partial u_2}{\partial x_3} + \frac{\partial u_3}{\partial x_2})$

## Implementation Files

### Core Implementation: `calculate_dissipation_integral.f90`

**Subroutine:** `Get_DissipationIntegral(QuadSc, ALPHA, total_dissipation, mfile)`

**Parameters:**
- `QuadSc`: TQuadScalar type containing velocity field data (`valU`, `valV`, `valW`)
- `ALPHA`: Particle indicator function array (0 in fluid, particle_id in particles)
- `total_dissipation`: Output - computed integral value
- `mfile`: File unit for logging/output

**Key Features:**
- Uses Q2 finite element basis functions for accurate velocity gradient computation
- Follows FeatFloWer's element loop + cubature point structure (same as `GetForces`)
- Implements interface element detection to exclude particle interiors
- Uses MPI reduction to sum contributions across all processors
- 3×3×3 Gaussian cubature for high accuracy integration

### Interface Functions: `dissipation_interface.f90`

**Simple Interface:** `Call_DissipationCalculation(mfile)`
- Easy-to-use wrapper that can be called from main transport routines
- Automatically uses current velocity field and particle configuration

**Effective Viscosity:** `Calculate_EffectiveViscosity(reference_shear_rate, effective_viscosity, mfile)`
- Calculates effective viscosity from dissipation rate: $\mu_{eff} = \frac{2\Phi}{\dot{\gamma}^2}$
- Where $\Phi$ is dissipation rate per unit volume and $\dot{\gamma}$ is imposed shear rate

## Integration into FeatFloWer

### In Main Transport Routine

Add to `Transport_q2p1_UxyzP_fc_ext` after velocity solve:

```fortran
! After velocity correction step
CALL Velocity_Correction()

! Calculate strain rate dissipation (optional)
IF (calculate_dissipation) THEN
  CALL Call_DissipationCalculation(mfile)
END IF

! Continue with particle updates
call fbm_updateForces(...)
```

### Compilation

Add to your Makefile:
```makefile
DISSIPATION_OBJS = calculate_dissipation_integral.o dissipation_interface.o

your_application: $(DISSIPATION_OBJS) $(OTHER_OBJS)
	$(FC) -o $@ $^ $(LDFLAGS)
```

## Computational Details

### Element Classification
- **Fluid elements**: All nodes have `ALPHA = 0`, full contribution
- **Particle elements**: All nodes have `ALPHA ≠ 0`, zero contribution (rigid body motion)
- **Interface elements**: Mixed nodes, simplified partial contribution

### Numerical Integration
- Uses same 27-point 3×3×3 Gaussian cubature as other FeatFloWer routines
- Exact for polynomials up to degree 5 (adequate for Q2 velocity gradients)
- Physical coordinate transformation via isoparametric mapping

### Performance Considerations
- Computational cost similar to `GetForces` calculation
- Scales linearly with number of elements
- MPI parallel with automatic load balancing

## Applications

### 1. Effective Viscosity Calculation
For suspension rheology:
```fortran
REAL*8 :: shear_rate = 1.0D0  ! Applied shear rate
REAL*8 :: mu_eff

CALL Calculate_EffectiveViscosity(shear_rate, mu_eff, mfile)
```

### 2. Energy Balance Verification
The dissipation integral should balance energy input:
$$\text{Power input} = \mu_{fluid} \int_{\Pi_o} \mathbf{D(u)} : \mathbf{D(u)} \, dx$$

### 3. Convergence Studies
Monitor dissipation integral during:
- Mesh refinement studies
- Time step convergence analysis
- Steady-state convergence

## Validation

### Test Cases
1. **Couette Flow**: Analytical solution $\dot{\gamma} = V/H$, uniform dissipation
2. **Poiseuille Flow**: Parabolic profile with known dissipation rate
3. **Stokes Flow Around Sphere**: Compare with analytical drag calculations

### Expected Results
- Should converge with mesh refinement (Q2 elements → 4th order accuracy)
- Should scale correctly with viscosity and shear rate
- Should remain positive definite (rate of strain tensor properties)

## Further Development

### Enhancements Needed
1. **Complete Jacobian Implementation**: Currently simplified
2. **Precise Interface Integration**: More sophisticated fluid domain detection
3. **Volume Fraction Weighting**: For partially filled interface elements
4. **Output Integration**: Connect to existing FeatFloWer output systems

### Potential Extensions
1. **Stress Power**: $\int \boldsymbol{\sigma} : \mathbf{D(u)} \, dx$ for non-Newtonian flows
2. **Component-wise Analysis**: Separate normal and shear contributions
3. **Particle Contribution**: Include lubrication and contact dissipation

## File Structure
```
FeatFloWer/
├── calculate_dissipation_integral.f90    # Core implementation
├── dissipation_interface.f90             # Easy-to-use interfaces  
└── docs/md_docs/
    └── strain_rate_dissipation_calculation.md  # This documentation
```

## References
- FeatFloWer documentation on finite element implementation
- GetForces subroutine structure (`source/src_quadLS/QuadSc_force.f90:25`)
- Element routines E013 (`source/Elements/e013.f`)
- Cubature routines CB3H (`extern/libraries/feat3d/src/cb3h.f`)