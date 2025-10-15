# Velocity Evaluation at Element Midpoints

## Mathematical Background

For Q2 finite elements, the velocity at any point within an element is computed using:

$$\mathbf{u}(\xi_1, \xi_2, \xi_3) = \sum_{j=1}^{27} \mathbf{u}_j \phi_j(\xi_1, \xi_2, \xi_3)$$

Where:
- $\mathbf{u}_j$ are the nodal velocity vectors (27 nodes for Q2 hexahedral elements)
- $\phi_j(\xi_1, \xi_2, \xi_3)$ are the Q2 basis functions
- $(\xi_1, \xi_2, \xi_3)$ are reference coordinates in $[-1,1]^3$

For element midpoints, we evaluate at the element center: $(\xi_1, \xi_2, \xi_3) = (0, 0, 0)$.

The physical coordinates of element centers are computed using trilinear mapping:
$$\mathbf{x}_{center} = \frac{1}{8} \sum_{k=1}^{8} \mathbf{x}_k$$

Where $\mathbf{x}_k$ are the coordinates of the 8 element vertices.

## Implementation Files

### Core Implementation: `evaluate_velocity_midpoint.f90`

#### Main Subroutine: `Evaluate_Velocity_Midpoint`

```fortran
SUBROUTINE Evaluate_Velocity_Midpoint(QuadSc, velocity_midpoints, element_centers, mfile)
```

**Parameters:**
- `QuadSc`: TQuadScalar type containing velocity field data
- `velocity_midpoints(3, NEL)`: Output array of velocity components at element centers
- `element_centers(3, NEL)`: Output array of physical coordinates of element centers
- `mfile`: File unit for logging

**Features:**
- Evaluates Q2 basis functions at element center $(0,0,0)$
- Computes both velocity and position for each element
- Handles proper array dimensioning and error checking
- Compatible with FeatFloWer's element numbering and data structures

#### Specialized Functions:

1. **`Evaluate_Single_Element_Velocity`** - Evaluate velocity at any point within a specific element
2. **`Write_Velocity_Midpoints_VTK`** - Export results to VTK format for visualization

### Interface Module: `velocity_midpoint_interface.f90`

#### Simple Interface: `Call_VelocityMidpointEvaluation`

```fortran
SUBROUTINE Call_VelocityMidpointEvaluation(QuadSc, output_vtk, mfile)
```

Easy-to-use wrapper that:
- Automatically handles memory allocation
- Provides optional VTK output for ParaView visualization
- Reports velocity statistics
- Manages cleanup

#### Analysis Functions:

1. **`Report_VelocityMidpoint_Statistics`** - Compute min/max/average/RMS statistics
2. **`Find_Elements_With_HighVelocity`** - Identify elements exceeding velocity threshold
3. **`Evaluate_Velocity_Along_Line`** - Create velocity profiles along specified lines

## Usage Examples

### Basic Usage

```fortran
USE VelocityMidpoint_Interface

! After velocity solve in main transport routine
CALL Call_VelocityMidpointEvaluation(QuadSc, .TRUE., mfile)
```

### Advanced Usage with Analysis

```fortran
! Allocate arrays
REAL*8, ALLOCATABLE :: vel_mid(:,:), centers(:,:)
INTEGER :: high_vel_elements(1000), num_found

ALLOCATE(vel_mid(3, NEL), centers(3, NEL))

! Evaluate velocities
CALL Evaluate_Velocity_Midpoint(QuadSc, vel_mid, centers, mfile)

! Find high-velocity regions
CALL Find_Elements_With_HighVelocity(QuadSc, 5.0D0, high_vel_elements, num_found, mfile)

! Create VTK output for visualization
CALL Write_Velocity_Midpoints_VTK(vel_mid, centers, 'flow_analysis', mfile)
```

### Integration with Main Transport Routine

Add to `Transport_q2p1_UxyzP_fc_ext`:

```fortran
! After velocity correction
CALL Velocity_Correction()

! Optional velocity analysis
IF (analyze_velocity_field) THEN
  CALL Call_VelocityMidpointEvaluation(QuadSc, .TRUE., mfile)
END IF

! Continue with particle updates
call fbm_updateForces(...)
```

## Output Formats

### Console Output
```
Evaluating velocity at element midpoints...
Number of elements: 125000

=== Element-Centered Velocity Statistics ===
Number of elements: 125000

U-component:
  Min/Max/Avg:  -2.3456E+00   3.4567E+00   1.2345E-02
  RMS:           1.4567E+00

V-component:
  Min/Max/Avg:  -1.8765E+00   2.9876E+00  -3.4567E-03
  RMS:           9.8765E-01

W-component:
  Min/Max/Avg:  -4.5678E-01   5.6789E-01   2.3456E-04
  RMS:           2.3456E-01

Magnitude:
  Min/Max/Avg:   1.2345E-03   4.5678E+00   1.8765E+00
  RMS:           1.9876E+00
=============================================
```

### VTK Output

Generated VTK files contain:
- Element center coordinates as points
- Velocity vectors at each point
- Velocity magnitude as scalar field
- Compatible with ParaView, VisIt, and other VTK-based tools

## Applications

### 1. Flow Visualization
- Create streamlines starting from element centers
- Identify recirculation zones and stagnation points
- Visualize velocity magnitude distributions

### 2. Quality Assessment
- Check for unphysical velocity values
- Identify elements requiring mesh refinement
- Monitor convergence during iterative solves

### 3. Post-Processing Analysis
- Extract velocity profiles along specific lines
- Calculate flow rates through cross-sections
- Analyze velocity gradients and shear rates

### 4. Validation and Verification
- Compare with analytical solutions for simple geometries
- Verify mass conservation using element-centered data
- Create reference solutions for method development

## Performance Considerations

### Memory Usage
- Arrays scale as `O(NEL)` where NEL is number of elements
- Minimal memory footprint for velocity evaluation
- Optional arrays can be deallocated immediately after use

### Computational Cost
- Linear scaling with number of elements
- Basis function evaluation is cached for efficiency
- Typically < 1% of total simulation time

### Parallel Efficiency
- Element-wise operations are embarrassingly parallel
- No communication required between processors
- Scales perfectly with number of MPI processes

## Technical Details

### Q2 Basis Functions
- 27-node triquadratic hexahedral elements
- Uses `E013A` subroutine for basis function evaluation
- Center point evaluation: $\phi_j(0,0,0)$

### Element Geometry
- 8-vertex hexahedral elements for geometry definition
- Trilinear isoparametric mapping for coordinate transformation
- Element centers computed as average of vertex coordinates

### Error Handling
- Array dimension validation
- Element number bounds checking
- Reference coordinate range warnings

## Validation Test Cases

### 1. Uniform Flow
- **Setup**: Constant velocity field $\mathbf{u} = (u_0, 0, 0)$
- **Expected**: All element centers should have velocity $(u_0, 0, 0)$
- **Tolerance**: Machine precision

### 2. Linear Velocity Profile
- **Setup**: Couette flow $u(y) = \gamma y$
- **Expected**: Element center velocities match analytical profile
- **Tolerance**: Q2 interpolation accuracy

### 3. Quadratic Profile
- **Setup**: Poiseuille flow $u(y) = u_{max}(1 - (y/h)^2)$
- **Expected**: Exact representation by Q2 elements
- **Tolerance**: Machine precision

## Extensions and Future Development

### Enhancements
1. **Robust Element Search** - For arbitrary point evaluation
2. **Higher-Order Geometry** - Q2 isoparametric mapping
3. **Gradient Evaluation** - Velocity derivatives at element centers
4. **Parallel I/O** - Efficient VTK output for large meshes

### Additional Functionality
1. **Interpolation to Regular Grids** - For uniform data analysis
2. **Particle Tracking** - Using element-centered velocity data
3. **Flux Calculations** - Through element faces
4. **Vorticity Computation** - From velocity gradients

## File Structure
```
FeatFloWer/
├── evaluate_velocity_midpoint.f90      # Core implementation
├── velocity_midpoint_interface.f90     # Easy-to-use interfaces
└── docs/md_docs/
    └── velocity_midpoint_evaluation.md # This documentation
```

## References
- FeatFloWer Q2 finite element documentation
- Element routines E013/E013A (`source/Elements/e013.f`)
- VTK file format specification
- ParaView user documentation for visualization