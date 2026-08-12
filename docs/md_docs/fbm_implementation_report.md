# Fictitious Boundary Method Implementation Report

This note maps the fictitious boundary method (FBM) concepts from
`fbm_fictitious_boundary_method.pdf` to the current implementation.

The short answer is that the codebase contains a real FBM implementation in
the sense of the Turek/Wan/Rivkind paper: it keeps a fixed Eulerian finite
element mesh, repeatedly classifies mesh degrees of freedom as fluid or
fictitious solid, and enforces internal Dirichlet velocity constraints by
filtering solution, defect, and matrix data during the iterative solve.

The alpha-field hydrodynamic force path is also real, but it is a
resolved-particle extension around the same FBM geometry machinery. It computes
forces from a volume integral using the gradient of a discrete 0/1 particle
indicator. That force-evaluation method is not the central method described in
`fbm_fictitious_boundary_method.pdf`, whose focus is implicit treatment of
Dirichlet conditions by filtering.

## Implementation Map

The central runtime path is:

1. `Transport_q2p1_UxyzP_fc_ext`
   - `source/src_quadLS/QuadSc_main.f90:332`

2. Geometry/alpha update before the fluid solve
   - `source/src_quadLS/QuadSc_main.f90:346`

3. Momentum solve with boundary/FBM filtering
   - `source/src_quadLS/QuadSc_main.f90:364`
   - `source/src_quadLS/QuadSc_boundary.f90:538`

4. Pressure correction
   - `source/src_quadLS/QuadSc_main.f90:512`

5. Hydrodynamic force computation
   - `source/src_quadLS/QuadSc_main.f90:572`
   - `source/src_fbm/fbm_main.f90:273`

6. Particle/PE update
   - `source/src_quadLS/QuadSc_main.f90:599`
   - `source/src_fbm/fbm_main.f90:1705`

## Alpha Field

The alpha field is implemented as `FictKNPR`, plus a 64-bit identity companion
in the PE path:

- `FictKNPR(i) != 0` means Q2 degree of freedom `i` lies inside a fictitious
  solid region.
- In PE mode, `FictKNPR(i)` is effectively solid/fluid, while
  `FictKNPR_uint64(i)` stores which particle owns the DOF.
- The point classification is performed for vertices, edge midpoints, face
  centers, and element centers in `QuadScalar_FictKnpr`.
  - `source/src_quadLS/QuadSc_boundary.f90:60`

For PE particles, the geometric point-in-solid test calls C++ PE geometry:

- Fortran wrapper:
  - `source/src_fbm/fbm_main.f90:550`
- C++ point containment:
  - `libs/pe/src/interface/object_queries.cpp:404`
- Accelerated hash-grid variant:
  - `libs/pe/src/interface/object_queries.cpp:539`

So alpha is not solved as a PDE. It is rebuilt geometrically every step by
asking whether each finite element DOF is inside a particle or body.

## Implicit Treatment

In the paper, "implicit treatment" means the fictitious boundary is handled as
an internal Dirichlet constraint inside the linear/nonlinear solve, instead of
by adding an explicit penalty force, body force, or viscosity-density blockage.

In this code, the fictitious boundary DOFs are marked in `FictKNPR`, then
treated like constrained velocity DOFs:

- `FictKNPR` is built by geometric point-in-body queries in
  `QuadScalar_FictKnpr`.
  - `source/src_quadLS/QuadSc_boundary.f90:60`
- The actual point classifier is called for Q2 vertices, edges, faces, and
  element centers.
  - `source/src_quadLS/QuadSc_boundary.f90:108`
- The matrix rows/couplings for fictitious DOFs are filtered in
  `Boundary_QuadScalar_Mat`.
  - `source/src_quadLS/QuadSc_boundary.f90:668`
- For the coupled 9-block non-Newtonian/stress matrix, the same idea is in
  `Boundary_QuadScalar_Mat_9`.
  - `source/src_quadLS/QuadSc_boundary.f90:724`

The important block is:

- `source/src_quadLS/QuadSc_boundary.f90:695`

If `FictKNPR(I) != 0`, off-diagonal matrix entries for that DOF are zeroed. In
the 9-block version, cross-component couplings are also zeroed:

- `source/src_quadLS/QuadSc_boundary.f90:762`

That is the implicit Dirichlet mechanism: fictitious boundary conditions are
embedded into the algebraic solve.

## Iterative Filtering

The iterative filtering shows up as repeated calls that project velocity and
defect data back onto the admissible boundary-condition subspace during the
solve.

The main flow is in `Transport_q2p1_UxyzP_fc_ext`:

- Before solving, geometry/alpha is refreshed.
  - `source/src_quadLS/QuadSc_main.f90:346`
- The RHS/defect is filtered.
  - `source/src_quadLS/QuadSc_main.f90:385`
- The current velocity vector is filtered.
  - `source/src_quadLS/QuadSc_main.f90:393`
- After assembling the implicit matrix/defect, the defect is filtered again.
  - `source/src_quadLS/QuadSc_main.f90:405`
- Inside the nonlinear iteration, after each solve/reassembly cycle, the defect
  is filtered again.
  - `source/src_quadLS/QuadSc_main.f90:459`

The two filtering operators are:

- Defect filter:
  - `source/src_quadLS/QuadSc_boundary.f90:538`
- For fictitious DOFs:
  - `source/src_quadLS/QuadSc_boundary.f90:548`
- Value filter:
  - `source/src_quadLS/QuadSc_boundary.f90:578`
- For fictitious DOFs it calls `fbm_velBCUpdate`, which imposes the internal
  boundary velocity:
  - `source/src_quadLS/QuadSc_boundary.f90:604`

For moving particles, that imposed value is rigid-body motion,
`U + omega x (x - X)`:

- `source/src_fbm/fbm_main.f90:1188`

## Inside The Iterative Solvers

The multigrid solver receives the constrained-DOF masks:

- `source/src_quadLS/QuadSc_def.f90:2767`

The smoother then skips or specially treats constrained DOFs through `KNPR`.
For example, the scalar SOR smoother is:

- `source/src_quadLS/QuadSc_solver.f:301`

It only performs the normal SOR update when `KNPR(IEQ) == 0`:

- `source/src_quadLS/QuadSc_solver.f:315`

The vector velocity smoother does the same for coupled velocity components:

- `source/src_quadLS/QuadSc_solver.f:23`

This is the solver-side realization of the paper's iterative filtering idea:
the constrained subspace is preserved through smoothing and multigrid cycles.

## Dirichlet FBM Enforcement

This is the part closest to `fbm_fictitious_boundary_method.pdf`.

The paper's key idea is not a smeared material model, but an internal-boundary
filtering method: DOFs belonging to fictitious boundaries are treated like
constrained Dirichlet values during iterative solution.

That is exactly what happens here:

- `Boundary_QuadScalar_Def` zeroes the defect/residual on fictitious-body DOFs.
  - `source/src_quadLS/QuadSc_boundary.f90:548`
- `Boundary_QuadScalar_Val` overwrites velocity values at fictitious-body DOFs.
  - `source/src_quadLS/QuadSc_boundary.f90:602`
- The imposed value is rigid-body velocity:
  `U + omega x (x - X)`.
  - `source/src_fbm/fbm_main.f90:1188`
- The PE-aware version retrieves local/remote particle state and matches by
  64-bit system ID.
  - `source/src_fbm/fbm_main.f90:1232`

That is a faithful fixed-grid fictitious-boundary implementation: the particle
surface is not mesh-fitted, but the velocity inside the particle region is
constrained to the body motion.

## Hydrodynamic Force Computation

The alpha-force implementation has two main variants.

Legacy/local particle-index version:

- `source/src_quadLS/QuadSc_force.f90:2`

PE-aware version:

- `source/src_quadLS/QuadSc_force_extension.f90:874`

The algorithm is:

1. Loop over particles.
2. Loop over Q2 hexahedral elements.
3. Skip elements where all 27 Q2 DOFs are inside or all are outside the current
   particle.
   - `source/src_quadLS/QuadSc_force.f90:103`
4. At Gauss points, interpolate:
   - velocity and velocity gradients,
   - pressure,
   - alpha and alpha gradients.
5. Build `-grad(alpha)`.
   - `source/src_quadLS/QuadSc_force.f90:480`
6. Compute Newtonian Cauchy stress action, `sigma * (-grad alpha)`, with the
   symmetric viscous gradient.
   - `source/src_quadLS/QuadSc_force.f90:492`
7. Integrate force and torque.
   - `source/src_quadLS/QuadSc_force.f90:500`
8. Send forces to PE.
   - `source/src_quadLS/QuadSc_force_extension.f90:454`
   - `source/src_particles/dem_query.f90:416`
   - `libs/pe/src/interface/object_queries.cpp:1165`

This matches the formula documented in
`docs/md_docs/hydrodynamic_force_computation.md`:

```text
F_i = - integral sigma · grad(alpha_i) dV
```

## Why This Is Called A Sharp Interface Method

It is called sharp because the body/fluid distinction is binary at the
discrete DOF level:

- inside body: alpha/`FictKNPR` is nonzero,
- outside body: alpha/`FictKNPR` is zero.

There is no smooth phase-field transition, no diffuse color-function
transport, and no viscosity-density ramp across several cells. The solid
boundary is detected by exact geometric containment queries, then imposed as
an internal Dirichlet boundary.

That said, it is not body-fitted. The physical surface can cut through mesh
cells, and the finite element representation sees the interface through the
Q2 interpolation of a 0/1 nodal indicator. So "sharp" here means sharp
material/constraint classification on a fixed grid, not exact
boundary-conforming quadrature.

## Assessment

This is a genuine FBM implementation for fixed-grid incompressible flow with
moving internal boundaries. The strongest evidence is the repeated pattern:

- rebuild fictitious geometry/alpha,
- filter velocity values and defects at internal-body DOFs,
- solve Navier-Stokes on the unchanged mesh,
- compute resolved hydrodynamic force from stress and `grad(alpha)`,
- update particle dynamics,
- rebuild geometry for the next step.

The main caveat is documentation wording: `fbm_fictitious_boundary_method.pdf`
is primarily about implicit Dirichlet boundary treatment and multigrid
filtering. The alpha-gradient force formula in
`hydrodynamic_force_computation.md` belongs to the resolved-particle
force-evaluation layer built on top of that FBM machinery. The docs should
explicitly distinguish those two layers when cross-referencing the paper and
implementation.

