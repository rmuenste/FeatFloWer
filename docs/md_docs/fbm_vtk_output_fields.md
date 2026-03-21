# FBM Per-Cell VTK Output Fields

This document describes the quantities computed per mesh cell $T$ and written
to the `fbm_forces.pvtu` VTK output file.

## Background

The Fictitious Boundary Method (FBM) computes hydrodynamic forces on an
immersed rigid body by integrating the stress tensor against the gradient of an
indicator function $\alpha$ over the volume mesh.  The indicator field
$\alpha(\mathbf{x})$ equals 1 at vertices inside the rigid body and 0
elsewhere.  Its gradient $\nabla\alpha$ is non-zero only in cells that
intersect the body surface, so the volume integral approximates the surface
integral of the traction vector.

## Stress tensor

The discrete Cauchy stress tensor used throughout is the **symmetric** form:

$$
\sigma_h = -p\,\mathbf{I} + \mu\bigl(\nabla\mathbf{u} + (\nabla\mathbf{u})^\top\bigr)
$$

where $p$ is the pressure, $\mu$ the dynamic viscosity, and $\mathbf{u}$ the
velocity field.

> **Note:** An earlier version of the code used the non-symmetric approximation
> $\sigma_h \approx -p\,\mathbf{I} + \mu\,\nabla\mathbf{u}$, which
> underestimates drag.  The symmetric form is now active in both the global
> force summation (`GetForcesCyl`) and the per-cell output.

## Output fields

### PointData

| Field | Formula | Description |
|-------|---------|-------------|
| **Alpha** | $\alpha_i \in \{0,1\}$ | FBM indicator: 1 if vertex $i$ is inside the rigid body, 0 otherwise. Derived from `FictKNPR`. |

### CellData

All cell fields are volume integrals over the cell $T$ evaluated with
numerical cubature (3D hexahedral Gauss rule, `ICUB`).  Let
$\omega_q$ denote the cubature weight at Gauss point $q$ (including the
Jacobian determinant), and let
$\mathbf{n} = -\nabla\alpha = (-\partial_x\alpha,\;-\partial_y\alpha,\;-\partial_z\alpha)$
be the approximate inward surface normal.

#### CellForceX / CellForceY / CellForceZ

Per-cell hydrodynamic force integral, scaled by `ForceScale`:

$$
F_i(T) = s_i \sum_q \bigl[\sigma_h \cdot \mathbf{n}\bigr]_i\;\omega_q
$$

Written out component-wise at each Gauss point:

$$
f_1 = -p\,n_1 + \mu\bigl[2\,u_{1,x}\,n_1 + (u_{1,y}+u_{2,x})\,n_2 + (u_{1,z}+u_{3,x})\,n_3\bigr]
$$

$$
f_2 = -p\,n_2 + \mu\bigl[(u_{2,x}+u_{1,y})\,n_1 + 2\,u_{2,y}\,n_2 + (u_{2,z}+u_{3,y})\,n_3\bigr]
$$

$$
f_3 = -p\,n_3 + \mu\bigl[(u_{3,x}+u_{1,z})\,n_1 + (u_{3,y}+u_{2,z})\,n_2 + 2\,u_{3,z}\,n_3\bigr]
$$

where $u_{i,j} = \partial u_i / \partial x_j$ and $s_i$ is the
`Properties%ForceScale` factor for component $i$.

#### ScaledForceX / ScaledForceY / ScaledForceZ

Surface traction approximation — force divided by a characteristic surface
area:

$$
\tilde{F}_i(T) = \frac{F_i(T)}{|T|^{2/3}}
$$

where $|T| = \sum_q \omega_q$ is the cell volume.  The scaling
$|T|^{2/3}$ converts from a volume integral to a quantity with units of
force per unit area, providing a local approximation of the surface
traction on the immersed boundary.

#### CellNormalX / CellNormalY / CellNormalZ

Integrated approximate surface normal per cell:

$$
N_i(T) = \sum_q n_i\;\omega_q = -\sum_q \frac{\partial\alpha}{\partial x_i}\;\omega_q
$$

This is non-zero only in cells that intersect the rigid body surface.
The direction approximates the outward surface normal of the body.

#### SymGradX / SymGradY / SymGradZ

Viscous part of the traction without viscosity coefficient and without
pressure — the symmetric velocity gradient contracted with the normal:

$$
G_i(T) = \sum_q \bigl[(\nabla\mathbf{u}+(\nabla\mathbf{u})^\top)\cdot\mathbf{n}\bigr]_i\;\omega_q
$$

This isolates the deformation-rate contribution to the surface traction,
independent of $\mu$ and $p$.  Useful for comparing the relative
magnitudes of viscous and pressure forces.

#### VolScale

The cell-volume-based scaling factor:

$$
V_s(T) = |T|^{2/3}
$$

where $|T| = \sum_q \omega_q$.  This is the denominator used in the
ScaledForce fields.  Output separately so the user can apply alternative
scaling in post-processing.

## Relationships between fields

The full per-cell force can be reconstructed from the components:

$$
F_i(T) = s_i\Bigl[-p\;N_i(T) + \mu\;G_i(T)\Bigr]
$$

where the pressure contribution uses the CellNormal and the viscous
contribution uses the SymGrad field (before ForceScale is applied).

## File format

The output is a parallel VTK UnstructuredGrid (`.pvtu` + `.vtu` pieces).
The mesh uses Q1 hexahedral cells (VTK type 12) with 8 vertices per
element. Each MPI rank writes one `.vtu` piece file; rank 1 writes the
`.pvtu` master file referencing all pieces.
