# FEM Machinery Notes

This document collects recurring implementation details for the finite element building blocks that appear across FeatFloWer. It focuses on the Q2 hexahedral element implementation `E013`, because most FBM/FEM coupling routines rely on it through the external `ELE` symbol.

## Element Evaluation (`E013`)

* Location: `source/Elements/e013.f`
* Element: tri-quadratic (27-node) hexahedron used for the Q2 velocity field.
* Interface: `SUBROUTINE E013(XI1, XI2, XI3, IPAR)`
  * `XI1`, `XI2`, `XI3` – reference coordinates in `[-1, 1]` where the element should be evaluated.
  * `IPAR` selects the action:
    * `0` – compute values/derivatives at the supplied `(XI1, XI2, XI3)`.
    * `-1` – return the element number (`13`) so the caller can query the correct `NDFL` value.
    * `-2` – precompute basis data on all cubature points and store them in the internal `DHELP` cache.
    * `-3` – reuse the cached cubature entry referenced by `/CUB/` to avoid recomputation.

The routine stores its results in `COMMON /ELEM/`, so the caller must ensure that the Jacobian (`DJAC`), determinant (`DETJ`), vertex coordinates (`DX/DY/DZ`), and vertex indices (`KVE`) already describe the physical element before calling.

## `BDER` Switches

`BDER(1:10)` is the flag array inside `/ELEM/` that tells `E013` which quantities to assemble:

| Index | Meaning                            | Stored in `DBAS(1, :, index)`    |
|-------|------------------------------------|----------------------------------|
| 1     | Basis function value               | `\varphi_i`                      |
| 2     | Physical derivative in x-direction | `\partial \varphi_i / \partial x` |
| 3     | Physical derivative in y-direction | `\partial \varphi_i / \partial y` |
| 4     | Physical derivative in z-direction | `\partial \varphi_i / \partial z` |
| 5–10  | Not supported for this element     | Trigger error `IER = -131`       |

Best practice: set only the entries you need to `.TRUE.` before calling `ELE`, because the routine checks the flags and skips unnecessary work. Most velocity interpolation tasks simply enable `BDER(1)` when only function values are required.

## Mapping Reference to Physical Space

`E013` calls the helper `E013A` to obtain reference-space values and first derivatives with respect to `\xi_1`, `\xi_2`, `\xi_3`. The helper stores the raw polynomials in the temporary array `DHELP(NNBAS, 4, *)`. Afterward, `E013` multiplies the derivative entries by the inverse Jacobian constructed from `DJAC` and `DETJ` and writes the transformed quantities into `DBAS`:

```fortran
DBAS(1, idfl, 2) = (1/DETJ) * ( ... linear combination of DJAC ... )
```

Because this transformation happens inside `E013`, callers receive derivatives already expressed in physical coordinates and can directly use them for gradients, stresses, or post-processing.

## Cached Cubature Support

When loops repeatedly evaluate the element at Gauss points, the recommended sequence is:

1. Initialize the cubature rule (`CB3H`, etc.) so `/CUB/` contains `DXI` (reference points) and `DOMEGA` (weights).
2. Call `E013(0D0, 0D0, 0D0, -2)` once to populate `DHELP(:, :, icub)` for every cubature index.
3. During integration, set `ICUBP` in `/CUB/` and call `E013(0D0, 0D0, 0D0, -3)` to reuse the precomputed entry instead of recomputing basis polynomials.

This pattern explains why flux computations toggle `IPAR = -3` inside element loops when visiting each Gauss point.

## Practical Checklist Before Calling `ELE`

1. Run `SETLEV` so `mg_mesh%level(ilev)` is active and mesh data are loaded.
2. Use `NDFGL` to fill `KDFG`, `KDFL`, and set `IDFL = NDFL(IELTYP)`.
3. Copy the hexahedron vertex coordinates into `DX/DY/DZ` and vertex indices into `KVE`.
4. Build the Jacobian `DJAC` and determinant `DETJ` for the current reference point.
5. Configure `BDER` according to the required quantities.
6. Call `ELE(xi1, xi2, xi3, ipar)` and read results from `DBAS`.

### Common Cubature Choices

`CB3H` defines several tensor-product rules on the unit cube. In FBM diagnostics the default is `ICUB = 9`, which corresponds to the 3×3×3 Gauss rule:

- `NCUBP = 27` points located at `\xi = \{ -0.7745966692, 0, 0.7745966692 \}` along each axis (tensor product).
- `DXI(icubp,1:3)` stores the reference coordinates, `DOMEGA(icubp)` the product of the 1D weights.
- Once `CB3H(9)` has been called, `/CUB/` carries those coordinates and weights for all later element loops.

Pick other `ICUB` values when lower-order integration suffices (e.g., `ICUB=7` for 2×2×2 Gauss). The same caching steps apply regardless of the selected rule.

### Local/Global DoF Mapping (`NDFGL`)

`NDFGL` (`extern/libraries/feat3d/src/ndfgl.f`) provides the consistent mapping between global degrees of freedom and the local numbering expected by `DBAS`:

- Inputs: element index `IEL`, element type `IELTYP` (e.g., 13 for Q2), connectivity arrays `KVERT`, `KEDGE`, `KAREA`, and the desired sorting behavior via `IPAR`.
- Outputs: `KDFG(:)` filled with the global DoF indices that belong to the current element; if `IPAR=1`, `KDFL(:)` records the associated local indices.
- For Q2 elements (`IELTYP=13`), `NDFGL` first collects the 8 vertex nodes, then the 12 edge midpoints, then the 6 face midpoints, optionally sorts them (`NGLS`), and finally appends the element-center DoF (`NVT+NET+NAT+IEL`) to reach 27 entries.
- The combination `KDFG(il)`/`KDFL(il)` is exactly what the `ELE` routines expect when looping over `IDFL` shape functions.

When new elements are introduced, `IELTYP` must map to the correct branch in `NDFGL` so all higher-level operators (force integration, Reynolds sampling, etc.) continue to receive consistent local/global mappings.

### Particle Indicator Arrays (`FictKNPR`)

`QuadScalar_FictKnpr` (`source/src_quadLS/QuadSc_boundary.f90`) fills two related arrays that FBM diagnostics consume:

- `FictKNPR(:)` — integer flag per Q2 DoF, `0` for fluid nodes and `1` for solid nodes. This is the simple “alpha” field used by routines that only need to distinguish inside/outside.
- `FictKNPR_uint64(:,:)` — byte representation of the PE particle `SystemID` (unsigned 64-bit). The helper `longIdMatch(ig, particle%bytes)` reconstructs the ID in Fortran so individual particles can be identified even when the physics engine distributes IDs across domains.

When classifying interface elements, prefer `longIdMatch` to verify whether a DoF belongs to a specific particle (see `source/src_quadLS/QuadSc_force_extension.f90:205-214`). This avoids ambiguity when multiple particles share the same fluid indicator value.

### Pressure DoFs and Scaling Factors

- **Pressure representation:** `source/src_quadLS/QuadSc_force.f90:420-445` documents the discontinuous P1 layout used in Q2/P1 setups. Each element stores four coefficients—constant plus three linear slopes relative to the element centroid—so the pressure at any Gauss point `(XX,YY,ZZ)` is reconstructed as

  ```fortran
  JJ = 4*(IEL-1) + 1
  Press = P(JJ) + (XX-DX0)*P(JJ+1) + (YY-DY0)*P(JJ+2) + (ZZ-DZ0)*P(JJ+3)
  ```

  The field is element-local and intentionally discontinuous across boundaries.

- **Force scaling factors:** The `factors(1:6)` array passed into `GetForces` multiplies the accumulated hydrodynamic forces/torques before storing them in `myFBM%ParticleNew(IP)%ResistanceForce/TorqueForce` (`source/src_quadLS/QuadSc_force.f90:561-582`). This legacy hook compensates for symmetry planes: e.g., with a half-domain symmetric in x, set `factors(1) = 2.0` to double the computed x-force. Leave entries at 1.0 when no symmetry scaling is needed.

Once these steps are satisfied, higher-level routines (Reynolds number samplers, hydrodynamic force integrators, etc.) can rely on `DBAS` to provide consistent Q2 values and derivatives.
