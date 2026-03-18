# Cylinder Attraction Smoother

## Overview

The cylinder attraction smoother is a mesh node redistribution technique for the `q2p1_fac3d` application (DFG flow-around-cylinder benchmark). It concentrates mesh nodes near the cylinder surface to improve the boundary layer resolution of the Fictitious Boundary Method (FBM).

The smoother operates during mesh initialization in `app_init.f90` and is activated by the parameter `bFAC3D_CylUmbrellaWeight`.

## Algorithm

The attraction works in two phases:

### Phase 1: Radial attraction toward the cylinder surface

Each mesh node is classified by its signed distance `dist` to the cylinder barrel surface (negative = inside, positive = outside). Nodes near the caps are skipped — only the radial (barrel) direction is treated.

For each node, the nearest point on the cylinder barrel is computed by radial projection:

```
nearX = CylCenter(1) + CylRadius * dx / rxy
nearY = CylCenter(2) + CylRadius * dy / rxy
```

The node is then moved toward this nearest point:

```
x_new = x_old + omega * alpha * (nearX - x_old)
y_new = y_old + omega * alpha * (nearY - y_old)
```

where `omega = 0.2` is a global relaxation factor, and `alpha` is a distance-dependent attraction weight.

#### Exterior nodes (dist > 0)

Two distance bands control the attraction strength:

| Band | Range | Alpha profile |
|------|-------|---------------|
| Strong | `[0, dBand1)` = `[0, 0.04)` | `0.6 + 0.4 * dist/dBand1` (linear, 0.6 to 1.0) |
| Fade-out | `[dBand1, dBand2]` = `[0.04, 0.15]` | `((dBand2 - dist) / (dBand2 - dBand1))^2` (quadratic decay) |
| None | `> dBand2` = `> 0.15` | 0 (no movement) |

The strong band pulls nodes close to the surface firmly. The quadratic fade-out prevents abrupt transitions that would create enlarged elements at the boundary of the attraction zone.

#### Interior nodes (dist < 0)

Interior nodes receive a gentler outward push toward the surface:

| Band | Range | Alpha profile |
|------|-------|---------------|
| Push | `[0, dBandIn)` = `[0, 0.04)` | `dAlphaIn * (1 - |dist|/dBandIn)` with `dAlphaIn = 0.45` (linear fade) |
| None | `> dBandIn` | 0 (no movement) |

The interior attraction is weaker than the exterior to avoid mesh destruction inside the cylinder body. The effective maximum movement factor is `dAlphaIn * omega = 0.45 * 0.2 = 0.09` per iteration.

### Phase 2: Post-attraction Laplacian (umbrella) smoothing

After 8 attraction iterations, 2 iterations of Laplacian umbrella smoothing (`UmbrellaSmoother_STRCT`) are applied to heal the transition zone between attracted and unaffected regions. This prevents the enlarged elements that would otherwise appear at the fade-out boundary.

The number of post-attraction smoothing iterations was tuned experimentally: too many iterations (e.g. 10) undo the attraction effect, too few leave visible transition artifacts.

## Code flow in app_init.f90

```
1. Cylinder attraction (8 iterations, omega=0.2)
   - CylinderAttraction(nvt, dcorvg, 0.2)
   - ParametrizeBndryPoints_STRCT (re-project boundary nodes after each iteration)

2. Post-attraction umbrella smoothing (2 iterations)
   - UmbrellaSmoother_STRCT(0.0, 1)

3. ProlongateCoordinates (propagate fine-level coords to extended level)

4. Final projection to NLMAX+1
   - ParametrizeBndryPoints_STRCT(mg_mesh, nlmax+1)

5. VTK debug output
   - Writes initial_mesh.pvtu with CylinderDistance and AttractionWeight fields
```

## Key subroutines

All three subroutines are defined at the end of `applications/q2p1_fac3d/app_init.f90`:

- **`CylinderAttraction(nVtx, dcorvg, dOmega)`** — moves nodes radially toward the barrel surface
- **`ComputeCylDist(nVtx, dcoor, dist)`** — computes signed distance to the cylinder surface
- **`ComputeCylAttractionWeight(nVtx, dist, weight)`** — computes the attraction alpha for VTK visualization (mirrors the alpha logic in `CylinderAttraction`; negative values indicate interior nodes)

The cylinder geometry parameters (`dFAC3D_CylCenter`, `dFAC3D_CylRadius`, `dFAC3D_CylLength`) are read from `var_QuadScalar`.

## VTK debug output

The initial mesh is written to `_vtk/initial_mesh.pvtu` with two scalar fields:

- **CylinderDistance**: signed distance to the cylinder surface (negative = inside)
- **AttractionWeight**: the alpha value used for attraction (negative = interior push)

The VTK output routines (`Output_VTK_mesh_piece_with_2fields`, `Output_VTK_mesh_main_with_2fields`) are defined in `source/OutputProfiles.f90`.

## Known issues and fixes

### ILEV global variable corruption by umbrella smoothers

#### Problem

The umbrella smoother routines modify the global/common variable `ILEV` (from `def_FEAT`) as a side effect:

- `UmbrellaSmoother_STRCT` (in `source/Umbrella.f90`) sets `ILEV = NLMAX` internally
- `us_UmbrellaSmoother` (in `source/src_mesh/umbrella_smoother.f90`) sets `ILEV = mgMesh%nlmin` on exit

The VTK output routines in `OutputProfiles.f90` use `KNVT(ILEV)` and `KNEL(ILEV)` to determine array sizes for the VTK file header. If `ILEV` does not match the level of the data being written, the VTK files contain incorrect array sizes and ParaView crashes when loading them.

In the specific case of `q2p1_fac3d`: the initialization DO loops leave `ILEV = NLMAX + 1` (Fortran DO loop exit semantics). The VTK output writes data from `mg_mesh%level(mg_mesh%maxlevel)`, which corresponds to `KNVT(NLMAX + 1)`. When the post-attraction umbrella smoothing runs between the initialization loops and the VTK output, `ILEV` gets changed to `NLMAX`, causing a size mismatch.

#### Fix applied

Both umbrella smoother routines now save `ILEV` on entry and restore it on exit:

```fortran
! In UmbrellaSmoother_STRCT (source/Umbrella.f90):
integer :: ILEV_save
...
ILEV_save = ILEV    ! save before any level changes
...
ILEV = ILEV_save    ! restore on exit

! Same pattern in us_UmbrellaSmoother (source/src_mesh/umbrella_smoother.f90)
```

**Important:** The restore must NOT call `SETLEV(2)` because the saved `ILEV` value may be outside the valid range `[NLMIN, NLMAX]` that `SETLEV` enforces (e.g. `NLMAX + 1` from a DO loop exit). The `ILEV` variable is only restored for array indexing purposes (`KNVT(ILEV)`, `KNEL(ILEV)`) — the internal state managed by `SETLEV` was already correct before the smoother was called and does not need to be re-established.

#### Status

This fix resolves the observed VTK corruption. However, the broader issue — that `ILEV` is a global variable used both as a loop counter and as implicit state by output routines — remains a latent source of bugs. Any routine that modifies `ILEV` without restoring it can silently corrupt downstream code. When working with umbrella smoothers or other routines that call `SETLEV`, verify that `ILEV` has the expected value before any code that depends on it.

## Parameter tuning history

The attraction parameters were tuned iteratively through visual inspection in ParaView:

1. **Initial attempt** (30 iterations, omega=0.4): too aggressive, nodes collapsed onto the surface
2. **Reduced** (8 iterations, omega=0.2): mesh no longer destroyed
3. **Alpha profile redesign**: changed from uniform alpha to distance-banded profile with `0.6 + 0.4*dist/dBand1` near surface and quadratic fade further out
4. **Added post-attraction smoothing**: 2 umbrella iterations to heal the transition zone (10 was too many, 5 still too many)
5. **Interior attraction**: added with conservative parameters (`dAlphaIn=0.45`, `dBandIn=0.04`) after initial attempt with `ABS(dist)` destroyed the mesh
