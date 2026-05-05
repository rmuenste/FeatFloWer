# Modular attraction smoother proposal (post-`ccf1ec55f67098dfba991e11872b35d9d2ea1cc9`)

## Current state in branch

The current cylinder-attraction implementation is split across application and core mesh-smoother code:

- `applications/q2p1_fac3d/app_init.f90` hard-enables cylinder weighting and hardcodes geometry defaults (`bFAC3D_CylUmbrellaWeight`, center/radius/length), runs attraction loops with fixed iteration counts, and defines cylinder-specific routines (`ComputeCylDist`, `CylinderAttraction`, `ComputeCylAttractionWeight`).
- `source/src_mesh/umbrella_smoother.f90` also contains cylinder-specific behavior (`ComputeFAC3DCylDistance`) and branch logic in `GetWeight`/distance setup for the weighted umbrella variant.
- `source/src_quadLS/QuadSc_var.f90` stores global flags and geometry values, but only for the FAC3D cylinder case.

This means the feature is not yet reusable for other applications/geometries and has duplicated cylinder-distance logic.

## Target architecture for a switchable feature usable by all applications

## 1) Create a dedicated module for attraction smoothing

Introduce a new module (suggested location: `source/src_mesh/attraction_smoother.f90`) with three layers:

1. **Configuration/state** (`tAttractionConfig`)
2. **Distance provider interface** (application/pluggable geometry callback)
3. **Execution kernel** (generic node attraction + optional post-smoothing)

### Suggested API surface

- `as_init_from_simpar(...)`:
  - Reads/takes all attraction-related parameters
  - Validates ranges and derives defaults
- `as_register_distance_provider(proc)`:
  - Registers a callback `distance(nVtx, dcorvg, dist, meta)`
  - Callback is application-defined (cylinder, sphere, wall STL, signed-distance field, ...)
- `as_apply(mg_mesh, ilev, cfg, distance_proc)`:
  - Performs `nAttractSteps`
  - Computes distance per step
  - Maps distance to attraction strength via configurable kernel
  - Applies displacement with configurable `omega`
  - Reprojects boundary points (`ParametrizeBndryPoints_STRCT`)
- `as_post_smooth(mg_mesh, cfg)`:
  - Optional umbrella cleanup (`nPostUmbrellaSteps`)

With this API, applications only provide geometry distance behavior and parameter values; smoother mechanics remain centralized.

## 2) Distance provider abstraction

The key abstraction is distance computation plus a small provider-specific
movement kernel:

- **Current special case:** `dist = max(r_xy - radius, z_gap)` for finite cylinder barrel/caps.
- **Generalized requirement:** the smoother can ask a provider for `dist(i)` and
  then apply the matching displacement rule for that provider.

The tempting interface was to make every provider return a signed distance plus
a target projection point. That would make the movement formula generic, but it
is the wrong default contract here: if a mesh generator already knows robust
target projection points for the relevant grid nodes, it usually no longer
needs this attraction smoother.

The implemented direction is therefore pragmatic:

- Use signed distance as the shared quantity for weighting, diagnostics, and
  umbrella distance weighting.
- Keep the actual attraction displacement provider-specific. The current
  cylinder provider uses radial barrel attraction and deliberately skips cap
  projection.
- Add new providers by adding their distance/movement kernel in the attraction
  module, while keeping app initialization and umbrella smoothing independent of
  the geometry details.

## 3) Parameter model to expose via `SimPar@...`

The modular component should collect all currently hardcoded tuning knobs and geometry selection flags.

### Core enable/select parameters

- `SimPar@AttractionEnable` (Yes/No, default No)
- `SimPar@AttractionMode` (`distance_weight_only`, `node_attraction`, `both`)
- `SimPar@AttractionDistanceProvider` (e.g. `cylinder`, `sphere`, `stl`, `user`)

### Iteration/scheduling parameters

- `SimPar@AttractionSteps` (currently fixed `8` in FAC3D)
- `SimPar@PostAttractionUmbrellaSteps` (currently fixed `2`)
- `SimPar@AttractionOmega` (currently fixed `0.2`)
- `SimPar@AttractionApplyOnLevel` (`NLMAX` by default)

### Exterior/interior weighting bands (currently hardcoded)

- `SimPar@AttractionBandOuterStrong` (currently `dBand1 = 0.04`)
- `SimPar@AttractionBandOuterMax` (currently `dBand2 = 0.15`)
- `SimPar@AttractionBandInner` (currently `dBandIn = 0.04`)
- `SimPar@AttractionInnerMaxAlpha` (currently `dAlphaIn = 0.45`)
- `SimPar@AttractionOuterMinAlpha` (currently `0.6` in near band)
- `SimPar@AttractionOuterMaxAlpha` (currently `1.0` in near band)
- `SimPar@AttractionOuterFarExponent` (currently quadratic fade-out `^2`)

### Distance-weighted umbrella integration

- `SimPar@UmbrellaDistanceWeightEnable` (replaces FAC3D-specific boolean naming)
- `SimPar@UmbrellaWeightCap` (currently `25`)
- `SimPar@UmbrellaWeightPower` (currently `2.3`)
- `SimPar@UmbrellaKernelProfile` (currently piecewise `KernelFunction`)

### Geometry-specific parameter block (cylinder provider)

- `SimPar@AttractionCylinderCenter = x y z`
- `SimPar@AttractionCylinderRadius = r`
- `SimPar@AttractionCylinderLength = L`
- Optional: `SimPar@AttractionCylinderUseCaps = Yes/No`

## 4) Refactor plan (incremental, low risk)

1. **Extract cylinder routines** from `q2p1_fac3d/app_init.f90` and `umbrella_smoother.f90` into the new attraction module.
2. **Replace hardcoded constants** with config fields read from parser defaults.
3. **Introduce generic names** in `var_QuadScalar` (keep FAC3D aliases temporarily for backward compatibility).
4. **Hook parser** (`source/src_util/param_parser.f90`) to populate the new config from `SimPar@Attraction...` keys.
5. **Provide one shared call site** after initial umbrella smoothing in each application template:
   - `CALL AttractionSmoother_ApplyIfEnabled(...)`
6. **Migrate FAC3D** to new keys while preserving old cylinder keys as deprecated aliases for one transition period.

## 5) Backward compatibility strategy

To avoid breaking existing FAC3D setups:

- Keep reading old variable names/keys (if present), but print deprecation warnings.
- If both old and new keys are present, new keys win.
- Keep old defaults matching current behavior when `AttractionEnable=Yes` and `AttractionDistanceProvider=cylinder`.

## 6) Why this solves the modularization request

- **Optional/switchable:** global enable + mode flags.
- **Reusable across all applications:** shared module + one optional call in app init flow.
- **User-supplied distance:** provider callback contract.
- **Configurable tuning:** all attraction/weight constants become `SimPar` parameters instead of hardcoded literals.
