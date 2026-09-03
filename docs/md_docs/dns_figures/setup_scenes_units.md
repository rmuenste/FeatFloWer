# Unit systems of the setup-scene renders

Minimal parameter sets behind `d11_hasimoto_setup.pvsm`,
`d31_beetstra_p020_setup.pvsm`, and `d31_beetstra_p010_re30_setup.pvsm`,
and the derivation of the dimensionless readings. All three runs are
nondimensional by construction: the fields in the state files carry
**code units**, and only the dimensionless groups below are physical.

## Common frame

Periodic unit cell, gravity off, flow driven by a constant body force in z:

| quantity | symbol | value (all scenes) |
|---|---|---|
| cell edge | L | 1 |
| fluid density | ρ | 1 |
| cell volume | V = L³ | 1 |
| driving | f (body force / volume, +z) | per scene below |

With ρ = 1 the dynamic and kinematic viscosities coincide numerically
(μ = ρν = ν). The natural velocity unit is ν/L, so a raw velocity value
u from the state file has the dimensionless form

    û = u · L / ν        (a Reynolds number based on the cell edge)

In steady state the total drag on the solid phase balances the driving:
F_total = f · V_cell (Hasimoto's own convention, eq. 2.14), and the
comparison velocity is the **superficial** (whole-cell average) velocity
U_sup — column 2 of each run's `bulk_flow.log` (Hasimoto eqs. 2.5/4.8;
Beetstra eq. 1).

## Scene 1 — Hasimoto cell (`d11_l4`)

| quantity | value |
|---|---|
| sphere radius | a = 1/6 (centered, D/h = 24 at L4) |
| viscosity | ν = μ = 1 |
| forcing | f = 10⁻² |
| measured U_sup (t = 4) | 1.779 × 10⁻³ |

Derived:

- Velocity numbers on the colorbar ≈ their Reynolds number:
  Re_L = U·L/ν ≈ 2 × 10⁻³ → deep Stokes flow.
- Drag coefficient (the deliverable):
  K = F / (6π μ a U_sup) = 10⁻² / (6π · 1 · (1/6) · 1.779×10⁻³) ≈ **1.79**
  vs. Hasimoto's analytic K = 1.8322 at this φ = (4/3)π a³ ≈ 0.0194
  (the ≈ −2.3 % residual is the documented O(h) interface widening at
  D/h = 24; see datasheet rows d11_l5_postfix / d11_rh_collapse).

## Scene 2 — Beetstra Stokes array (`d31_p020_s1`)

| quantity | value |
|---|---|
| spheres | N = 16 parents, r = 0.143971 (d = 0.287941), RSA seed 1 |
| solid fraction | φ = N·(4/3)πr³ / V = 0.200 |
| viscosity | ν = μ = 1 |
| forcing | f = 10⁻² |
| measured U_sup (t = 4) | 4.479 × 10⁻⁵ |

Derived (Beetstra 2007 conventions, datasheet row d31_convention_pinned):

- Re = ρ U_sup d / μ ≈ 1.3 × 10⁻⁵ → Stokes regime.
- Per-sphere fluid force: F_total/N = f·V/N = 6.25 × 10⁻⁴.
- Drag part (their force split, p. 490): F_d = (1−φ)·F_total/N = 5.0 × 10⁻⁴.
- Normalized drag (the deliverable):
  F* = F_d / (3π μ d U_sup) = 5.0×10⁻⁴ / (3π · 0.287941 · 4.479×10⁻⁵)
  ≈ **4.11** (this seed; ensemble 4.37 ± 0.25 vs. Beetstra eq. 6 = 4.19).

## Scene 3 — Beetstra finite-Re array (`d31_p010_s1_re30`)

| quantity | value |
|---|---|
| spheres | N = 16 parents, r = 0.114270 (d = 0.228539), RSA seed 1 |
| solid fraction | φ = 0.100 |
| viscosity | ν = μ = 10⁻² (the knob that buys finite Re at O(1) velocity) |
| forcing | f = 1.75074 |
| measured U_sup (t = 4) | 1.170 |

Derived:

- Re = ρ U_sup d / μ = 1.170 · 0.228539 / 0.01 ≈ **26.7** (the "Re ≈ 27"
  tier of the datasheet, Beetstra eq. 2 — superficial velocity).
- F_d = (1−φ)·f·V/N = 0.9 · 0.109421 = 9.848 × 10⁻².
- F* = F_d / (3π μ d U_sup) ≈ **3.91** (this seed; ensemble at this tier
  3.63 ± 0.39 — moderate-Re microstructure dispersion, row
  d31_corner_discriminator).

Because ν differs by 100× between scenes 2 and 3, their colorbar values
are **not on a common axis**; normalize each field by its own U_sup to
compare the pictures side by side.

## Mapping to SI, if ever needed

Pick a target fluid (ν_SI) and cell size (L_SI); then any code velocity
converts as

    u_SI = (u · L/ν)_code · ν_SI / L_SI

All dimensionless deliverables (K, F*, Re, φ) are unchanged by
construction. Example: water (ν_SI = 10⁻⁶ m²/s) in a 1 cm cell maps the
Hasimoto U_sup to 1.78×10⁻³ · 10⁻⁶/10⁻² ≈ 0.18 µm/s.

## Notes

- `particles.xyz` in the array run directories lists parents **plus
  periodic image spheres** (34 and 30 entries respectively); only the
  16 parents enter φ and the force normalization. Image spheres exist
  because the FBM indicator is not periodic (row d31_periodic_indicator).
- Solver output fields carry no unit metadata; `Pressure_V` and
  `Velocity` are code units per the tables above.
