# ten Cate (2002) E1–E4 reference trajectories

Digitized experimental (PIV) curves from A. ten Cate et al., *Phys.
Fluids* 14 (2002) 4012–4025 (`ten_cate_piv.pdf`, repo root), provided by
R. Münster 2026-07-31. This restores the `tc-ref/` fileset referenced by
`pipemesh_v1/handoff_euler_lagrange_drag_validation_ten_cate.md` and is
the physics reference for the DNS campaign's D0.4/D1 gates
(`dns-validation-campaign-plan.md`).

## Files and format

Two space-separated columns, no header; irregular time stamps (digitized
from the paper's figures).

- `case_E<n>_h.csv` — column 1: time [s]; column 2: **gap height h/d**
  (sphere-bottom-to-wall distance normalized by d = 0.015 m; starts at
  (0.1275 − 0.0075)/0.015 = 8.0 at release).
- `ref_E<n>.dat` — column 1: time [s]; column 2: **settling velocity v_z
  [m/s]** (negative = downward).

## Case parameters (paper Table I)

| Case | ρ_f [kg/m³] | μ [Pa·s] | u_∞ [m/s] | Re_∞ |
|---|---|---|---|---|
| E1 | 970 | 0.373 | 0.038 | 1.5 |
| E2 | 965 | 0.212 | 0.060 | 4.1 |
| E3 | 962 | 0.113 | 0.091 | 11.6 |
| E4 | 960 | 0.058 | 0.128 | 31.9 |

Sphere: d = 15 mm, ρ_p = 1120 kg/m³, release z = 0.1275 m; box
100×100×160 mm.

## Peak settling velocities in these curves (sanity-checked 2026-07-31)

| Case | v_min [m/s] | \|v_min\|/u_∞ |
|---|---|---|
| E1 | −0.0372 | 0.979 |
| E2 | −0.0594 | 0.989 |
| E3 | −0.0866 | 0.952 |
| E4 | −0.1230 | 0.961 |

Consistent with the confined-box ratio ≈ 0.955 quoted in the EL-campaign
notes. **Gate note:** DNS comparisons must target these confined-PIV
curves, not the unbounded u_∞ — e.g. the documented serial E4 result
−0.1329 m/s is +3.8% vs u_∞ = 0.128 but **+8.1% vs the PIV peak
−0.1230**, which is the honest baseline for the D1 resolution study.
