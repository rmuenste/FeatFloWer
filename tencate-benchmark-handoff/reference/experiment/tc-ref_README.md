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

## Peak settling velocities: digitized curves vs the paper's printed ratios

The paper's Table II (p. 4018) ends with the **experimentally obtained**
u_max/u_∞ ratios — printed numbers, more authoritative than any
digitization. Audit (2026-08-01):

| Case | Table II u_max/u_∞ | → v_peak [m/s] | digitized-curve v_min | digitized ratio | curve error |
|---|---|---|---|---|---|
| E1 | 0.947 | −0.03599 | −0.0372 | 0.979 | **+3.4%** |
| E2 | 0.953 | −0.05718 | −0.0594 | 0.990 | **+3.9%** |
| E3 | 0.959 | −0.08727 | −0.0866 | 0.952 | −0.7% |
| E4 | 0.955 | −0.12224 | −0.1230 | 0.961 | +0.6% |

The digitized E3/E4 peaks are consistent with print (≲0.7%), but the
digitized **E1/E2 peaks are ~3.5–4% too fast** — also unphysical in
trend (wall retardation should strengthen, not weaken, as Re drops).

**Gate note:** peak-velocity gates use the Table II printed values above
as primary reference; the digitized curves remain the reference for
*shape* (v(t), h/d trajectories, timing) with the E1/E2 amplitude caveat.
Never gate against unbounded u_∞ — e.g. the documented serial E4 result
−0.1329 m/s is +3.8% vs u_∞ = 0.128 but +8.7% vs Table II's −0.12224.
