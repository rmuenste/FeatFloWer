# Di Felice drag: two conventions in use (found 2026-10-14)

## The forms

- FeatFloWer `ELDragModel = difelice` (campaign form, `EL_DRAG_CLOSURE`):
  F = ½ ρ C_D(Re) A |u_r|² ε^(−χ), Re = ρ ε d |u_r|/μ, u_r the interstitial
  slip from the kernel sample; Stokes limit 2.88 π μ d |u_r| ε^(−χ−1).
- Standard form (Kafui et al. 2002, Zhu et al. 2007 review, CFDEM
  `DiFeliceDrag.C`: `pow(voidfraction, 2−Xi)`): the same times ε², i.e.
  ε^(2−χ); Stokes limit 2.88 π μ d |u_r| ε^(1−χ). Available since commit
  40d33611 as `ELDragModel = difelice_std` (opt-in; default unchanged and
  bitwise, twin job 149716).

Ratio campaign/standard = ε^(−2): +5% at ε = 0.975, +16% at 0.927, +56% at 0.8.

## Origin (liggghts-public-34, CFDEM survey §14.4, from the Di Felice 1994 abstract; full text not checked)

Di Felice defines F = F₀ ε^(−χ) with F₀ "the drag force on an unhindered
particle subject to the same volumetric flux of fluid", i.e. F₀ is evaluated
at the SUPERFICIAL slip U = ε u_r (and Re = ρ d U/μ, which equals our Re).
Rewriting with the interstitial slip gives |U|² = ε² |u_r|², hence the ε².
The campaign form keeps the superficial Re but the interstitial slip in F:
it mixes the two definitions. The code comment that justified ε^(−χ) as
"the original Di Felice form" was therefore mistaken.

Consistency check (finite-Re bookkeeping, C_D ~ Re^(−a)): with the
pressure share and the return flow present, the standard form gives
n = χ + 1 at Re → 0 (4.7) and (χ+1)/2 at high Re (2.35), i.e. it
reproduces Richardson–Zaki's 4.65 / 2.39 — presumably how χ(Re) was
calibrated. The campaign form gives 6.7 and 3.96, overshooting RZ.

| exponent n at Re ≈ 80 (χ 3.10, a 0.46) | both | pressure only | return flow only | neither |
|---|---:|---:|---:|---:|
| campaign form | 3.96 | 2.96 | 3.31 | 2.31 |
| standard form | 2.66 | 1.66 | 2.01 | 1.01 |

Measured column (campaign form, neither): 1.72; DNS 4.57; unbounded RZ 3.09.

## Consequences

- Every campaign result at ε < 1 is convention-sensitive: V1b (ε_eff 0.975:
  the standard form lowers the drag 5%, u_t rises ~3% at Re 30–80, which
  removes part of the "self-voidage bias"), V3 slopes, V4 pipe, kroupa_*
  (ε 0.7–0.95: tens of percent on the drag). None of these has been re-run.
- Prediction for the standard-form column u_t (job 149717): 0.357 → ≈ 0.369
  (ε_seen 0.988 at N→1 … the lone-particle ε_eff 0.972 gives drag ×0.945,
  u ∝ F^(1/(2−a)) → +3.7%), close to the Dallavalle fixed point 0.369.
- Recommendation: make `difelice_std` the reference form for all new work
  and re-run the rheology cases under it on a clean worktree, recording the
  deltas as datasheet rows (never retune). The campaign form stays as
  `difelice` for bitwise reproduction of the record.

Sources: Di Felice (1994) Int. J. Multiphase Flow 20, 153 (abstract);
CFDEM forum thread "Calculation of Di Felice's drag force / source code";
LIGGGHTS-PUBLIC md_docs/CFDEM_SURVEY.md §14.4 (commit 619f6788).
