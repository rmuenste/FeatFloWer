# Oblate spheroid in plane Couette flow with inertia — literature digest for the DNS validation case

Sources (both in `literature/`):

- **DA00** = Ding & Aidun, *J. Fluid Mech.* 423 (2000) 317–344, `print_ding_aidun_2000.pdf` (page numbers below are journal pages; PDF page = journal page − 316).
- **DG25** = Di Giusto, Bergougnoux & Guazzelli, *J. Fluid Mech.* 1017 (2025) A41, `orientation-of-flat-bodies-of-revolution-in-shear-flows-at-low-reynolds-number.pdf` (page numbers A41-n).

Figure values were read off 200–400 dpi renders with pixel calibration against the printed axis ticks; stated uncertainties are the reading uncertainty only. Tags: `[DA00 §6 p.338]`, `[DG25 fig 10 p.A41-14]`, etc. "not stated" = the paper is silent.

Notation reminder: DA00's `Re = ρ_f G d²/μ` with `d = 2b` (major *diameter*); DG25's `Re_p = ρ_f a² γ̇/μ` with `a` the major *radius*. Hence **Re_p(DG25) = Re(DA00)/4** for the same particle and shear rate. DA00's Re_c = 81 is Re_p ≈ 20.

---

## A. Ding & Aidun 2000

### A.1 Frame, governing model and dimensionless groups (§1, §2)

- Frame: `(x', y', z')` fixed in space, shear flow `(u, v, w) = (0, 0, −G y')`, so **z' is the flow direction, y' the gradient direction, x' the vorticity direction**. `χ` is the angle from the `y'`-axis to the body `y`-axis (the `b` semi-axis). `[DA00 §1 p.318–319, fig 1 caption p.319]`
- Jeffery solution used as the Re = 0 reference (ellipsoid with one principal axis `x` kept parallel to vorticity, rotating about it):
  `χ̇/G = (b² cos²χ + c² sin²χ)/(b² + c²)` — eq (3); period `GT = 2π(b² + c²)/(bc)` — eq (4). For `b = 2c`: `GT = 5π = 15.71`, `max χ̇/G = 0.8` (at χ = 0), `min χ̇/G = 0.2` (at χ = π/2). `[DA00 §1 p.318–319]`
- Non-dimensional equations (5)–(7): `Re(∂u/∂t + u·∇u) = −∇p + ∇²u`, `∇·u = 0`, and the particle equation is the **single scalar rotation equation** `α q Re d²χ/dt² = N`, "where N is the sum of torques on the particle about the x-axis". Scales: length `d`, velocity `Gd`, time `1/G`, torque `μ G dⁿ` (n = 2 in 2-D, n = 3 in 3-D). `Re = ρ_f G d²/μ`, `α = ρ_s/ρ_f`, `q = I/(ρ_s d^{n+2})` with `I` the moment of inertia about the x-axis; "d is an appropriate particle length scale, such as 2b for an ellipse". `[DA00 §2 p.320]`
- Numerics: lattice-Boltzmann (single-relaxation-time BGK, `c_s² = 1/3`, `ν = (2τ−1)/6`), particle advanced by 4th-order Runge–Kutta, no-slip via bounce-back with boundary at mid-link; the geometric radius (not a hydrodynamic radius) is used, justified by a viscosity-independence test at Re = 3.192 for a cylinder with r = 7.99 lattice units in a 278 × 32 domain. `[DA00 §2 p.321]`
- Translational degrees of freedom: not stated (only the rotation equation (7) is written; whether the particle centre was held or free is not stated anywhere).

### A.2 The 3-D ellipsoid case (§6, p.338–340; figs 20–22)

| Item | Value | Tag |
|---|---|---|
| Geometry | "an ellipsoid with semi-axes **a = b = 2c (oblate spheroid)**"; `a` is along `x` (vorticity), `b` and `c` in the shear plane (from fig 1 / eq (1) conventions). The **symmetry axis is the short `c`-axis**, which lies in the flow–gradient plane and rotates in it (tumbling configuration). | `[DA00 §6 p.338; §1 eq (1)]` |
| Density ratio | "a neutrally buoyant three-dimensional particle" → α = 1; fig 21 caption "α = 1". | `[DA00 §6 p.338; fig 21 p.339]` |
| Rotation constraint | Quote: **"The axis of rotation is always the x-axis."** No statement whether this was imposed or observed; but the dynamical model (7) has only the scalar `χ` about `x`, so effectively a **1-DOF rotation model** (consistent with §8: "as long as the axis of rotation of the particle is perpendicular to the plane of shear"). Free 3-D rotation (out-of-plane tilt / log-rolling drift) is **not** in the model. | `[DA00 §6 p.338; §2 eq (7) p.320; §8 p.343]` |
| Re definition (3-D) | "the particle Reynolds number, as defined here, is given by **4Gb²/ν**" — i.e. `d = 2b`. | `[DA00 §6 p.338]` |
| Initial orientation | not stated for the 3-D runs. Fig 21 curves start at the *maximum* rate at Gt = 0 (read ≈ 0.85–0.9 by the printed axis), which corresponds to χ = 0 (b-axis along the gradient direction, cf. eq (3)). The 2-D map of fig 16 states `χ_initial = 0`. | `[DA00 fig 21 p.339; fig 16 caption p.336]` |
| Re values run | "the range of Reynolds numbers 5 to 90 is presented in figures 20 to 22": Re = 5 (fig 20), 50, 70, 90 (fig 21); the fig 22 symbols are at Re = 50, 60, 65, 70, 75 (read, ±0.5). | `[DA00 §6 p.338; figs 20–22]` |
| Domain | "The computational domain for three-dimensional simulations is **40 × 200 × 80 lattice nodes**". Axis assignment **not stated**. Inference (flagged): the 2-D runs are quoted as `L × H` (flow × gradient) with `H = 320 = 10 b` for `b = 32`; if the 3-D 80-node direction is the gradient (`H`) with the 16 × 16 × 8 particle (b = 8) then `H/b = 10` as in 2-D, `L = 200 = 25 b`, and the vorticity extent is `40 = 5 a`; with the 8 × 8 × 4 particle (b = 4): `H/b = 20`, `L = 50 b`, vorticity extent `10 a`. | `[DA00 §6 p.340; §4 p.327]` |
| Periodic vs walls | not stated for 3-D. For the cylinder: "a cylinder suspends freely in between two parallel plates moving in opposite directions" (walls at ±H/2, counter-moving); streamwise/vorticity boundary conditions never stated. | `[DA00 §3 p.323]` |
| Wall speed / shear generation | Counter-moving plates (above); wall speed never given explicitly (2-D: `G = 1/2048` in lattice units with `H = 320` ⇒ `U_wall = ±G H/2 = ±0.078`). For 3-D: not stated. | `[DA00 §4 p.327; §3 p.323]` |
| Particle resolution | "a combination of **8 × 8 × 4 or 16 × 16 × 8** lattice nodes is used to discretize the ellipsoid" (which Re used which: not stated; whether these are full axes 2a × 2b × 2c or semi-axes: not stated — full axes is the natural reading, giving b = 4 or 8 lattice units). | `[DA00 §6 p.340]` |
| Cost | Re = 70, Gt = 0→100, 51 200 steps, 16 SP2 processors, 3.5 h. (⇒ `G Δt = 1.95 × 10⁻³`.) | `[DA00 §6 p.340]` |
| Re_c and fit | "The transition from time-periodic to steady state in the three-dimensional case, considered here, occurs at **Re_c = 81**." "Application of the scaling law, presented in (12), with **C = 200** provides a good fit to the rotation period of the ellipsoid from Re = 50 to Re = 81". | `[DA00 §6 p.340; fig 22 caption]` |

Equations behind the fit `[DA00 §4 p.327–329]`:
- eq (9): `χ̇_min/G ∼ (Re_c − Re)` (linear vanishing of the minimum rate; 2-D data in fig 9 with Re_c ≃ 29).
- eq (10)–(11): `χ̇/G ∼ ε + A(χ − χ₀)²` near the slowest angle χ₀, `T ∼ ∫ dχ/χ̇` ⇒ `T ∼ ε^{−1/2}`.
- eq (12): **`GT = C (Re_c − Re)^{−1/2}`**; 2-D ellipse (b = 2c, α = 1): C = 100, Re_c = 29 (fig 10); 3-D oblate spheroid: C = 200, Re_c = 81 (fig 22); along the 2-D path α = 4 − Re/15: C = 106, Re_c = 34.5 (fig 19). Generalisation `T ∝ |p − p_c|^{−1/2}` for any parameter p (saddle-node bifurcation). `[DA00 §5 p.337–338]`

**Figure 22 read-off (GT vs Re, oblate spheroid, α = 1)** `[DA00 fig 22 p.340]` — five "+" symbols; axis calibrated on the printed ticks (Re 40–90, GT 30–105). Reading uncertainty ≈ ±0.5 in Re, **±1.0 in GT** (marker half-width ≈ 0.8 GT units; the 400 dpi scan is slightly blurred).

| Re (read) | GT (read) | eq (12) with C = 200, Re_c = 81 | Jeffery GT = 5π |
|---|---|---|---|
| 50 | 36.7 | 35.9 | 15.71 |
| 60 | 44.7 | 43.6 | |
| 65 | 50.6 | 50.0 | |
| 70 | 59.5 | 60.3 | |
| 75 | 82.3 | 81.6 | |

(The plotted fit curve is drawn from Re ≈ 45 to ≈ 77 and ends at GT ≈ 96.) All five symbols lie within ≈ 1 GT unit of the fit; the fit is *not* a free two-parameter regression of these points — Re_c = 81 was stated as the observed transition and C = 200 chosen to match.

**Figure 21 read-off (χ̇/G vs Gt, α = 1; Re = 0 analytical, 50, 70, 90)** `[DA00 fig 21 p.339]`. Line styles as printed: (a) Re = 0 solid, (b) Re = 50 **dotted**, (c) Re = 70 **dashed**, (d) Re = 90 solid, flat. The rate signal has two peaks per full rotation (χ → χ + π symmetry), so **GT = 2 × peak spacing**.

| Curve | max χ̇/G (read) | min χ̇/G (read) | peak spacing (Gt) → GT | waveform |
|---|---|---|---|---|
| (a) Re = 0 | 0.70 ± 0.01 | 0.10 ± 0.01 | 7.9 → 15.8 (Jeffery 15.71) | symmetric Jeffery pulses |
| (b) Re = 50 | 0.68 ± 0.02 | 0.030–0.035 | 18.3 ± 0.5 → 36.6 (fig 22: 36.7) | narrow pulses (width ≈ 6 Gt) separated by long flat valleys near 0.03 |
| (c) Re = 70 | 0.66 ± 0.02 | 0.011–0.017 | 29.5 ± 0.5 → 59 (fig 22: 59.5) | same, valleys near 0.01, ≈ 24 Gt long |
| (d) Re = 90 | starts ≈ 0.85–0.9 at Gt = 0 | settles to 0.00 ± 0.01 by Gt ≈ 10 after a small undershoot (≈ −0.01 at Gt ≈ 5) | — | monotone arrest, no further rotation to Gt = 80 |

Caveat on the absolute rate levels: by the printed ordinate the Re = 0 curve oscillates between 0.10 and 0.70, whereas Jeffery eq (3) for b = 2c gives 0.20–0.80 (and fig 20, same paper, shows the Re = 0 curve correctly between 0.2 and 0.8). Either the (a) curve or the fig 21 tick labels are offset by ≈ 0.1; the Re = 90 curve sits at exactly 0.0 by the printed labels, so the labels are probably right for (b)–(d) and the (a) curve is mis-drawn, but this cannot be settled from the paper. **Use periods (fig 22) rather than absolute rate extrema (fig 21) as the quantitative gate; treat the fig 21 minima as ±0.1 systematic.** Peak *widths* and the ratio of pulse width to period are robust.

**2-D vs 3-D bound (conjecture)** `[DA00 §6 p.338]`: "A three-dimensional ellipsoid of a major semi-axis b can be constructed from infinite slices of elliptical cylinders with infinitesimal thickness and major semi-axis in the range of 0 to b. Hence we conjecture that the rate of rotation of the ellipsoid at some Reynolds number, say Re, should be in between the rate of rotation of the two-dimensional cases at Reynolds number 0 and Re." Verified at Re = 5 (fig 20: 3-D curve lies between the 2-D Re = 0 and Re = 5 curves); consequence: 3-D deviation from Jeffery is smaller than 2-D at equal Re, and Re_c(3-D) = 81 > Re_c(2-D) = 29. `[DA00 §6 p.338–340]`

### A.3 2-D elliptical cylinder — what transfers (§4, §5)

Set-up for §4 `[DA00 §4 p.327]`: grid 1600 × 320 (L × H), ellipse b = 32, c = 16 (b = 2c), G = 1/2048, α = 1, viscosity varied to set Re = 5, 10, 15, 20, 24, 26, 28, 30, 40, 50. Confinement therefore **H = 10 b = 5 (2b), L = 50 b**. Periodic rotation for Re ≤ 28, stationary for Re ≥ 30; Re_c ≃ 29 `[fig 8, fig 9, eq (9)]`. Minimum rate vs Re (fig 9, read ±0.001): Re 15 → 0.055, 20 → 0.032, 24 → 0.017, 26 → 0.011, 28 → 0.0037 (text: ε = 0.00367 on 1600 × 320, 0.00376 on 3200 × 640; GT = 82.4 vs 88.2 — "T is very sensitive to the minimum value of the rate of rotation near the point of transition") `[DA00 §4 p.333–335, fig 15]`.

**Confinement dependence of Re_c (ellipse): not stated** — no H/b variation was run for the ellipse (only the resolution doubling at fixed H/b = 10, fig 15). §7 asserts the *exponent* −1/2 is universal "independent of any details of the system, such as the confinement of the channel, the density ratio, and the aspect ratio", but gives no Re_c(H/b) data `[DA00 §7 p.342]`. Confinement effects are quantified only for the circular cylinder: plateau χ̇/G = 0.48 (H/r = 8, Re ≤ 3) and 0.42 (H/r = 4, Re ≤ 9) vs 0.5 unbounded; high-Re power law exponent −0.28 at H/r = 4 instead of −1/2 `[DA00 §3 p.322–326, fig 2]`.

**Density-ratio dependence of Re_c (fig 16, χ_initial = 0, b = 2c)** `[DA00 §5 p.335–336, fig 16]`: "the critical Reynolds number increases with α". Symbols read off fig 16 (Re ±0.5, α ±0.05); + = periodic, ○ = stationary:

| α | periodic (+) at Re | stationary (○) at Re | ⇒ Re_c(α) |
|---|---|---|---|
| ≈ 0.25 | 10, 28 | 30 | 28–30 |
| 1.0 | 0.08(?), 5, 10, 15, 20, 24, 26, 28 | 30, 40, 50 | ≃ 29 (text) |
| 1.6 | — | 36 | < 36 |
| 2.0 | 30 | 40, 50 | 30–40 (dashed curve ≈ 31) |
| 3.0 | 15, 50 | — | > 50 |
| 4.0 | 50 | — | > 50 |
| path α = 4 − Re/15 | (15, 3.0), (18, 2.8), (21, 2.6), (24, 2.4), (27, 2.2), (30, 2.0) all + | — | Re_c = 34.5 along the path (fig 19) |

Dashed boundary Re_c(α): vertical at Re ≈ 28.7 for α ≲ 1, bending to α ≈ 2.0 at Re ≈ 31, α ≈ 2.5 at Re ≈ 37, α ≈ 2.8 at Re = 50. Limit **Re_m = lim_{α→0} Re_c ≈ 28.5–29** (arrow in fig 16). Text: at Re = 50, α_c lies between 2 and 3 (α ≥ 3 rotates forever, α ≤ 2 arrests) `[DA00 §5 p.335]`. At Re = 10 the minimum rate is 0.40, 0.29, 0.09, 0.085 for α = 100, 20, 1, 0.25 (extrapolated 0.083 for α → 0); heavier particles have longer transients and smaller rate fluctuations `[DA00 §5 p.335, fig 17]`.

**Steady orientation after arrest (2-D)** `[DA00 §4 p.329–331]`: χ₀ (slowest angle) = π/2 at Re = 0, ≃ 0.476π at Re = 1, ≃ 0.45π at Re = 28; above Re_c the stable stationary angle decreases with Re, "when Re = 50 the stable angle is about 0.37π" (66.6°, measured from the gradient axis y' to the b-axis, i.e. the long axis is ≈ 23° from the flow direction). Fixed-ellipse torque curves (fig 12): at Re = 40 two zeros, stable χ_a (smaller angle) and unstable χ_b (larger); if released at χ = 80° at Re = 40 "it will initially rotate in a clockwise manner before reaching the stationary orientation" (fig 11c). Fig 8(d) (Re = 30) and fig 21(d) (3-D, Re = 90) show the rate undershooting slightly below zero before settling — the approach to the fixed point is monotone in angle after that, no oscillation. Phase portraits (fig 14): saddle-node at Re_c; at Re = 100 the trajectory converges to the node. 3-D stable angle: **not stated**.

---

## B. Di Giusto, Bergougnoux & Guazzelli 2025

### B.1 Set-up (§2, p.A41-2 – A41-9)

- Shear cell: tank 500 × 40 × 90 mm; transparent Mylar belt driven over two cylinders; belt inner faces **L_y = 27 mm** apart; operative volume 140 mm (flow x) × 30 mm (vorticity/gravity z) × 27 mm (gradient y); `[DG25 §2.2 p.A41-4]`. Confinement `κ = 2a/L_y` = 0.17–0.77 (table 1). Constant shear rate through the depth verified by PIV: γ̇ = 3.18 ± 0.1 s⁻¹ (μ = 0.02 Pa s) to 3.63 ± 0.01 s⁻¹ (μ = 0.8 Pa s) `[DG25 §2.2 p.A41-5]`.
- Fluid: water + Ucon oil (viscosity 0.02–0.8 Pa s) + citric acid (density match to ±4 kg m⁻³); neutrally buoyant, so **Re_p = ρ_f a² γ̇/μ = St** (particle and fluid inertia indistinguishable) `[DG25 §2.1–2.2 p.A41-3, A41-5]`.
- Bodies (table 1, p.A41-4): disks D002–D004 (r = 0.026–0.044), **one oblate spheroid ELL06 (r = 0.561 ± 0.002, a = 2.291 mm, half-thickness 1.286 mm, κ = 0.170, stereolithography, ρ ≈ 1200 kg m⁻³)**, circular rings R009–R05 (r = 0.087–0.452), triangular rings TR003–TR04 (r = 0.031–0.397), L- and T-section rings RL01, RT01. Aspect ratio `r = (thickness 2ℓ)/(diameter 2a)`.
- Re_p range: 0.031–4.9 overall (table 2, p.A41-6). **ELL06: Re_p = 0.031 (5 runs, 5.4 periods), 0.474 (11 runs, 6.6 periods), 0.731 (5 runs, 6.2 periods), 1.029 (6 runs, 8.1 periods)** — all with a finite `t_run/T`, i.e. **rotating at every Re_p tested** `[DG25 table 2]`.
- Orientation measurement: two synchronized cameras (top: flow–gradient plane; side: flow–vorticity plane), 7.5 fps, ≈ 20 px/mm; particle tracked by Watershed; 3-D orientation vector **n** inferred by a two-stream LeNet-5 CNN trained on Blender renders of the .stl (3200 synthetic image pairs; residual < 10 %) `[DG25 §2.3, App. A–B]`. Initial orientation "essentially random"; runs 5–10 per Re_p `[DG25 §2.3.1 p.A41-5]`. Period from the Fourier transform of n-components; rotating vs aligned distinguished by a power-spectrum peak criterion `[DG25 §2.3.3 p.A41-9]`. Migration/lift not measured (runs too short) `[§2.3.1]`.

### B.2 Theories and simulations they compare with

Named in the paper: Jeffery (1922); Bretherton (1962) (equivalent aspect ratio); Subramanian & Koch (2005, 2006); **Einarsson, Candelier, Lundell, Angilella & Mehlig (2015a PRE, 2015b Phys. Fluids)**; **Dabade, Marath & Subramanian (2016)**; Marath & Subramanian (2017, 2018); **Rosén, Einarsson, Nordmark, Aidun, Lundell & Mehlig (2015a PRE 92 063022)** (lattice-Boltzmann + steady-state simulations, bifurcation map); Rosén, Lundell & Aidun (2014); Rosén, Do-Quang, Aidun & Lundell (2015b); **Ding & Aidun (2000)**; Zettner & Yoda (2001); Singh, Koch & Stroock (2013); Borker, Stroock & Koch (2018); Harris & Pittman (1975). `[DG25 §1 p.A41-1–2; §3.2–3.3; refs]`. Curves in figs 7–9 are Jeffery (solid) and Einarsson et al. 2015b (dashed).

### B.3 Orbit-selection result for oblate bodies

- Statement of the theory as used `[DG25 §1 p.A41-2]`: "Prolate bodies are attracted towards the tumbling orbit, while oblate bodies are carried towards the spinning or the tumbling orbit … In the latter case of oblate bodies, a bifurcation between stable and unstable tumbling is predicted at an aspect ratio ≈ 0.14 by asymptotic theories (Einarsson et al. 2015b; Dabade et al. 2016) and is shown to survive up to particle Reynolds numbers of the order of 5 in simulations (Rosén et al. 2015a). However, this bifurcation is unexpectedly not observed experimentally (Di Giusto et al. 2024)."
- Precise threshold `[DG25 §3.3 p.A41-14]`: "This bifurcation was predicted to occur at a critical aspect ratio of **0.137** in the asymptotic theories (Einarsson et al. 2015b; Dabade et al. 2016) … extended in the time-resolved lattice Boltzmann simulations of Rosén et al. (2015a) and shown to survive up to Re_p = 5, even with a confinement of κ = 0.2."
- Which orbit is stable `[DG25 §3.2 p.A41-12]`: "the bifurcation towards **a single stable spinning orbit above a critical aspect ratio of approximately 0.14**, predicted by asymptotic theory of Einarsson et al. (2015b) as well as by the numerical simulations of Rosén et al. (2015a) at Re_p = 1−5". Terminology: **spinning = symmetry axis along vorticity (log-rolling); tumbling = symmetry axis rotating in the flow–gradient plane** `[DG25 §3.2 and fig 9 caption: "alignment … in the plane of shear (third column) or in the vorticity direction (fourth column)"]`.
  - So for the theory: **r < 0.137 (flat): bistable — spinning and tumbling both attracting, selection by initial orientation; r > 0.137 (e.g. r = 0.5): tumbling unstable, spinning (log-rolling) the unique attractor** at small Re_p; Rosén et al. 2015a keep this up to Re_p ≈ 5 (brown dashed line in fig 10 runs from r_eq = 0.137 at Re_p ≤ 0.5, bending to r_eq ≈ 0.25 at Re_p ≈ 6, where it meets the tumbling→fixed-orientation points).
- Fig 10 (p.A41-14), bifurcation map Re_p vs r_eq: solid line `Re_p,c = 15 r` (oblate slender limit, Rosén et al. 2015a) for tumbling → stable fixed orientation; brown dashed line stable/unstable tumbling; cyan octagons = Rosén et al. 2015a LB/steady-state transitions to a fixed point at approximately (r_eq, Re_p) = (0.10, 2.5), (0.13, 3.5), (0.15, 4.5), (0.20, 5.5), (0.25, 8), (0.33, 13), **(0.5, 35)** (read ±15 %; prolate branch at r_eq 4–10, Re_p 20–35). *Caveat: how DG25 converted Rosén's Re definition to Re_p is not stated; taken at face value the r_eq = 0.5 octagon (Re_p ≈ 35, i.e. DA00-Re ≈ 140) sits well above DA00's Re_c = 81 (Re_p ≈ 20).*
- What DG25 actually observed:
  - Circular ring R05 (r = 0.45, r_eq = 0.33) `[fig 7 p.A41-11; §3.2 p.A41-12]`: at Re_p = 0.06 Jeffery orbits with no drift; at Re_p = 0.15 "well described by the asymptotic theory of Einarsson et al. (2015b) and a slow drift toward the spinning orbit may be discernable"; at Re_p = 1.19 "the ring is predicted to drift into a spinning orbit, but in the experiments, the ring remains in a tumbling orbit". "…does not appear to be clearly observed for a circular ring … may necessitate a considerably longer observational timeframe, i.e. over more than 10 rotational periods".
  - Triangular ring TR008 (r = 0.09) `[fig 8]`: drift at Re_p = 0.15 in agreement with Einarsson et al. 2015b; period increases to Re_p = 1.07; **permanent alignment in the plane of shear at Re_p = 4.90** (text says 4.09 once, caption 4.90).
  - Disk D003 (r = 0.03, r_eq = 0.07) `[fig 9 p.A41-13]`: at Re_p = 0.52 drift to **either spinning or tumbling depending on initial orientation** (both limiting orbits predicted by Einarsson et al. 2015a); at Re_p = 2.06 **alignment either in the plane of shear or along vorticity**, again set by the initial orientation.
  - **Oblate spheroid ELL06 (r = 0.561): rotating at all four Re_p (0.03–1.03); no orbit-drift time series shown and it is absent from the bifurcation map (fig 10 legend lists only rings and disks).** Its only appearances are the period plots: fig 5 (letter "O", Re_p ≤ 0.5, period closer to the Harris & Pittman correlation `r_eq = 1.14 r/r^0.156` than to Jeffery's `r_eq = r`) and fig 6 (T_J/T ≈ 1.0 up to Re_p ≈ 0.7, dropping to ≈ 0.75–0.8 at Re_p = 1.03 — read from the grey symbols, ±0.1; the exact values are in the JFM notebook `figure_6/Figure_6.ipynb`, not fetchable from here). Spheroid orbit-drift observations are in the companion paper Di Giusto, Bergougnoux, Marchioli & Guazzelli 2024, JFM 979 A42 (not in `literature/`), which DG25 summarises as: drift "towards two attracting limiting orbits, namely the spinning orbit or the tumbling orbit, contingent upon their initial orientation, as previously observed for oblate spheroids and disks (Di Giusto et al. 2024)" and "a clear bifurcation between stable and unstable tumbling is not observed" `[DG25 §4 p.A41-16]`.
- Dependence on Re_p / St: only Re_p is varied (St = Re_p). "For Re_p < 0.8, inertia has a negligible effect on the period … As Re_p ≳ 1, the period exhibits a noticeable increase and, in the case of flat disks and rings with triangular cross-sections, may even reach infinity." `[DG25 §3.1 p.A41-9–10, fig 6]`. Possible confinement influence (κ ≈ 0.2–0.7) acknowledged, not quantified `[§3.2 p.A41-10; §4]`.

### B.4 Permanent alignment and the critical scaling (quantitative gates available)

- Which shapes/Re: alignment observed for **flat disks (D002 at Re_p = 3.16, 4.84; D003 at 2.06) and triangular rings (TR003 at 0.87, 1.25; TR005 at 1.39, 1.76; TR008 at 4.90; TR01 at 4.27)** — entries with "–" in table 2 / full symbols in fig 10. **No alignment in the Stokes regime** for any ring (contrary to Borker et al. 2018), "even for the lowest aspect ratio of r ≈ 0.03" `[DG25 §3.1 p.A41-9; §4 p.A41-16]`. Circular rings, L/T rings and **the spheroid never aligned** in the tested range.
- Period divergence `[DG25 fig 11 p.A41-15; §3.3]`: `T γ̇ = (67 ± 13) (Re_p,c − Re_p)^{−1/2}` with **Re_p,c = 15 r_eq** taken from the Rosén et al. 2015a slender-oblate asymptote (not fitted); valid for `Re_p,c − Re_p ≲ 3`; "in agreement with the prediction of Ding & Aidun (2000)". (Compare DA00's 3-D fit `GT = 200 (Re − Re_c)^{−1/2}` in DA00 units; converting DA00 to Re_p: `GT = 200 (4(Re_p,c − Re_p))^{−1/2} = 100 (Re_p,c − Re_p)^{−1/2}` with Re_p,c = 20.25 — same exponent, prefactor 100 vs 67 ± 13 for much flatter bodies.)
- Alignment angle `[DG25 fig 12; §3.3]`: flow–shear-plane angle `φ_a = π/2 + r_eq + (30 r_eq)^{1/2} (Re_p,c − Re_p)^{1/2}/15` (Rosén et al. 2015a, r_eq → 0); data "qualitatively follow", deviation 90° − φ_a between ≈ −2° and −15° for r_eq = 0.03–0.11.
- Quantitative orbit-drift rate / orbit-constant decay C(t): **not stated** (no drift rate, no C(t) fit is given; only n(t) traces vs the Einarsson et al. 2015b prediction in figs 7–9, with the data in the JFM notebooks).

---

## C. Synthesis for the case designer (aspect-ratio-2 oblate spheroid, a = b = 2c, r = 0.5)

**(i) Is the DA00 in-plane tumbling orbit an attractor for a free body at Re ≈ 5–80?** Not according to the small-inertia theory as relayed by DG25: for oblate spheroids with r > 0.137 the tumbling orbit is *unstable* and the spinning/log-rolling orbit (symmetry axis along vorticity) is the sole attractor at Re_p ≪ 1, and Rosén et al. 2015a extend that up to Re_p ≈ 5 (DA00-Re ≈ 20) `[DG25 §1, §3.2, §3.3, fig 10]`. DA00 never tested this: their model is a 1-DOF rotation about the vorticity axis (eq (7)) and they simply state "The axis of rotation is always the x-axis" `[DA00 §2, §6]`, so their Re_c = 81 and the fig 22 periods are conditional on that constraint (or on the exact mirror symmetry of the initial condition, which a DNS preserves only up to round-off). DG25's own experiments did not see the theoretical drift clearly (rings, disks and, per their 2024 paper, spheroids kept tumbling over the ≤ 10 periods observed), so the drift is slow, but the DNS gate must not rely on it being absent. What happens to the free-body orbit selection between Re_p ≈ 5 and 20 is **not covered by either paper** (Rosén et al. 2015a/2015b would be the reference). Therefore: **constrain the rotation to the vorticity axis (or enforce the mirror symmetry) for the DA00 comparison, and run the free-body case separately as a second experiment** whose expected small-Re outcome is drift toward log-rolling.

**(ii) Cleanest single gate from each paper.** DA00: the **period GT at fixed Re on the fig 22 curve**, e.g. GT = 36.7 ± 1 at Re = 50 (Jeffery 15.7) and 59.5 ± 1 at Re = 70, with the fit `GT = 200 (81 − Re)^{−1/2}`; the period is measured from peak spacing of χ̇ (two peaks per period) and is insensitive to the fig 21 ordinate ambiguity; the secondary gate is arrest (χ̇ → 0, no rotation) at Re = 90. DG25: for a spheroid of r ≈ 0.5 the only directly comparable data are ELL06's period ratio **T_J/T ≈ 1.0 for Re_p ≤ 0.73 and ≈ 0.75–0.8 at Re_p = 1.03 (fig 6; = DA00-Re 4.1)** — i.e. a free DNS at Re_p ≈ 1 should show the Jeffery period lengthened by ≈ 20–30 % and still rotate, with any drift toward log-rolling too slow to complete within ≈ 8 periods; the DG25 critical-scaling fit (67 ± 13 prefactor, Re_p,c = 15 r_eq) is for flat disks/rings and is *not* a gate for r = 0.5.
