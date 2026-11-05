# CASES — E1 … E4, one section per case

Conventions: all four cases share the sphere (d = 15 mm, ρ_p = 1120 kg/m³),
the box (100 × 100 × 160 mm) and the release (bottom apex 120 mm above the
bottom, i.e. centre z = 0.1275 m, h/d = 8, from rest) [paper p.4013]. Only the
fluid changes. "Table I/II" = ten Cate 2002 p.4013/p.4018. Derived numbers are
computed from the Table I values (formulas in `reference/reference_values.csv`).
FeatFloWer values: `reference/featflower/ff_peaks_summary.csv` (re-derived from
the certified logs; identical to the datasheet rows named). Errors are vs the
Table II experimental peak (ratio × u_∞); + means faster than experiment.
"Time to wall" is the FeatFloWer touchdown time, first step with h/d < 0.05
(the criterion of `tools/compare_tencate.py`); "gap = 2h" is the first step with
gap ≤ 2 h_min, below which the force is sub-grid (guide §5). Wall-clock times
are the `Overall time` of the run's `_data/Statistics.txt`, 32 MPI processes on
one node (nx-class nodes, guide §10).

Resolution per level (all cases, quarter-box mesh): L2 h = 1.314 mm, D/h = 11.4;
L3 h = 0.627 mm, D/h = 23.9; L4 h = 0.306 mm, D/h = 49.1 (`DNS_RESOLUTION` log lines).

Common recommendation on output: log particle position, velocity and
hydrodynamic force at every time step (FeatFloWer: `particle_force.log`); the
time of peak is a plateau argmin and the approach phase lasts only a few
hundred steps, so coarser cadence loses both. Volume output is optional and
large (guide/stage tool: once-refined-mesh VTK is ~8× the file size).

---

## E1 — Re 1.5 (the Stokes-like case, known-hard)

Parameters [Table I]: ρ_f = 970 kg/m³, μ_f = 0.373 Pa·s, u_∞ = 0.038 m/s,
Re = 1.5, St = 0.19. Derived: ν = 3.845e-4 m²/s, ρ_p/ρ_f = 1.1546,
Re_recomputed = 1.48, τ_p = ρ_p d²/(18 μ) = 0.0375 s, τ_ν = d²/ν = 0.585 s
(paper p.4019: 0.59 s), τ_adv = d/u_∞ = 0.395 s (paper: 0.39 s), buoyant weight
2.6004e-3 N.

Experiment: u_max/u_∞ = 0.947 → u_max = 0.03599 m/s [Table II]. Digitised
velocity curve: 26 samples, t = 0.053–4.113 s, minimum −0.03720 m/s at
t = 1.145 s (3.4 % too fast — row `tc_ref_audit`); digitised trajectory
49 samples, t = 0.013–4.092 s, ends at h/d = 0.118. Paper's own LBM S1:
0.894 (−5.6 % vs experiment) [Table II].

FeatFloWer (rows `e1_l3`, `e1_l3_dt0p5`, `e1_l4`, `e1_l3_dt0p5_sync` numbers
in `e4_l3_dt_ladder_sync`):

| run | D/h | dt [ms] | u_peak [m/s] | ratio u_peak/u_∞ | err vs Table II | t_peak [s] | gap = 2h [s] | time to wall [s] | rest (|u|<1e-4) [s] | wall clock |
|---|---|---|---|---|---|---|---|---|---|---|
| L3 | 23.9 | 1.0 | −0.03484 | 0.917 | −3.19 % | 1.896 | 3.777 | 3.852 | 4.170 | 3221 s (54 min) |
| L3 | 23.9 | 0.5 | −0.03499 | 0.921 | −2.76 % | 1.889 | 3.763 | 3.839 | 4.027 | 6199 s |
| L4 | 49.1 | 1.0 | −0.03410 | 0.897 | −5.23 % | 1.979 | 3.968 | 3.938 | 4.264 | 25242 s (7.0 h, 60 GB) |

No L2 run exists for E1 (guide §2 table "—").

Known residual: refinement moves E1 *away* from the experiment (−3.2 → −5.2 %)
and toward ten Cate's own LBM (S1 = 0.894; FeatFloWer L4 = 0.897, +0.3 %;
Richardson h→0 estimate 0.886–0.899 brackets 0.894 — rows `e1_l4`,
`sim_v_sim`, `err_decomposition`). The paper itself: "At the lowest Reynolds
number, a systematic underprediction of the velocity ratio of approximately 5 %
is found … the deviation is independent of the resolution" [p.4022–4023]. The
campaign verdict (row `sim_v_sim`, guide §2): the 2 %-vs-experiment gate is
unreachable by simulation at Re 1.5; the honest reference band includes the
paper's S-series. A radius-smearing explanation was tested and refuted (row
`e1_l4`). Treat an E1 result in the 0.89–0.92 band as consistent with both
published simulations; a result at 0.947 would be remarkable and should be
examined for a compensating error.

Lubrication-ON variant (row `d22_g3_tencate`, run `e1_l3_g3def`): identical to
the base run down to gap = 2h, then 15 % slower at gap = 1h and 12 % at 0.5h;
time from gap = 2h to rest 334 ms vs 301 ms base; both settle to a numerical
equilibrium gap of 0.13–0.14 mm (FBM/PE contact equilibrium, not physics).

Recommended run: ≥ 4.3 s of physical time (4300 steps at 1 ms; PIV window ends
4.11 s) — the sphere needs ≈ 3.85 s to reach the wall and ≈ 4.2 s to rest at
L3. Noise floor at the plateau: rms of detrended F_z 0.57 % (L3) / 0.59 % (L4)
of the mean force (row `d12_plateau_noise`).

## E2 — Re 4.1

Parameters [Table I]: ρ_f = 965 kg/m³, μ_f = 0.212 Pa·s, u_∞ = 0.060 m/s,
Re = 4.1, St = 0.53. Derived: ν = 2.197e-4 m²/s, ρ_p/ρ_f = 1.1606,
Re_recomputed = 4.10, τ_p = 0.0660 s, τ_ν = 1.024 s, τ_adv = 0.250 s, buoyant
weight 2.6870e-3 N.

Experiment: u_max/u_∞ = 0.953 → 0.05718 m/s [Table II]. Digitised velocity
curve: 27 samples, t = 0.068–2.501 s, minimum −0.05936 m/s at t = 1.831 s
(3.9 % too fast — row `tc_ref_audit`, `e2_l3`); trajectory 31 samples,
t = 0–2.467 s, ends h/d = 0.052. LBM S2: 0.950.

FeatFloWer (rows `e2_l2`, `e2_l3`, `e2_l4`):

| run | D/h | dt [ms] | u_peak [m/s] | ratio | err vs Table II | t_peak [s] | gap = 2h [s] | time to wall [s] | rest [s] | wall clock |
|---|---|---|---|---|---|---|---|---|---|---|
| L2 | 11.4 | 1.0 | −0.05923 | 0.987 | +3.59 % | 1.453 | 2.185 | 2.259 | 2.289 | 429 s |
| L3 | 23.9 | 1.0 | −0.05752 | 0.959 | +0.60 % | 1.505 | 2.299 | 2.329 | 2.469 | 2071 s (35 min) |
| L4 | 49.1 | 1.0 | −0.05623 | 0.937 | −1.66 % | 1.519 | 2.395 | 2.383 | 2.533 | 16228 s (4.5 h, 60 GB) |

Row `e2_l3` is the row that exposed the digitisation defect: against the
digitised minimum the L3 result read −3.09 % (FAIL), against the printed Table
II ratio +0.59 % (PASS). vs LBM S2: FeatFloWer L3 0.959 is +0.9 %.

Lubrication-ON variant (`e2_l3_g3def`, row `d22_g3_tencate`): u(1h) 16 %
slower, u(0.5h) 17 %; 2h-to-rest 140 ms vs 124 ms.

Recommended run: ≥ 2.7 s (2700 steps at 1 ms; PIV ends 2.50 s; wall at
≈ 2.33 s, rest ≈ 2.47 s at L3). Plateau force noise 1.58 / 0.84 / 1.12 % at
L2/L3/L4 (row `d12_plateau_noise`; the L4 value is a mild non-monotonicity).

## E3 — Re 11.6

Parameters [Table I]: ρ_f = 962 kg/m³, μ_f = 0.113 Pa·s, u_∞ = 0.091 m/s,
Re = 11.6, St = 1.50. Derived: ν = 1.175e-4 m²/s, ρ_p/ρ_f = 1.1642,
Re_recomputed = 11.62, τ_p = 0.1239 s, τ_ν = 1.915 s, τ_adv = 0.165 s, buoyant
weight 2.7390e-3 N.

Experiment: u_max/u_∞ = 0.959 → 0.08727 m/s [Table II] — the highest ratio of
the four; the paper attributes the maximum at E3 to decreasing wall hindrance
with Re, and the lower E4 ratio to the sphere still accelerating when it
reaches the bottom [p.4020–4021]. Digitised velocity curve: 22 samples,
t = 0.015–1.635 s, minimum −0.08665 m/s at 1.130 s (−0.7 %, clean);
trajectory 30 samples, t = 0–1.612 s, ends h/d = 0.039. LBM S3: 0.955.

FeatFloWer (rows `e3_l2`, `e3_l3`, `e3_l4`):

| run | D/h | dt [ms] | u_peak [m/s] | ratio | err vs Table II | t_peak [s] | gap = 2h [s] | time to wall [s] | rest [s] | wall clock |
|---|---|---|---|---|---|---|---|---|---|---|
| L2 | 11.4 | 1.0 | −0.09023 | 0.992 | +3.39 % | 1.166 | 1.480 | 1.510 | 1.633 | 286 s |
| L3 | 23.9 | 1.0 | −0.08753 | 0.962 | +0.30 % | 1.180 | 1.546 | 1.556 | 1.580 | 1369 s (23 min) |
| L4 | 49.1 | 1.0 | −0.08546 | 0.939 | −2.07 % | 1.205 | 1.595 | 1.591 | 1.619 | 11131 s (3.1 h, 60 GB) |

Row `e3_l3` also records a timing lag of +0.050 s of the FeatFloWer curve
relative to the digitised PIV curve (plateau-argmin/time-origin caveat, guide
§7). vs LBM S3: L3 0.962 is +0.7 %.

Recommended run: ≥ 1.8 s (1800 steps; PIV ends 1.63 s; wall ≈ 1.56 s).
Plateau force noise 2.44 / 1.66 / 1.02 % at L2/L3/L4.

## E4 — Re 31.9 (the campaign's workhorse / twin-gate case)

Parameters [Table I]: ρ_f = 960 kg/m³, μ_f = 0.058 Pa·s, u_∞ = 0.128 m/s,
Re = 31.9, St = 4.13. Derived: ν = 6.042e-5 m²/s, ρ_p/ρ_f = 1.1667,
Re_recomputed = 31.78 (print rounding), τ_p = 0.2414 s, τ_ν = 3.724 s (paper
p.4020: 3.72 s), τ_adv = 0.117 s (paper: 0.12 s), buoyant weight 2.7737e-3 N.
The sphere reaches the bottom before the wake is developed (τ_ν ≫ time to
bottom) and "hardly decelerates prior to contact" [p.4020].

Experiment: u_max/u_∞ = 0.955 → 0.12224 m/s [Table II]. Digitised velocity
curve: 16 samples, t = 0.060–1.198 s, minimum −0.12303 m/s at 0.806 s
(+0.6 %, clean); trajectory 31 samples, t = −0.007–1.191 s, ends h/d = 0.026.
Time to bottom "approximately 1.3 seconds" [p.4020] vs curve end 1.19 s. LBM
S4: 0.947 (and S6 at r_0 = 8 lu: 0.947; S8 at 2 lu: 0.921) [Table II].

FeatFloWer (rows `e4_ladder_l2/l3/l4`, `e4_l3_dt_ladder_sync`, `e4_l4_dt0p5`,
`dt_stability_refuted`):

| run | D/h | dt [ms] | u_peak [m/s] | ratio | err vs Table II | t_peak [s] | gap = 2h [s] | time to wall [s] | rest [s] | wall clock |
|---|---|---|---|---|---|---|---|---|---|---|
| L2 | 11.4 | 1.0 | −0.12763 | 0.997 | +4.41 % | 0.938 | 1.096 | 1.113 | 1.123 | 210 s |
| L3 | 23.9 | 1.0 | −0.12323 | 0.963 | +0.81 % | 0.956 | 1.145 | 1.150 | 1.160 | 949 s (16 min) |
| L3 | 23.9 | 0.5 | −0.12396 | 0.968 | +1.40 % | 0.959 | 1.139 | 1.145 | 1.154 | 1971 s |
| L3 | 23.9 | 0.25 | −0.12461 | 0.974 | +1.94 % | 0.947 | 1.135 | 1.140 | 1.189 | 3827 s |
| L4 | 49.1 | 1.0 | −0.12033 | 0.940 | −1.56 % | 0.980 | 1.177 | 1.175 | 1.186 | 7295 s (2.0 h, 60 GB) |
| L4 | 49.1 | 0.5 | −0.12133 | 0.948 | −0.74 % | 0.968 | 1.169 | 1.167 | 1.177 | 14528 s (4.0 h) |

The dt = 0.25 ms run is the one that refuted the "stability floor" (row
`dt_stability_refuted`: stable, 0 sawtooth warnings, once PE and CFD used the
same dt). The spatial ladder at fixed dt = 1 ms is non-monotone (+4.4 / +0.8 /
−1.6 %) because spatial (+) and temporal (−) errors have opposite signs (row
`e4_ladder_l4`); the fit gives T(1 ms) ≈ −1.5/−1.1 pp and S(L4) ≈ −0.0…−0.5 pp
(row `e4_l3_dt_ladder_sync`). Finest configuration L4/0.5 ms: 0.948 vs LBM S4
0.947 (+0.03 %, row `sim_v_sim`). Velocity-curve RMS vs digitised PIV: 16.4 %
of peak as-is, 4.8 % after a −40 ms time-origin shift (row `e4_ladder_l3`).
Fluid CFL max 0.24 (L3) / 0.51 (L4) at 1 ms.

Historical note (row `d0_visc`, `tc-ref_README.md`): the pre-campaign E4 deck
carried μ = 0.053 instead of 0.058 Pa·s (−8.6 %) and produced −0.1329 m/s,
+8.7 % vs Table II — the error that triggered the recertification.

Recommended run: ≥ 1.3 s (1300 steps at 1 ms; PIV ends 1.20 s; wall ≈ 1.15 s).
Plateau force noise 4.13 / 2.35 / 0.92 % at L2/L3/L4 (row `d12_plateau_noise`;
controlled constant-velocity ladder 3.84 / 1.79 / 0.88 %, row
`d12_prescribed_ladder`).

## Observed in all cases (from the extracted curves, not in any row)

- The final resting gap is resolution dependent: ≈ 0.60–0.62 mm (≈ 0.46 h)
  at L2, 0.13–0.16 mm (≈ 0.21–0.25 h) at L3, 0.0 at L4 (sphere on the PE
  contact plane, part of the weight carried by the contact: final F_z < buoyant
  weight at L4). This is the numerical FBM/contact equilibrium, not a physical
  observable — do not compare resting gaps between codes.
- Hydrodynamic force at rest = buoyant weight to 4 digits at L2/L3 (e.g. E1
  L3 2.6007e-3 vs 2.6004e-3 N), a cheap consistency check for any code.
