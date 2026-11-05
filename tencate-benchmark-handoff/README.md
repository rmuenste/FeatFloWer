# ten Cate sphere-sedimentation benchmark (E1–E4) — hand-off for an external CFD library

Assembled 2026-10-02 from the FeatFloWer DNS validation campaign (branch
`feature/dns-validation`). Self-contained: everything cited below is either in
this folder or is a page of the paper PDF in `literature/`. Every number carries
its source in brackets; `row <id>` means a row of
`reference/datasheet_rows_tencate.csv` (verbatim copy of the campaign
datasheet rows), `guide §n` means a section of
`campaign_docs/dns_practitioners_guide.md`, `paper p.N` means the ten Cate
2002 PDF page. Where sources disagree this file says so (section 9).

Internal hand-off only: the PDFs in `literature/` are copyrighted and are not
tracked in the FeatFloWer repository for that reason — do not redistribute.

## 1. Purpose

Give another code the complete, traceable definition of the benchmark, the
experimental reference values, the FeatFloWer results it should expect to
reproduce or beat, and the list of traps the campaign fell into, so a second
implementation can be compared against (a) the published experiment and (b)
FeatFloWer's certified numbers without access to the FeatFloWer repository.

## 2. The benchmark in one page

Source: A. ten Cate, C. H. Nieuwstad, J. J. Derksen, H. E. A. Van den Akker,
"Particle imaging velocimetry experiments and lattice-Boltzmann simulations on
a single sphere settling under gravity", *Phys. Fluids* 14 (11), 4012–4025
(2002), doi 10.1063/1.1512918 [paper front page, p.4012].

- Container: depth × width × height = 100 × 100 × 160 mm, closed box, filled
  with silicon oil; the free surface is at the top [paper p.4013, Sec. II].
- Sphere: precision Nylon bearing, d_p = 15 mm, ρ_p = 1120 kg/m³ [p.4013].
- Release: the sphere hangs at the capillary tip of a Pasteur pipette held by
  vacuum, 120 mm from the bottom of the tank (bottom apex of the sphere at
  120 mm, i.e. sphere centre at z = 127.5 mm, initial gap h/d_p = 8), and is
  released from rest by opening an electronic valve [p.4013; the h/d = 8 start
  is also the first point of every digitised trajectory,
  `reference/experiment/case_E*_h.csv`].
- Four fluids (E1–E4), chosen so the sphere Reynolds number spans 1.5–31.9 and
  the Stokes number 0.19–4.13 [Table I, p.4013]. Re is based on the terminal
  velocity u_∞ of the sphere in an *infinite* medium, computed from the Abraham
  (1970) drag correlation C_d = 24/9.06² · (9.06/√Re + 1)² [p.4013, eq. (1)];
  St = (1/9) Re ρ_p/ρ_f [p.4012, p.4020].
- Measured: sphere trajectory (gap height h/d_p vs t, at pixel accuracy from
  the coloured top of the sphere) and settling velocity vs t, from release to
  rest at the bottom, with cross-correlation PIV of the fluid at 60–248 Hz
  [Table I; p.4014; Fig. 5]. All measurements were done twice; the duplicate
  trajectories "practically coincide" [p.4014, Fig. 5(c,d)].
- Stated measurement accuracy: PIV displacement accuracy ≈ 0.1 pixel, i.e. a
  relative error of 2 % for the highest velocities and ≈ 17 % for displacements
  below 0.5 pixel; sphere position determined at pixel accuracy [p.4014].
- No rebound at bottom impact for St = 0.19–4.13 [abstract p.4012; p.4013].
- Time to reach the bottom at Re = 31.9: "approximately 1.3 seconds" [p.4020];
  the digitised E4 trajectory ends at t = 1.19 s (h/d = 0.026)
  [`reference/experiment/case_E4_h.csv`] — see section 9.

## 3. The four cases

Table I of the paper, verbatim [paper p.4013, "TABLE I. Setup of the
sedimentation experiments."]:

| Case | ρ_f [kg/m³] | μ_f [N s/m²] (printed as 10⁻³) | u_∞ [m/s] | Re [-] | St [-] | Camera frequency [s⁻¹] | Resolution |
|---|---|---|---|---|---|---|---|
| E1 | 970 | 373 ×10⁻³ = 0.373 | 0.038 | 1.5 | 0.19 | 60 | low |
| E2 | 965 | 212 ×10⁻³ = 0.212 | 0.060 | 4.1 | 0.53 | 100 | low |
| E3 | 962 | 113 ×10⁻³ = 0.113 | 0.091 | 11.6 | 1.50 | 170 | high |
| E4 | 960 | 58 ×10⁻³ = 0.058 | 0.128 | 31.9 | 4.13 | 248 | high |

Note on the viscosity column: the paper's column header is "μ_f [Ns/m²]" and
the entries are 373, 212, 113, 58; the campaign, the staging tool
(`tools/stage_tencate_case.py`) and `reference/experiment/tc-ref_README.md`
all read these as mPa·s (0.373 … 0.058 Pa·s). Re recomputed from
ρ_f u_∞ d_p/μ_f with these values gives 1.48 / 4.10 / 11.62 / 31.78
(`reference/reference_values.csv`, column `Re_recomputed`), consistent with
the printed 1.5 / 4.1 / 11.6 / 31.9 to rounding — the reading is right.

Experimental peak velocity ratios, Table II bottom rows, verbatim [paper
p.4018, "TABLE II … At the bottom of the table, the experimentally obtained
velocity ratio is included for comparison."]:

| Case | Re | u_max/u_∞ (experiment) | ⇒ u_max = ratio × u_∞ [m/s] |
|---|---|---|---|
| E1 | 1.5 | 0.947 | 0.03599 |
| E2 | 4.1 | 0.953 | 0.05718 |
| E3 | 11.6 | 0.959 | 0.08727 |
| E4 | 31.9 | 0.955 | 0.12224 |

(The products are the campaign's primary peak references; row `tc_ref_audit`,
`e2_l3`, guide §7.) The paper's own base LBM simulations S1–S4 give
u_max/u_∞ = 0.894 / 0.950 / 0.955 / 0.947 [Table II p.4018]; the paper states
"The maximum sedimentation velocity predicted by the simulations is generally
within 1 % of the experimental result, except at the lowest Reynolds number,
where the difference is approximately 5 %. Increased resolution does not result
in an improvement." [p.4021].

## 4. What to compare

1. **Peak settling velocity** u_max against ratio × u_∞ from Table II (table
   above). This is the gated quantity in the campaign (2 % tolerance for
   E2–E4; E1 is a special case — section 6).
2. **Time of peak / shape of u(t)** against the digitised PIV velocity curves
   `reference/experiment/ref_E*.dat`. Shape only: the digitised E1/E2 peak
   amplitudes are +3.4 % / +3.9 % faster than the printed Table II ratios
   (digitisation error, row `tc_ref_audit`, `tc-ref_README.md`); the time of
   peak is a plateau argmin and not a sharp observable (row `e4_ladder_l3`:
   velocity-curve RMS vs PIV 16.4 % of peak unshifted, 4.8 % after a −40 ms
   time-origin shift; guide §7 "~40 ms PIV time-origin caveat").
3. **Trajectory h/d(t)** against `reference/experiment/case_E*_h.csv`.
4. **Approach to the wall** (deceleration from gap ≈ 1 d to rest): record,
   do not gate. Below a gap of ≈ 2 mesh cells the resolved force is no longer
   a converged DNS result in any immersed-boundary code (guide §5; the paper
   itself needed a sub-grid lubrication force below one lattice spacing,
   p.4017 eq. (10), p.4023 Fig. 13). See PITFALLS.md.
5. **Rest**: the sphere must come to rest at the bottom without rebound
   (St ≤ 4.13, paper p.4012).

## 5. How FeatFloWer ran it (certified configuration)

All facts from the certified rundir decks copied to
`reference/featflower/decks/<run>/` (`q2p1_param.dat`, `example.json`,
`cube.json`, `job.sbatch`) and the run logs, unless a guide/row is cited.

- Method: Q2/P1 finite elements (triquadratic velocity, discontinuous linear
  pressure), fictitious boundary method (FBM) for the sphere, geometric
  multigrid, MPI domain decomposition; rigid-body side = the `pe` engine in
  **serial mode** (every rank integrates the same single body; forces reduced
  through the CFD MPI layer), constraint solver `HardContactAndFluid`
  (guide §1; row `d0_solver`). Application `q2p1_bench_sedimentation`
  [`job.sbatch`]. Lubrication add-on **OFF** in all base runs (the two
  `*_g3def` runs are the lubrication-ON variants, row `d22_g3_tencate`).
- Time scheme Crank–Nicolson (`SimPar@TimeScheme = CN`), dt = 1.0 ms base
  (`SimPar@TimeStep = 0.0010d0`), with synced dt-ladder runs at 0.5 and
  0.25 ms (deck `TimeStep` == json `stepsize_`, row `pe_stepsize_mismatch`).
- **Domain: a quarter box**, not the full 100 × 100 × 160 mm box. The
  certified decks point at `_adc/benchSym/bench.prj`
  (`SimPar@ProjectFile`), the mesh `mesh/quarterbox_benchSym/mesh.tri`:
  x ∈ [−0.05, 0], y ∈ [−0.05, 0], z ∈ [0, 0.16] m, 876 coarse hexahedra, 1239
  vertices [mesh file header and coordinates]. Boundary tags
  [`mesh/quarterbox_benchSym/*.par`]: `x.par` = `Symmetry100` on x = 0,
  `y.par` = `Symmetry010` on y = 0 (the two symmetry planes through the sphere
  axis), `xwall.par`/`ywall.par` = `Wall` on x = −0.05 / y = −0.05,
  `bot.par` = `Wall` on z = 0, `top.par` = `Outflow` on z = 0.16 with the deck
  setting `SimPar@NoOutflow = Yes`. The sphere therefore sits on the box axis
  and only one quarter of the flow is computed; the hydrodynamic force is
  reconstructed with `Prop@ForceScale = 0d0,0d0,4d0,0d0,0d0,0d0` (z-force
  × 4, lateral forces and torques zeroed — the symmetry assumption; see
  PITFALLS.md P12). The full-box kit `mesh/fullbox_ten_cate_mesh_v1/` exists
  in the repository but was **not** used for any certified row (section 9).
- Sphere: radius 0.0075 m, centre (0, 0, 0.1275) m at t = 0, at rest
  [`cube.json`; first line of every `particle_force.log`: px = py = 0,
  pz = 0.1275, v = 0], ρ_p = 1120 kg/m³ [`example.json`
  `particleDensity_`]. A PE plane at z = 0 is the bottom wall for the contact
  model [`cube.json`]. Gravity (0, 0, −9.81) m/s² acts on the body in the PE
  json only; the deck has `Prop@Gravity = 0d0,0d0,0d0` (fluid gravity off:
  the PE side seeds the net weight, so the deck value must be zero or
  buoyancy is counted twice — memory note translated in PITFALLS.md P13).
- Fluid: `Prop@Density = 970d0,30d0` and `Prop@Viscosity = 373d-3,1d0` for
  E1 (second entries are a second, unused phase), i.e. the deck carries the
  Table I values ρ_f and μ_f; the PE json repeats them as `fluidDensity_`,
  `fluidViscosity_` (both sides must agree — row `d0_visc`). On which
  viscosity (dynamic vs kinematic) the slot means, see section 9 / PITFALLS P1.
- Mesh levels and resolution (h = dvol^(1/3) of the finest Q2 element,
  `DNS_RESOLUTION` lines of the run logs; element counts = 876 · 8^(L−1)):

  | Level | h_min [m] | D/h | finest-level elements (quarter box) | DOFs inside the sphere (`dofs_per_particle` at t = 1 ms; runs e4_l2 / e1_l3 / e1_l4) |
  |---|---|---|---|---|
  | L2 | 1.31425429e-3 | 11.413 | 7 008 | 1 190 |
  | L3 | 6.27151451e-4 | 23.918 | 56 064 | 8 450 |
  | L4 | 3.05512350e-4 | 49.098 | 448 512 | 63 728 |

  The coarse quarter-box mesh is graded in x and y (cells of ≈ 1.9–3.5 mm
  next to the axis, 14 mm at the far walls) and uniform in z (36 layers of
  4.44 mm) — see `mesh/MESH.md`. "D/h counts Q2 elements; nodal spacing is
  h/2 — halve literature values quoted for LBM/IBM node counts when
  comparing" (guide §2).
- Parallel layout: 32 MPI processes (1 master + 31 subdomains,
  `_mesh/NEWFAC/sub0001/GRID0001..0031.tri`, recursive partitioning of the
  876-element coarse mesh, ≈ 28 coarse elements per subdomain), one node,
  25 GB (L2/L3) or 60 GB (L4) [`job.sbatch` of each run; guide §10].
- Run length: 4300 / 2700 / 1800 / 1300 steps at dt = 1 ms for E1 / E2 / E3 /
  E4 (`MaxNumStep`), chosen to cover the PIV time window plus margin
  (`tools/stage_tencate_case.py`, `CASES` table: PIV ends 4.11 / 2.50 / 1.63 /
  1.20 s).
- Output: particle state every step (`particle_force.log`: time, force,
  torque, position, velocity; `SED_BENCH_VEL/POS` and `DNS_PART_STATE` lines
  in stdout). VTK off for gate runs.
- Diagnostics observed in the certified logs (`reference/featflower/ff_peaks_summary.csv`
  and run logs): fluid CFL max 0.056 (E1 L3) … 0.51 (E4 L4) at dt = 1 ms; 0
  `DNS_SAWTOOTH_WARNING` lines in every run; hydrodynamic force at rest equals
  the buoyant weight (ρ_p−ρ_f)·(π/6)d³·g to 4 digits at L2/L3 (e.g. E1 L3:
  2.6007e-3 N vs 2.6004e-3 N).

## 6. Certified results

Peak settling velocity vs the Table II experimental peak (ratio × u_∞),
dt = 1.0 ms (guide §2 table; rows `e*_l*`; all 15 peaks re-derived from the
raw logs in `reference/featflower/ff_peaks_summary.csv` agree with the rows
to the last printed digit):

| D/h (level) | E1, Re 1.5 | E2, Re 4.1 | E3, Re 11.6 | E4, Re 31.9 |
|---|---|---|---|---|
| 11.4 (L2) | — | +3.6 % (−0.05923) | +3.4 % (−0.09023) | +4.4 % (−0.12763) |
| 23.9 (L3) | −3.2 % (−0.03484) | +0.6 % (−0.05752) | +0.3 % (−0.08753) | +0.8 % (−0.12323) |
| 49.1 (L4) | −5.2 % (−0.03410) | −1.7 % (−0.05623) | −2.1 % (−0.08546) | −1.6 % (−0.12033) |

(values in m/s; sign of the error: + = FeatFloWer faster than experiment.)
Same numbers as ratios u_peak/u_∞ (experiment 0.947 / 0.953 / 0.959 / 0.955):
L2 — / 0.987 / 0.992 / 0.997; L3 0.917 / 0.959 / 0.962 / 0.963;
L4 0.897 / 0.937 / 0.939 / 0.940 [`ff_peaks_summary.csv`].

Time step ladder, synced PE/CFD step (rows `e4_l3_dt_ladder_sync`,
`e4_l4_dt0p5`, `e1_l3_dt0p5`; guide §3): E4 L3 +0.81 / +1.40 / +1.94 % at
dt = 1.0 / 0.5 / 0.25 ms (all stable, 0 sawtooth warnings); E4 L4 −1.56 →
−0.74 % at 1.0 → 0.5 ms; E1 L3 −3.19 → −2.76 %. Additive error model
err = S(L) + T(dt) (`tools/tencate_error_decomposition.py`, row
`e4_l3_dt_ladder_sync`): T(1 ms) = −1.5/−1.1 pp for E4 and −0.8/−0.6 pp for E1
(p = 1 / p = 2 readings, order not pinned), S(L4) = −0.0 … −0.5 pp. The good
dt = 1 ms numbers at L3 "ride on partial cancellation" of a positive spatial and
a negative temporal error (guide §3; row `e4_ladder_l4`: "L3 agreement partly
fortuitous").

Cross-case structure (row `d13_matrix`, guide §2): E2–E4 ladders are identical
within ~0.5 pp despite 8× in Re; the L3→L4 shift is uniform (−2.1 … −2.4 pp)
across all four cases — the spatial error is interface-representation
dominated, not flow-regime dependent.

Simulation vs simulation (row `sim_v_sim`): FeatFloWer's finest configurations
give u_max/u_∞ = 0.897 (E1 L4) / 0.959 (E2 L3) / 0.962 (E3 L3) / 0.947 (E4 L4,
dt 0.5 ms) against ten Cate's own LBM S1–S4 = 0.894 / 0.950 / 0.955 / 0.947:
+0.3 / +0.9 / +0.7 / +0.03 %. The E1 gap to the *experiment* (≈ −5 %) is in
the 2002 paper too (S1 = 0.894 vs 0.947 = −5.6 %), so "the E1 2 %-vs-experiment
gate is unreachable by simulation; the reference band must include the paper's
own sim" (row `sim_v_sim`; guide §2 "Low Re … do not chase it with
resolution").

Timing (from the extracted curves, `ff_peaks_summary.csv`, plateau argmin):
t_peak = 1.90 / 1.51 / 1.18 / 0.96 s (E1–E4, L3); touchdown (h/d < 0.05)
3.85 / 2.33 / 1.56 / 1.15 s; gap = 2h reached at 3.78 / 2.30 / 1.55 / 1.15 s.
Digitised PIV: argmin at 1.15 / 1.83 / 1.13 / 0.81 s (broad plateaus — not a
sharp observable), trajectories end at 4.09 / 2.47 / 1.61 / 1.19 s.

## 7. Reference discipline (why the ladder compares against printed values)

Guide §7 and rows `tc_ref_audit`, `e2_l3`: peak gates use the **printed
Table II ratios × Table I u_∞**, not the digitised curve minima. The audit
(2026-08-01) found the digitised E1/E2 minima at 0.979 / 0.990 of u_∞ against
printed 0.947 / 0.953 (+3.4 / +3.9 %, "also unphysical in trend: wall
retardation should strengthen, not weaken, as Re drops"); E3/E4 digitised
within 0.7 %. Consequently the earliest E4 rows (`e4_ladder_l2/l3/l4`) quote
errors against the digitised E4 peak −0.12303 m/s (L3: +0.16 %) while all
later rows and the guide table quote against Table II −0.12224 m/s (L3:
+0.81 %). Both are in the datasheet extract; use Table II. Never gate against
the unbounded u_∞ itself (the box slows the sphere to ≈ 95 % of u_∞, paper
p.4020); the campaign's pre-recertification serial result −0.1329 m/s "is
+3.8 % vs u_∞ = 0.128 but +8.7 % vs Table II" (`tc-ref_README.md`) and was
traced to a deck viscosity typo (row `d0_visc`, PITFALLS P1).

## 8. Folder map

- `README.md` (this file), `CASES.md`, `PITFALLS.md`, `LITERATURE.md`,
  `MANIFEST.md`.
- `reference/experiment/` — digitised PIV curves (zip + unpacked) and their
  README; `reference/reference_values.csv` — Table I/II numbers + derived
  quantities with source column; `reference/featflower/` — FeatFloWer settling
  curves per run (CSV), the raw `particle_force.log` files, the decks, and
  `ff_peaks_summary.csv`; `reference/datasheet_rows_tencate.csv` — the cited
  datasheet rows verbatim; `reference/figures/`; `reference/site_lubrication_export/`.
- `mesh/` — both coarse meshes, the generator, `MESH.md`.
- `tools/` — staging, comparison and analysis scripts (FeatFloWer-specific;
  read for the comparison *protocol*).
- `campaign_docs/` — the full practitioner's guide and plan extracts.
- `literature/` — the three PDFs.

## 9. Source conflicts and gaps found while assembling this folder

1. **Quarter box vs full box.** The campaign plan (§D1.3, `campaign_docs/plan_extracts.md`)
   specifies the "Full-box mesh (`ten_cate_mesh_v1`)" for the E1–E4 matrix,
   and §11 item 5 promises "one full-box vs quarter-box twin row in D1". Every
   certified rundir deck instead uses the quarter-box `benchSym` mesh
   (`SimPar@ProjectFile = "_adc/benchSym/bench.prj"`, `ForceScale 0,0,4`),
   and no full-box twin row exists in the datasheet (grep for
   quarter/benchSym/symmetry in the datasheet returns no ten Cate row). The
   certified numbers in section 6 are quarter-box, symmetry-plane results.
   The full-box kit shipped here (`mesh/fullbox_ten_cate_mesh_v1`) is a 6 × 6 × 9
   coarse grid (16.7 × 16.7 × 17.8 mm cells) that would need two more
   refinement levels than the quarter box to reach the same D/h (MESH.md).
2. **Viscosity slot.** A FeatFloWer source comment (`QuadSc_main.f90:364`,
   quoted in PITFALLS P1) and the parameter parser treat `Prop@Viscosity` as
   *kinematic* (the `Prop@DynVisc` input is divided by density to fill the
   same slot). The certified DNS decks nevertheless put the Table I *dynamic*
   viscosity (0.373 … 0.058 Pa·s) into `Prop@Viscosity`, the staging tool
   writes it that way, row `d0_visc` compares the slot directly with μ, and the
   FBM force routine receives the slot value without any density factor. The
   resolved-DNS results are physically right only if the slot is consumed as μ
   on that path (ν = 0.058 m²/s would put E4 at Re ≈ 0.03). Conclusion
   recorded here: in this code the same input is dynamic viscosity for the
   resolved Navier–Stokes/FBM path and kinematic for the point-particle (EL)
   closures; the decks are consistent with the DNS path. Not independently
   re-derived from the matrix assembly code in this hand-off.
3. **E4 peak error quoted two ways** (+0.16 % vs digitised, +0.81 % vs
   Table II) — a reference change, explained in section 7.
4. **Time to bottom, E4**: paper text "approximately 1.3 seconds" (p.4020) vs
   digitised trajectory end 1.19 s and FeatFloWer touchdown 1.15 s (L3).
5. **Re of E4**: Table I prints 31.9; ρ_f u_∞ d/μ_f with the printed values
   gives 31.78 (rounding of u_∞ or μ_f in print).
6. **`example.json` leftovers**: `benchStartPosition_ = [1.0, 0.01, 0.1275]`
   and a `domainBoundary_` block referencing `atc_boundary_param_zero.obj`
   (file not present in the rundir) are template residue; the operative
   sphere position is `cube.json` (0, 0, 0.1275), confirmed by the logs.
7. **Pre-sync dt runs**: `q2p1_dns_rundir_e4_l3_dt0p7` (row `e4_l3_dt0p7`)
   and the un-suffixed `*_dt0p5` / `*_dt0p25` rundirs ran the PE integrator at
   1 ms while the CFD used the deck dt (row `pe_stepsize_mismatch`); they are
   excluded from `reference/featflower/`. Only `*_sync` dt runs are included.
8. **Experimental uncertainty of u_max**: the paper gives PIV accuracies
   (2 % / 17 %) and pixel-accurate sphere positions but no error bar on the
   Table II ratios; the duplicate-run statement (p.4014) is the only
   reproducibility evidence.
9. **Digitisation source figure**: `tc-ref_README.md` says the curves were
   digitised "from the paper's figures" without naming the figure (Fig. 5
   shows all four cases; Fig. 8 shows E1/E4 with simulations).
