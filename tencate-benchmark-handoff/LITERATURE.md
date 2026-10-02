# LITERATURE

Internal hand-off only: the PDFs below are copyrighted; the FeatFloWer
repository keeps them untracked for that reason (`literature/README.md` of
the repo, campaign plan §13). Do not redistribute this folder outside the two
groups.

## Primary reference (the benchmark)

**A. ten Cate, C. H. Nieuwstad, J. J. Derksen, H. E. A. Van den Akker**,
"Particle imaging velocimetry experiments and lattice-Boltzmann simulations on
a single sphere settling under gravity", *Physics of Fluids* **14** (11),
4012–4025 (November 2002). doi 10.1063/1.1512918. Received 26 December 2001,
accepted 19 August 2002, published 3 October 2002. Kramers Laboratorium voor
Fysische Technologie, Delft University of Technology.
Verified from the PDF front page and the running heads
(`literature/ten_cate_nieuwstad_derksen_vandenakker_2002_phys_fluids_14_4012.pdf`,
copy of `ten_cate_piv.pdf` at the repo root). Note: the hand-off brief spelled the second author "Nijstad"; the PDF and the
campaign plan §13 print **Nieuwstad**.
Used for: the case definitions (Table I, p.4013), the experimental and LBM
peak ratios (Table II, p.4018), the hydrodynamic-radius sensitivity (Table
III, p.4022), the measurement accuracy statement (p.4014), the time scales
(p.4019–4020), the lubrication force eq. (10) and Fig. 13 (p.4017, p.4023),
and the digitised trajectories (`reference/experiment/`, from its figures).
Page map of the PDF: 4012 abstract/intro; 4013 setup + Table I; 4014 PIV
accuracy, duplicates; 4015–4017 LBM method, internal-fluid correction eq. (8),
calibration eq. (9), lubrication eq. (10); 4018 Table II; 4019–4021 results,
time scales, Fig. 8 comparison; 4022 Table III, radius sensitivity; 4023
lubrication results, Fig. 13; 4024–4025 conclusions, references.

**F. Abraham**, "Functional dependence of drag coefficient of a sphere on
Reynolds number", *Phys. Fluids* **13**, 2194 (1970) — cited by ten Cate as
ref. 17 (p.4024) as the C_d correlation defining u_∞ and hence Re. Not in the
folder.

## Resolution reference

**M. Uhlmann, J. Dušek**, "The motion of a single heavy sphere in ambient
fluid: A benchmark for interface-resolved particulate flow simulations with
significant relative velocities", *International Journal of Multiphase Flow*
**59** (2014) 221–243. Verified from the PDF front page
(`literature/uhlmann_dusek_2014.pdf`). Spectral/spectral-element reference
data for a sphere at density ratio 1.5, Galileo number 144–250 (Re 185–365),
plus immersed-boundary runs at D/Δx = 15 … 48 that quantify the error vs
resolution per regime (abstract; Sec. on IBM results). Used by the campaign
as the independent published resolution-requirement study (plan §D1.4,
"documented resolution requirements (D/h up to ~24), directly comparable to our
guideline study"); the optional D1.4 cross-check was not executed. Remember the
node-count convention when comparing: FeatFloWer's D/h counts Q2 elements
(nodal spacing h/2), guide §2.

## Near-wall reference

**H. Brenner**, "The slow motion of a sphere through a viscous fluid towards a
plane surface", *Chem. Eng. Sci.* **16** (1961) 242–251
(`literature/brenner_1961.pdf`; citation per the repo's `literature/README.md`
and ten Cate's ref. 2, p.4024). Exact Stokes-flow resistance of a sphere
approaching a plane; the campaign's crossover rule (guide §5, rows
`d21_prescribed`, `d22_g2_brenner`) and `reference/figures/d21_brenner_crossover.png`
are measured against it.

## Cited but not copied

- **P. Causin, J.-F. Gerbeau, F. Nobile**, "Added-mass effect in the design of
  partitioned algorithms for fluid–structure problems", *Comput. Methods Appl.
  Mech. Eng.* **194** (2005) 4506–4527 — the theory originally invoked for the
  dt "stability floor" (row `dt_stability`); the attribution was refuted (row
  `dt_stability_refuted`, PITFALLS P2). In the repo as `literature/causin_2005.pdf`.
- **H. Hasimoto**, "On the periodic fundamental solutions of the Stokes
  equations and their application to viscous flow past a cubic array of
  spheres", *J. Fluid Mech.* **5** (1959) 317–328 — ten Cate's radius
  calibration (eq. (9), p.4017) and the campaign's effective-radius study
  (PITFALLS P11). In the repo as `literature/hasimoto_1959.pdf`.
- Method lineage of the FeatFloWer solver (plan §13, citations confirmed
  there): D. Wan, S. Turek, *Int. J. Numer. Meth. Fluids* 51 (2006) 531–566
  (FBM volume-integral force); R. Münster, O. Mierka, S. Turek, *Int. J.
  Numer. Meth. Fluids* 69 (2012) 294–313 (3D FEM-FBM).

## Campaign documents (FeatFloWer repository, branch `feature/dns-validation`)

- `docs/md_docs/dns_practitioners_guide.md` v2.5 (2026-10-02) — copied in full
  to `campaign_docs/`. §1 prerequisites, §2 resolution table, §3 dt, §4 noise,
  §5 near-contact, §6 contact parameters, §7 reference discipline, §10 costs.
- `docs/md_docs/dns_validation_datasheet.csv` — the campaign ledger; the
  ten Cate-relevant rows are copied verbatim to
  `reference/datasheet_rows_tencate.csv`.
- `dns-validation-campaign-plan.md` — §D0.4 (truth chain), §D1.3 (E1–E4
  matrix), §11 (risks), §13 (literature) extracted to
  `campaign_docs/plan_extracts.md`.
- `tc-ref/README.md` — digitised-data README, copied to
  `reference/experiment/tc-ref_README.md`.
- `pipemesh_v1/handoff_euler_lagrange_drag_validation_ten_cate.md` — an
  earlier (EL-era) hand-off on the same benchmark; superseded on the DNS side
  by this folder (its "~2 % at level 3" claim predates the recertification).
  Not copied.
- Owner memory notes (not copied; translated in PITFALLS.md):
  `ff-viscosity-convention`, `dns-dt-stability-floor`,
  `dns-pe-solver-selection`, `ff-deck-staging-pitfalls`,
  `enable-lubrication-legacy-kroupa`.
