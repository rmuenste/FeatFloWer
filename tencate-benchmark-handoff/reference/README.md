# reference/ — provenance of every data file

## experiment/
- `tenCateData.zip` — byte copy of `tenCateData.zip` at the FeatFloWer repo root
  (4.7 kB, dated 2026-07-31). Its unpacked content is identical (md5) to the
  repo's `tc-ref/` directory and to the website source data
  `ff-redesign/scripts/source-data/sedimentation/`.
- `case_E1_h.csv … case_E4_h.csv` — digitised experimental trajectories:
  column 1 time [s], column 2 gap height h/d (bottom apex of the sphere to the
  bottom wall, over d = 0.015 m). Space separated, no header, irregular time
  stamps. 49 / 31 / 30 / 31 samples.
- `ref_E1.dat … ref_E4.dat` — digitised experimental settling velocity:
  column 1 time [s], column 2 v_z [m/s] (negative = downward). 26 / 27 / 22 /
  16 samples.
- `tc-ref_README.md` — the repo's own README for these files (format, the
  digitisation audit, the gate note). The zip carries no README of its own;
  the digitised curves were provided by R. Münster on 2026-07-31 from the
  paper's figures (figure not named).
- Known defect: the E1/E2 velocity minima are +3.4 % / +3.9 % faster than the
  paper's printed Table II ratios (row `tc_ref_audit`). Use the curves for
  shape, the printed ratios for the peak.

## reference_values.csv
Table I and Table II values verbatim plus derived quantities; the `source`
column states the formula for every derived row.

## featflower/
- `ff_<case>_L<level>_dt<ms>.csv` — settling curve of one certified
  FeatFloWer run, every CFD step: `t_s, z_center_m, gap_m, h_over_d, u_z_mps,
  F_z_N`. Extracted (script `extract_ff_curves.py` in this folder's
  hand-off scratch, formula in the header comments) from the run's
  `particle_force.log`. Sampling = the CFD time step (1.0 / 0.5 / 0.25 ms).
  F_z is the hydrodynamic force on the whole sphere (quarter-domain value ×4
  by `ForceScale`); at rest it equals the buoyant weight at L2/L3.
  `*_lubON.csv` are the lubrication-add-on-ON variants (row `d22_g3_tencate`);
  the base runs have the add-on off.
- `ff_peaks_summary.csv` — per run: peak velocity and time, error vs the
  Table II peak, ratio to u_∞, time at which gap = 2 h_min, touchdown time
  (h/d < 0.05, same criterion as `tools/compare_tencate.py`), time of rest
  (|u| < 1e-4 m/s after the peak), final gap and force, buoyant weight, and
  the peak value recorded in the datasheet row for cross-check (all 15 agree).
- `raw_logs/particle_force_<run>.log` — the untouched source logs (header
  `# time ip fx fy fz tx ty tz px py pz vx vy vz`). Origin
  `q2p1_dns_rundir_<run>/particle_force.log` in the FeatFloWer repo root.
- `decks/<run>/` — `q2p1_param.dat` (CFD deck), `example.json` (PE/particle
  config), `cube.json` (bodies: sphere + bottom plane), `job.sbatch`;
  `decks/e1_l3/start/` the two PE boilerplate files the application reads.
  Run naming: `e<case>_l<level>[_dt0p5_sync|_dt0p25_sync][_g3def]`.
- Pre-sync dt runs (PE at 1 ms while CFD ran at the deck dt — row
  `pe_stepsize_mismatch`) are deliberately NOT included.

## datasheet_rows_tencate.csv
Verbatim rows of `docs/md_docs/dns_validation_datasheet.csv` cited in this
hand-off (header + 49 rows). Columns:
`suite,case,quantity,expected,expected_source,measured,rel_error,tolerance,verdict`.
Row ids are the `case` column.

## figures/
- `e1_l3_vs_piv.png` — from `q2p1_dns_rundir_e1_l3/` (dated 2026-08-01, one
  minute after the run's log; layout = `tools/compare_tencate.py --plot`:
  v_z(t) and h/d(t), DNS line vs PIV dots).
- `d22_g3_tencate.png` — `docs/md_docs/dns_figures/`, produced by
  `tools/d22_g3_tencate_analysis.py`: E1/E2 bottom approach, FBM without and
  with the lubrication add-on vs digitised PIV (paper Fig. 13 analogue).
- `d21_brenner_crossover.png` — `docs/md_docs/dns_figures/`: resolved FBM
  wall-approach force vs Brenner's exact solution, the gap ≲ 2h crossover
  (guide §5, row `d21_prescribed`).

## site_lubrication_export/
Copied read-only from `ff-redesign/scripts/source-data/sedimentation/lubrication/`
(the benchmark website's source data):
- `approach_E1_base.csv`, `approach_E2_base.csv` — `time,u,gap` of the
  certified E1/E2 L3 runs in the approach window (t = 3.0–4.3 s / 1.9–2.7 s).
  Cross-checked against the raw logs: at t = 3.000 s E1 gives
  u = −3.356687e-02, gap = 2.029036e-02 in both. Redundant with
  `ff_E1_L3_dt1p0.csv` but kept as the website's provenance.
- `approach_E1_lub.csv`, `approach_E2_lub.csv` — same window from the
  lubrication-ON runs, with `f_lub` and `n_pairs` columns.
- `cases.csv` — window bounds, h_min = 6.27151451e-4 m, clamp factor 2.
- `brenner_bands.csv` — row `d22_g2_brenner`/`d22_g2b_deficit` numbers:
  deviation from Brenner in the 1h–2h and sub-1h gap bands for FBM alone
  (−14.8 % / −25.9 %), FBM + full analytic resistance (+74.6 % / +67.4 %,
  double counting), FBM + deficit form (+7.2 % / +24.0 %).
