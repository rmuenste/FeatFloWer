# DNS validation campaign — expected-vs-actual datasheet

Ledger of every quantitative claim in the DNS (FBM) validation campaign.
Same conventions as `el_validation_datasheet.md`: one row per claim;
verdicts PASS / RECORDED / FAIL / RESOLVED; `expected_source` cites a
RUNBOOK or literature; `tolerance` may be a non-numeric gate. Machine
copy: `dns_validation_datasheet.csv`. Campaign design:
`dns-validation-campaign-plan.md` (repo root).

| suite | case | quantity | expected | expected_source | measured | rel_error | tolerance | verdict |
|---|---|---|---|---|---|---|---|---|
| d0_infra | d0_smoke_twin | 10-step SED_BENCH_VEL trajectory new ForceScale binary vs old SED_BENCH binary | byte-identical | d0_smoke_twin RUNBOOK | 10/10 lines identical | 0 | byte_identity | PASS |
| d0_infra | d0_smoke_twin | old-binary rerun vs 1300-step reference log head (determinism) | byte-identical | d0_smoke_twin RUNBOOK | 10/10 lines identical | 0 | byte_identity | PASS |
| d0_infra | d0_solver | FBM dynamics under pe master default solver (HardContactEulerLagrange) | nonzero settling trajectory | ten Cate E4 physics | force and velocity exactly 0 no error | n/a | physics_finding | RECORDED |
| d0_infra | d0_smoke_twin | parallel-PE new vs old binary 10-step trajectory | match to printed precision | d0_smoke_twin RUNBOOK | matches all steps | 0 | printed_precision | PASS |
| d0_infra | d0_smoke_twin | DNS_RESOLUTION global inside-DOF count parallel vs serial | same-order (interface duplication only) | d0_smoke_twin RUNBOOK | 1108 vs 1190 (6.9% duplication gap) | 0.069 | same_order_gate | PASS |
| d0_infra | d0_hmin | mode-consistent global h_min on identical mesh | fine-mesh h_min equal to serial | mesh identity (byte-identical .tri) | was min-over-ranks of local MAX h (inverted negation trick); dvol itself correct | 4x | physics_finding | RESOLVED |
| d0_infra | d0_hmin | h_min and D_over_h after reduction fix, serial vs parallel | identical | byte-identical mesh | 1.31425429e-3 / 11.413 in both modes | 0 | exact_match | PASS |
| d0_infra | d0_visc | E4 deck CFD viscosity vs paper and PE side | 0.058 (ten Cate Table I; example.json fluidViscosity_) | ten_cate_piv.pdf; example.json | deck had 53d-3 (-8.6%); Schiller-Naumann predicts +3.6% u_t from this alone (measured overshoot vs PIV: +8.1%) | 0.086 | physics_finding | RESOLVED |
| d1_metrology | e4_ladder_l2 | E4 peak settling velocity vs confined PIV (corrected deck, D/h=11.4) | -0.12303 m/s | tc-ref/ref_E4.dat | -0.12763 m/s | +0.0374 | level_ladder | RECORDED |
| d1_metrology | e4_ladder_l3 | E4 peak settling velocity vs confined PIV (corrected deck, D/h=23.9) | -0.12303 m/s | tc-ref/ref_E4.dat | -0.12323 m/s | +0.0016 | 2pct | PASS |
| d1_metrology | e4_ladder_l3 | E4 velocity-curve RMS vs PIV after -40ms origin shift | small residual | tc-ref/ref_E4.dat | 4.8% of peak (16.4% unshifted; broad-plateau argmin artifact + ~40ms origin offset) | 0.048 | shape_documentation | RECORDED |
| d1_metrology | e4_ladder_l4 | E4 peak settling velocity vs confined PIV (corrected deck, D/h=49.1, dt=1.0ms) | -0.12303 m/s | tc-ref/ref_E4.dat | -0.12033 m/s | -0.0219 | level_ladder | RECORDED |
| d1_metrology | e4_ladder_l4 | spatial ladder monotonicity at fixed dt=1.0ms (L2/L3/L4 = +3.74/+0.16/-2.19%) | monotone approach to PIV | level-ladder convention | non-monotone; spatial and temporal errors have opposite signs, L3 agreement partly fortuitous; motivates D1 dt ladder at L4 | n/a | physics_finding | RECORDED |
| d1_metrology | e3_l3 | E3 (Re=11.6) peak settling velocity vs confined PIV (D/h=23.9, dt=1.0ms) | -0.08665 m/s | tc-ref/ref_E3.dat | -0.08753 m/s; timing lag +0.050 s | +0.0102 | 2pct | PASS |
| d1_metrology | e2_l3 | E2 (Re=4.1) peak settling velocity vs confined PIV (D/h=23.9, dt=1.0ms) | -0.05936 m/s | tc-ref/ref_E2.dat | -0.05752 m/s; digitized E2 ratio v_PIV/u_inf=0.990 is off-trend (E3 0.952, E4 0.961; retardation should grow as Re drops); vs 0.955*u_inf expectation run is +0.38%; dt/2 rerun would bound temporal term | -0.0309 | 2pct | FAIL |
