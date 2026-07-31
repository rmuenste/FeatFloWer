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
