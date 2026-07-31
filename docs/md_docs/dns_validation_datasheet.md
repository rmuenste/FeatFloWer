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
| d0_infra | d0_dvol | parallel-PE dvol at DOF-matched level | fine-mesh h_min equal to serial | mesh identity (byte-identical .tri) | h_min 6.980e-3 = coarse h vs serial 1.745e-3 | 4x | physics_finding | RECORDED |
