# DNS validation campaign — expected-vs-actual datasheet

Ledger of every quantitative claim in the DNS (FBM) validation campaign.
Same conventions as `el_validation_datasheet.md`: one row per claim;
verdicts PASS / RECORDED / FAIL / RESOLVED; `expected_source` cites a
RUNBOOK or literature; `tolerance` may be a non-numeric gate. Machine
copy: `dns_validation_datasheet.csv`. Campaign design:
`dns-validation-campaign-plan.md` (repo root).

| suite | case | quantity | expected | expected_source | measured | rel_error | tolerance | verdict |
|---|---|---|---|---|---|---|---|---|
| d0_infra | d0_smoke_twin | 10-step SED_BENCH_VEL trajectory, new ForceScale binary vs old SED_BENCH binary | byte-identical | d0_smoke_twin RUNBOOK | 10/10 lines identical | 0 | byte_identity | PASS |
| d0_infra | d0_smoke_twin | old-binary rerun vs 1300-step reference log head (determinism) | byte-identical | d0_smoke_twin RUNBOOK | 10/10 lines identical | 0 | byte_identity | PASS |
| d0_infra | d0_solver | FBM dynamics under pe master default solver (HardContactEulerLagrange) | nonzero settling trajectory | ten Cate E4 physics | force and velocity exactly 0, no error | n/a | physics_finding | RECORDED |
