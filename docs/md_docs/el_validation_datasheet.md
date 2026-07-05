# Euler-Lagrange Validation Expected-vs-Actual Datasheet

Source branch: `feature/euler-lagrange-phase1`

| suite | case | quantity | expected | expected_source | measured | rel_error | tolerance | verdict |
| --- | --- | --- | ---: | --- | ---: | ---: | --- | --- |
| tier2 | terminal_velocity | v_z | -3.556e-02 | Stokes/DiFelice terminal-velocity balance in `applications/q2p1_el_pipeflow/tier2_cases/terminal_velocity/RUNBOOK.md` | -3.556384e-02 | 1.080e-04 | 2.0e-03 | PASS |
| tier2 | straddling_conservation | rel_error | 0 | Conservation reference in tier2 straddling RUNBOOK and baseline | 0.0 | 0.0_abs | 1.0e-10 | PASS |
| tier2 | straddling_conservation | residual | 0 | Feedback conservation reference in tier2 straddling RUNBOOK and baseline | 5.170037e-26 | 5.170037e-26_abs | 1.0e-10 | PASS |
| tier2 | momentum_conservation | drift_rel | 0 | Total momentum conservation reference in tier2 momentum RUNBOOK and baseline | 4.290711e-07 | 4.290711e-07_abs | 1.0e-05 | PASS |
| tier2 | momentum_conservation_semi | drift_rel | 0 | Semi-implicit total momentum drift characterized as first-order about 3*dt | 1.510406e-03 | 1.510406e-03_abs | 3.1e-03 | PASS |
| tier2 | momentum_conservation_long | drift_rel | 0 | Long-run total momentum conservation reference in tier2 momentum_long RUNBOOK and baseline | 1.000013e-05 | 1.000013e-05_abs | 2.0e-05 | PASS |
| tier2 | saffman_lift | drag_x | 2.356194e-05 | Analytic Stokes drag budget in tier2 saffman_lift baseline | 2.356194e-05 | 2.081e-07 | 1.0e-09 | PASS |
| tier2 | saffman_lift | lift_z | 1.276770e-06 | Analytic Saffman lift budget in tier2 saffman_lift baseline | 1.276770e-06 | 3.091e-07 | 1.0e-09 | PASS |
| validation | pipe_hp_check | u_mean_level_2 | 0.03125 | Hagen-Poiseuille `u_mean = f*R^2/(8*mu)` in pipe_hp_check RUNBOOK | 0.0313647 | 3.670e-03 | 3.0e-02 | PASS |
| validation | pipe_hp_check | u_mean_level_3 | 0.03125 | Hagen-Poiseuille `u_mean = f*R^2/(8*mu)` in pipe_hp_check RUNBOOK | 0.0312310 | 6.080e-04 | 3.0e-02 | PASS |
| validation | v3_ss_frozen | neutral_lift_on_rmean_t1.3 | transient inward kick then static | v3_ss_frozen RUNBOOK W4.1 plumbing expectation | 0.4995905 | NA | 0.02 | PASS |
| validation | v3_ss_frozen | neutral_lift_off_rmean_t1.3 | 0.5000000 | v3_ss_frozen RUNBOOK radial-freeze expectation with `ELLiftModel=none` | 0.5000000 | 0.0 | 0.02 | PASS |
| validation | v3_ss_frozen | dense_lift_on_rmean_t1.3 | steady inward migration | v3_ss_frozen RUNBOOK sustained-migration proof | 0.4980670 | NA | sign_and_monotonic | PASS |
| validation | v3_ss_frozen | dense_outer_drdt_vs_saffman | 1.72e-03 | Pure-Saffman estimate in v3_ss_frozen RUNBOOK | 1.30e-03 | 2.442e-01 | Mei_trim_ratio_0.73_to_0.76 | PASS |
| validation | v3_ss_frozen | dense_inner_drdt_vs_saffman | 5.72e-04 | Pure-Saffman estimate in v3_ss_frozen RUNBOOK | 4.18e-04 | 2.692e-01 | Mei_trim_ratio_0.73_to_0.76 | PASS |
| stage0 | seeding | achieved_phi_seed_12345 | 0.05 | `docs/md_docs/el_stage0_acceptance.md` random seeding reproducibility | 0.049997399884330371 | 5.200e-05 | 2.0e-02 | PASS |
| stage0 | seeding | achieved_phi_seed_99991 | 0.05 | `docs/md_docs/el_stage0_acceptance.md` random seeding reproducibility | 0.049997399884330371 | 5.200e-05 | 2.0e-02 | PASS |
| stage0 | seeding | byte_identity_seed_12345 | clean_diff | `docs/md_docs/el_stage0_acceptance.md` PE 1x1x3 vs 3x3x3 comparison | clean_diff | 0 | byte_identity | PASS |
| stage0 | seeding | byte_identity_seed_99991 | clean_diff | `docs/md_docs/el_stage0_acceptance.md` PE 1x1x3 vs 3x3x3 comparison | clean_diff | 0 | byte_identity | PASS |
| stage0 | file_mode | terminal_series_diff | clean_diff | `docs/md_docs/el_stage0_acceptance.md` Part-A terminal vs file-mode rerun | clean_diff | 0 | byte_identity | PASS |
| validation | v1b_tencate_settling | u_t_E1_oneway | 3.82e-02 | Di Felice fixed point (eps=1) in v1b_tencate_settling RUNBOOK; ten Cate u_inf 0.038 | 3.4547e-02 | 9.6e-02 | eps_eff_0.975_corrected_2pct | PASS |
| validation | v1b_tencate_settling | u_t_E2_oneway | 5.86e-02 | Di Felice fixed point (eps=1); ten Cate u_inf 0.060 | 5.3611e-02 | 8.5e-02 | eps_eff_0.975_corrected_2pct | PASS |
| validation | v1b_tencate_settling | u_t_E3_oneway | 8.66e-02 | Di Felice fixed point (eps=1); ten Cate u_inf 0.091 | 8.0134e-02 | 7.5e-02 | eps_eff_0.975_corrected_2pct | PASS |
| validation | v1b_tencate_settling | u_t_E4_oneway | 1.197e-01 | Di Felice fixed point (eps=1); ten Cate u_inf 0.128 | 1.11905e-01 | 6.5e-02 | eps_eff_0.975_corrected_2pct | PASS |
| validation | v1b_tencate_settling | u_t_E4_twoway | ~0.123 | ten Cate confined 0.955*u_inf; feedback co-flow raises u_t (see RUNBOOK finding 3) | 1.18638e-01 | 3.5e-02 | informative | RECORDED |
| validation | v3_ss_frozen | inertial_outer_dsdt_Rc30 | -1.72e-04 | Matas 2004 fig.14 Rc=30 table at s=0.9 via ghat*1.055e-5 (eps=1) | -1.147e-04 | 1.0e-01 | eps_selfvoidage_0.89_ratio | PASS |
| validation | v3_ss_frozen | inertial_inner_dsdt_Rc30 | 1.79e-05 | Matas 2004 fig.14 Rc=30 table at s=0.1 via ghat*1.055e-5 (eps=1) | 1.59e-05 | 1.1e-01 | eps_selfvoidage_0.89_ratio | PASS |
| validation | v3_ss_frozen | inertial_slope_ratio_outer_inner | 7.12 | Matas table ghat(0.9)/ghat(0.1) profile-shape check | 7.20 | 1.1e-02 | 5.0e-02 | PASS |
| infrastructure | pe_z_periodic_wrap | volume_conservation_through_wrap | 0 | decomposePeriodicZ3D verification in v3_ss_frozen RUNBOOK (2000 steps; 2 wraps) | 2.5e-15 | 2.5e-15_abs | 1.0e-10 | PASS |

For rows with zero expected value, the `rel_error` entry is an absolute error marker because relative error is undefined.
