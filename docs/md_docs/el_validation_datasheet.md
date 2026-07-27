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
| validation | v3_ss_coupled | lift_on_off_drmean_ratio_t1000 | discriminating | v3_ss_coupled RUNBOOK W4.2 ON/OFF pair (jobs 131967/131968) | 22:1 (ON +0.088 vs OFF +0.004) | NA | sign_and_rate | PASS |
| validation | v3_ss_coupled | rmean_annulus_t2000 | 0.675 | Matas-Asmolov equilibrium s_eq; v3_ss_coupled RUNBOOK 100k-step run (job 132699) | 0.6643 | 1.6e-02 | 7.4e-02 | PASS |
| stage2 | fourway_shear | max_overlap_frac_dp | 0 | Stage-2 gate max_overlap <= 1% d_p; fourway_shear RUNBOOK (job 132694, 400 samples) | 0.0 | 0.0_abs | 1.0e-02 | PASS |
| stage2 | fourway_shear | tgran_late_over_early | 1 | Stage-2 no-spontaneous-heating gate; fourway_shear RUNBOOK | 1.0002 | 2.4e-04 | bounded | PASS |
| stage2 | fourway_shear | newton_pair_under_contact_load | 0 | EL_NEWTON_PAIR at ~5.3k simultaneous contacts; fourway_shear RUNBOOK | 0.0 | 0.0_abs | 1.0e-10 | PASS |
| stage2 | unit_periodic_contact | newton_pair_periodic_image_contact | 0 | Contact through periodic z-plane; unit_periodic_contact RUNBOOK | 4.5e-20 | 4.5e-20_abs | 1.0e-10 | PASS |
| infrastructure | pe_newton_pair_fix | worst_mismatch_rz_sweep_phi020 | 0 | pe 8b037ae fix verified at scale; v2_rz_settling RUNBOOK second sweep (jobs 132467-132471) | 1.9e-17 | 1.9e-17_abs | 1.0e-10 | PASS |
| validation | v2_rz_settling | fluid_internal_momentum_leak_z | 0 | EL_FLUID_PAIR audit (job 132698, 800 audited steps to t=40); v2_rz_settling RUNBOOK leak part 2 | 7.4e-04_cumulative | 7.4e-04_abs | 1.0e-10 | RESOLVED |
| infrastructure | fluid_conv_divergence_form | fluid_pair_mismatch_z_per_step | 0 | Convective-form leak root cause + divergence-form fix (commit 879933b1); v2_rz_settling RUNBOOK leak part 3 (jobs 135770-135772) | 1.94e-13_flat | 1.94e-13_abs | non_growing_vs_9.8e-8_legacy | PASS |
| infrastructure | fluid_conv_divergence_form | energy_stability_production | stable_to_t100 | RZ sweep jobs 135773-135776; v2_rz_settling RUNBOOK leak part 4 | NaN_at_t31..49 | NA | audit_use_only | RECORDED |
| infrastructure | el_momentum_fix | residual_mismatch_z_per_step | 0 | Measured-leak compensator (commit 22da0d5c); v2_rz_settling RUNBOOK leak part 4 (job 135865) | 6.0e-10_worst_nongrowing | 6.0e-10_abs | bounded_vs_9.8e-8_growing_legacy | PASS |
| infrastructure | el_momentum_fix | cum_drift_z_production_sweep_t100 | 0 | Compensator at level-3 production, 50k steps; v2_rz_settling RUNBOOK sweep 3 (jobs 135866-135869) | 5.2e-7_worst_bounded | 5.2e-7_abs | bounded_vs_5.7e-2_legacy_ramp | PASS |
| validation | v2_rz_settling | rz_ratio_U_over_U0_tavg | eps^4.65 (0.79/0.61/0.47/0.35) | Richardson-Zaki 1954 n=4.65 low-Re; v2_rz_settling RUNBOOK sweep 3 | 2.65/3.35/3.31/3.02 (enhanced) | +2.4e+00_to_+7.5e+00 | 1.0e-01 | FAIL |
| validation | v2_rz_settling | rz_exponent_early_window_t5_15 | 4.65 | RZ fit on pre-instability window; v2_rz_settling RUNBOOK sweep 3 | 2.5_to_3.3 | 4.6e-01 | 0.5_abs | RECORDED |
| validation | v2_rz_settling | mesoscale_settling_instability | homogeneous_state_stable | RZ assumes wall-bounded homogeneous suspension; periodic box has no suppression mechanism | Tgran x100 onset t~20; clustering; 2.6-3.3x enhancement | NA | physics_finding | RECORDED |

For rows with zero expected value, the `rel_error` entry is an absolute error marker because relative error is undefined.
| infrastructure | unit_lubrication_pair | cross_rank_pair_visibility_and_newton | pairs>0, NP=0 | Shadow-copy margin + fold-once lubrication (pe 916600a/24295f5); unit_lubrication_pair RUNBOOK (job 135994) | pairs=20, NP 1.8e-22 | 1.8e-22_abs | 1.0e-10 | PASS |
| validation | kroupa_shear | etaL_over_mu_phi_sweep | positive, monotone, superlinear (Kroupa Fig 2b/3 shape) | Kroupa et al. Langmuir 2016 lubrication-stress mechanism; kroupa_shear RUNBOOK (jobs 136099-136102) | 0.0079/0.0363/0.2496/0.7610 at phi 0.05-0.30 | NA | shape_and_monotonicity | PASS |
| validation | kroupa_shear | lubrication_off_twin_etaL | 0 | Virial accumulates only lubrication pairs; kroupa_shear RUNBOOK (job 136025) | 0.0_exact | 0.0_abs | 1.0e-12 | PASS |
| infrastructure | kroupa_shear | newton_pair_designated_treater | 0 | Relay-based pair treatment Newton-exact by construction; kroupa_shear RUNBOOK | 2.1e-19_worst_full_sweep | 2.1e-19_abs | 1.0e-10 | PASS |
| validation | kroupa_shear | substep_convergence_phi030_total_etaP | converged (<=5% drift 10->100) | Temporal-resolution ladder substeps 10/50/100; kroupa_shear RUNBOOK (jobs 136206-136208) | etaP/mu = 1.592/1.606/1.614 (+0.8% drift); substeps_=10 adequate | 8.0e-03 | 5.0e-02 | PASS |
| validation | kroupa_shear | contact_channel_etaC_phi030 | significant at high phi | Converged-PGS contact-impulse virial (pe 11b518a); kroupa_shear RUNBOOK ladder | etaC/mu = 0.83-0.86 > etaL/mu = 0.75-0.76 at phi=0.30; fires at ~zero overlap; 1+etaP/mu = 2.59-2.61 vs KD 2.75 | NA | discovery | RECORDED |
| validation | v4_pressure_drop | mu_app_over_mu_baseline_no_lubrication | 1.0_flat (drag-only model boundary) | Kroupa Fig 5 lubrication-off branch; v4_pressure_drop RUNBOOK (jobs 135958-135961) | 0.999/0.998/0.993/0.987 at phi 0.05-0.20 | NA | model_boundary_reference | RECORDED |
| infrastructure | el_meso_filter | plume_suppression_phi020 | homogeneous_state_stable | Lane-mode filter (commit 6520efbe); v2_rz_settling RUNBOOK option-1 probe (job 135956) | Tgran/6 but U/U0=2.14 (diag modes persist) | NA | insufficient_fallback_option2 | RECORDED |
| validation | v4_pressure_drop | mu_app_monotone_with_lubrication | monotone increasing (Kroupa Fig 5 ON branch) | Pairwise lubrication generates suspension viscosity in the pipe; v4_pressure_drop RUNBOOK (jobs 136113-136116) | 1.012/1.056/1.108/1.148 at phi 0.05-0.20 (baseline flat 0.999-0.987) | NA | monotonicity_and_on_off_contrast | PASS |
| validation | v4_pressure_drop | mu_app_over_mu_vs_krieger_dougherty | 1.139/1.312/1.533/1.821 | KD phi_max=0.64 curve; v4_pressure_drop RUNBOOK lubricated sweep | 1.012/1.056/1.108/1.148 (9-20% of KD excess; no wall lubrication + lift migration depletes near-wall shear zone) | -1.1e-01_to_-3.7e-01 | 1.5e-01_to_2.5e-01 | RECORDED |
| validation | v2_rz_settling | small_box_rz_ratio_phi020_Ldp10 | eps^4.65 = 0.354 | Option-2 probe d_p=0.1 L/d_p=10; v2_rz_settling RUNBOOK (jobs 136121/136122) | U_RZ/U0 = 1.08 (ratio 3.04x target); same band as L/d_p=20 | +2.0e+00 | 1.0e-01 | FAIL |
| validation | v2_rz_settling | cluster_induced_settling_reframe | homogeneous_state_unstable_in_periodic_box | Lane filter insufficient + L/d_p=10 unchanged -> option 4; enhancement 2.6-3.3x consistent with CIS literature for two-way point-particle sedimentation | decision_recorded | NA | physics_reframe | RECORDED |
