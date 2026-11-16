# s05_ergun_cfdem_twin - FeatFloWer + pe side of CFDEM case 05 (Ergun packed bed)

Files here: this RUNBOOK (chronological; later sections correct earlier ones where they say so),
decks `q2p1_param_E1.dat` (frozen periodic bed, spec deck), `q2p1_param_E2_adv.dat` /
`_E2_cons.dat` / `_E2_cons_pf_fb.dat` (column: advective form, conservative form, conservative form
with pressure force and feedback), `example_E1.json` / `example_E2.json` (pe), the settled column
bed `particles_E2_settled.xyz` with its provenance `particles_E2.xyz.json`, `tools/`
(post-processing and predictions), `runs_versions.txt` (tree, pe and binary of every run).
Not committed: meshes (scale_ff_mesh.py --scale 0.0277 on pipe_hp_check PIPEZ27; clean E1 bed
z scale 0.0267) and the E1 beds (make_bed_xyz.py of the peer spec; commands below). Run directories:
/data/warehouse17/rmuenste/runs/ergun05/<run>. Peer spec: LIGGGHTS-PUBLIC
md_docs/case_specs/05_ergun/SPEC.md. Summary of the outcome: sections of 2026-11-13.

## Stage E1: frozen periodic bed (2026-11-11)

Spec: LIGGGHTS-PUBLIC/md_docs/case_specs/05_ergun/SPEC.md sections 6.2/6.3, gates G-E1-0 and G5.
Binary: FeatFloWer-va/build-va q2p1_el_pipeflow (369d512b; EL sources identical to master 010f430b,
master differs only in solution_io / app_initialization / fbm_main / param_parser), pe df8c8a6.
Staging: `base/` = scaled PIPEZ27 mesh (scale_ff_mesh.py --scale 0.0277), spec deck
E1_frozen_periodic unchanged, bed ergun_tallbed/particles_E1.xyz (N 36 708), benchRadius_ =
radius_pe 4.93148e-4. Runs: 28 ranks, level 3 (82 944 elements), tools/el_slurm_submit_np.sh.

| Run | job | change vs base | wall | result |
|---|---|---|---|---|
| E1_smoke | 159316 | 10 steps, audit 1 | 1 min 23 s | G-E1-0 PASS: set-up completes, up = 0 exactly, drag 0.1643 N at step 10, no clip, no NaN |
| E1_fb200 | 159318 | f_b 200 | 12 min 35 s | uf_super_z 6.0447e-3 (51 samples), drag_z 6.5129e-2 N |
| E1_fb500 | 159319 | f_b 500 (deck as shipped) | 12 min 24 s | uf_super_z 1.49037e-2, drag_z 1.62824e-1 N |
| E1_fb1000 | 159320 | f_b 1000 | 15 s | FAIL: NaN in step 1 (start-up: no drag in the first step, the fluid accelerates freely) |
| E1_fb1000_dt5e-4 | 159325 | f_b 1000, dt 5e-4, 200 steps | 25 min | uf_super_z 2.95492e-2 (101 samples), drag_z 3.25648e-1 N |
| E1_fb500_outMAX | 159322 | OutputLevel MAX (fields incl. eps) | 12 min 28 s | all EL streams bitwise vs E1_fb500 |

Verified:
- Radius: eps deposition and drag both use the pe radius (benchRadius_): deposited volume
  1.844087e-5 m^3 = 36 708 V_p(d 0.98630 mm) to 3e-15; `EL_DRAG_CLOSURE` takes the diameter from the
  particle state. The fluid sees solid fraction 0.5532, eps_mean 0.446752 (not 0.4243).
- Mesh: V 3.33320e-5 m^3, cross-section 6.01661e-4 m^2 = 0.99839 of pi 0.01385^2. Clipped cells 0
  (eps min 0.294, floor 0.2).
- Force balance: drag_z / (rho f_b V) = 0.97697 / 0.97698 / 0.97698 at the three rungs; the missing
  2.3 % is wall friction on the no-slip hull. The balance is n F against rho f_b (unweighted
  momentum; `ELVolumeAveraged none`).
- Steady from about step 20; TAVG over the last half of each run.
- Bed structure seen by the fluid (E1_fb500_outMAX, tools/e1_series.py): the periodic slab has a
  nearly empty layer at the seam (five z layers, z = 0.3 and 53.6-55.1 mm, layer eps up to 0.89: the cut at 60.4 mm reaches the loose top
  of the 60.5 mm bed, as the profile in particles.xyz.json already shows, phi 0.48 and 0.09 in its
  last two bins); median layer eps 0.436; core
  (r < 12 mm) eps 0.4356 carrying 1.027 x the mean flux, wall ring eps 0.4791 carrying 0.921 x (no
  bypass: the no-slip wall and the kernel-smeared ring outweigh the looser packing).

G5 (machinery check), measured U_sup against three references (tools/e1_predict.py, e1_series.py;
coded difelice_std with d = 0.98630 mm):

| f_b | measured | spec table (d 1 mm, eps 0.451) | uniform bed at eps_mean 0.4468, d_pe | z-series on the actual eps(z), drag share 0.977 |
|---|---|---|---|---|
| 200 | 0.0060447 | 0.00707 (0.855) | 0.006578 (0.919) | 0.005966 (1.013) |
| 500 | 0.0149037 | 0.01743 (0.855) | 0.016216 (0.919) | 0.014711 (1.013) |
| 1000 | 0.0295492 | 0.03455 (0.855) | 0.032150 (0.919) | 0.029174 (1.013) |

- Linearity: measured / closure is constant to 1e-3 over the three rungs: PASS (2 %).
- Level: FAIL against the uniform-bed column (-8.1 %; gate 5 %); +1.3 % against the closure applied
  in series to the measured eps(z) with the measured drag share. The shortfall is the bed (seam gap
  raises eps_mean above the packed value; resistance is convex in eps), not the machinery.
- Ratio to the consistent model (n F = eps rho f_b) at the same uniform reference: 2.21 (0.007334 at
  f_b 500): the size of the missing eps weighting.

Inferred / open:
- The Q2 velocity has element-scale structure: mean u_z from level-3 vertices 0.03036, from all Q2
  nodes 0.03241, solver quadrature 0.03310. Field splits above are therefore relative numbers.
- The f_b 1000 start-up failure at dt 1e-3 needs a ramp or the smaller step; not investigated.
- Record committed 2026-11-14 on feature/el-ergun-column (this directory); the run directories stay
  in /data/warehouse17/rmuenste/runs/ergun05.

## Clean E1 bed (2026-11-12)

`bed_clean/` = make_bed_xyz.py --mode periodic --z0 0.005 --Lp 0.0534 (slab ends 2.1 mm below the
bed top): N 35 651, 408 dropped at the seam, pe solid fraction 0.5565. Mesh `mesh_e1c`
(scale 0.0277, zscale 0.0267), PeriodicLength 0.0534, OutputLevel MAX. Fluid sees eps_mean 0.442558,
V 3.212867e-5 m^3, 0 clipped cells, drag share 0.9770.

| f_b | job | measured U_sup | uniform at eps_mean, d_pe | same x drag share | z-series on eps(z) |
|---|---|---|---|---|---|
| 200 | 159816 | 0.0059315 | 0.006310 (0.940) | 0.006167 (0.962) | 0.005900 (1.005) |
| 500 | 159817 | 0.0146247 | 0.015556 (0.940) | 0.015203 (0.962) | 0.014549 (1.005) |
| 1000 (dt 5e-4) | 159818 | 0.0289952 | 0.030844 (0.940) | 0.030142 (0.962) | 0.028853 (1.005) |

A periodic splice of a random packing still leaves a loose seam (four z layers, layer eps up to
0.79; median 0.436). G5 against the uniform column: -6.0 % (was -8.1 %), -3.8 % after the measured
wall friction, +0.5 % against the bed-matched series reference. Linear to 1e-3.

## E2 column (2026-11-12): worktree FeatFloWer-ergun, feature/el-ergun-column 188bcead, pe feature/el-cylinder-wall 43e8fb5

Patches: P1 inlet (boundary type 975, SimPar@ELInletVelocityStart/End/RampTime, echo line, stop
when the mesh carries Inflow975 without the keys), P2 pe cylinder wall (cylinderWallEnabled_,
cylinderWallRadius_, contract test). Default path: ten certified EL configurations bitwise
(jobs 159838-159847, runs/coupling_p4p5/ergun_gates). Build build-el (build-va recipe).

| Run | job | what | result |
|---|---|---|---|
| E2_dry | 159827 | 40 steps, inlet 1e-12 | G0 echo present, wall reported, 10 000 contained, max overlap 0.12 % d, no NaN; bed still settling (top particle -0.3 d in 0.02 s): G-E2-0 height criterion not met on the raw bed |
| E2_settle | 159836 | 200 steps dry | mean z 7.4808 -> 7.4618 mm, drift 1 um per 0.01 s at the end, max overlap 0.08 % d; positions at 0.1 s = `bed_e2_settled/particles.xyz` (core solid fraction 0.583 with r_pe, H about 15.4 mm top / 14.9 mm equivalent) |
| E2_ramp0.4_pfNo | 159853 | ramp 0.002 -> 0.02 over 0.4 s, 800 steps, continuity (advective, uniform sink), no pressure force | 1 h 52 min; drag slope 3784 Pa s/m; fluid dp inlet-outlet 77.2 Pa = drag/A 76.6 Pa at U_in 0.02; no fluidisation (hydro force 0.46 of bed weight) |
| E2_ramp0.4_pfYes | 159854 | + ELPressureForce | drag slope 3781, drag + pressure force 5892 Pa s/m; 0.73 of bed weight at the end |
| E2_ramp0.4_pfYes_fbYes | 159855 | + ELPressureFeedback | 3710 / 8738 Pa s/m; reaches 0.997 of bed weight at U_in 0.02 (bed starts to lift, mean w 0.85 mm/s) |
| E2_ramp0.4_pfNo_none | 159875 | ELVolumeAveraged none | 1 h 51 min; u about U_in everywhere (div u = 0), superficial flux in the bed 0.36 U_in (vertex read); drag slope 2331 Pa s/m = 0.93 of n F(U) H 2490-2510 (slip = U, not U/eps); fluid dp 49.9 Pa vs drag/A 47.2 Pa at U 0.02; 0.28 of bed weight |
| E2_ramp0.4_pfNo_formC | 159876 | conservative form + BDF2 + corrector + alpha | REFUSED at start-up (abort 907): the O-grid pipe mesh has curved Q2 hull edges, the flux kernel's layout check accepts trilinear cells only (recorded limitation of form C). The guard worked; form C cannot run on this mesh until the curved-Q2 follow-up is done |

Finding (tools/e2_zprofile.py, frame at U_in 0.02, vertex-based so about 8 % low in absolute terms):
the superficial flux through the bed is only about 0.58-0.60 of U_in. Two causes:
1. The inlet prescribes the fluid velocity u = U_in on a plane where the kernel porosity is 0.59
   (the bed rests on the inlet), and u is the interstitial velocity under the continuity stage, so
   the superficial inflow is eps_inlet U_in. (CFDEM's lesson 5 in FeatFloWer form.)
2. The compatibility sink of the continuity stage is designed for closed or periodic boxes; with an
   open outlet it adds a uniform source over the whole column: above the bed u_z rises linearly from
   0.71 to 0.94 U_in.
Relative to the flux actually passing the bed the drag slope is 3784 / 0.60 = 6300 Pa s/m against
n F(U/eps) H = 6040-6100 for the settled bed (unweighted momentum), i.e. the machinery is consistent
within the vertex-read uncertainty; CFDEM R2 (improved set-up) is 12 129 Pa s/m, the consistent
closure 14 480.

## E2 with the inlet and sink fixes and the curved-cell extension (2026-11-13)

Code (worktree FeatFloWer-ergun, feature/el-ergun-column, nothing pushed):
- 6cc26f37 inlet fix: SimPar@ELInletSuperficial; the type-975 Dirichlet value is U_in / eps_h at
  the inlet nodes (eps_h = transport porosity, built in the transfer pass), so the superficial
  inflow is U_in. Sink fix: ELVACompatibility = none (no compatibility sink), required with
  NoOutflow = No under the continuity stage and refused otherwise.
- d3c7a172 (worktree FeatFloWer-formc-curved, feature/el-va-formc-curved), merged as b8baaf2f: the
  form C layout check accepts curved Q2 cells (mid-node nearest to its own reference position of
  the 19, within 0.2 h); unit tests on a curved mesh (cancellation, telescoping, moment identities,
  volume to rounding).
- Gates: ten certified EL configurations bitwise with both fixes (jobs 159962-159971,
  runs/coupling_p4p5/ergun_gates2); negative test 159973 refused as intended; curved extension
  T10a4 / VB_L2 / tier-2 bitwise (159974-159976). DNS twins NOT rerun (QuadSc_user.f90 and
  QuadSc_boundary.f90 are shared files): open before any merge.
- Binaries: build-el (pre-merge, job 159972), build-el2 (merged, job 159977).

Same ramp as above (0.002 -> 0.02 m/s over 0.4 s, dt 5e-4, 800 steps, settled bed, no pressure
force), ELInletSuperficial = Yes, ELVACompatibility = none. No NaN, no clipping, bed static.

| Run | job | continuity form | drag slope Pa s/m | flux leaving / U_in | hydro force / bed weight at U 0.02 |
|---|---|---|---|---|---|
| E2f_ramp0.4_pfNo | 159972 | advective (BE, no corrector) | 6351 | 1.195 (solver, domain mean) ; vertex read 1.02 in the bed, 1.21 above | 0.768 |
| E2f_ramp0.4_pfNo_formC | 159977 | conservative + BDF2 + 1 corrector | 5554 | 1.002 (EL_VA_FLUXMOM uf_eps_z 0.020021-0.020041 at U_in 0.02) ; vertex read 0.90 in the bed, 0.996 above | 0.675 |

Readings:
1. Both fixes work: the bed now receives the prescribed superficial flux (0.6 U_in before).
2. The advective form is not mass conservative at the porosity jump of the bed top: the flux rises
   from the bed value to 1.21 U_in within two element layers and stays flat to the outlet. With an
   open outlet nothing removes that error (in the closed s04 box the sink hid it). Its slope 6351 is
   therefore the slope at about 1.14 x the nominal flux (6351 / 5554 = 1.144; vertex flux ratio in
   the bed 1.02 / 0.90 = 1.13).
3. Form C on its first open-boundary run: flux = U_in to 0.1-0.2 % (solver quadrature), boundary
   flux residual Phi 5.5e-9 m3/s against an inflow of 1.2e-5 (4.6e-4), unsinked residual L1 2e-8,
   id_z 6.7e-7. The vertex-based profile reads 0.90 inside the bed and 0.996 above: the known
   element-scale Q2 structure (vertex mean below quadrature where eps varies), not a flux loss.
4. Form C slope 5554 Pa s/m = 0.91-0.92 of the uniform closure n F(U/eps) H (6040-6100), the same
   ratio as E1 measured/uniform (0.919, d_pe): column and periodic bed agree. Fluid dp across the
   bed 110.2 Pa at U_in 0.01912 (5764 Pa s/m incl. wall friction and end layers).
5. Against CFDEM R2 (12 129 Pa s/m) the E2 slope is 0.458; the gap is the momentum formulation
   (unweighted momentum, drag on slip U/eps without the pressure-gradient share), not flux or
   continuity. Consistent closure 14 480. G-E2 level gate vs CFDEM: FAIL (recorded); linearity
   and flux balance: PASS with form C.

Open: pressure-force / feedback rungs with the fixes and form C; porosity-weighted momentum stage;
DNS twins for the branch; pe PR for feature/el-cylinder-wall; home of case 05 in the repository.

## Why CFDEM R2 is 0.838 of the consistent closure, and the like-for-like FF prediction (2026-11-13)

Tool: `tools/bed_structure_predict.py <xyz> --d D [--eps-core E --H H]` (spheres split into sample
points, z layers of 1.065 mm = CFDEM dz, core ring r < 12 mm and wall ring; Di Felice chi(Re) at
U 0.005; per-bin gradient n F(U/eps) [unweighted] or n F / eps [consistent]; layers in series, the
two rings of a layer in parallel). It reproduces the peer's uniform reference exactly (14 472).

| Bed | uniform reference | z series | + rings parallel | measured |
|---|---|---|---|---|
| CFDEM bed_column.xyz, d 1 mm, eps 0.416, H 14.88 mm, consistent | 14 472 | 13 254 (0.916) | 12 806 (0.885) | 12 129 (0.838) |
| FF bed_e2_settled, d_pe 0.99390 mm, eps 0.4176, H 14.65 mm, unweighted | 5 896 | 5 488 (0.931) | 5 372 (0.911) | 5 554 form C (0.942; 1.034 of the structured value) |
| same, consistent | 14 116 | 12 990 (0.920) | 12 561 (0.890) | rungs 159990 / 159991 running |

- The 16 % belongs to the bed and the reference, not to CFDEM: floor and top layers have eps 0.47
  and 0.50 against 0.416 in the core (-8 %), the wall ring has 0.45 against 0.40 (-3 %), and the
  rest (about -5 %) is cell-scale eps scatter (sd 0.025 per layer in the peer's run; the closure
  has d ln G / d eps = -10.6, so cells in parallel give 1 - 0.5 * 130 * sd^2 = 0.96). The peer's
  layer check (SPEC 12.2: measured 11.91 against own-Ksl series 11.94 layer units) says the same.
- Correction of this RUNBOOK: the E2 particles have d_pe 0.99390 mm (benchRadius_ 4.96948e-4, bed
  eps 0.4176), not the E1 value 0.98630. "6040-6100" above used eps 0.416-0.417 with H 14.9; the
  like-for-like uniform value is 5 896 and form C is 0.942 of it (not 0.91-0.92), 1.034 of the
  structured prediction (the 2.5 d kernel smooths the cell-scale scatter and part of the wall ring).
- Missing piece in FF: the fluid sees -grad p, model A needs -eps grad p. For a fixed bed
  ELPressureForce + ELPressureFeedback supplies exactly that (particle -V grad p, reaction
  +phi grad p deposited on the fluid => eps |grad p| = n F). Pre-fix evidence: 3710 -> 8738 =
  2.355, predicted (1/eps) x (0.890/0.911) = 2.34. Not implemented and not needed here (Re 0.1,
  fixed bed): eps in the inertia and viscous terms.
- Forecast for 159991 (pf + feedback, form C): 12 561 x 1.034 = about 13 000 Pa s/m (1.07 x CFDEM).

## E2 pressure-force and feedback rungs with form C (2026-11-13)

Clones of E2f_ramp0.4_pfNo_formC with only the named keys changed; binary build-el2 (b8baaf2f).
No NaN, no clipping. Slopes over U in (0.003, 0.0075).

| Run | job | keys | drag slope | drag + pressure force slope | hydro force / bed weight | fluid dp at U_in 0.01912 |
|---|---|---|---|---|---|---|
| E2f_ramp0.4_pfNo_formC | 159977 | - | 5554 | 5554 | 0.675 at U 0.02 | 110.2 Pa |
| E2f_ramp0.4_pfYes_formC | 159990 | ELPressureForce | 5553 | 8643 | 0.998 at U 0.02 (bed starts to lift, w 1.2 mm/s) | 104.9 Pa |
| E2f_ramp0.4_pfYes_fbYes_formC | 159991 | + ELPressureFeedback | 5450 | **12 686** | reaches 1.0 at U 0.01375 (t 0.261), max 1.019 | 170.5 Pa (bed weight/A 166.5) |

- Feedback rung against CFDEM R2 (12 129 Pa s/m): **1.046**. Against the structured consistent
  prediction for this bed (12 561): 1.010; forecast was 13 000. Ratio to the drag-only rung 2.284.
- U_mf: force first reaches the bed weight at 0.01375 m/s; the fixed-bed line reaches it at
  166.53 / 12 686 = 0.0131 m/s (CFDEM R2 line: 170.85 / 12 129 = 0.0141; -7 %). Plateau 1.02 x
  bed weight; fluid dp across the column 170.5 Pa = 1.024 x bed weight / A.
- After U_mf the bed lifts from the floor (eps 0.47 in the lowest 4 mm at t 0.38, H 15.44 ->
  15.92 mm, mean w 4 mm/s at the end); the form C boundary flux residual grows to 6e-7 m3/s and the
  mean eps u_z drops to 0.93 U_in: the bed is moving, particle volume flux is part of the balance.
- Inconsistent set (pf without feedback): predicted n F (1 + phi) = 5553 x 1.582 = 8787,
  measured 8643 (0.984).
- G6 (spec section 8: R2 +-10 % slope, +-15 % U_mf, plateau 1 +- 0.03): values inside all three
  tolerances with pressure force + feedback + form C. The spec ties a formal PASS to "the momentum
  stage"; here the -eps grad p part is supplied by the feedback key, inertia and viscous terms
  remain unweighted. Verdict wording left to the owner and the peer.

## DNS twins on feature/el-ergun-column (2026-11-13)

Tree 46fb4c55 (= b8baaf2f + CRLF restoration of QuadSc_boundary.f90, content identical), pe 43e8fb5,
new build tree FeatFloWer-ergun/build-dns (certified 3-step recipe, gcc 13.2, HardContactAndFluid,
serial PE, FBM acceleration; build_twin_dns.sh). Rundirs in the DNS clone:
q2p1_dns_rundir_twin_ergun_{e4l3,dkt}, inputs identical to the master twins.

| Twin | job | result vs certified head |
|---|---|---|
| e4_l3 | 160029 | particle_force.log ca98e7af225b, bulk_flow.log 89abce88b675, seven log streams: bitwise |
| DKT | 160030 | particle_force.log 82411f64f661, bulk_flow.log cfb8b42d0969, five log streams: bitwise |

The EL gates (ergun_gates2) ran on b8baaf2f binaries; 46fb4c55 differs only in line endings.

## Momentum models B / Bfull / A on the E2 column (2026-11-15)

Branch feature/el-momentum-models, final commit 08f7d420 (binary FeatFloWer-mommodels/build-el-e,
versions.txt dirty 0), pe 43e8fb5. Deck `q2p1_param_E2_model.dat` = `q2p1_param_E2_cons_pf_fb.dat`
with the eight coupling keys replaced by the single key `SimPar@ELMomentumModel` (the model
implies them; shown value Bfull, the B and A runs differ in this one value). Same bed, ramp,
mesh and partition as the feedback rung (job 159991). 800 steps each, no NaN, no clipping. Slopes
over U in (0.003, 0.0075), bed weight / A = 166.53 Pa, CFDEM R2 slope 12 129 Pa s/m.

| Run | job | wall time | drag slope | total slope (drag + pressure force) | / CFDEM R2 | / structured prediction 12 561 | force reaches bed weight at | plateau max | fluid dp at the end |
|---|---|---|---|---|---|---|---|---|---|
| E2f pf + fb, no model (reference rung) | 159991 | - | 5450 | 12 686 | 1.046 | 1.010 | 0.01375 m/s | 1.019 | 170.5 Pa |
| E2n_ramp0.4_B | 160559 | 1 h 51 | 13 069 | 13 069 | 1.077 | 1.040 | 0.01449 m/s | 1.014 | 170.4 Pa |
| E2n_ramp0.4_Bfull | 160560 | 2 h 57 | 5412 | 13 301 | 1.097 | 1.059 | 0.01390 m/s | 1.018 | 170.2 Pa |
| E2n_ramp0.4_A | 160561 | 2 h 57 | 5465 | 12 731 | 1.050 | 1.014 | 0.01370 m/s | 1.017 | 170.2 Pa |

- All three models carry the bed: slope within +5 ... +10 % of CFDEM R2 (spec G6 +-10 %), force
  reaches the bed weight at 0.0137 ... 0.0145 m/s (CFDEM R2 line 0.0141; -3 ... +3 %, spec +-15 %),
  plateau 1.014 ... 1.018 x bed weight (spec 1 +- 0.03), fluid dp across the column 1.022 ... 1.023
  x bed weight / A.
- Spread between the models in the fitted slope: 4.5 % (A lowest, Bfull highest). The fit window
  lies in the ramp and the lines have different intercepts (B +0.1, Bfull -5.7, A +2.6 Pa); total
  dp at fixed U: 65.3 / 61.4 / 66.4 Pa at U 0.005 and 131.5 / 124.8 / 129.5 Pa at U 0.01
  (B / Bfull / A). Chord slopes 0.005 -> 0.01: 13 242 / 12 676 / 12 624. The models agree in the
  steady fixed-bed limit only in the continuum; here B takes drag / eps with the local eps, Bfull
  and A add V_p grad p from the sampled pressure gradient. The spread is not broken down further.
- Viscous particle force (Bfull, A; EL_MM_FORCE): total z -1.7e-5 N at the end, negligible against
  drag and pressure force in this case.
- Flux: superficial fluid flux / U_in 1.003 at U 0.005 and 1.002 at U 0.01 in all three. After
  the bed starts to move the domain mean of eps u_z falls below U_in (particle volume flux):
  B 0.970, Bfull 0.970, A 0.921 at U 0.02 (reference rung 0.934). **Open:** with models B and
  Bfull the audit values of predictor and corrector passes agree to 5e-4 at the end; with model A
  and with the no-model rung they alternate (1.840e-2 / 2.102e-2 / 1.843e-2 m/s at step 800 for A),
  i.e. the source corrector does not converge on the moving bed where the residual divisor is eps
  (none, A), and does where it is 1 (B, Bfull). The fixed-bed part of the curve is not affected
  (values agree to 3 digits up to U 0.011). Not investigated further; a moving-bed statement with
  model A needs more correctors or the divisor looked at.
- Runs of the code before the two reviews (E2m_ramp0.4_{B,Bfull,A}, jobs 160459 / 160464 / 160465):
  slopes 13 064 / 13 301 / 12 730, i.e. the review fixes changed nothing here beyond 0.04 %.
- Cost: Bfull and A 1.6 x the wall time of B (per-step viscous integration and, for A, the
  operator terms); every model rebuilds the weighted mass matrix once per step.
- G6 verdict: the momentum stage now exists, all three models are inside the three spec
  tolerances. Owner verdict 2026-11-15: PASS (the three model rows of the datasheet; the rows of the
  feedback rung without a model, job 159991, stay RECORDED).

## Model A with the weighted pressure matrix (2026-11-15)

Follow-up of the open point above (the continuity iteration of model A alternates on the moving
bed). Branch feature/el-va-weighted-pressure, commits bf2fd6b0 (key) and f06d70d6 (review
fixes); binary FeatFloWer-pweight/build-el-d (f06d70d6, dirty 0). Mechanism, key and reviews:
`docs/md_docs/el_momentum_models_2026-11-15.md` section 7.

**Probes** (model A, constant-rate ramp from U 0.016 so that the bed moves within 120 steps;
first-round binary bf2fd6b0; values at step 120):

| Pressure matrix | correctors | job | net boundary flux of eps u after the passes, m3/s | residual L1 at the end of the step | mean eps u_z over the passes, m/s | wall time |
|---|---|---|---|---|---|---|
| start-up (divisor eps) | 1 | 160576 | -1.9e-7 / +3.4e-7 / -2.2e-7 | 2.3e-7 | 0.01593 / 0.01658 / 0.01589 | 26:37 |
| start-up (divisor eps) | 4 | 160577 | +2.8e-7 / -2.1e-7 / +8.5e-8 / -2.6e-8 / +6.9e-9 | 7.0e-9 | settles at 0.01615 | 29:31 |
| weighted (divisor 1) | 1 | 160578 | 7.8e-15 / 2.4e-17 | 3.5e-9 | 0.016155 / 0.016149 / 0.016149 | 28:47 |
| weighted (divisor 1) | 2 | 160579 | 7.4e-15 / 1.9e-17 / -1.3e-18 | 3.8e-10 | 0.016154 | (node shared) |

The start-up matrix over-corrects the through-flow with alternating sign (factors about -0.7,
-0.4, -0.3 per pass): it converges, with four correctors. The weighted matrix needs one. Drag and
pressure force at step 120 agree to 0.4 % between all four.

**Full ramp**, deck `q2p1_param_E2_model_A_weighted.dat` (= the model deck with
`SimPar@ELMMPressureOperator = weighted`, one corrector), job 160618, 800 steps, 3 h 02 min,
no NaN:

| Run | job | drag slope | total slope | / CFDEM R2 | force reaches bed weight at | plateau max | eps u_z / U_in at U 0.02 | boundary flux at step 800, m3/s | residual L1 at step 800 |
|---|---|---|---|---|---|---|---|---|---|
| E2n_ramp0.4_A (start-up matrix) | 160561 | 5465 | 12 731 | 1.050 | 0.01370 m/s | 1.017 | 0.921 | -8.1e-7 (alternating) | 9.5e-7 |
| E2p_ramp0.4_A_w2 (weighted) | 160618 | 5463 | 12 718 | 1.049 | 0.01370 m/s | 1.016 | 0.970 | 5.0e-17 | 4.6e-9 |

- The fixed-bed slope, the onset and the plateau are unchanged (slope -0.1 %); the earlier model A
  numbers stand.
- On the moving bed the flux now agrees with models B and Bfull (0.970 each) and the three audit
  values of a step agree to 3e-5 (1.93901e-2 / 1.93903e-2 / 1.93899e-2 m/s at step 800). The open
  point is closed for model A with the key; without a model the start-up matrix with more
  correctors remains the way.
- Cost: wall time 3 h 02 against 2 h 57 (+3 %); the solver statistics give 0.3 s per step for the
  matrix rebuild. A first copy of the run (job 160615) ran at 21 s per step on node nx-02
  against 13.6 s elsewhere and hit its 4 h limit at step 680; the cause is the node (the time is
  outside the solver statistics, same deck and binary), its numbers up to there agree.
## Reruns with the dynamic-viscosity default (2026-11-15/16)

`SimPar@ELFluidViscosity` defaults to `dynamic` since 5a882280 (branch feature/el-fluid-viscosity;
case s07_pipe_viscosity): the fluid's viscous operator carries mu = rho nu = 1.5e-3 Pa s instead
of nu = 1.5e-4 (rho_f = 10 here). Same decks as the recorded runs, no key added; binary
FeatFloWer-visc/build-el-g (5a882280, dirty 0). 800 steps each, no NaN. Every earlier number of
this RUNBOOK belongs to the `kinematic` convention and is reproduced with that key.

| Run (E2v_ramp0.4_...) | job | total slope, Pa s/m | before (kinematic) | change | / CFDEM R2 | force reaches bed weight at | plateau max | fluid dp at U_in 0.01912 |
|---|---|---|---|---|---|---|---|---|
| pfNo_formC (drag only) | 160701 | 5577 | 5554 | +0.4 % | 0.460 | not reached (0.677 at U 0.02) | - | 111.1 Pa (110.2) |
| pfYes_formC (+ pressure force) | 160702 | 8682 (drag 5570) | 8643 | +0.5 % | 0.716 | not reached (0.997 at U 0.02) | - | 105.4 Pa (104.9) |
| pfYes_fbYes_formC (+ feedback) | 160703 | 12 791 (drag 5478) | 12 686 | +0.8 % | 1.055 | 0.01366 m/s (0.01375) | 1.017 (1.019) | 170.7 Pa (170.5) |
| B | 160696 | 13 103 | 13 069 | +0.3 % | 1.080 | 0.01469 m/s (0.01449) | 1.014 (1.014) | 170.8 Pa (170.4) |
| Bfull | 160697 | 13 378 (drag 5436) | 13 301 | +0.6 % | **1.103** | 0.01384 m/s (0.01390) | 1.019 (1.018) | 170.5 Pa (170.2) |
| A (start-up matrix, 1 corrector) | 160698 | 12 782 (drag 5478) | 12 731 | +0.4 % | 1.054 | 0.01363 m/s (0.01370) | 1.017 (1.017) | 170.5 Pa (170.2) |

- The bed pressure drop is drag-dominated: the ten times larger fluid viscosity raises every
  slope by 0.3 ... 0.8 % (wall shear and the velocity profile in and above the bed).
- **Spec G6 with the corrected operator:** B and A stay inside the +-10 % slope tolerance; Bfull
  is at +10.3 % and no longer inside it (it was at +9.7 %). Onset (+-15 %: -3 ... +4 %) and
  plateau (1 +- 0.03) are met by all three. The datasheet rows of the model gates are restated
  below for this convention; the owner's PASS of 2026-11-15 was given on the kinematic numbers.
- Particle viscous force (Bfull / A): -1.6e-4 N at the end (was -1.7e-5 N; factor 9.3 for the
  factor 10 in the coefficient), still 0.2 % of the bed weight (0.100 N).
- Model A on the moving bed with the start-up pressure matrix still alternates (flux 0.920 U_in
  at U 0.02 against 0.970 for B and Bfull): unchanged by the viscosity; the weighted pressure
  matrix (branch feature/el-va-weighted-pressure) is the cure, its rerun with this default is
  recorded on the integration branch.
- Slope definitions (noted by the artifact work, 2026-11-15): the FeatFloWer slopes here are
  force-based with an intercept, CFDEM's 12 129 is the inlet-outlet pressure difference fitted
  through the origin; CFDEM force-based is 12 009, on which basis all ratios are 1 % higher.

## Model A with the weighted pressure matrix and the dynamic-viscosity default (2026-11-16)

Both changes together, on the integration branch integration/el-pweight-visc (merge of
feature/el-va-weighted-pressure and feature/el-fluid-viscosity; binary FeatFloWer-int/build-el-b,
ed573db0, dirty 0; the later merges of this branch add records only). Deck
`q2p1_param_E2_model_A_weighted.dat` unchanged (no viscosity key: default `dynamic`). 10-step
probe job 160706, full ramp job 160707, 800 steps, 3 h 01 min, no NaN.

| Run | job | viscosity | pressure matrix | drag slope | total slope | / CFDEM R2 | force reaches bed weight at | plateau max | eps u_z / U_in at U 0.02 | boundary flux at step 800, m3/s |
|---|---|---|---|---|---|---|---|---|---|---|
| E2n_ramp0.4_A | 160561 | kinematic | start-up | 5465 | 12 731 | 1.050 | 0.01370 m/s | 1.017 | 0.921 | -8.1e-7 (alternating) |
| E2p_ramp0.4_A_w2 | 160618 | kinematic | weighted | 5463 | 12 718 | 1.049 | 0.01370 m/s | 1.016 | 0.970 | 5.0e-17 |
| E2v_ramp0.4_A | 160698 | dynamic | start-up | 5478 | 12 782 | 1.054 | 0.01363 m/s | 1.017 | 0.920 | (alternating) |
| E2v_ramp0.4_A_w | 160707 | dynamic | weighted | 5476 | 12 768 | 1.053 | 0.01366 m/s | 1.017 | 0.969 | 3.7e-17 |

- The two changes act independently: the viscosity raises the slope by 0.4 % with either matrix,
  the weighted matrix lowers it by 0.1 % with either viscosity and removes the alternation on the
  moving bed (residual 4.4e-9 at step 800, audit values of a step 1.93913e-2 / 1.93815e-2 /
  1.93811e-2 m/s).
- Model A in its intended configuration (dynamic viscosity, weighted matrix): slope 1.053 of
  CFDEM R2, onset 0.01366 m/s (-3 %), plateau 1.017, fluid dp 170.4 Pa: inside spec G6.
