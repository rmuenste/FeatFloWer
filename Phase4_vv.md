 Phase 4: Deep Verification & Validation of the Euler-Lagrange Scheme

 Context

 Phases 1-3 built the unresolved Euler-Lagrange (CFD-DEM) coupling: kernel
 particle-mesh transfer (P1), one-way particle forces with a semi-implicit PE drag
 update (P2), and explicit dilute two-way fluid feedback (P3). The tests written so
 far are shallow — they check plumbing (sign, Σb_a = −F, no-op-when-off,
 exact tstep scaling) on a synthetic per-rank mesh. They do not verify the
 scheme: no convergence order, no analytic-solution comparison, no global
 conservation law, and (by the harness's own admission) no real partitioned MPI
 path. A real bug already slipped through this gap (the Phase-3 feedback
 double-count on shared DOFs).

 Phase 4 adds no new physics. It builds a verification suite that gives true
 verification: measured convergence order, analytic-benchmark agreement, and exact
 conservation laws. Two outcomes matter as much as "tests pass": (a) the suite will
 quantify the scheme's conservation error and so expose the known Phase-3
 explicit-vs-implicit drag inconsistency (issue D), turning a hand-wave into a
 number; (b) it finally exercises the real partitioned path, closing the standing
 Phase-1 distributed-straddling TODO.

 Two tiers

 Tier 1 — standalone deep tests (ff_quadLS_app, synthetic mesh, no PE, fast,
 deterministic, CI-friendly). Verify the discrete operators and closures:
 convergence order and analytic exactness.

 Tier 2 — runtime end-to-end tests (q2p1_el_pipeflow + PE, real partitioned
 mesh, via the tools/featflower_test YAML harness). Verify the coupled scheme:
 global conservation, terminal velocity, lift, and the distributed path. References
 are analytic (drag law, conserved momentum, Saffman), not prior output — this
 is verification, not regression.

 ---
 Tier 1: standalone convergence / analytic tests

 Reuse the parameterizable structured mesh already in
 source/src_el/tests/test_el_transfer.f90 (build_structured_mesh(mesh, ncell),
 build_dof_coordinates) — sweep ncell ∈ {4,8,16,32} for refinement studies. The
 operators under test are EL_INTEGRATE_PARTICLE / EL_DEPOSIT_PARTICLE
 (source/src_el/el_quadrature.f90, with the EL_SAMPLE_* offsets) and the
 closures in source/src_el/el_forces.f90.

 New executable source/src_el/tests/test_el_convergence.f90 (a small MMS/EOC
 harness; no existing MMS infrastructure, build from scratch), registered like
 the others in cmake/modules/ProjectFiles.cmake (add_test
 el-convergence-serial + a helper to compute observed order = log2(e_h/e_{h/2})).

 1. G2P interpolation order (MMS). Impose a smooth manufactured field
 (trigonometric u = sin(2πx)…, and a separate polynomial field) on the Q2
 DOFs; sample at a particle via EL_INTEGRATE_PARTICLE. Targets: (a) constant
 and linear/quadratic fields reproduced exactly (extends current spot checks
 to a proper assertion); (b) for the smooth field, kernel-weighted value error
 → 0 at the expected order as h→0 at fixed kernel/h ratio; (c) velocity-
 gradient and pressure-gradient sampling order (EL_SAMPLE_GRAD_U/_P).
 2. Kernel-width study. Hold mesh fixed, vary el_kernel_width_factor; verify
 the documented behavior (smoothing vs resolution trade-off) and that the
 normalization stays exact (sample(1)>0, partition-of-unity).
 3. P2G spreading order & moments. Deposit a force; assert zeroth moment exact
 (Σb_a = −F, already covered — fold in), and first moment
 Σ_a b_a⊗x_a = −F⊗x_p to O(h^p) (new — verifies the spread is centered, not
 just conservative).
 4. Gather/scatter adjoint consistency. Verify EL_DEPOSIT is the weighted
 transpose of EL_INTEGRATE for the same kernel/particle (structural property
 that guarantees discrete conservation; catches asymmetry regressions).
 5. Force-law analytic checks (extend test_el_kernel_forces.f90): Di Felice /
 Stokes / Schiller-Naumann drag vs literature C_d(Re) across a Re sweep;
 Saffman-Mei lift vs the analytic value in a prescribed linear shear
 (magnitude and direction); the ε^(2−χ) voidage exponent sensitivity
 documented in el_forces.f90.
 6. Semi-implicit drag ODE (extend libs/pe/tests/interface/pe_el_semi_implicit_drag_test.cpp):
 already verifies exact-exponential convergence and stiff stability; add the
 F_other≠0 equilibrium and a measured first-order temporal EOC.

 Tier 2: runtime end-to-end tests

 New conserved-quantity diagnostic (new module
 source/src_el/el_diagnostics.f90, or a routine in el_transfer.f90) that each
 step computes and prints keyword lines for the harness to parse — mirror the
 existing EL step … print (el_transfer.f90, gated on showid /
 el_write_diagnostics) and the SED_BENCH_VEL pattern the sedimentation YAML
 already consumes:

 - Global fluid momentum ∫ρu dV: element-cubature loop reusing the E013 /
 domega·|detJ| pattern from EL_INTEGRATE_PARTICLE and mesh%level%dvol;
 MPI_Allreduce over MPI_COMM_SUBS. (No pre-assembled Q2 mass matrix is
 exposed — integrate directly.)
 - Global particle momentum Σ m_i v_i: from getAllParticles
 (source/src_particles/dem_query.f90) + MPI_Allreduce.
 - Print e.g. EL_MOMENTUM time= … fx fy fz px py pz tx ty tz, EL_TERMINAL_VEL ….

 Hook the call into the EL step (after EL_ADVANCE_PARTICLES) or
 postprocessing_app (source/postprocessing/solution_io.f90:1416).

 Case inputs live next to the app (applications/q2p1_el_pipeflow/_data/): reuse
 example.json (set gravity_, particleDensity_, fluidDensity_,
 benchRadius_) + q2p1_param.dat. Each case gets a tools/featflower_test
 YAML definition modeled on q2p1_bench_sedimentation_pe_serial.yaml
 (keyword_columns parser, fortran_d_or_e, compare: tolerance) with an
 analytic baseline.

 7. Global momentum conservation — the Newton-pair test (highest value).
 Periodic box (3D periodicity is hard-set in
 applications/q2p1_dns_drag/app_init.f90), no gravity, no body force, two-way
 ON, particle(s) given an initial slip. Assert M_fluid + M_particles constant
 over many steps. This is the deepest true test of the coupling and will
 measure the Phase-3 explicit/implicit drag drift; set the tolerance from the
 expected O(Δt) mismatch and flag if it exceeds a physical bound (motivates the
 future "spread the implicit reaction" fix).
 8. Single-particle terminal velocity. Gravity on, one particle; assert
 terminal velocity equals the drag-law balance (net submerged weight
 (ρ_p−ρ_f)Vg = drag) and the transient follows u(t)=u_t(1−e^{−t/τ_p}).
 Directly mirrors the sedimentation-benchmark metric pattern.
 9. Single particle in linear shear — Saffman lift. Prescribe a shear inflow;
 assert cross-stream lift / migration matches Saffman-Mei. Live-scheme analog of
 Tier-1 #5.
 10. Distributed-straddling conservation (closes the Phase-1 TODO and validates
 the Phase-3 fix end-to-end). Real partitioned mesh, particle straddling rank
 boundaries; assert global volume Σ|K|α_p = ΣV_i and feedback Σb = −ΣF
 across the partition (the real-MPI version the synthetic harness cannot do).
 This is the production-ordering regression for capture-before-E013Sum3 that
 test_partitioned_feedback_no_double_count explicitly cannot cover.

 ---
 Critical files

 Reuse: source/src_el/el_quadrature.f90 (operators + EL_SAMPLE_*),
 source/src_el/el_forces.f90 (closures), source/src_el/tests/test_el_transfer.f90
 (build_structured_mesh, build_dof_coordinates, check),
 source/src_particles/dem_query.f90 (getAllParticles),
 tools/featflower_test/testcases/definitions/q2p1_bench_sedimentation_pe_serial.yaml
 (YAML pattern), applications/q2p1_el_frozen_trace/example.json (case config),
 cmake/modules/ProjectFiles.cmake (add_test).

 New: source/src_el/tests/test_el_convergence.f90 (Tier-1 MMS/EOC);
 source/src_el/el_diagnostics.f90 (global momentum/terminal-velocity printouts);
 extensions to test_el_kernel_forces.f90 and pe_el_semi_implicit_drag_test.cpp;
 per-case _data/ inputs under applications/q2p1_el_pipeflow/; new
 tools/featflower_test/testcases/definitions/*.yaml + baselines/*.yaml for the
 runtime cases; add_test entries in ProjectFiles.cmake.

 Sequencing

 1. Tier-1 convergence harness (test_el_convergence.f90) + force-law/lift and
 semi-implicit extensions — fast, no PE, immediate deep value.
 2. el_diagnostics.f90 momentum/terminal-velocity printouts + unit coverage.
 3. Runtime test #7 (momentum conservation) — the single most valuable true test;
 quantifies the conservation error.
 4. Runtime #8 (terminal velocity), #9 (Saffman), #10 (straddling).

 Verification (how to test the tests)

 - Tier 1: ninja test_el_convergence test_el_kernel_forces then
 ctest -R el- in build/ (local USE_PE=OFF Zen5 build already configured);
 observed convergence orders must match expected (≈3 for Q2 value, ≈2 for
 gradients) within a tolerance band; force-law/Saffman within a few %.
 - Tier 1 PE math: g++ … pe_el_semi_implicit_drag_test.cpp (header-only) or
 via the PE build's pe-el-semi-implicit-drag CTest.
 - Tier 2: build with USE_PE=ON and run through tools/featflower_test
 (validate → run → compare → report) per
 docs/md_docs/featflower_test_usage_guide.md; momentum-conservation drift,
 terminal velocity, and Saffman lift compared to analytic baselines within
 tolerance (likely on the remote, since PE-coupled runs are heavier).
