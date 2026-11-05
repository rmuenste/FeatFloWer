# Plan: Finish Phase-4 Tier-2 V&V (WP-A) + pay down EL physics debt (WP-B)

## Context

The Euler-Lagrange coupling on `feature/euler-lagrange-phase1` has Phases 1–3 done and
the Tier-2 momentum-conservation test green (drift_rel 4.3e-7, np=28). Before the
validation campaign (Richardson–Zaki, Segré–Silberberg, pipe pressure drop) can produce
trustworthy numbers, two work packages remain:

- **WP-A**: actually run the three remaining Tier-2 verification cases + a long-horizon
  momentum run. Exploration found concrete blockers: wrong `mpi_ranks` in the terminal
  YAML; saffman/straddling reference a nonexistent mesh; the saffman case is physically
  degenerate (no shear anywhere ⇒ lift ≡ 0 ⇒ verifies nothing).
- **WP-B**: three physics-debt items that would invalidate validation results:
  (1) Di Felice drag convention is provisional and inconsistent with the sampled
  (interstitial) velocity; (2) the lift "wall correction" is a placeholder (ε_f used as
  wall factor); (3) drag feedback to the fluid is explicit-only — unstable for stiff
  drag at φ→20%.

User decisions (defaults chosen after no response, revisitable):
**Saffman case → frozen analytic linear-shear field. Di Felice → plan-doc §1.2 spec
exactly.**

Work environment: primary checkout `/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer`;
runs use `/home/user/rmuenste/nobackup/code/FF-EL/FeatFloWer/build-el-phase2-pe-gcc14`
(USE_PE=ON). Run recipe precedent:
`applications/q2p1_el_pipeflow/tier2_cases/momentum_conservation/RUNBOOK.md`
(mpirun resolved from binary rpath, `--wdir`, np = subdomains+1). One commit per item.

---

> **STATUS: WP-A1–A4 are DONE and committed** (44b3df0b terminal, c202ccf5
> straddling, f5ab6a16 momentum-long, 999d72a5 saffman) with passing harness runs
> (terminal vz=−0.03556384; straddling residual 5.2e-26, record_ranks=2;
> momentum-long drift_rel=1.000013e-5; saffman lift_z=1.27677e-6). The A-sections
> below are kept for reference only. Remaining work = WP-B (B1 → B2 → B3) plus the
> WP-A follow-ups in review amendment 4.

## WP-A1 — Terminal velocity: fix ranks, run, record

- `tools/featflower_test/testcases/definitions/q2p1_el_pipeflow_tier2_terminal.yaml`:
  `mpi_ranks: 28 → 4` (and slurm ntasks). QBOX9Z3 mesh = 3 subs; PE MPI_Aborts unless
  processesX·Y·Z == commsize (`example.json` has processesZ_=3).
- Stage via cmake target `q2p1_el_pipeflow_tier2_terminal_velocity_stage`, run np=4
  per RUNBOOK pattern, 120 steps.
- Accept: last `EL_TERMINAL_VEL` vz within ±2e-3 of committed analytic baseline
  −3.556e-2 (u_t = (2/9)(ρp−ρf)gr²/μ). Record measured value in commit message.

## WP-A2 — Straddling conservation: re-point at QBOX9, fix particle, run

- `tier2_cases/straddling_conservation/q2p1_param.dat`: `MeshFolder="QBOX9"`,
  `SubMeshNumber=27`, project file `quiescent_box_9/file.prj`; copy `Periodic=Yes`,
  `NoOutflow=Yes` from the momentum case (QBOX9 faces are Periodic-tagged).
- `example.json`: `processesX_/Y_/Z_ = 3/3/3` (np stays 28, YAML already correct).
- `particles.xyz`: ONE particle at `0.32733 0.5273 0.6191` — center clearly owned
  (6.7e-3 off the x=1/3 rank boundary, avoiding PE ownsPoint face ambiguity) while the
  kernel support [0.3148, 0.3398] (half-width 2.5·r = 0.0125, r=0.005) straddles the
  boundary. PE setup instantiates only particles.xyz[0] — keep single sphere.
- `applications/q2p1_el_pipeflow/CMakeLists.txt:65-112`: add per-case mesh-source
  override so straddling (and A3/A4 below) stage the mesh from
  `tier2_cases/momentum_conservation/mesh` (mesh copied only `if(EXISTS .../mesh)`
  today). Add a short `MESH.md` in the case dir noting the shared mesh.
- Run np=28, 5 steps. Accept: `EL_VOLUME_CONSERVATION rel_error < 1e-10` and
  `EL_FEEDBACK_CONSERVATION residual < 1e-10` every step; confirm in the log that ≥2
  ranks hold records for the particle (halo straddling actually exercised).

## WP-A3 — Long-horizon momentum case (opt-in, separate dir)

- New `tier2_cases/momentum_conservation_long/` = copy of momentum_conservation with
  `MaxNumStep=10000`, `MaxSimTime=5.0001` (dt 5e-4 — BOTH caps must be raised),
  `ELMomentumAuditFreq=100`; example.json `timesteps_=10000` for consistency. Mesh
  shared via the A2 CMake override; add to the foreach case list.
- New definition `q2p1_el_pipeflow_tier2_momentum_long.yaml` (`enabled: false`,
  slurm time 02:00:00) + baseline (0.0/0.0): `EL_MOMENTUM_ELEMINT drift_rel` col 16
  tol 1e-5, `EL_FEEDBACK_CONSERVATION residual` col 12 tol 1e-10, occurrence last.
- RUNBOOK.md with acceptance greps: final drift_rel < 1e-5 AND
  `grep -c 'void-fraction clipping'` = 0 (absence can't be a YAML metric).
- Run once (~25+ min, np 28); record final drift_rel in commit message.

## WP-A4 — Saffman case: frozen linear-shear carrier + analytic baseline

Mechanism (q2p1_el_frozen_trace is checkpoint/dump-driven — not reusable cheaply):
an analytic prescribed-carrier mode inside q2p1_el_pipeflow. A linear field is exact
in Q2, so sampled u_f and ∇u at the particle are exact.

- `source/src_el/el_config.f90`: `el_prescribed_field = 'none' | 'linear_shear'`,
  `el_shear_rate` (G); validate + print. Field convention: `u = (G·(z−z_c), 0, 0)`,
  z_c = domain z-midpoint.
- `applications/q2p1_el_pipeflow/app_init.f90`: parse `ELPrescribedField`,
  `ELShearRate` (pattern: the ELMomentumAuditFreq CASE added recently).
- `source/src_quadLS/QuadSc_main.f90` `Transport_q2p1_UxyzP_el`: helper
  `EL_IMPOSE_PRESCRIBED_FIELD()` (workers; pin ILEV=NLMAX + SETLEV per the /TRIAD/
  note; Q2 DOF coords via `SetUp_myQ2Coor`/`myQ2Coor`) writes the field into
  `QuadSc%valU/V/W`; when active, call it before the pre-advance
  `EL_PARTICLE_MESH_PASS` and SKIP `Transport_q2p1_UxyzP_fluid_core` + momentum
  diagnostics. Default `'none'` = bit-identical behavior; routine only used by
  q2p1_el_pipeflow.
- Case files (`tier2_cases/saffman_lift/`): mesh → QBOX9 (CMake override);
  `MaxNumStep=5`, `TimeStep=1e-3`; `ELPrescribedField=linear_shear`, `ELShearRate=1.0`;
  `ELDragModel=stokes`, `ELLiftModel=saffman_mei` (flipped to `saffman` in B2);
  `ELApplyParticleForces=No` (particle stays static ⇒ stationary configuration, every
  step identical), `ELApplyFluidFeedback=No`, `ELPressureForce=No`;
  particles.xyz: one particle at `0.5 0.5 0.75`; example.json 3/3/3, v0=0, g=0.
- Analytic baseline (r=0.005, ρ_f=1, ν=1e-3 ⇒ μ=1e-3, G=1, z_c=0.5 ⇒ slip=(0.25,0,0),
  ω=(0,1,0)):
  - `lift_z = 1.615·d²·√(ρμ|ω|)·0.25 = 1.27678e-6`
  - `drag_x = 3πμd·0.25 = 2.35619e-5` (Stokes cross-check)
- YAML: parse `EL_FORCE_BUDGET`, columns **`{drag_x: 4, lift_z: 14}`** (0-based token
  count of `EL_FORCE_BUDGET step= N drag= 3v pressure= 3v lift= 3v ...`; NOTE the Plan
  agent's 5/15 indices were off by one — 4/14 verified against the WRITE statement;
  existing lift_y:13 was correct but tests the wrong component for this geometry),
  occurrence last, tol 1e-9 each; baseline file gets the two values above.
- Verify: run np=28; all 5 steps print identical budgets (static config); harness
  `featflower-test validate` + compare.

---

> **Review amendments (accepted 2026-07-03):**
> 1. B3 is renamed **"fluid-matrix implicit drag sink/source"** — particle-side
>    semi-implicit integration + issue-D-consistent deposition ALREADY exist
>    (EL_SEMI_IMPLICIT_DRAG_IMPULSE, el_transfer.f90:~207); B3 adds only the missing
>    fluid-operator reaction term. Scope unchanged.
> 2. B1 commit message + rewritten code comment MUST state the convention resolution
>    explicitly: the code comment's "drop to 1−χ" advice belongs to the Zhou
>    set-A/Schiller-Naumann family; we instead adopt the complete Di Felice closure
>    as one package (Dallavalle-type C_D, ε in Re_p, ε^(−χ), interstitial slip) —
>    consistent with the comment's own "Di Felice original" bullet. Richardson–Zaki
>    (V2) remains the empirical arbiter.
> 3. B2 hard rule: source the exact Zeng et al. (2009) constants BEFORE coding.
>    **RESOLVED 2026-07-03: constants extracted from the full paper
>    (zeng_2009.pdf at repo root — Zeng, Najjar, Balachandar & Fischer,
>    Phys. Fluids 21, 033302 (2009)); see the sourced-correlation block in WP-B2.**
>    Structural discovery: Zeng's composite is a STANDALONE lift closure valid for
>    all wall distances within 1<Re<200 (decays to the unbounded-shear value far
>    from the wall, including the lift sign change at Re_cr≈60–200 depending on L)
>    — NOT a multiplicative factor →1. The original "factor →1 far from wall"
>    test spec is withdrawn; WP-B2's model structure is revised accordingly.
> 4. WP-A follow-ups to ride along with WP-B: widen long-momentum tolerance
>    1.05e-5 → 2e-5 (2× measured; document the presumed linear solver-noise
>    accumulation) + RUNBOOK step to archive the workdir simulation log (shared
>    harness workdir gets overwritten by the next case); derive `zc` in
>    EL_IMPOSE_PRESCRIBED_FIELD from mesh z-extents (natural fit in B2's
>    EL_SET_DOMAIN_BOX work); gate the EL_HALO_RECORDS allreduces on
>    el_write_diagnostics.

## WP-B1 — Di Felice drag: implement plan-doc §1.2 exactly

`source/src_el/el_forces.f90` `EL_DRAG_CLOSURE` ('difelice' branch only; stokes/SN
untouched):
- `Re_p = MAX(1e-12, ρ·ε_f·d·|slip|/μ)` (ε now included; keep ε-free Re for other
  models), `C_D = (0.63 + 4.8/√Re_p)²` (new local expression — do NOT change
  `EL_SPHERE_CD`), χ unchanged in form but on the new Re basis, voidage factor
  **ε^(−χ)** replacing ε^(2−χ) (sampled u_f is interstitial —
  el_quadrature.f90:116-141).
- Zero-slip drag_B fallback: slip→0 limit is `2.88π·μ·d·ε^(−χ−1)` (derivation in
  comment; replaces `3πμd·ε^(2−χ)`).
- Replace the provisional comment block (~95-137) with the final convention note; document
  that C_D(Re→0) = 23.04/Re ≈ 4% below Stokes' 24/Re — intrinsic to the correlation.

Tests (`source/src_el/tests/test_el_kernel_forces.f90`):
- STOP 40 (difelice≈stokes within 1%) is unsatisfiable with the new C_D → replace by
  (a) difelice vs its OWN closed form at Re≈1e-8, ε=1, tol 1e-12; (b) ratio
  |F_difelice|/|F_stokes| = (0.63√Re+4.8)²/24 ≈ 0.96025 at Re=1e-6 within 1e-3
  (pins the 4% offset as a tripwire).
- STOP 52 (Re sweep): expected force = new closed form. STOP 53 (ε^(2−χ) ratio,
  line 125): replace by full closed-form evaluation at ε=0.5 (ε changes Re_p→χ→C_D).
  Keep STOP 47 (monotonicity — still holds for ε^(−χ)). Add drag_B fallback
  continuity check near slip=0.
- No Tier-2 baseline uses 'difelice' — no case impact.
- Verify: `ctest -R el-kernel-forces` + `ctest -R 'el-transfer|el-convergence'`.

## WP-B2 — Real lift: Mei (1992) factor + Zeng (2009) wall factor + wall distance

- New `source/src_el/el_geometry.f90` (register in BOTH source lists in
  `cmake/modules/ProjectFiles.cmake`, ~lines 280 and 307): module-cached
  `el_domain_box(6)` + `EL_SET_DOMAIN_BOX(...)` +
  `EL_WALL_DISTANCE(position)` — analytic axis-aligned-box min-face distance; returns
  HUGE when unset (⇒ wall factor 1, one-time warning). Comment names the future
  general-geometry backend: CGAL `getdistanceid`
  (FullC0ntact/inshape3dcore/fortrancppinterface/cppinterface.h:213-221, usage pattern
  fbm_main.f90:2038-2043).
- Wire-up: `applications/q2p1_el_pipeflow/app_init.f90` calls `EL_SET_DOMAIN_BOX`
  right after `get_global_domain_extents` (line ~441).
- `el_forces.f90` `EL_LIFT_CLOSURE`:
  - `EL_MEI_LIFT_FACTOR(re_p, beta)` (pure, literature-cited: Saffman 1965 JFM 22:385;
    Mei 1992 IJMF 18:145): Re_p≤40 → `(1−0.3314√β)e^(−Re_p/10) + 0.3314√β`;
    Re_p>40 → `0.0524√(β·Re_p)`; β = clamp(0.5·Re_shear/Re_p, 0.005, 0.4) documented.
  - **Sourced Zeng correlation** (Zeng, Najjar, Balachandar & Fischer, Phys. Fluids
    21, 033302 (2009) — full paper at repo root `zeng_2009.pdf`; cite eq. numbers in
    code comments). Nondimensionalization: L = wall distance of the particle CENTER
    in diameters (contact L=0.5), δ = L − 1/2, Re = |u_slip|·d/ν (paper: ambient
    velocity at particle center; stationary particle ⇒ slip); lift coefficient
    normalized as `F_L = C_Ls · ½·ρ_f·|u_slip|²·(π d²/4)`, directed wall-normal,
    positive away from the wall:
    - Eq. (19), contact: `C_Ls,w = 3.663 / (Re² + 0.1173)^0.22`
      (Re→0 limit 5.87 = Leighton & Acrivos; large-Re 3.663·Re^−0.44).
    - Eq. (28), composite, valid 1 < Re < 200, any L > 1/2:
      `C_Ls = C_Ls,w · exp(−0.5·δ·(Re/250)^(4/3)) · [exp(α_sL·δ^β_sL) − λ_sL]`
    - Eq. (29): `α_sL = −exp(−0.3 + 0.025·Re)`, `β_sL = 0.8 + 0.01·Re`,
      `λ_sL = (1 − exp(−δ))·(Re/250)^(5/2)`.
    - Reduces to C_Ls,w at δ→0; approximates unbounded-shear finite-Re lift at
      large L (paper Fig. 21 matches Kurose–Komori at L=4), including the lift
      SIGN CHANGE at Re_cr (paper Table VI: L=0.75→198, L=1→126, L=2→75, L=4→59).
    - Sourced test points: Eq.19 Re→0 ⇒ 5.87; paper Table I DNS at L=0.505:
      C_L = 2.653 (Re=2), 1.305 (Re=10), 0.3384 (Re=200) — Eq. 19 reproduces
      these within ~5% (assert accordingly); Eq.28 δ→0 ≡ Eq.19 to 1e-12.
    - (Noted, out of scope: the paper also provides wall-corrected DRAG
      correlations, Eqs. 11/13/14 — candidate future `ELDragWallCorrection`.)
  - Model semantics (el_config validation list gains `'saffman'`):
    `saffman` = classic (today's formula, what Tier-2 pins); `saffman_mei` = classic ×
    Mei factor; wall-aware model = Zeng composite per the STRUCTURE DECISION below.
    Delete the `MIN(1,MAX(0,eps_f))` placeholder. Wall distance flows via the module
    cache — no signature changes in el_transfer.
  - **STRUCTURE DECISION — option A chosen by default (user unavailable; revisit
    before B2 lands):** `saffman_mei_wall` (name kept for config stability) =
    Re-based switch: Zeng composite (slip-based Re, wall-normal direction via
    EL_WALL_DISTANCE + wall-normal from the box face that realizes the min
    distance) for Re_p ≥ 2; Saffman–Mei for Re_p ≤ 0.5; linear-in-log(Re) blend of
    the two force vectors for 0.5 < Re_p < 2 (blend documented in the comment as an
    implementation choice, not literature). Wall distance unset (no
    EL_SET_DOMAIN_BOX call) ⇒ fall back to saffman_mei with a one-time warning.
    Outside Zeng validity (Re_p > 200) clamp Re to 200 and warn-once (document).
- Same commit: flip `tier2_cases/saffman_lift` `ELLiftModel → saffman`; baseline
  invariant (1.27678e-6); rerun case to confirm.
- Tests: `saffman` reproduces current STOP 48 value; Mei factor tabulated points
  (Re→0 ⇒ 1; branch continuity at Re=40, 0.0524√40≈0.3314; Re=100,β=0.1 ⇒ 0.16572).
  Zeng tests (replace STOP 49): Eq.19 Re→0 ⇒ 5.87 within 0.1%; Eq.19 vs paper
  Table I DNS (L=0.505: 2.653@Re=2, 1.305@Re=10, 0.3384@Re=200) within 6%;
  Eq.28 δ→0 ≡ Eq.19 to 1e-12; sign change bracket at L=4 (C_Ls>0 at Re=50,
  C_Ls<0 at Re=70, paper Table VI Re_cr≈59); blend continuity at Re_p=0.5 and 2
  (force vector continuous to 1e-12); unset-box fallback = saffman_mei exactly;
  switchability of all four model names.

## WP-B3 — Fluid-matrix implicit drag sink/source (`ELDragCoupling = explicit | semi_implicit`)

(Particle-side semi-implicit drag + issue-D-consistent explicit deposition already
exist; this item adds only the fluid-operator reaction term — see review amendment 1.)

Default `explicit`; every new path a strict no-op unless
`el_apply_fluid_feedback AND semi_implicit AND field allocated`. Shared-solver safety
is the top risk — single-IF gating, pattern `bNS_Stabilization`
(QuadSc_def.f90:2328-2337).

- **Config** (`el_config.f90`): `el_drag_coupling` string + derived logical
  `el_drag_semi_implicit`; validation: semi_implicit REQUIRES el_apply_fluid_feedback
  (STOP at startup otherwise). Parse in both app_init.f90 files.
- **Fields** (`el_fields.f90`): add `drag_B_ml(:)` (Q2-nodal scalar, distributed/
  unsummed accumulator) + `drag_B_source(:)` (pre-advance snapshot), mirroring the
  force_rhs / fluid_feedback_source pair; allocate/zero/capture/release alongside;
  restart_version 1→2 with the new field.
- **Deposit** (`el_quadrature.f90` `EL_DEPOSIT_PARTICLE`): optional scalar `drag_b`
  arg, spread with the same basis_integral/normalization weights (lines ~335-339).
- **Transfer** (`el_transfer.f90`): widen owned_result/record_result 4→5 comps
  (slot 5 = drag_B; halo reduce/broadcast are ncomp-generic — no halo change). Split,
  derived in a comment: reaction on fluid = −B·u_f + B·U_p − F_lift; the −B·u_f part
  goes implicit, so under semi_implicit the pre-advance pass sets
  `feedback_used = lift − drag_B·U_p` (deposit negates ⇒ fluid RHS gets B·U_p − lift).
  PE-side arming unchanged. `EL_FEEDBACK_CONSERVATION` stays a deposit-fidelity check
  in both modes; print a `mode=` tag; comment that the Newton-pair assertion moves to
  the Tier-2 momentum drift metric (new O(Δt) pairing mismatch by design).
- **Matrix** (`QuadSc_def.f90` `Matdef_general_QuadScalar`):
  - Insertion A (idef=−1, after the DO ILEV loop ends ~line 2325, pointers at NLMAX):
    gated loop `A11/22/33mat(qMat%LdA(I)) += thstep·drag_B_source(I)` — covers all
    Newtonian/non-Newtonian branches at once; thstep=tstep·θ there. Stored
    distributed ⇒ summed once by E013UVWMAT (called 2728/2733 after BC rows). BC
    routines untouched (diagonal kept, defect zeroed).
  - Insertion B (idef=1, old-time pass, thstep=tstep·(1−θ)):
    `defU/V/W −= thstep·drag_B_source·valU/V/W` (vanishes for BE).
  - **Level decision: fine-level-only B** (coarse MG ops are preconditioner-only;
    converged solution defined by fine matrix+defect). Documented consequence:
    possible MG-rate degradation for stiff B — record MG iteration counts in the
    acceptance run; named follow-up = restrict drag_B_source via existing
    E013_Restriction. No new mg_ array needed.
- **Output**: `OutputProfiles.f90` add `ELDragB` scalar PointData export next to the
  existing EL exports (~line 2259; mirror the pvtu/header occurrence).
- **Tests**: test_el_transfer — deposit with known drag_B, assert ΣB_ml = B_p (serial
  + MPI variants); semi_implicit deposited force_rhs = B·U_p − lift. Config
  accept/reject in test_el_kernel_forces.
- **Re-verification protocol (acceptance gate)**:
  1. Explicit regression: rerun momentum_conservation (np 28) — EL_MOMENTUM_ELEMINT /
     EL_FEEDBACK_CONSERVATION sequences must match pre-B3 log (any diff = bug).
  2. Semi-implicit run: same case + `ELDragCoupling=semi_implicit`; record final
     drift_rel, clip count (must be 0), MG iteration counts vs explicit.
  3. Provisional gate drift_rel < 1e-3; then pin the measured value in a new opt-in
     case dir `momentum_conservation_semi/` + YAML (tolerance = 2× measured), and a
     short `docs/md_docs/el_semi_implicit_drag_note.md` (audit semantics + θ placement
     + measured numbers).
  4. Default stays `explicit` until step 3 numbers are reviewed with the user.
  5. One non-EL smoke (existing FBM/DNS featflower test or q2p1_fc_ext run) to prove
     no shared-path regression.

---

## Sequencing & risks

Order: A1 → A2 → A3 → A4 → B1 → B2 → B3 (A1-A3 config-only, unblock verification now;
A4 independent new code; B1/B2 pure el_forces layer — Tier-2 terminal/straddling/
momentum pin `stokes` drag so B1 can't invalidate them; saffman pins `saffman` classic
so B2 leaves its baseline invariant; B3 last, gated by the momentum rerun).

- B3 touches shared solver code → gated single-IF, explicit-mode byte-identity check,
  revert = one commit (fields/keys inert).
- Saffman 1e-9 tolerance assumes exact Q2 linear field + interior particle; if cubature
  noise exceeds it, widen to 1e-8 and record measured.
- Straddling ownsPoint ambiguity: if 0.32733 still lands wrong, shift to 0.325 (support
  still straddles).
- Zeng constants must come from the paper at implementation time — flagged, not
  invented.
- /TRIAD//ILEV: no new NDFGL sites; new field-imposition helper pins level per
  docs/md_docs/el_triad_ndfgl_common_block_note.md; deposits stay inside the guarded
  pass (EL_ASSERT_KDFG_IN_BOUNDS).

## Verification summary

- Unit: `ctest -R 'el-kernel-forces|el-transfer|el-convergence'` in the gcc14 build.
- Tier-2 runs per case RUNBOOKs (np 4 / 28), metrics via
  `featflower-test validate|run <definition.yaml>`; acceptance values above.
- B3: three-step momentum re-verification protocol; explicit mode must be
  bit-identical, semi_implicit drift measured and pinned.
