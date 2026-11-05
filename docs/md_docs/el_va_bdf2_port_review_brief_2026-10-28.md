# Review brief: deck-switchable BDF2 continuity source + projection correctors (FeatFloWer, Fortran 90 + MPI)

Repository (git worktree, branch feature/el-volume-averaged):
  /data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va
The change under review is the single commit 14d478a5 (parent 4fcf14ca). See it with
  git -C /data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va show 14d478a5
Files touched (absolute paths under that worktree):
  source/src_el/el_config.f90              (two new deck keys, validation, echo, EL_VA_SOURCE_CORRECTORS_ENABLED)
  source/src_el/el_fields.f90              (deps_f_dt_prev, rate_valid, prev_rate_valid; restart version 3)
  source/src_el/el_transfer.f90            (rate_valid = .TRUE. where the genuine backward difference is formed)
  source/src_el/el_volume_averaged.f90     (BDF2 rate, pass labels, EL_VA_TSRC / EL_VA_RESIDUAL logging)
  source/src_quadLS/QuadSc_el_va.f90       (NEW include of module Transport_Q2P1: EL_VA_Projection_Correctors)
  source/src_quadLS/QuadSc_main.f90        (one guarded call in Transport_q2p1_UxyzP_fluid_core, icorr=0 on the predictor, include line)
  applications/q2p1_el_pipeflow/app_init.f90   (parsing of SimPar@ELVASourceTimeScheme / SimPar@ELVASourceCorrectors)
  source/src_el/tests/test_el_kernel_forces.f90 (restart round-trip test extended)
  docs/md_docs/el_volume_averaged_stage.md, validation_cases/s04_voidfraction_cfdem_twin/q2p1_param_VB_bdf2.dat (docs, sample deck)

## What the code does
FeatFloWer is a parallel Q2/P1 finite-element Navier-Stokes solver (projection scheme: momentum
predictor, pressure Poisson solve, Velocity_Correction / Pressure_Correction, see
source/src_quadLS/QuadSc_corrections.f90). Rank 0 is a CFD master that owns no mesh data and no
Euler-Lagrange (EL) fields but takes part in every collective solve (Solve_General_LinScalar).
The EL volume-averaged continuity stage (deck key SimPar@ELVolumeAveraged = continuity) imposes
div u = g, g = -(d eps/dt + u.grad eps)/eps, by adding G/dt to the pressure right-hand side in
EL_ADD_CONTINUITY_SOURCE (called in the fluid core after Matdef_General_LinScalar(...,1)).

The commit ports a scheme that was validated on a test branch:
1. SimPar@ELVASourceTimeScheme = BE | BDF2 (default BE): with BDF2 the rate in g is
   1.5*(d eps/dt)^{n+1} - 0.5*(d eps/dt)^n. The previous-interval rate is stored in
   el_field_data%deps_f_dt_prev, copied from deps_f_dt in EL_BEGIN_FIELD_UPDATE when the
   history advances (that routine is called once per time step with advance_history=.TRUE.
   in the pre-advance pass, and again with advance_history=.FALSE. in the post-advance pass;
   the genuine rate deps_f_dt = (eps^{n+1}-eps^n)/dt is formed in EL_FINALIZE_VOID_FRACTION,
   el_transfer.f90, in the ELSE IF (dt>0) branch). rate_valid / prev_rate_valid guard the
   first step and restarts without history (fall back to BE).
2. SimPar@ELVASourceCorrectors = n (default 0): after the projection, EL_VA_Projection_Correctors
   runs n passes: save valP as the base, zero valP, reassemble (1/k)B^T u with the corrected
   velocity, add the recomputed source, copy to rhsP, collective solve, Velocity_Correction,
   Pressure_Correction. P^n is kept in a local copy so that afterwards LinSc%valP_old = P^n
   again and LinSc%P_new = 1.5 P^{n+1} - 0.5 P^n is rebuilt from the final P^{n+1}
   (P_new is read by the force/torque evaluation, QuadSc_force_torque_calc.f90).
   The guard in the fluid core is EL_VA_SOURCE_CORRECTORS_ENABLED() (deck keys only), NOT
   EL_VOLUME_AVERAGED_ACTIVE() (which also tests field allocation and is therefore .FALSE.
   on the master): a guard on allocation deadlocked the collective solve on the test branch.
3. Logging: EL_VA_TSRC line always on audit steps; the ' corr= k' tag on EL_VA_SOURCE /
   EL_VA_MOMENT / EL_VA_TSRC and the EL_VA_RESIDUAL lines only when correctors > 0.
   The log compare tool used for the bitwise gates ignores extra trailing tokens and
   streams absent from the reference, so the default-key logs must be byte-identical
   for the common streams (gates running now: tier-2 twin and VB_L2 at the defaults).

The deck is applications/..../_data/q2p1_param.dat, parsed line by line in app_init.f90
(SimPar@Key = value); EL keys are validated in el_config.f90 (EL_VALIDATE_CONFIG or similar)
and echoed to the protocol file.

For reference only, the original test-branch implementation (not for merge) is in the worktree
/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va-fixb, commits 4c0dd4f3 and 9fb6157c
(it also carried a third key, ELVASourceEpsLevel, deliberately NOT ported).

## What the owner wants checked (their words)
"make sure its implementation is modular and can be easily switched for example in the
q2p1_param.dat file ... double check it for modularity, that it does not have unwanted side
effects. That it is truly easily switchable by the user over the param file."

Specifically:
- Is behaviour at the default keys provably identical to the pre-commit code (same floating
  point operations in the same order, same log lines for the common streams, no extra
  collective calls, no extra allocations that change anything)? Point at any path where
  the defaults could change results or logs.
- Is the switch truly deck-only: does any part of the new behaviour depend on compile-time
  options, on build flags, on a key other than the two documented ones, or on the ordering
  of keys in the deck? Are invalid values rejected with a clear message? Are the keys
  rejected or silently ignored when ELVolumeAveraged is not 'continuity'?
- Side effects of the corrector loop on solver state: valP_old, P_new, rhsP, defP, the
  timing statistics myStat, the Protocol_LinScalar output, anything else in the
  Transport_Q2P1 module scope that Velocity_Correction / Pressure_Correction /
  Solve_General_LinScalar touch (e.g. the GAMMA > 0 AddDiffPrec branch in
  Pressure_Correction uses rhsP). Does any routine after the fluid core read defP or rhsP?
- MPI safety: every rank must reach the collective solve the same number of times; the
  master (myid == 0) must not touch EL fields or LinSc arrays it does not own; the
  Allreduces in EL_ADD_CONTINUITY_SOURCE and EL_VA_LOG_CONSTRAINT_RESIDUAL are over
  MPI_COMM_SUBS (workers only) and must be reached by every worker the same number of
  times, including the audit-step-gated ones.
- BDF2 history correctness: is deps_f_dt_prev the rate of the PREVIOUS interval at the moment
  the source is evaluated (post-advance pass of step n+1)? Are rate_valid /
  prev_rate_valid set and reset in the right places, including after a restart read and
  after EL_FINALIZE/EL_INITIALIZE? Is the step-1 fallback to BE correct?
- Restart format: version 3 writes extra data; older readers break on new files (is that
  acceptable / documented?), new reader on old files is handled.
- Modularity: is the corrector loop isolated enough (one include file, one call site, one
  guard function), is the fluid-core change minimal, are there leftovers (unused USE
  items, unused variables, public symbols that need not be public)?
- Anything else: Fortran pitfalls (OPTIONAL argument handling, MERGE with character
  constants of different length, SAVE semantics, allocatable locals, integer/real kinds),
  compiler-warning noise, naming consistency with the rest of src_el.

## Output
Findings ordered by severity (blocker / should fix / nit). For each: file:line, what is wrong,
why it matters, a concrete fix. Then a short verdict on the three owner questions
(modular? side effects at defaults? truly deck-switchable?). Say explicitly which claims you
verified by reading code and which you could not verify.
Do not modify any files in either worktree. Read-only review.
