# Review brief: ELVACompatibility = alpha (fix A) port + restructuring of the EL continuity-source module (FeatFloWer, Fortran 90 + MPI)

Repository (git worktree, branch feature/el-volume-averaged):
  /data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va
Commit under review: cef8eeaa (parent 57fa2e1f). See it with
  git -C /data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va show cef8eeaa
Files:
  source/src_el/el_volume_averaged.f90     (restructured module: named reduction slots, EL_VA_SINK_MODE / EL_VA_SINK, EL_VA_LOG; new mode 'alpha')
  source/src_el/el_config.f90              (key SimPar@ELVACompatibility = uniform|alpha, validation, echo)
  applications/q2p1_el_pipeflow/app_init.f90 (parser case)
  docs/md_docs/el_volume_averaged_stage.md, validation_cases/s04_voidfraction_cfdem_twin/q2p1_param_VB_bdf2_alpha.dat
For context, the previous two commits on the same module (14d478a5, cafa91bc: BDF2 source + projection
correctors, reviewed two days ago by you and a second model) and the corrector include
source/src_quadLS/QuadSc_el_va.f90, which calls EL_ADD_CONTINUITY_SOURCE with icorr = 1..n and -1.

## What the code does
FeatFloWer is a parallel Q2/P1 finite-element Navier-Stokes solver. Rank 0 is a CFD master without mesh
or EL fields that joins every collective pressure solve; the workers reduce over MPI_COMM_SUBS. The EL
volume-averaged continuity stage (SimPar@ELVolumeAveraged = continuity) imposes div u = g with
g = -(d eps/dt + u.grad eps)/eps through the pressure right-hand side (EL_ADD_CONTINUITY_SOURCE). In a
closed/periodic box the pressure matrix has a constant nullspace, so the volume integral int_g of the raw
source must be removed. Until now this was done uniformly (g_e - gbar, gbar = int_g/V). The new mode
'alpha' removes it through the sink s_e = int_g alpha_p(e)/W, W = sum alpha_p|e| (sink where the
particles are; uniform fallback when W < 1e-30). It was validated on a test branch (8aeb26c8): it closes
the level-3 return-flow deficit of the s04 twin (G5 0.9926 -> 1.0021) and lifts level 2 from 1.0039 to
1.0055, matching the peer's predictions to the printed digits.

The commit ALSO restructures the module. Before: one routine with local(1..13) magic-index accumulators,
every audit line written twice (file + stdout) with duplicated format strings, the removal inline. After:
named parameters for three reduction blocks (C_* core block = the campaign 8-entry Allreduce, must stay
bitwise; T_* temporal block reduced on audit steps only; A_* alpha block reduced in mode alpha only),
EL_VA_SINK_MODE (deck -> mode) and EL_VA_SINK (mode -> per-element sink; the expressions are kept
verbatim from the validated runs: uniform 'gbar', alpha 'int_g*alpha_e/w_alpha'), EL_VA_LOG (one line to
both units from rank showid), and four labelled sections in EL_ADD_CONTINUITY_SOURCE.

Bitwise requirements (gates running now, you cannot see them; reason from the code):
- default keys: identical floating-point result and identical common log streams to the pre-commit code
  (tier-2 twin and VB_L2 default gates against certified logs);
- BDF2 keys on, uniform: identical to the T8b log of commit cafa91bc;
- alpha keys on: identical to the T9a/T9b logs of the test branch 8aeb26c8 (whose implementation used
  local(13) accumulators and 'gc = gk(iel) - global(1)*ak(iel)/w_alpha').

## What the owner wants (their words)
"Move fix A to the VA branch, same review and gates. [...] The implementation of new features should be
clean, maintainable, extensible code. We want a sustainable codebase."
and, from the previous round: modular, no unwanted side effects at the defaults, truly switchable from
q2p1_param.dat.

## Please check
1. Sustainability: is the restructured module clean, readable and extensible? Would adding a third sink
   mode or a third time scheme be a local change? Are names, comments and the section layout right? Is
   there anything that a maintainer would trip over (hidden coupling between EL_VA_SINK_MODE() called
   twice - once at the start and once in the audit block - and the run-time fallback to uniform that
   rebinds the local 'mode'; the audit block's test on the deck mode vs the effective mode; the
   'alpha' local reused for two purposes; LOG_LEN; MERGE with padded literals)? Anything you would do
   differently for a sustainable codebase: say so concretely, with the code you would write.
2. Bitwise identity at the defaults: same operations in the same order for the raw source, the core
   accumulators, gbar and the subtraction g_e - gbar (now via EL_VA_SINK returning gbar); the same single
   8-entry Allreduce; no new collective at the defaults; log bytes identical when the line is first
   written to a CHARACTER buffer and then printed with '(A)' and TRIM (watch for trailing blanks or
   truncation: LOG_LEN = 320 vs the longest line).
3. Bitwise identity in the alpha mode vs the test-branch expression (int_g*alpha_e/w_alpha evaluated
   left to right; the accumulators alpha*vol and zc*alpha*vol; the corrected-integral sum).
4. MPI: every worker reaches the same Allreduces (core always; tsrc on audit steps; alpha block when the
   deck says alpha; the corrected-integral reduce on audit steps when the deck says alpha - note the
   run-time fallback to uniform must not change which collectives are reached, because W is a global
   quantity, identical on all workers after the reduce; verify this reasoning).
5. Deck switch: validation and echo in el_config.f90, parser in app_init.f90, rejection outside the
   continuity stage, lower-casing, defaults; the corrector passes also get the alpha removal (they call
   the same routine).
6. Leftovers: unused variables, imports, PUBLIC symbols; compiler-warning noise; naming consistency with
   the rest of src_el; the header comment's accuracy.

## Output
Findings ordered by severity (blocker / should fix / nit), each with file:line, what is wrong, why it
matters, a concrete fix (code where useful). Then a short verdict on: sustainable (clean, maintainable,
extensible)? side effects at defaults? truly deck-switchable? Say what you verified by reading and what
you could not. Do not modify any files. Read-only review.
