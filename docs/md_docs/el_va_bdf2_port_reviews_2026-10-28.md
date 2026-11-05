# Reviews of the BDF2 source / projection-corrector port (commit 14d478a5, 2026-10-28)

Both reviews were read-only and ran on 14d478a5 before the fixes in cafa91bc. Brief: `el_va_bdf2_port_review_brief_2026-10-28.md` (same text for both). Verbatim final messages; the Fable text is the agent's report as returned.

## Codex (gpt-6.1-sol, high)

I would fix the direct-solver lifetime issue before merging. The change is reasonably modular and switchable from the startup deck, but strict default bitwise equivalence is not provable from this implementation.

1. **Blocker — correctors can reuse a freed pressure-solver factorization.**  
   [QuadSc_el_va.f90:73](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_quadLS/QuadSc_el_va.f90:73)

   With `MatrixRenewal` containing `C2` or `C3`, the predictor’s `MG_Solver` frees the master’s UMFPACK/PARDISO factorization on exit ([QuadSc_mg.f90:135](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_quadLS/QuadSc_mg.f90:135)). The corrector immediately calls the solver again without the matrix setup that creates that factorization.

   This affects pressure coarse-solver types 2/4 and, when compiled, 9. PARDISO explicitly aborts when called without factorization; UMFPACK receives an invalidated numeric handle. The supplied sample uses `C1`, so its gates do not cover this combination.

   **Fix:** retain the factorization through the predictor and all correctors, then perform the existing cleanup once. A smaller interim fix is to reject the unsupported combination explicitly.

2. **Should fix — the default source reduction changes the floating-point execution contract.**  
   [el_volume_averaged.f90:161](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_el/el_volume_averaged.f90:161)

   The worker Allreduce changes unconditionally from 8 doubles to 10. MPI may select a different reduction algorithm or segmentation based on message length, changing summation order. `global(1)` and `global(2)` determine `gbar`, so this can affect the pressure RHS and subsequent solution, as well as common log streams.

   At BE, the element-level source arithmetic otherwise retains the old operand order: `d_cont` simply receives `deps_f_dt`. The unchanged arithmetic does not establish bitwise equivalence across the changed collective.

   **Fix:** preserve the original eight-entry reduction for default production values. BE’s temporal L1 is already `global(3)`; its signed temporal diagnostic could be derived from existing source/advection totals, with that derivation documented. Alternatively, explicitly relax the strict collective-equivalence requirement for diagnostics.

   Other unconditional default changes are real, although I found no direct numerical dependency on them:

   - An extra `REAL*8(nel)` allocation and history copy in `EL_FIELDS`.
   - Version-3 restart output.
   - Two configuration echo lines and the new `EL_VA_TSRC` stream.
   - Additional diagnostic arithmetic and associated timing overhead.

3. **Should fix — BDF2 silently assumes equal timesteps.**  
   [el_volume_averaged.f90:81](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_el/el_volume_averaged.f90:81)

   `1.5*d_current - 0.5*d_previous` is endpoint BDF2 only for equal adjacent intervals. No previous interval length is stored or checked.

   Timestep changes are reachable: `postprocessing_app` calls `SetSimulationTimeStep`, and `Reload_SimPar` can change `TimeStep` while synchronizing PE’s timestep. Existing valid rate flags survive such changes, so the next source is labelled BDF2 despite using incorrect coefficients.

   **Fix:** store the previous interval length and use variable-step coefficients, or invalidate BDF2 history whenever the timestep changes and use BE for that transition. Document the chosen behavior.

4. **Should fix — rate-validity invariants are incomplete around restart and non-difference finalization.**  
   [el_transfer.f90:475](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_el/el_transfer.f90:475), [el_fields.f90:202](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_el/el_fields.f90:202)

   Only the genuine-difference branch assigns `rate_valid`. The initial-history and `dt<=0` branches overwrite `deps_f_dt` with zero without clearing an existing true flag. Also, version 3 restores the rate flags but does not restore `history_valid`; a fresh read can therefore produce valid rate flags alongside invalid epsilon history. The BDF2 guard checks only `prev_rate_valid`.

   This does **not** invalidate the normal positive-timestep pre/post sequence I traced: the post pass replaces the temporary pre-pass rate before source evaluation. It does leave the exported field/restart API capable of promoting a synthetic zero rate into valid previous history.

   **Fix:** clear `rate_valid` whenever no genuine difference is formed, restore or derive `history_valid` consistently on restart, and require current as well as previous rate validity for BDF2. Extend the restart test to continue through field updates and source evaluation; the present test only checks stored arrays and flags.

5. **Nit — malformed deck values lack reliable, key-specific rejection.**  
   [app_init.f90:803](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/applications/q2p1_el_pipeflow/app_init.f90:803)

   Unknown scheme strings and out-of-range integers receive clear validation messages. Malformed integers instead encounter an unchecked list-directed `READ`. List-directed null input or `/` can also leave the previous/default value unchanged, passing validation.

   **Fix:** parse into temporaries with `IOSTAT`, reject missing/null values, and report the key and supplied value.

6. **Nit — corrector observability and exit-state comments overstate what is preserved.**  
   [QuadSc_el_va.f90:31](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_quadLS/QuadSc_el_va.f90:31), [QuadSc_el_va.f90:87](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_quadLS/QuadSc_el_va.f90:87)

   - On non-audit steps, MG can leave `defP` holding a solver residual. On audit steps, the final assembly overwrites it with the final constraint assembly. `rhsP` retains the last solve’s RHS, including solver boundary treatment. They do not generally both hold “the last pass’s assembly.”
   - The final audit assembly and its three reductions are outside the added timing brackets.
   - `Protocol_LinScalar` does not print `cTitle`; its title-output code is commented out. `"Pressure-Poisson corrector"` therefore does not distinguish those protocol blocks.

   **Fix:** document the actual work-array state, account for the final diagnostic cost, and emit an explicit corrector/pass label when enabled.

By code inspection, the following parts are sound:

- **History ordering:** the pre pass copies the completed previous-interval rate before recomputing fields; the post pass preserves that copy and forms the current interval’s rate. Fresh step 1 uses BE; step 2 can use BDF2. `EL_FINALIZE` followed by initialization resets both rate flags. `EL_INITIALIZE` alone remains an “ensure allocation” operation when dimensions already match.
- **Pressure state:** `p_n` preserves the predictor’s `P^n`; each pass bases pressure accumulation on the preceding corrected pressure. Exit restores `valP_old=P^n` and rebuilds `P_new` from final pressure before force evaluation.
- **`GAMMA>0`:** `Pressure_Correction` consumes the current pass’s `rhsP` before the final diagnostic assembly. I found no stale predictor-RHS use introduced here.
- **Downstream readers:** I found no reader of pressure `rhsP` or `defP` after the fluid core in this application’s normal continuation/postprocessing path. The next pressure assembly overwrites them.
- **Statistics:** each extra solve updates `tMGP` and `iLinP`; the wrapper adds assembly and correction times. Solver output parameters describe the last solve.
- **Fortran details:** optional-argument use is safe here; `MERGE('BDF2','BE  ',...)` has matching character lengths; `p_n` has ordinary local allocatable lifetime. I found no obvious new unused import/local problem.

MPI participation is correct for the normal, nonempty, consistently allocated worker partitions: every rank enters exactly `n` additional pressure solves, while the new direct field/array accesses exclude the master. Every worker enters the source reduction each pass; audit steps add two residual reductions per corrector and three reductions for the final assembly.

There remains a defensive MPI hazard: rank-local early returns before collectives. The source’s empty/undersized-field returns predate this commit; the new residual routine adds an undersized-array return before its Allreduces ([el_volume_averaged.f90:217](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_el/el_volume_averaged.f90:217)). Empty workers should contribute zero, and inconsistent array sizes should cause a coordinated failure rather than a local return.

The version-3 reader accepts versions 1–2 and disables their rate history. Older readers reject version 3. That compatibility break can be acceptable, especially since these routines currently have only test callers, but the documentation should explicitly state that new files cannot be read by older binaries.

**Owner verdict:**

- **Modular?** Yes, broadly: one include, one guarded call site, and a compact configuration guard. Solver-factorization lifetime needs integration beyond that wrapper.
- **Side effects at defaults?** No additional pressure solves, and common source/moment text formatting is preserved if values match. Strict numerical/collective identity is unproven because of the enlarged Allreduce; memory, restart format, diagnostics, and timing also change.
- **Truly deck-switchable?** Yes at startup in an existing supported `q2p1_el_pipeflow` build, without a new feature flag or rebuild. Validation occurs after parsing, so distinct-key ordering is immaterial. Nondefault controls are rejected outside `continuity`; explicit `BE`/`0` are accepted. Existing prerequisites still apply: the application requires PE, serial PE is rejected, and `ELPrescribedField != none` bypasses the fluid solve. Runtime `Reload_SimPar` uses the generic parser and does not reload these new keys.

All conclusions above come from reading commit `14d478a5`, its parent, and the surrounding call paths. I did not build, run MPI tests, or verify the pending bitwise gates, numerical reproduction, or current compiler warnings. No files were modified.
## Claude Fable 5.1 (xhigh)

**Scope note.** Reviewed the commit as asked. The worktree at HEAD = 14d478a5 carried uncommitted edits to 5 files (el_fields.f90, el_transfer.f90, el_volume_averaged.f90, test_el_kernel_forces.f90, QuadSc_el_va.f90) that already address findings F1-F4 below; marked where so. Everything else is against the committed code. No build, run, or state-changing git command was used.

### Findings (by severity)

No blockers on the gated decks. One deck-reachable crash path exists (F1).

**F1 (high) - corrector loop can solve with a freed coarse factorization; purely deck-triggered.**
`source/src_quadLS/QuadSc_el_va.f90:79-104` runs a second `Solve_General_LinScalar` per step. `source/src_quadLS/QuadSc_mg.f90:135-144`: at the end of every pressure `MG_Solver` call, on the master, `IF (cVariable=="Pressure" .AND. myid==0 .AND. myMatrixRenewal%C>=2) CALL myUmfPack_Free()` (or `myPardiso_Free()`). `source/UmfpackSolver.f90:9-32`: `myUmfPack_Solve` reuses the stored `numeric` handle, no lazy refactorization; refactorization only happens in `QuadSc_solver_coarse.f90:110/113/194` at operator regeneration. Pressure coarse types that use it: 2 and 4 (UMFPACK) and 9 (PARDISO). So `SimPar@MatrixRenewal = ...C2|C3` + `Pres@MGCrsSolverType = 2|4|9` + `ELVASourceCorrectors > 0` calls `umf4sol` on a freed handle in the corrector. Every EL deck in the repo uses `C1` + type 4, so no gated run hits it, but it is exactly the "unwanted side effect reachable from the param file" the owner asked about. Fix: reject the combination with a message (uniform on all ranks), or skip the free while correctors are on. *Worktree already adds the ERROR STOP guard inside `EL_VA_Projection_Correctors` (fails on the first fluid step; acceptable, or move to app_init after `EL_VALIDATE_CONFIG` where `LinSc%prm` is visible so it fails at startup).*

**F2 (medium) - default-path identity of the source Allreduce is empirical, not provable.**
`MPI_Allreduce(local, global, 10, ...)` replaces the previous count 8. Exact identity of elements 1-8 rests on the MPI library, not on the code. *Worktree splits it into the original count-8 call plus an audit-gated 2-element reduce; endorse.*

**F3 (medium) - BDF2 coefficients assume equal steps; `rate_valid` is never cleared.**
`EL_VA_CONTINUITY_RATE` hard-codes 1.5/-0.5; `el_transfer.f90` only ever sets `rate_valid = .TRUE.`. Not a live bug: the EL loop uses a constant step (`q2p1_el_pipeflow.f90:30-31`, `dt = tstep`; `IADTIM` is parsed but unused on the EL path), and `dt<=0` cannot occur. But the module routine is generic and silently wrong under variable dt. *Worktree stores `rate_dt`/`rate_dt_prev`, uses a = dt/(dt+dt_prev) (exactly 0.5 for equal steps), and clears the flags on non-genuine branches; endorse.*

Nits:
- F4 `EL_VA_LOG_CONSTRAINT_RESIDUAL`: rank-local `IF (SIZE(defp).LT.4*nel) RETURN` ahead of two Allreduces. Unreachable (`LinSc%ndof = 4*nel`) but a deadlock trap. Same pre-existing pattern in `EL_ADD_CONTINUITY_SOURCE`. *Worktree fixes the new one.*
- F5 `tools/el_stream_compare.py` TAGS lacks `EL_VA_TSRC` and `EL_VA_RESIDUAL`, so the keys-on reproduction gate never compared the new streams. Add both.
- F6 `MERGE('BDF2','BE  ',...)` without TRIM gives `scheme= BE  ` with trailing blanks. Cosmetic.
- F7 `el_config.f90:343` mis-indented echo line.
- F8 `Protocol_LinScalar` ignores `cTitle` (title WRITEs commented out): the corrector's pressure block is indistinguishable from the predictor's in the protocol. Pre-existing; optional.
- F9 Unit test covers only the restart round-trip. `EL_BEGIN_FIELD_UPDATE`'s prev-capture semantics are MPI-free and cheap to test.
- F10 Docs: production restarts never call `EL_READ_RESTART` (unit test is the only caller), so a `StartingProc` restart always begins with one BE step; worth one sentence.

### Verified fine
- Step ordering: pass with `-itns` -> `advance_history=.TRUE.`, `EL_BEGIN_FIELD_UPDATE` precedes `EL_FINALIZE_VOID_FRACTION` in each pass; `EL_ADVANCE_PARTICLES`; pass with `+itns`; fluid core. So `deps_f_dt_prev` = the post-advance rate of step n at the moment the step n+1 source is evaluated; the pre-advance FINALIZE's ~0 value is overwritten before any consumer. Step 1 falls back to BE, step 2 onward BDF2. Correct.
- FP identity at defaults: `d_cont = deps_f_dt(iel)` exactly; `g` and `local(6)` same expressions; `corr_tag=''` so the trailing `A` descriptor writes zero bytes; default logs gain only the two config echo lines and `EL_VA_TSRC` on audit steps. No collectives added at defaults (the guard is a pure function of two deck variables, identical on all ranks).
- Corrector side effects: `valP_old` restored to P^n; `P_new` rebuilt over `SIZE(P_new)=ndof=4*nel=lMat%nu`, consistent with `Pressure_Correction`; `defP`/`rhsP` are work arrays (`Resdfk_General_LinScalar` reads them but only for `Tracer`; MG uses them as its own buffers; no reader in `postprocessing_app`/`handle_statistics`). `myStat`, `LinSc%prm%MGprmOut`, `iLinP` count the extra solves (expected). The audit-gated final assembly cannot change results.
- GAMMA>0 branch uses `rhsP` = the corrector's residual RHS, so the rotational term telescopes across passes. Not reachable from the EL deck (`myGDATNEW` has no GAMMA key).
- MPI: master skips every worker block exactly as the predictor does and joins `Solve_General_LinScalar` and `Protocol_LinScalar`. Worker Allreduces per pass: 1 in the source, +2 on audit steps, gated on values identical on all workers. `EL_VOLUME_AVERAGED_ACTIVE` is uniform on workers, `.FALSE.` on the master, so the deck-only guard is the right one.
- Deck-only switch: both keys parsed in `myGDATNEW` on every rank; `EL_VALIDATE_CONFIG` runs after the file is closed, so key order is irrelevant; values lower-cased; invalid value -> `STOP 1` with a message on all ranks; non-default values without `continuity` -> STOP. No compile flag beyond the pre-existing parser gate. A non-numeric integer gives a Fortran runtime I/O error rather than a message, same as every other integer key.
- Restart v3: old binary on a v3 file -> clean `ERROR STOP`; new reader on v1/v2 -> flags cleared.
- Fortran: keyword `icorr=0` with an explicit module interface OK; `int_g_out` zeroed before the early returns; MERGE lengths equal; `corr_tag` LEN 16 >= 9 needed; host association in the include identical to the fluid core's use; local `USE` beside the module-level `USE` is legal; CMake needs no change (include fragments are not listed in `cmake/modules/ProjectFiles.cmake`).
- Port fidelity vs 9fb6157c: same operation order in the loop; differences (local `p_n` instead of the module SAVE array, audit gating of the residual log and final assembly, dropped `ELVASourceEpsLevel`, predictor passing `icorr=0` unconditionally) are all equivalent.

### Verdict
- Modular? Yes. One include file, one guarded call site, one guard function; the fluid-core diff is three lines plus a USE; `src_el` changes are self-contained.
- Side effects at defaults? None in arithmetic or control flow found by reading. The one point not provable from code is the Allreduce count change (F2; fixed in the worktree). Logs gain only new streams.
- Truly deck-switchable? Yes: two keys, validated on every rank with clear messages, order-independent, no build-flag dependence, with one exception: the `C2/C3` + direct-coarse-solver combination (F1) must be rejected before this is safe for arbitrary user decks.

Not verified: the bitwise gate jobs (logs not inspected); the MPI library's elementwise behaviour for count 8 vs 10; the default value of `GAMMA`; the s04 RUNBOOK numbers; compilation of the include file.
