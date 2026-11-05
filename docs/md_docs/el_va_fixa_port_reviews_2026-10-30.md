# Reviews of the ELVACompatibility (fix A) port and module restructuring (commit cef8eeaa, 2026-10-30)

Both reviews were read-only and ran on cef8eeaa before the fixes in f93eb55e. Brief: `el_va_fixa_port_review_brief_2026-10-30.md` (same text for both). The owner's emphasis for this round: clean, maintainable, extensible code.

## Codex (gpt-6.1-sol, high)

The restructuring is substantially clearer, and I found no source-level change to the intended numerical expressions. There is one blocker against the stated **identical alpha log-stream requirement**.

1. **Blocker — alpha audit records changed order.**  
   [el_volume_averaged.f90:307](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_el/el_volume_averaged.f90:307)

   In `8aeb26c8`, each audit group is:

   ```text
   COMPAT → SOURCE → MOMENT → TSRC
   ```

   In `cef8eeaa`, it is:

   ```text
   SOURCE → MOMENT → TSRC → COMPAT
   ```

   Thus the chronological alpha streams cannot be byte-identical, even if every number matches. Comparisons performed separately by record type could conceal this.

   **Fix:** move the complete alpha audit block—including its worker Allreduce—before the `EL_VA_SOURCE` write. Keep the reduction outside any `showid` guard. This restores reference ordering without changing arithmetic or collective order relative to other collectives.

2. **Should fix — fallback reports a misleading moment shift.**  
   [el_volume_averaged.f90:270](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_el/el_volume_averaged.f90:270), [el_volume_averaged.f90:314](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_el/el_volume_averaged.f90:314)

   When `W < 1e-30`, the effective sink becomes uniform, but `z_alpha` remains zero. The logged shift therefore becomes `-z_centre*int_g`, although alpha weighting is no longer applied. This matters whenever the fallback coincides with a nonzero raw integral.

   This behavior is **inherited from the test branch**, rather than introduced by the restructuring.

   **Fix:** retain zero as the undefined `z_alpha` convention, but report zero shift for the uniform fallback:

   ```fortran
   mshift = 0.0d0
   IF (mode.EQ.SINK_ALPHA) &
     mshift = (z_alpha-EL_DOMAIN_Z_CENTER())*core_g(C_INT_G)
   ```

   Use `mshift` in the existing format. Active-alpha output remains unchanged; fallback output intentionally changes.

3. **Nit — make requested versus effective mode explicit.**  
   [el_volume_averaged.f90:201](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_el/el_volume_averaged.f90:201), [el_volume_averaged.f90:307](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_el/el_volume_averaged.f90:307)

   The current code is correct, but its safety depends on understanding that `mode` changes meaning after fallback, while the second `EL_VA_SINK_MODE()` call recovers the original request. That is an avoidable maintenance trap.

   **Fix:** capture both meanings once:

   ```fortran
   INTEGER :: requested_mode, effective_mode

   requested_mode = EL_VA_SINK_MODE()
   effective_mode = requested_mode

   ! Allocation, weight accumulation and reductions:
   IF (requested_mode.EQ.SINK_ALPHA) ...

   ! After the weight reduction:
   IF (requested_mode.EQ.SINK_ALPHA) THEN
     IF (w_alpha.LT.1.0d-30) effective_mode = SINK_UNIFORM
   END IF

   ! Element sink:
   g = gk(iel) - EL_VA_SINK(effective_mode, gbar, &
                          core_g(C_INT_G), alpha, w_alpha)

   ! Corrected-integral reduction and compatibility audit:
   IF (requested_mode.EQ.SINK_ALPHA) ...
   ```

   Use `effective_mode` for the logged mode as well. This preserves the present collective behavior and makes the distinction visible.

4. **Nit — extension documentation overpromises locality.**  
   [el_volume_averaged_stage.md:374](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/docs/md_docs/el_volume_averaged_stage.md:374), [el_volume_averaged.f90:48](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_el/el_volume_averaged.f90:48)

   A third sink requires changes beyond `EL_VA_SINK` and configuration: mode decoding, weight storage/accumulation/reduction, fallback handling, and labels. The current `MERGE` would label any non-alpha third mode as `uniform`.

   Likewise, a third time scheme needs attention to eligibility/history, coefficient calculation, rate evaluation, and the binary `use_bdf2` log label.

   **Fix now:** describe those extension points accurately. For example:

   > Additional schemes require updates to selection, history eligibility, coefficients, rate evaluation and labels. Additional sinks require mode decoding, weights/reductions, fallback handling, sink evaluation and labels.

   When introducing another sink, replace the binary label expression with an explicit `EL_VA_SINK_NAME(mode)` selector. A larger abstraction is unnecessary for this port.

5. **Nit — the opening compatibility description is stale.**  
   [el_volume_averaged.f90:26](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_el/el_volume_averaged.f90:26)

   It still says the volume mean is always removed through one Allreduce. The later text correctly describes configurable removal and additional reductions.

   **Fix:** replace that passage with:

   ```fortran
   ! Compatibility: the raw source integral is removed through the configured
   ! sink distribution; uniform mean removal is the default (see below).
   ! Reductions use the worker communicator; the CFD master does not participate.
   ```

   Also wrap the new 118-character validation message at [el_config.f90:283](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va/source/src_el/el_config.f90:283) to follow the repository’s 100-character convention.

What I verified by reading:

- **Defaults and BDF2/uniform:** vertex averaging, epsilon clamp, rate evaluation, advection, raw source, each core accumulator’s expression and element summation order, volume safeguard, and `gbar` division are preserved. Moving alpha extraction before the first core accumulators introduces no dependency change. The uniform helper returns `gbar`, preserving subtraction followed by multiplication by volume and division by `dt`. The parent module also matches `cafa91bc`.
- **Collectives:** the core remains one eight-double `MPI_SUM` Allreduce. The temporal block remains two doubles on audit steps. Uniform adds no collective. Alpha adds its two-double weight reduction and audit-only scalar corrected-integral reduction.
- **Alpha arithmetic:** `int_g*alpha_e/w_alpha`, `alpha*vol`, `zc*alpha*vol`, and the corrected-integral accumulation retain the reference expression structure and element order. No reciprocal precomputation or regrouping was introduced.
- **Fallback and MPI:** the weight reduction occurs before fallback. `W` is globally reduced, so workers select the same effective mode. The final alpha audit reduction tests the deck request, so fallback retains participation and diagnostics. This assumes the existing entry preconditions hold; the unchanged rank-local early returns still require consistent allocation and valid partitions.
- **Logging:** existing formats and arguments are preserved. Current records end in a numeric field or trimmed tag/name, so `TRIM(line)` removes buffer padding without removing meaningful trailing bytes. A conservative maximum is **191 characters** with 32-bit integers; `LOG_LEN=320` is ample. Both padded `MERGE` expressions have equal-length alternatives and are valid.
- **Deck switch:** the parser assigns the key, validation lowercases its value, accepts only `uniform|alpha`, and rejects alpha outside `continuity`. Omission defaults to uniform. Both output destinations echo it. The sample differs from the BDF2 deck only by the alpha key.
- **Correctors and leftovers:** predictor, corrector and final diagnostic assemblies share the sink path. I found no unused declarations/imports or unnecessary new public symbols. Reusing `alpha` across the two loops is harmless; naming the corrected source `g_corr` would improve readability.

**Verdict:** sustainable and materially more maintainable, with the small cleanups above. No numerical default side effect is apparent from the source, although extra local diagnostic work exists and the configuration echo adds a line. It is truly deck-switchable, including correctors. Strict alpha stream identity currently fails on record ordering.

I compared the parent, reference implementation, configuration, parser, callers and sample deck. I did not compile, run MPI cases, inspect the running gates, or certify compiler-generated bitwise results. No files were modified.
## Claude Fable 5.1 (xhigh)

Read-only review of cef8eeaa; nothing built or run. Line numbers from the post-commit files.

### Blockers
None.

### Should fix

**S1. A collective is buried inside the audit-lines section** (`el_volume_averaged.f90:307-310`). The `int_gc` Allreduce sits under `IF (audit) THEN ... IF (EL_VA_SINK_MODE().EQ.SINK_ALPHA)` between three `EL_VA_LOG` calls, in a block whose heading is "audit lines". Every other dual-write block in `src_el` (`el_diagnostics.f90:309,580,655,730`) is wrapped in `IF (myid.EQ.showid)`; the first maintainer who applies that house pattern to this block deadlocks the workers in alpha mode. Today it is correct (all workers evaluate `audit` and the deck mode identically), but the layout invites the bug. Hoist it to the reductions, directly after the corrected-source loop; the collective sequence and the bytes are unchanged.

**S2. Deck mode and effective mode share one variable** (`el_volume_averaged.f90:201, 210, 250, 263, 271-277, 307, 315`). `mode` is the deck mode until line 273 rebinds it to the run-time fallback; afterwards the audit block re-asks the deck via a second `EL_VA_SINK_MODE()` call while the label uses the rebound `mode`. The semantics are right (collectives keyed on the deck, label on the effective mode, exactly as 8aeb26c8), but a reader has to reconstruct which `mode` each test means. Two names and one decision point: `mode_deck = EL_VA_SINK_MODE()`, `mode = EL_VA_EFFECTIVE_SINK_MODE(mode_deck, w_alpha)` with the fallback rule in that function (W is global, so every worker takes the same branch and no collective selection changes).

**S3. A third sink mode is not yet a local change** (`el_volume_averaged.f90:88-94, 139-156, 250-254, 283-285, 315`; `el_config.f90:274-279`). A hypothetical weighted mode needs the name in the EL_CONFIG CASE list, the mapping in EL_VA_SINK_MODE, a new CASE plus a new argument in EL_VA_SINK, new accumulators and a per-element array in the element loop, the placeholder dance in the corrected loop (`alpha = 0; IF (ALLOCATED(ak)) alpha = ak(iel)`, `alpha` reused as a dummy), and the two-way MERGE label that cannot take a third name. Six places in two files; the docs line "one function case plus a config entry" is only true for the time scheme. The alpha-specific plumbing is really "weighted sink" plumbing: one name table (deck value = label, validated against the same table), `EL_VA_SINK_WEIGHT(mode, alpha_e)` (one CASE per weighted mode), a generic `wk(:)` weight array with `W_SUM`/`W_ZSUM` accumulators, and `EL_VA_SINK` with `CASE (SINK_UNIFORM) gbar; CASE DEFAULT int_g*wk/w_total` (expressions verbatim). With this, a new mode is one CASE plus one table entry.

**S4. The alpha gate cannot see the new stream, and the line order differs from the test branch** (`tools/el_stream_compare.py:7-10`; `el_volume_averaged.f90:291-318` vs 8aeb26c8 lines 225-260). TAGS has no EL_VA_COMPAT, so "alpha keys identical to the T9a/T9b logs" is a statement about SOURCE/MOMENT/TSRC/RESIDUAL only; add the tag. 8aeb26c8 wrote EL_VA_COMPAT before EL_VA_SOURCE; this commit writes it after EL_VA_TSRC. Per-tag comparison is unaffected, a raw diff of the two logs is not byte-identical even when every stream is; say so in the record. (The int_g_corr field is a ~1e-21 cancellation residual; report it, do not gate on it.)

### Nits
- Header lines 26-31 still say "the volume mean of g is therefore removed"; now mode-dependent. The "Collectives" paragraph omits the alpha-and-audit int_gc reduce.
- LOG_LEN = 320: ample (longest line EL_VA_COMPAT about 190 bytes with ' corr= -1'); comment which line sets it. Overflow of an internal WRITE is a runtime error, not silent truncation.
- EL_VA_LOG beside EL_VA_LOG_CONSTRAINT_RESIDUAL: the prefix now means both "write a line" and "compute and log"; EL_VA_WRITE_LINE matches the EL_WRITE_* verbs of src_el. The dual-write pattern is hand-copied five times in el_diagnostics.f90; the helper belongs in a shared place in a follow-up, not in this commit.
- Padded-literal MERGE only ever holds two names; the mode one goes away with the table.
- EL_VA_SINK, EL_VA_SINK_MODE can be PURE.
- The residual line is formatted on every worker and discarded on all but showid; negligible.
- Each corrector pass re-reduces W (2 doubles) although alpha_p is unchanged between passes; harmless, keeps the routine pass-agnostic.
- Pre-existing: the rank-local early returns (nel.LE.0, SIZE checks) precede the collectives; EL_LOWERCASE exists in both EL_GEOMETRY and EL_CONFIG.

### Checklist
1. Sustainability: clear improvement over local(1..13) and the twelve duplicated WRITEs; slot names and the four sections read well; a third time scheme is one CASE plus config. A third sink mode is not yet local (S3); the single mode variable and the collective in the log block are the two things a maintainer will trip over (S1, S2).
2. Bitwise at the defaults, verified by reading: raw source, gk and the eight core accumulators are textually the same statements in the same order; gbar as before; g = gk - gbar then defp + g*vol/dt is the same operation sequence as (gk-gbar)*vol/dt. Build is gfortran 14.3 -O2 -march=znver1 (FMA contraction active): the add's operand is a quotient in both versions, so no contraction either way; the new in-loop statements share no product with an existing one that could change a contraction decision. Collectives at the defaults: one 8-entry Allreduce plus the 2-entry tsrc on audit steps, exactly as before; the alpha reduces are unreachable with uniform. Logs: no item ends in a blank, so TRIM(line) strips only buffer padding; '(A)' to the same units in the same order gives identical bytes. BDF2 + correctors + uniform takes exactly the pre-commit path.
3. Alpha mode vs 8aeb26c8: accumulators alpha*vol, zc*alpha*vol; fallback test w_alpha.LT.1d-30; z_alpha = sum/w_alpha; sink int_g*alpha_e/w_alpha (same association as global(1)*ak(iel)/w_alpha); int_gc + g*vol as local(13) + gc*vol; COMPAT fields identical; same fallback label. Collective order core, tsrc, alpha, int_gc in both. Only the line order differs (S4).
4. MPI: audit is a function of istep, the deck mode of a global config string; both identical on all workers; the fallback is taken on the reduced w_alpha and selects no collective. Master returns at the top; corrector and final passes call the routine inside IF (myid.ne.0) on every worker. Verified; the only issue is layout (S1).
5. Deck switch: default 'uniform', lower-cased, validated, rejected with the other keys when ELVolumeAveraged /= continuity, echoed to protocol and terminal; app_init.f90 is the only EL-key parser; EL_VALIDATE_CONFIG runs on all ranks. Predictor, correctors and the final pass all use the alpha sink. Sample deck differs by the one key.
6. Leftovers: no unused variables, imports or PUBLIC symbols; no module cycle; -Wall should be quiet.

### Verdict
Sustainable: yes, clearly better than the parent, and the three should-fix items are small local edits that would make the module genuinely one-CASE extensible for sink modes. Side effects at the defaults: none found by reading. Truly deck-switchable: yes. Not verified: any run; the gates confirm the FMA/vectorisation reasoning; the COMPAT stream is gated only once the tag is added.
