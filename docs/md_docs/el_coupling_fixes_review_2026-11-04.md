# Review of the P4/P5 coupling fixes (commit 5404e693, 2026-11-04)

Read-only review by Codex (gpt-6-astra, high) of the first fix commit; brief: `el_coupling_fixes_review_brief_2026-11-04.md`. The follow-up commit applies its should-fix items.

P4 and P5 are correctly repaired at the feedback-construction level. I found no blocker introduced by these two arithmetic changes, but the diagnostics and several claims need correction.

1. **Should fix — the existing fluid momentum audit omits semi-implicit drag, and its compensator can alter the physics.**  
   [el_diagnostics.f90:570](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va-coupling/source/src_el/el_diagnostics.f90:570)

   `expected = dt*(global_src + ext_force*global_vol)` counts the captured explicit source but omits the drag operator’s contribution. In semi-implicit mode, the expected coupling impulse must also include
   \[
   -\Delta t\sum_i D_i[(1-\theta)u_i^n+\theta u_i^*],
   \]
   where `D = drag_B_source` and \(u^*\) is the momentum-stage velocity before projection.

   Consequently, `EL_FLUID_PAIR` reports legitimate drag as a momentum mismatch. With `ELMomentumFix=Yes`, lines 572–574 feed that mismatch into a compensating acceleration, tending to cancel legitimate drag.

   **Concrete fix:** include both drag contributions in the audit using the captured coefficient and the corresponding old and momentum-stage velocities. Until then, reject `ELMomentumFix=Yes` with semi-implicit drag. This is **pre-existing**, but directly affects interpretation and use of the repaired P4 path.

2. **Should fix — `EL_FLUID_SOURCE` does not always describe an applied source; the post-pass line remains a deposition audit.**  
   [el_transfer.f90:428](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va-coupling/source/src_el/el_transfer.f90:428), [el_transfer.f90:445](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va-coupling/source/src_el/el_transfer.f90:445)

   The new line prints in one-way mode, although the fluid never consumes that source. It also prints with a prescribed fluid field, when the fluid solve is skipped entirely. In semi-implicit mode, its `rhs` is only the captured explicit force vector—not the complete fluid coupling force or timestep-scaled impulse.

   The post-pass `EL_FEEDBACK_CONSERVATION` remains mathematically valid as a **spread/deposition conservation check**. For frozen semi-implicit particles it reports the full old-velocity force, whereas the captured source contains the compensation. Its small residual therefore does not certify the applied fluid reaction.

   **Concrete fix:** label the new quantities `phase=pre`, `rhs_explicit`, and `applied=yes/no`, with application status determined by both feedback enablement and whether the fluid solve runs. Label/document the post-pass line as `phase=post`, `diagnostic_only`, and `deposition_residual`. Report the full applied coupling impulse after the momentum solve. Do not change the post-pass force merely to make the two lines agree.

3. **Should fix — “No = frozen particles” is not the implemented behavior.**  
   [el_config.f90:25](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va-coupling/source/src_el/el_config.f90:25)

   `ELApplyParticleForces=No` skips force assignment and hydro-state arming. It does not skip `EL_ADVANCE_PARTICLES`: [QuadSc_main.f90:394](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va-coupling/source/src_quadLS/QuadSc_main.f90:394) still advances PE, whose position integrator updates positions from velocity.

   Thus a mobile particle with nonzero initial velocity can drift; contacts can also affect its motion. This matters for interpreting the “frozen” control and configuring E1.

   **Concrete fix:** describe the key as disabling application of the evaluated hydrodynamic forces. Document the separate PE configuration needed to hold particles stationary, and verify that configuration in the frozen-bed controls. P4’s algebra remains correct for the sampled particle velocity.

4. **Nit — comments overstate conservation and bitwise log preservation.**  
   [el_transfer.f90:279](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va-coupling/source/src_el/el_transfer.f90:279), [el_config.f90:50](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va-coupling/source/src_el/el_config.f90:50)

   The stated equality
   `-drag_force_eff - f_explicit = B*u_f - drag_force_eff - f_explicit - B*u_new`
   requires `u_new=u_f`. Generally, the split introduces an additional drag correction. That behavior predates this patch, but the updated comment retains the false equality.

   Also, “keeps every certified log bitwise” is false literally: the new audit line appears under default diagnostic settings at the audit cadence. The pre/post difference is not specific to pressure feedback; effective drag and changed particle positions already produce differences.

   **Concrete fix:** describe the semi-implicit term as a linearization around the sampled carrier velocity, state its residual correction, and limit preservation claims to existing numerical outputs, excluding the new line.

**Verified signs, scaling, and pass selection**

[EL_DEPOSIT_PARTICLE:596](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va-coupling/source/src_el/el_quadrature.f90:596) deposits **minus** `feedback_used` with normalized basis weights. It deposits `drag_B` with the same weights. The pre-pass captures both arrays before parallel summation; resetting the post-pass working arrays does not overwrite those captures.

[EL_APPLY_FLUID_FEEDBACK_SOURCE:131](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va-coupling/source/src_el/el_fields.f90:131) adds **`+dt*fluid_feedback_source`**, and its caller passes `tstep`. There is no additional density division or sign reversal there.

One correction to the brief: [QuadSc_def.f90:2516](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va-coupling/source/src_quadLS/QuadSc_def.f90:2516) is an **old-velocity RHS subtraction**, not another diagonal addition. The implementation uses:

- RHS: `-tstep*(1-theta)*D*u_old`;
- matrix: `+tstep*theta*D`;
- captured explicit source: `+tstep*source`.

For P4, writing \(E=L+P\) when pressure feedback is enabled,
\[
feedback\_used=B(u_f-u_p)+E-Bu_f=-Bu_p+E.
\]
Its negative deposit supplies \(Bu_p-E\), completing the theta-weighted drag operator. At `theta=1`, this is the backward-Euler split in the brief. **The compensation belongs only in the captured pre-pass.**

**Pressure and PE consistency**

[EL_SEMI_IMPLICIT_DRAG_IMPULSE:439](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va-coupling/source/src_el/el_forces.f90:439) computes
\[
D_{\rm eff}=\frac{m(u_{\rm end}-u_{\rm start})}{\Delta t}-f_{\rm other},
\qquad f_{\rm other}=P+L+G.
\]

I also read PE’s hydro-state setter, substep loop, and reference-velocity update. They implement the corresponding recurrence. Pressure already shapes PE’s velocity through `f_other`; adding it there again would double-count it. P5 correctly adds it to the **fluid reaction** through `f_explicit`.

The pressure channel now has equal-and-opposite force accounting, up to floating-point and deposition error, for advanced particles with active fluid feedback. This does **not** establish exact total Newton-pair conservation in semi-implicit fluid mode: its additional velocity-dependent drag correction remains. Frozen particles exchange momentum with their constraint/support rather than gaining the corresponding particle momentum.

I found no remaining pressure-dropping rebuild in an actively solved, two-way explicit or semi-implicit path. Remaining ineffective combinations are:

| Configuration | Effect |
|---|---|
| `ELApplyFluidFeedback=No` | No reaction reaches the fluid. |
| `ELPressureForce=No` | Pressure force is zero, so pressure feedback adds zero. |
| Prescribed fluid field | Fluid solve is skipped. |
| Nonnegative transfer pass | Diagnostic deposit only; no source capture. |

Warnings for the first three combinations when pressure feedback is explicitly requested would prevent confusion; they need not all be fatal. The existing rejection of semi-implicit drag without fluid feedback is appropriate.

**Maintainability and determinism**

`f_explicit` is usable, though `feedback_nondrag` would distinguish it more clearly from PE’s `f_other`. A shared helper for selecting non-drag reactions would reduce the duplication that caused P5. I would keep pass selection explicit rather than force all three assignments into one opaque function. Any refactor must preserve addition order; avoid reconstructing non-drag feedback through subtraction.

No new MPI collectives, communicator changes, ownership changes, or reduction-order changes were introduced. The new log reuses existing reductions. Both new arithmetic paths are reachable; I found no new dead branch.

**Verdict:** P4 and P5 are fixed as scoped. With pressure feedback off and the frozen semi-implicit branch excluded, the four stated regression configurations retain the same source-level floating-point arithmetic; copying `lift` into `f_explicit` adds no arithmetic. That supports numerical equivalence, but does not prove compiled-run bitwise identity. Full logs necessarily change at the new audit cadence.

The patch is small and readable, but its diagnostic semantics and comments need cleanup. I verified the source paths and PE recurrence by reading, and `git diff --check` passed. I did not build, run MPI cases, or verify the pending regression results. No files were modified.