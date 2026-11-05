# Conservative continuity form: implementation reviews (2026-11-05)

Reviews of commit dc2221e0 (parent e7ae87fc) on feature/el-va-conservative-form, requested with the
owner's sustainability requirement and the design consultation's findings as the checklist. Brief,
then both answers verbatim. The fixes are in the commit after dc2221e0 (see
el_va_conservative_form_2026-11-05.md, section Reviews).

## Brief

```
# Review brief: conservative continuity form (commit dc2221e0 on feature/el-va-conservative-form)

Working directory: /data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va-formc (FeatFloWer, Fortran 90/2008,
gfortran 14.3 -O2 -march=znver1, OpenMPI 4.1.6, CMake/ninja). Review the commit dc2221e0 (parent e7ae87fc):
`git show dc2221e0 --stat`, `git diff e7ae87fc dc2221e0`. Do not modify any files.

## What was built (read docs/md_docs/el_va_conservative_form_2026-11-05.md first)

A deck-switchable third axis of the Euler-Lagrange volume-averaged continuity constraint,
`SimPar@ELVAContinuityForm = advective | conservative` (default advective = the campaign code):
div(eps u) = -d eps/dt in the weak form of all four P1 pressure test functions, with shared face fluxes
of a Q1-reconstructed continuous transport porosity, evaluated by a new MPI-free module
(source/src_el/el_va_flux.f90) and fed as the direct residual right-hand side -r_m/(eps_e dt) (row 1 minus
the compatibility sink) into the existing projection, which the existing corrector loop
(source/src_quadLS/QuadSc_el_va.f90) iterates to the fixed point. The design was consulted with GPT-6 Astra
before coding (docs/md_docs/el_va_conservative_form_design_review_2026-11-05.md: brief + answer verbatim);
the design note lists which recommendations were adopted. Files: el_va_flux.f90 (new), el_volume_averaged.f90,
el_config.f90, el_fields.f90, el_transfer.f90, QuadSc_el_va.f90, app_init.f90 (parser), ProjectFiles.cmake,
tests/test_el_va_flux.f90 (new unit test, passes), tools/el_stream_compare.py, docs, sample deck.

Owner requirements (hard): (1) the DEFAULT configuration must be bitwise identical to the parent commit -
no change in arithmetic, operation order, or MPI collective pattern on the advective path (the legacy element
loop and the legacy defp update are kept verbatim; new collectives only under the deck-keyed form); gates are
running on the cluster; (2) new code must be clean, maintainable, extensible (sustainable codebase); (3) MPI
master/worker discipline: rank 0 (master) holds no mesh/EL fields and joins only the collective pressure
solves; collectives are keyed on deck values and the step number, never on rank-local state; no rank-local
early return in front of a collective.

Smoke evidence (VB_L2 box, 11 ranks, 20 steps, job 157076): face flux = B divergence to 1e-17, boundary
flux and the three moment identities ~1e-20, unsinked residual 1.1e-3 -> 8.7e-6 -> 1.5e-7 over the
predictor + one corrector, eps-weighted sink integral -3e-11; reconstructed-eps superficial velocity
4.2418e-4 vs the exact 4.24115e-4. A bug found and fixed during the smoke: the reference-face normal must
carry sgn(det J) because FEAT accepts left-handed hexes (it only uses |det J|); the unit test now has a
mirrored cell.

## What to check (findings ordered by severity, each with file:line, the problem, why it matters, the fix)

1. Correctness of the discrete scheme as implemented against the design note: the residual r_m = F_m + d_cont Q_m,
   the right-hand side defp_m := -r_m/(eps_e dt) (row 1: - s_e |e|/dt), the sink bookkeeping gV_1 = D_1 - r_1/eps_e,
   the sign chain (defp = -(1/dt) B^T u at entry; Velocity_Correction gives D du = dt rhsP), the storage moments
   Q_m with dvol, the flux kernel (basis numbering vs E013A, outward normal with sgn(det J), weights), the
   reconstruction (vertex volume-weighted mean, E013Sum2 under the pinned level, Q1 embedding).
2. Bitwise-at-default by reading: is there ANY statement on the advective path whose arithmetic or order changed?
   Any new collective reached at the default? (Compare the legacy loop and the legacy defp block with the parent.)
3. MPI safety: the transfer-pass reconstruction (every worker, with or without particles), the flux block
   reductions (audit steps only), the eps-weighted sink reduce, the comm-scheme guard in QuadSc_el_va.f90
   (does it fire before any collective on all ranks?), the layout self-check ERROR STOP (rank-local; is that
   acceptable, or must it be made uniform?), allocation/size early returns.
4. Sustainability: module boundaries (EL_VA_FLUX MPI-free and stateless beyond the reference tables),
   the form dispatch, table-driven keys, naming, comments, the unit test's coverage; anything a maintainer
   would trip over; dead code; the interplay with the existing sink modes and the BDF2 rate.
5. Known limits to confirm or extend: closed impermeable boxes only (periodic/outflow: KNPRP rows,
   seam z-moment), no tolerance-based stopping, dvol single precision, el_eps_f_min clipping and relaxation
   make S nonzero, ELMesoFilter rejected, restart (transport eps rebuilt before first use?).
6. Performance: the kernel runs 3x per step at audit frequency 1 (predictor, corrector, final audit pass);
   the smoke took 1:45 vs 0:37 for 20 steps. Obvious cheap wins without changing results?

End with a one-paragraph verdict: ready for the gate ladder / needs fixes first, and what you did not verify.
```

## Codex (gpt-6-astra, reasoning effort high)

Reviewed `dc2221e0` against `e7ae87fc`, including both design notes and the projection call chain. No files were modified.

1. **High — conservative mode retains rank-local exits before worker collectives.**  
   [el_volume_averaged.f90:292](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va-formc/source/src_el/el_volume_averaged.f90:292) returns on allocation state, `nel <= 0`, and undersized arrays before the reductions. The predictor also has an allocation-dependent guard at [QuadSc_main.f90:801](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va-formc/source/src_quadLS/QuadSc_main.f90:801). A worker taking either exit can leave its peers in unmatched collectives. These are **inherited defects**, but the conservative path still violates the stated hard requirement.

   **Fix:** under the conservative deck key, enter from every worker, validate readiness collectively before accessing arrays, let empty workers contribute zeros, and abort coherently on inconsistent state. Preserve the existing advective path. Include transport-porosity allocation and size in this validation.

   The new rank-local `ERROR STOP`s at [el_va_flux.f90:151](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va-formc/source/src_el/el_va_flux.f90:151) and line 418 also lack explicit MPI-wide termination. OpenMPI normally terminates the remaining processes after an abnormal exit, but that is launcher behavior. Return an error status from the MPI-free module and handle it with coordinated validation or `MPI_Abort` in the caller.

2. **Medium — the communication-scheme guard runs too late.**  
   [QuadSc_el_va.f90:77](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va-formc/source/src_quadLS/QuadSc_el_va.f90:77) executes only after both transfer passes, conservative predictor assembly, and the first pressure solve. Reconstruction has already called `E013Sum2` at [el_transfer.f90:369](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va-formc/source/src_el/el_transfer.f90:369). For unsupported schemes, `E013Sum2` silently performs no exchange, producing partition-dependent porosity before rejection.

   **Fix:** validate the scheme on every solver rank during initialization, before the first transfer pass. The current guard is rank-uniform for the shared configuration, but it does **not** fire before communication or use of the reconstructed field.

3. **Medium — `m1rate` omits part of the implemented storage first moment.**  
   [el_volume_averaged.f90:575](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va-formc/source/src_el/el_volume_averaged.f90:575) accumulates only `zc*d_cont*vol`. The four-row equations instead imply the storage contribution
   \[
   \sum_e d_{\mathrm{cont},e}\bigl(x_{c,z,e}V_e+Q_{z,e}\bigr).
   \]
   `Q_z` is generally nonzero on warped trilinear cells. An independent calculation using the existing warped-cell test coordinates gives `Q_z/V ≈ 3.00069e-4`. Thus the diagnostic can disagree with a correctly converged flux solely because it reports an incomplete target.

   **Fix:** accumulate `d_cont*(xc(3)*vol + q(4))`. If the existing deposition-centre moment remains useful, report it separately. For incompatible storage, the flow identity also needs the sink first moment; a small global `sink_eps` alone does not measure that contribution.

4. **Low — the layout check imposes a geometry restriction and retains mesh-dependent state.**  
   [el_va_flux.f90:106](/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va-formc/source/src_el/el_va_flux.f90:106) tests whether higher-order coordinates equal vertex averages. That validates trilinear geometry as well as numbering: valid curved Q2 elements fail with a misleading “dof layout” error. Meanwhile, the process-global `layout_checked` flag at line 102 skips validation for subsequent meshes or levels.

   **Fix:** validate numbering against connectivity/FEAT mapping independently of coordinates. If trilinear geometry is an intentional restriction, name and document it explicitly. Keep validation state with the mesh lifecycle rather than inside the numerical module. Apart from this flag, the reference-table cache and pure kernels provide the intended stateless separation.

The discrete scheme itself checks out by inspection:

- `r_m = F_m + d_cont*Q_m`, `gV_1 = D_1-r_1/eps_e`, and all four direct `defp` updates match the design. The sink is applied only to row 1.
- `Build_BTMatP1` assembles positive divergence moments; `BT_Mul_U_mod` applies `-1/dt`; `Velocity_Correction` consequently gives `D Δu = dt rhsP`. The implemented signs are consistent.
- Storage moments consistently use `dvol`, with geometric moments normalized by geometric volume. The distinction is deliberate and correctly implemented.
- Basis numbering, cyclic tangent pairs, `sgn(det J)`, face weights summing to four, and volume weights summing to eight are consistent. Subtracting `eps_u` from the three linear face moments gives the required weak form.
- Reconstruction uses volume-weighted vertex contributions, pins global `ILEV` around `E013Sum2`, and embeds the resulting Q1 field into Q2. `/TRIAD/` is already pinned for the transfer pass. Workers without particles still execute reconstruction.
- Uniform/alpha compatibility and the globally determined zero-weight fallback are preserved. BE/BDF2 rate selection is shared with the advective form.

For the default requirement, I performed a byte comparison of the extracted legacy blocks: the **40-line raw-source loop and eight-line `defp` update loop are identical to the parent**. I found no changed ordering or arithmetic feeding the legacy numerical results, and no new collective reachable with the advective form. New form checks and unused diagnostic initializations execute, and configuration output gains a line. This supports numerical bitwise preservation by reading; compiler-generated behavior still requires the running gates.

Outside the early-exit issue, the MPI separation is correct: master avoids mesh reconstruction and source evaluation, joins the pressure solves, and follows the deck-controlled corrector count. Flux reductions and the eps-weighted sink reduction are gated by form and audit step; weighted-sink reductions use the requested deck mode.

The implementation’s boundaries and remaining limits are:

- **Closed impermeable domains:** confirmed. Pressure boundary handling zeros all four rows for `KNPRP=1`; open-boundary constraints therefore need separate treatment. Periodic coordinate moments require the seam contribution. These restrictions are documented rather than comprehensively enforced.
- **Convergence:** one corrector is an eligibility rule, not a convergence guarantee. Fixed counts are suitable for the ladder if acceptance uses final unsinked residuals for constant and linear rows.
- **Storage compatibility:** `dvol` carries single-precision rounding from `REAL(AAA)`. Clipping can create nonzero total storage rate; relaxation changes the rate and its relation to particle motion. At convergence, `sink_eps = -(S+Phi)`, not necessarily zero. The comment at line 427 saying it “must vanish” should be qualified.
- **Restart:** transport porosity is rebuilt in both transfer passes before the fluid solve reads it. Production restart does not currently restore EL rate history through the version-3 helper—the only callers are tests—so it starts with BE, as the existing documentation states.
- **Other combinations:** `ELMesoFilter` rejection is appropriate. `GAMMA>0` remains outside the initial certification scope because pressure correction applies its additional update on every pass.

The unit test provides useful kernel coverage, including mirrored orientation and warped-cell cancellation. It does not exercise production `E013Sum2`, differently oriented neighboring connectivity, source/sink assembly, or the complete projection iteration. Before broad certification, add focused coverage for a real partition containing a particle-free worker, four-row manufactured residuals with BE/BDF2, and deliberately incompatible storage. The reported smoke provides integration evidence, but does not replace those cases.

The clearest performance opportunities are:

- Skip the pre-advance transport reconstruction: its result is overwritten by the post-advance pass before any consumer reads it.
- Pass `audit` into conservative evaluation and avoid diagnostic-only sums, cube roots, unweighted velocity integrals, and divergence checks on other steps.
- Reuse scratch allocations; cache geometry and porosity quadrature values across predictor/corrector/audit evaluations with explicit invalidation.
- Evaluate face traces using their nine nonzero basis functions. Retain the geometry derivatives needed for the determinant sign; blindly restricting the entire Jacobian calculation to those nine nodes would be incorrect.

These can preserve numerical accumulation order. Preassembling velocity-to-flux coefficients offers further savings but can change rounding.

**Verdict:** needs fixes first to satisfy the stated MPI contract, particularly conservative-mode readiness handling and the late scheme guard; also correct `m1rate` before using warped meshes in the ladder. The closed-box discretization and default-path preservation look sound by inspection. I did not rebuild, rerun the Fortran tests or MPI smoke, inspect cluster gate results, or verify restarted-versus-continuous trajectories. Validation here comprised source tracing, byte comparison of the legacy loops, an independent warped-cell moment calculation, and a clean `git diff --check`.
## Claude Fable 5.1 (reasoning effort xhigh)

## Summary

No blockers. The default (advective) path is statement-for-statement identical to the parent
at source level, and every new collective is keyed on the deck form and the step number. The
conservative residual algebra, the sign chain, the flux kernel (numbering, face tables, signed
normal, weights, geometry map) and the reconstruction match the design note and the FEAT code
they rely on. The findings are hygiene and placement items that do not affect the gate ladder.

## Findings

### Blockers

None.

### Should-fix (none of these affects the gates; 1, 2, 3 are worth doing before merge)

**1. The communication-scheme guard sits after the call it guards.**
`source/src_quadLS/QuadSc_el_va.f90:77`
```
IF (EL_VA_CONSERVATIVE_FORM() .AND. iCommSwitch.NE.3 .AND. iCommSwitch.NE.4) THEN
```
is evaluated in `EL_VA_Projection_Correctors`, i.e. after both transfer passes of the first
step have already called `E013Sum2` (`source/src_el/el_transfer.f90:369`) and after the
predictor source has consumed the result. `E013Sum2` has no ELSE branch
(`source/src_quadLS/QuadSc_mpi.f90:2918-2928`: `if (iCommSwitch.eq.3) ... if (iCommSwitch.eq.4) ...`),
so for any other mode it is a silent no-op and the first predictor would run on an un-summed
reconstruction before the stop. Why it matters little today: `iCommSwitch` is a hard-coded
default (`source/src_quadLS/QuadSc_var.f90:32: INTEGER :: iCommSwitch=3`) and is assigned
nowhere else, so the guard cannot fire; no deadlock is possible either (the stop is uniform,
deck-keyed). Fix: put the same check inside the deck-keyed block of the transfer pass before
the `ALLOCATE(eps_num...)` at `el_transfer.f90:365` (one `USE var_QuadScalar, ONLY: iCommSwitch`),
or keep the corrector-loop check and add a comment that it is defensive only.

**2. Transport porosity is passed to `SIZE()` without an allocation check.**
`source/src_el/el_volume_averaged.f90:515-516`
```
CALL EL_VA_FLUX_MOMENTS(mesh, ilev, mesh%level(ilev)%dcorvg, &
  el_field_data%epsilon_transport_q2, val_u, val_v, val_w, ...)
```
and `source/src_el/el_va_flux.f90:416` evaluates `SIZE(eps_q2)`, which is undefined for an
unallocated allocatable. Today it is always allocated: `EL_VOLUME_AVERAGED_ACTIVE()` becomes true
only through `EL_ENSURE_FIELDS` (`el_transfer.f90:154`), the same pass reaches the reconstruction
block at line 363 with no rank-local RETURN in between (the only returns in the routine are
line 139 for the master and the PE_SERIAL_MODE abort), and both passes precede the fluid core
(`source/src_quadLS/QuadSc_main.f90:374-406` then `:435`). A maintainer who adds a RETURN in
that window or calls the source from a new site gets undefined behaviour instead of a message.
Fix: one line at the top of `EL_VA_CONSERVATIVE_RAW_SOURCE`:
`IF (.NOT.ALLOCATED(el_field_data%epsilon_transport_q2)) ERROR STOP 'EL_VA: transport porosity not built'`.

**3. The dof-layout self-check is rank-local, lazy and keeps hidden state.**
`source/src_el/el_va_flux.f90:420-423`
```
IF (.NOT.layout_checked) THEN
  CALL EL_VA_CHECK_DOF_LAYOUT(mesh, ilev, coor)
  layout_checked = .TRUE.
```
runs on the first flux evaluation, inside the fluid core of step 1, and `ERROR STOP`s on the
offending worker (`el_va_flux.f90:151`). Because it terminates the process rather than
returning, mpirun aborts the job; there is no desync hazard, which is acceptable under the
owner's rule (no rank-local *return* in front of a collective). But (a) the message appears on
whichever worker fails first, not on showid; (b) `layout_checked` is module SAVE state in a
module advertised as stateless beyond the tables (`el_va_flux.f90:37-38`), and it is never
reset; (c) the check is per-partition, so a partition without the offending element pattern
passes silently (not realistic for a global numbering convention, but the guarantee is weaker
than it reads). Fix: run the check once in the first transfer pass under the deck guard, reduce
`nbad` with `MPI_Allreduce(MPI_SUM, MPI_COMM_SUBS)` and stop uniformly with the message on
showid; make `EL_VA_CHECK_DOF_LAYOUT` return `nbad` instead of stopping, so it can also get a
negative unit test (see finding 8).

**4. Fifteen campaign-slot statements are duplicated.**
`source/src_el/el_volume_averaged.f90:544-562` (conservative) repeats `:347-365` (legacy): the
alpha lookup, the eight `core(...)` slots, the two `tsrc(...)` slots and the weight block. The
duplication is forced by the verbatim requirement on the legacy loop for this commit, and the
copy is exact (checked line by line). Fix later, not now: after the gate ladder, factor both
into one contained helper taking `(g, adv, eps, deps_dt, d_cont, alpha, zc, vol, dt)`, called
from both loops in the same statement order, and gate it bitwise (owner's "ports restructure,
keep arithmetic verbatim, bitwise gates").

**5. The pre-advance pass rebuilds the transport porosity that nobody reads.**
`source/src_el/el_transfer.f90:361-362` documents it ("The pre-advance pass rebuilds it too ...
overwritten by the post-advance pass before the fluid step reads it"). That is one `E013Sum2`
on two Q2 arrays plus the reconstruction per step for nothing. `istep` is uniform across ranks,
so `IF (EL_VA_CONSERVATIVE_FORM() .AND. istep.GE.0)` keeps the collective pattern deck/step-keyed.
Results unchanged. Optional.

### Nits

**6. Stale "longest audit line" comment.** `el_volume_averaged.f90:152`
`! longest audit line: EL_VA_COMPAT with ' corr= -1', about 190 characters` -- `EL_VA_FLUXMOM`
is now the longest at ~272 characters (still under `LOG_LEN = 320`; an overflow would be a
runtime error on the internal WRITE). Update the comment.

**7. Forward reference to a RUNBOOK section that does not exist.** `source/src_el/el_config.f90:115`
`certified for closed impermeable boxes (s04 RUNBOOK section 12)` -- the RUNBOOK ends at 11.2.
Point to `docs/md_docs/el_va_conservative_form_2026-11-05.md` until section 12 is written.

**8. Unit-test coverage gaps** (`source/src_el/tests/test_el_va_flux.f90`; all six tests are
sound and scale-aware, and test 1 pins the numbering to `E013A` at four interior points):
(a) the helper's edge/face tables (`el_test_helpers.f90:149-152`) are copies of the module's
(`el_va_flux.f90:84-87`), so the tests cannot detect an error both share -- only the runtime
layout check does (it passed on the real mesh in the smoke); a negative test that permutes one
`kedge` entry and asserts the layout check reports it would close this (needs finding 3);
(b) the per-element residual algebra of `EL_VA_CONSERVATIVE_RAW_SOURCE` (`Q_m`, `r_m`, `gV_1`,
`adv`) is untested; factoring lines 533-543 into a PURE function would allow a manufactured
compatibility negative test (consultation §9.6); (c) only trilinear geometry is exercised
(`BUILD_DOF_COORDINATES` places Q2 dofs at vertex means) -- fine for the exactness claims, but a
curved-Q2 shared-face test (displaced mid-face node shared by two cells) would certify the
"any conforming mesh" cancellation claim at `el_va_flux.f90:17-20`.

**9. `uf_p0_z` is not the `EL_MEAN_SLIP` integral.** `el_volume_averaged.f90:574`
`fsum(F_P0U_Z) + el_field_data%epsilon_f(iel)*vol_u(3,iel)` uses the 27-point Q2-geometry rule
of the kernel, while `EL_MEAN_SLIP` uses `EL_Q1_MAP`; identical on bricks only up to rounding.
The design note (line 73-75) should say "P0-weighted, Q2 rule" so nobody expects bitwise
equality with `EL_MEAN_SLIP`.

**10. Two centres, two predicates.** In `EL_VA_CONSERVATIVE_RAW_SOURCE` the vertex-mean `zc`
(line 519-524) feeds the campaign slots and the centre node `xc` (line 525) feeds the moment
identities; they coincide on straight-edged meshes (layout check, 1e-6 h) -- say so in a comment.
`EL_VA_FORM()` (`el_volume_averaged.f90:205`, table-driven) and `EL_VA_CONSERVATIVE_FORM()`
(`el_config.f90:454`, string compare) answer the same deck question; they agree because the
validation forces `continuity` for any non-advective form, but a third form will have to touch
both. Consider a single `EL_VA_FORM_ID()` in `EL_CONFIG`.

**11. Legacy blocks inside unlabeled ELSE branches at their old indentation.**
`el_volume_averaged.f90:326` and `:412`. Intentional (verbatim diff), but a one-line comment at
each `ELSE` ("campaign loop, verbatim; indentation kept for the diff") will stop a future
re-indentation commit from polluting the verbatim trail.

**12. Constant-row-only sink is undocumented as a choice.** Header `el_volume_averaged.f90:69`
states "row 1: - s_e|e|/dt" but neither the header nor the design note records the
consultation's caveat (on non-centroidal cells the constant-row sink is not the four-row weak
form of a spatially constant sink). One sentence.

## What was verified, with the key statements

**(a) Default path bitwise, by reading.**
- Legacy element loop: parent 259-298 vs current 327-366, `diff` identical. Legacy `defp` block:
  parent 319-326 vs current 413-420, identical. Reductions block: parent 300-316 vs current
  369-392 differs only by the inserted
  `IF (conservative .AND. audit) THEN ... MPI_Allreduce(fsum...) ... MPI_Allreduce(fmax...)`.
- Early returns (`:291-299`) unchanged from the parent. `conservative` is deck-only:
  `:304 conservative = (EL_VA_FORM().EQ.FORM_CONSERVATIVE)`. New allocations only under it
  (`:322 ALLOCATE(res(4,nel), epsk(nel))`). New collectives: `:379-383`, `:430-431`, all
  `conservative .AND. audit`. Audit lines `:462-480` under `conservative`.
- Other files at default: transfer block `el_transfer.f90:363 IF (EL_VA_CONSERVATIVE_FORM())`;
  corrector guard `QuadSc_el_va.f90:77` conservative only; `el_fields.f90` adds an unallocated
  component and a guarded DEALLOCATE; `el_config.f90` adds validation and one echo line
  `EL VA continuity form     = advective` to the protocol file (not an `EL_*` stream line, so
  `tools/el_stream_compare.py` gates are unaffected; a byte diff of the protocol would differ).
- Residual risk is compiler-level only (the legacy loop now sits inside an ELSE; gfortran
  defaults to `-ffp-contract=fast` on znver1, and contraction/scheduling decisions are not
  guaranteed invariant under surrounding-code changes). The running gates settle that.

**(b) Sign chain and residual algebra.**
- Entry: `QuadSc_proj.f` `BT_Mul_U_mod`: `LAX19(B1,...,DD,-1d0/DT,0D0)` then `-1d0/DT,1D0` for
  B2, B3 -> `defp = -(1/dt) B^T u`; `Matdef_General_LinScalar(...,1)` (`QuadSc_def.f90`) calls it.
  So `D_m = -dt*defp_m = int_e q_m div u^k`: `el_volume_averaged.f90:536 d(m) = -dt*defp(4*(iel-1)+m)`.
- P1 basis of B: `source/assemblies/QuadSc_BMatrix.f` (Build_BTMatP1) `DX0I = DCORVG(1,KDFG2(27))`,
  `XX = XX + DPP(1)*HBASI1` over the 27 Q2 nodes, `DBAS1(2) = XX-DX0I`. The kernel uses the same
  centre and map: `el_va_flux.f90:357 xc = x(:,NQ2)`, `:321 xq = xq + x(:,l)*t%fval(l,p,f)`,
  `:337 moments(2:4) + (xq - xc)*un`. Linear rows are therefore consistent with B.
- Correction: `QuadSc_corrections.f90` `B_Mul_U(..., -TSTEP, 0)` then
  `valU = valU - defU/MlRhoPmat` -> `u += dt Ml^-1 B p` -> `D u_new = D u + dt rhsP` with
  `rhsP = defP` (`QuadSc_el_va.f90:102`, `QuadSc_main.f90:810`).
- Residual and RHS: `:537 r(m) = flux(m,iel) + d_cont*q(m)` with `:533-534 q(1) = vol;
  q(2:4) = vol*geom(1:3,iel)/vgeom` (storage moments with dvol, design note line 48-49);
  `:540 gv = d(1) - r(1)/eps`, `:541 g = gv/vol`, `:543 adv = -(eps*gv + d_cont*vol)/vol`
  (so `g = -(d_cont+adv)/eps` keeps the campaign slot meaning);
  `:405 defp(4*(iel-1)+1) = -res(1,iel)/(epsk(iel)*dt) - s_e*vol/dt`,
  `:407 defp(4*(iel-1)+m) = -res(m,iel)/(epsk(iel)*dt)`. Hence `D u^{k+1}_1 = D_1 - r_1/eps - s_e|e| = gV_1 - s_e|e|`
  and at the fixed point `r_m = -delta_m1 eps_e |e| s_e`, i.e. `F_m = -d_cont Q_m - delta_m1 eps_e|e| s_e`
  as in the note. Row-1 compatibility: `sum_e rhsP_1 = -(1/dt)[sum r_1/eps + int_g]` with
  `int_g = sum(D_1 - r_1/eps)` gives `-(1/dt) sum_e D_1 = -(1/dt) oint u.n = 0` on the closed
  no-slip box (B's 3x3x3 rule is exact for the degree-(3,3,3) integrand on trilinear cells),
  the same condition the legacy path relies on; `divchk` 1e-17 in the smoke confirms face flux
  = B divergence per element. Sink bookkeeping: `:403 s_e = EL_VA_SINK(effective_mode, gbar,
  core_g(C_INT_G), w_e, w_total)` unchanged; `:410 sink_eps + epsk(iel)*s_e*vol`.
- `res` holds the unsinked `r` (`:539`), and the `F_RES1_*`, `F_RESL_L1` norms use it (`:566-569`).
  Identities: `:571 xc(m)*flux(1,iel) + flux(1+m,iel) + eps_u(m,iel)` = `oint_e x_m eps u.n`,
  telescoping to the boundary; `flux(2:4)` already carries `-eps_u` (`el_va_flux.f90:391`).

**(c) Flux kernel.**
- Numbering: `node_ijk` (`el_va_flux.f90:65-72`) checked against every `DHELP(l,1,·)` expression
  of `E013A` (`source/Elements/e013.f`, the file in `ProjectFiles.cmake:165`): e.g. node 9
  `Q8*(2-2X1^2)*X2(1-X2)*X3(1-X3)` = 1 at (0,-1,-1), node 21 `-Q8*(2-2X1^2)(2-2X2^2)X3(1-X3)` = 1
  at (0,0,-1), node 25 at (-1,0,0). All 27 agree; test 1 pins values and derivatives numerically.
- Dof gather `:282-291` = NDFGL layout for element 13 (`kvert`, `nvt+kedge`, `nvt+net+karea`,
  `nvt+net+nat+iel`), and B uses `KDFG2(27)` the same way. Coordinates: level `dcorvg` is a shared
  pointer to the finest level (`source/src_mesh/mesh_refine.f90:1257
  mgMesh%level(ilevel)%dcorvg => mgMesh%level(maxlevel)%dcorvg`), workers run with NLMAX+1 levels
  (`applications/q2p1_el_pipeflow/app_init.f90:201`), so `dcorvg(:,nvt+...)` on NLMAX are the Q2
  dof coordinates; the layout check verifies the midpoint property at runtime.
- Faces: `face_axis/face_coord` (`:75-76`) and `face_vertices` (`:86-87`) agree with the node
  coordinates of dofs 21-26; `face_ta/face_tb` (`:78-79`) satisfy `e_a x e_b = e_axis` for all six
  (e.g. face 2 axis 2: (3,1), `e3 x e1 = e2`). Edge table `:84-85` matches dofs 9-20.
- Normal: `:329-332` `nds = jac(:,a) x jac(:,b)` (components checked), then
  `face_sign(f)*SIGN(1.0d0, detj)*t%fw(p)*nds`. With `J e_a x J e_b = det J . J^{-T} e_axis` and
  `n dS = face_sign grad(xi_axis) |det J| dxi_a dxi_b` this is the outward `n dS` for either
  handedness; the mirrored-cell test (`test_el_va_flux.f90:233-276`) covers `det J < 0`.
- Weights: `fw = gw(i)*gw(j)` with `(5/9, 8/9, 5/9)`, points `+-0.7745966692414834`, sum 4 (face),
  `vw` sum 8 (volume), `w = vw*ABS(detj)` (`:385`). Exactness on trilinear cells: face integrand
  degree <= 5 per direction, 3-point Gauss exact to 5; the warped-cell `u = x` test exercises it.
- Geometry map: all 27 nodes in both the face and the volume Jacobian (`:317-324`, `:374-381`),
  the ILINT = 2 map of B.
- Reconstruction: `EL_VA_RECONSTRUCT_EPS_LOCAL` (`:444-465`) dvol-weighted vertex sums, Q2-sized
  arrays with zero non-vertex entries; `E013Sum2` under the pinned level
  (`el_transfer.f90:367-370 saved_mg_ilev = el_mg_ilev; el_mg_ilev = ilev; CALL E013Sum2(...)`);
  division `:371-372`; `EL_VA_EMBED_Q1_IN_Q2` (`:469-493`) edge/face/centre means. Bounded and
  linear-exact (test 6).

**(d) MPI safety.**
- Transfer pass: `E013Sum2` -> `E013SumSUPER` (mode 3) is worker-only point-to-point over
  `MGE013(ILEV)` (`QuadSc_mpi.f90`: `IF (myid.ne.MASTER)`, partners from `CommOrder`), so every
  worker must call it: the block is deck-keyed (`el_transfer.f90:363`), the master left at
  `:139`, and there is no rank-local return between `EL_ENSURE_FIELDS` (`:154`) and `:363`.
  Workers without particles take it (owned_particles(0), `:183`).
- Flux block (`el_volume_averaged.f90:379-383`) and eps-weighted sink (`:430-431`): keyed on
  `conservative .AND. audit`, placed after the campaign/temporal/weight reductions, before the
  audit section, which contains no collective. Header `:83-90` describes the order correctly.
- Corrector guard: uniform (deck + a compile-time constant), before any collective in that
  routine; see finding 1 for its position relative to `E013Sum2`.
- `EL_VA_FLUX_MOMENTS` size `ERROR STOP` (`el_va_flux.f90:418`) and the layout `ERROR STOP`
  (`:151`) are terminations, not returns; see finding 3. No new allocation-dependent return was
  added to `EL_ADD_CONTINUITY_SOURCE` (the pre-existing ones at `:297-299` are unchanged).
- Step order: pre-advance pass (`QuadSc_main.f90:374-380`, `-itns`), advance, post-advance pass
  (`:400-406`), fluid core (`:435`), so the transport porosity read by the predictor is
  `eps^{n+1}`; on restart the first mesh pass precedes the first fluid core the same way.

**(e) Sustainability.**
- `EL_VA_FLUX` is MPI-free, uses only `TYPES`, no COMMON, kernels PURE with the tables passed in
  (`:298`, `:348`); the driver (`:397`) is the only non-pure routine. Module state: `tables`
  (built once) and `layout_checked` (finding 3).
- Form dispatch is one `IF` outside the element loop (`:321`) and one at the `defp` block
  (`:396`), legacy code untouched; ids are table-driven (`FORM_*` <-> `el_va_form_names`,
  `EL_VA_FORM()` mirrors `EL_VA_SINK_MODE()`). Validation (`el_config.f90:317-343`) rejects
  unknown values, forms without `continuity`, `conservative` without a corrector, and the meso
  filter. `el_va_flux.f90` is in both source lists (`ProjectFiles.cmake:305`, `:335`), the test
  target and ctest entry exist (`:420-428`). Naming follows the consultation
  (`epsilon_transport_q2`, `flux_moments`, `EL_VA_FLUX`). No dead code found; `NFACE/NGPF/NGPV`
  are exported but unused outside the module (harmless).
- Interplay: `d_cont` comes from the unchanged `EL_VA_CONTINUITY_RATE` (BE/BDF2), the sink from
  the unchanged `EL_VA_SINK`/weights, so time scheme and sink modes compose with the form as
  intended; `S` and `sink_eps` expose clipping/relaxation incompatibility rather than hiding it.
- Sample deck differs from `q2p1_param_VB_bdf2_alpha.dat` by the single line
  `SimPar@ELVAContinuityForm = conservative`; `MatrixRenewal ...C1` with `Pres@MGCrsSolverType = 4`
  passes the corrector-loop restriction.

**(f) Known limits (confirmed as documented).** Closed impermeable boxes only (KNPRP rows,
periodic z-moment seam: design note lines 136-138); no tolerance stopping; dvol single precision
with `dvol_gap` logged (`:476`); `eps_e = MAX(epsilon_f, 1e-3)` is used consistently as the
preconditioning coefficient (`:526`, `:540`, `:405-407`, `:410`), separate from `el_eps_f_min`;
meso filter rejected (`el_config.f90:338-341`); restart: reconstruction rebuilt before first use
by the step order above, never stored.

**(g) Performance.** Kernel cost per element per call is roughly 6x9x27 face plus 27x27 volume
basis contributions (~35 kflop); three calls per step at audit frequency 1 should be well under
a second per call on the VB_L2 partition, so the 0:37 -> 1:45 smoke difference is unlikely to be
the kernel. Check first (consultation §8: time the blocks separately): the "Pressure-Poisson
corrector" iteration counts in the two protocols (the linear rows now carry nonzero targets),
and whether the baseline was the same corrector configuration. Cheap wins that leave the
conservative results bitwise unchanged: (i) in `EL_VA_ELEMENT_FACE_FLUX` only the nine nodes on
the face have nonzero values and tangential derivatives at the face points (the 1D factor is
`1 - 1 = 0` or `0.5*(-1)*0 = -0` exactly), so restricting the l-loop to them for `xq, uq, epsq`
and the two tangent columns adds the same nonzero terms in the same order; the sign of `det J`
is then taken once per element (volume loop) instead of per face point -- about 4x on the face
part; (ii) finding 5; (iii) `t%fw` is rewritten six times (trivial). None is needed for the
ladder; any change re-gates the conservative logs only.

## Verdict

Ready for the gate ladder. The default path is unchanged at source level (legacy loop and
legacy `defp` block mechanically identical to the parent, reductions and early returns
unchanged, all new work behind `conservative`), the conservative algebra and sign chain are
implemented as designed and consistent with the B matrices' own centre and geometry map, the
kernel's tables, signed normal and weights are correct, the reconstruction is bounded and summed
under the pinned level, and every new collective is deck/step-keyed with no rank-local return
in front of it. The should-fix items (guard position, one allocation check, the layout check's
placement/state, stale comments) are hygiene that can follow in a small commit after the gates
run; the duplication in finding 4 is a post-ladder refactor with its own bitwise gate.

## Not verified

- No build, no test run, no inspection of the smoke log (job 157076); the reported numbers were
  taken from the brief.
- Binary-level identity of the legacy loop under the new ELSE wrapping (FMA contraction and
  scheduling at `-O2 -march=znver1`): the cluster gates are the evidence.
- `E013Sum2Rec` (mode 4), periodic identifications in `MGE013`, and workers with zero elements.
- Solvability/convergence of the four-row RHS beyond the smoke residual drop (consultation §4).
- The effective `GAMMA` of the sample deck (no `GAMMA` line in it; default not traced).
- The restart path end to end (reasoned from the step order only).
- The FEAT claim that left-handed hexes are accepted (the `sgn(det J)` handling is correct
  either way).
- `tools/el_stream_compare.py` parsing of the two new tags against a real log.
- The test helper's `karea` convention against FEAT's; production is guarded by the runtime
  layout check, which passed in the smoke.
