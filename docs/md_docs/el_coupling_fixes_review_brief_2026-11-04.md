# Review brief: two Euler-Lagrange coupling fixes (FeatFloWer, Fortran 90 + MPI)

Repository (git worktree, branch fix/el-coupling-p4-p5, forked from feature/el-volume-averaged at e0164ff5):
  /data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va-coupling
Commit under review: 5404e693 (parent e0164ff5). See it with
  git -C /data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va-coupling show 5404e693
Files: source/src_el/el_transfer.f90 (the fixes, in EL_PARTICLE_MESH_PASS), source/src_el/el_config.f90 (comments only).
Read with it: source/src_el/el_forces.f90 (EL_COMPUTE_PARTICLE_FORCES builds feedback_force = drag + lift [+ pressure
under el_pressure_feedback]; EL_SEMI_IMPLICIT_DRAG_IMPULSE), source/src_el/el_fields.f90 (EL_CAPTURE_FLUID_FEEDBACK_SOURCE,
EL_APPLY_FLUID_FEEDBACK_SOURCE), source/src_el/el_quadrature.f90 (EL_DEPOSIT_PARTICLE: deposits -feedback_used into
force_rhs and drag_b into drag_B_ml), source/src_quadLS/QuadSc_def.f90 around lines 2417 and 2516 (thstep*drag_B_source
added to the velocity matrix diagonal when el_apply_fluid_feedback .AND. el_drag_semi_implicit), and the call sites of
the transfer pass in source/src_quadLS/QuadSc_main.f90 (Transport_q2p1_UxyzP_el: pre-advance pass with istep = -itns,
particle advance, post-advance pass with +itns, then the fluid core).

## The coupling, in short
Each time step runs the particle-mesh transfer pass twice. The pre-advance pass (advance_history = istep < 0) evaluates
the hydrodynamic forces, arms PE (the rigid-body engine) when ELApplyParticleForces = Yes, deposits the fluid feedback
force into force_rhs (EL_DEPOSIT_PARTICLE subtracts feedback_used) and the drag coefficient B into drag_B_ml, and
captures both for the fluid solve (EL_CAPTURE_FLUID_FEEDBACK_SOURCE: fluid_feedback_source = force_rhs, drag_B_source =
drag_B_ml). The post-advance pass (istep >= 0) repeats the evaluation at the new particle positions for diagnostics only
(EL_FEEDBACK_CONSERVATION, EL_FORCE_BUDGET, EL_MEAN_SLIP); nothing it deposits reaches the fluid. The fluid momentum
equation then gets the explicit source -feedback_used (EL_APPLY_FLUID_FEEDBACK_SOURCE) and, with ELDragCoupling =
semi_implicit and ELApplyFluidFeedback = Yes, the implicit term +B u^{n+1} on the matrix diagonal. The intended
semi-implicit split of the drag reaction -B(u^{n+1} - u_p) is: matrix +B u^{n+1}, explicit side +B u_p; since the
explicit feedback_force is the old-velocity drag B(u_f - u_p) (+ lift ...), the explicit deposit must be
feedback_used = feedback_force - B u_f (then -feedback_used = -B(u_f - u_p) + B u_f = +B u_p).

## The two defects (reported by a peer session reviewing case spec 05 E1, a frozen periodic packed bed) and the fixes
P5 (ELPressureFeedback not reaching the fluid): with ELApplyParticleForces = Yes the pre-advance pass rebuilds
feedback_used around the sub-cycled implicit drag impulse PE will apply (issue-D fix): feedback_used = drag_force_eff +
lift (- B u_f under semi_implicit). The pressure reaction that EL_COMPUTE_PARTICLE_FORCES adds to feedback_force under
el_pressure_feedback was dropped by that rebuild, while the post-advance diagnostic pass (feedback_used =
feedback_force) kept it, so EL_FEEDBACK_CONSERVATION reported a pressure contribution the fluid never received.
Fix: f_explicit = lift [+ pressure under el_pressure_feedback]; feedback_used = drag_force_eff + f_explicit [- B u_f].
Control (tier-2 periodic box, one particle, two-way, pressure force on, feedback No vs Yes, 60 steps, unfixed binary):
EL_MEAN_SLIP and EL_MOMENTUM_ELEMINT bitwise identical between No and Yes; only EL_FEEDBACK_CONSERVATION differed
(rhs_z 9.0587e-6 vs 9.1420e-6). The after-fix pair is running.

P4 (frozen particles feel the semi-implicit drag twice): with ELApplyParticleForces = No, ELDragCoupling = semi_implicit,
ELApplyFluidFeedback = Yes the matrix still gets +B u^{n+1} (drag_B is deposited and captured unconditionally) but the
pre-advance deposit was the full explicit feedback_force = B(u_f - u_p) + lift (the compensating - B u_f existed only
inside the ApplyParticleForces branch), so the fluid saw -B(2 u_f - u_p). Fix: in the pre-advance pass, when
el_drag_semi_implicit .AND. advance_history .AND. .NOT.el_apply_particle_forces, feedback_used = feedback_force - B u_f.
Control (same box, frozen particle with its initial slip velocity, explicit vs semi_implicit, unfixed binary):
uf_super_z at step 60 3.779e-5 (explicit) vs 3.420e-5 (semi-implicit). The after-fix pair is running.

New diagnostic: EL_FLUID_SOURCE, written in the pre-advance pass at the audit cadence (rhs = sum of the deposited
force_rhs, feedback = sum of feedback_used): the quantity the fluid actually receives, so a feedback gate cannot be
fooled again by the post-advance line.

Bitwise requirements (gates running: tier-2 twin one-way... no: the tier-2 momentum_conservation twin is two-way explicit
with particles advanced; the s04 VB_L2 case is one-way without feedback; the s01 settling twin two-way explicit; the
tier-2 semi-implicit two-way case with particles advanced): in all of these the deposited values must be unchanged
operation for operation, because (a) f_explicit = lift when el_pressure_feedback is .FALSE., so drag_force_eff +
f_explicit is the same expression as before, and (b) the new P4 branch is not entered when particles are advanced or
drag is explicit. Verify this reasoning; point at any path where a default-key run could change in floating point or
in its log lines (the EL_FLUID_SOURCE line is new and only written when el_write_diagnostics and audit cadence).

## Please check
1. Correctness of both fixes against the split described above; sign conventions (EL_DEPOSIT_PARTICLE subtracts
   feedback_used; QuadSc_def adds +thstep*B to the diagonal; EL_APPLY_FLUID_FEEDBACK_SOURCE adds the captured source
   with what sign and scaling?). Read those routines, do not take the brief's word.
2. Whether the P4 term belongs in the pre-advance pass only (the post-advance deposit is diagnostic; the matrix term
   uses the pre-advance drag_B capture) and whether the post-advance EL_FEEDBACK_CONSERVATION line is now misleading in
   the frozen semi-implicit case (it still reports the explicit deposit) - propose how the diagnostics should read.
3. Whether the semi-implicit impulse that PE applies (drag_force_eff from EL_SEMI_IMPLICIT_DRAG_IMPULSE) is consistent
   with adding the pressure reaction to f_explicit rather than to f_other (f_other = pressure + lift + grav_buoy shapes
   PE's implicit velocity; the fluid gets lift and pressure explicitly). Is the Newton pair (particle gains, fluid loses)
   exact for the pressure term now?
4. Any remaining path where ELPressureFeedback = Yes is silently ineffective (one-way mode, explicit vs semi-implicit,
   the final/diagnostic passes), and any other key combination with the same two-pass inconsistency.
5. Code quality for a sustainable codebase: naming (f_explicit), comment accuracy, whether the three feedback_used
   assignments should be one function, dead branches, and whether el_config.f90 should refuse or warn on any
   combination.
6. MPI and determinism: nothing collective changed; confirm.

## Output
Findings ordered by severity (blocker / should fix / nit), each with file:line, what is wrong, why it matters, a
concrete fix. Then a verdict: are P4 and P5 fixed; are default-key runs provably unchanged; is the code clean. State
what you verified by reading and what you could not. Do not modify any files. Read-only review.
