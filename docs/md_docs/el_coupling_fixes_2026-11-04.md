# EL coupling fixes P4 and P5 (branch fix/el-coupling-p4-p5, 2026-11-04)

Reported by the LIGGGHTS-PUBLIC session while preparing case spec 05 (ErgunTestMPI twin, stage E1 =
frozen periodic packed bed, `ELApplyParticleForces No`, `ELApplyFluidFeedback Yes`,
`ELDragCoupling semi_implicit`). Both confirmed by reading at e0164ff5 and by experiment; fixed in
5404e693, review fixes in df55a275 (`el_coupling_fixes_review_2026-11-04.md`, Codex gpt-6-astra high).

## The defects

**P5, pressure feedback never reached the fluid with particles advanced.** The pre-advance transfer pass
(the only one whose deposit the fluid receives, `EL_CAPTURE_FLUID_FEEDBACK_SOURCE`) rebuilds the feedback
around the sub-cycled implicit drag impulse PE applies (issue-D fix) and dropped the `ELPressureFeedback`
term; the post-advance diagnostic pass kept it, so `EL_FEEDBACK_CONSERVATION` reported a pressure share the
fluid never saw. The 2026-10-15 gate of the key (job 149907) read that line (+0.92 %) and was blind to it.
No recorded run used the key (the staged `vapf_*` / `contpf` reruns, never submitted, would have).

**P4, frozen particles felt the semi-implicit drag twice.** With `ELApplyParticleForces No` the fluid
matrix still carries the theta-weighted `+B u^{n+1}` (`QuadSc_def`, `drag_B_source`) but the explicit
deposit was the full old-velocity drag `B(u_f - u_p)`; the compensating `-B u_f` lived only inside the
particle-advance branch. The fluid saw `-B(2 u_f - u_p)`. No deck in the tree uses the combination; the
E1 deck does.

## The fixes (`source/src_el/el_transfer.f90`, pre-advance pass only)

- P5: `feedback_nondrag = lift [+ pressure under ELPressureFeedback]`; the rebuilt deposit is
  `drag_force_eff + feedback_nondrag [- B u_f]`. The pressure stays in `f_other` on the PE side (it shapes
  the implicit particle velocity) and enters the fluid once.
- P4: with semi-implicit drag and no particle advance the deposit is `feedback_force - B u_f`, completing
  the split `-B(u^{n+1} - u_p)` around the sampled `u_f`.
- New audit line `EL_FLUID_SOURCE phase= pre step= n mode= ... applied= yes|no rhs_explicit= ...
  feedback_explicit= ...`: the explicit source the fluid momentum equation receives (`applied` = fluid
  feedback on and no prescribed field). `EL_FEEDBACK_CONSERVATION` (post-advance) is a deposition
  residual of the diagnostic pass; its format is unchanged.
- Config: `ELPressureFeedback = Yes` is refused without `ELApplyFluidFeedback = Yes` and
  `ELPressureForce = Yes`; `ELMomentumFix = Yes` is refused with semi-implicit drag (the `EL_FLUID_PAIR`
  audit counts the captured explicit source only and would compensate the implicit drag operator as a
  mismatch; extending that audit is a separate item). `ELApplyParticleForces = No` does not freeze
  particles: PE advances them without the hydrodynamic forces (a particle at rest stays at rest, a moving
  one keeps its velocity); a bed is held by its PE configuration.

## Controls (tier-2 periodic box QBOX9, one particle with initial slip, 60 steps of 5e-4, 28 ranks,
`runs/coupling_p4p5/{before_e0164ff5,after_fix,after_fix2}`)

| Control | before (e0164ff5) | after (5404e693) | after (df55a275) |
|---|---|---|---|
| P5: pressure force on, feedback No vs Yes, particles advanced (jobs 156709/156710; 156715/156716; 156729/156730) | EL_MEAN_SLIP and EL_MOMENTUM_ELEMINT bitwise identical between No and Yes; only EL_FEEDBACK_CONSERVATION differs (rhs_z 9.0587e-6 vs 9.1420e-6) | fluid streams differ: uf_super_z at step 60 8.1455e-6 (No) vs 8.1997e-6 (Yes), EL_FLUID_SOURCE rhs_z 9.3158e-6 vs 9.4424e-6 | bitwise = 5404e693 |
| P4: particle forces off, feedback on, explicit vs semi_implicit (jobs 156711/156712; 156717/156718; 156731/156732) | uf_super_z at step 60 4.9635e-5 (explicit) vs 4.4400e-5 (semi, -10.5 %) | 4.9635e-5 vs 4.9651e-5 (+3e-4, the implicit correction) | bitwise = 5404e693 |
| key off / explicit paths vs before | - | P5 No-run and P4 explicit run bitwise vs e0164ff5 | bitwise |

## Gates at the defaults (every certified configuration; `tools/el_stream_compare.py`)

| Case | configuration | 5404e693 | df55a275 |
|---|---|---|---|
| tier-2 momentum_conservation twin vs `reference/run_g0_twin_va_2026-10-14.log` | two-way, explicit, advanced | 156720 bitwise | 156734 bitwise |
| tier-2 momentum_conservation_semi vs the e0164ff5 run (156713) | two-way, semi-implicit, advanced | 156719 bitwise | 156733 bitwise |
| s04 VB_L2 defaults vs the T4 log | one-way, no feedback | 156721 bitwise | 156735 bitwise |
| s04 T9a vs the T9a log | one-way, BDF2 + correctors + alpha | 156722 bitwise | 156736 bitwise |
| s01 G2_twoway_lev2 vs its log | two-way, explicit, advanced, 5 ranks | 156723 bitwise | 156737 bitwise |

The EL_FLUID_SOURCE stream is new (absent from every reference) and is compared by the tool from now on.

## Follow-ups
- Extend `EL_FLUID_PAIR` (el_diagnostics.f90) by the implicit drag operator, then lift the
  ELMomentumFix guard.
- Re-run the staged pressure-feedback rungs (`vapf_*`, `contpf`) with this binary: they are the first
  runs in which the key acts.
- Case spec 05 E1 can run on this branch; the peer's G-E1-0 smoke should also check the EL_FLUID_SOURCE
  `applied= yes` line.
