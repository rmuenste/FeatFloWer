# Euler–Lagrange Phase 2: outstanding test coverage

Reminder list of Phase 2 test items that were **specified but never coded**. Extracted
verbatim-in-substance from the Phase 2 annotated plan
(`source/src_quadLS/phase2-annotated.md`) before that planning document was removed —
the plan itself is superseded by the implementation and by
[`el_methods_inventory.md`](el_methods_inventory.md), but this checklist was still live.

Status at extraction: **14 items done, 9 outstanding.** The completed items are
registered tests and are not repeated here; see `test_el_kernel_forces` (STOP codes
40–49), `pe-el-semi-implicit-drag`, and the `el-transfer-*` CTest family.

Nothing below is a known defect. These are gaps in *verification*, not observed failures.

## Force closures (Fortran)

- [ ] **Di Felice limit form at ε_f = 1, moderate Re.** Verify the voidage factor is
      exactly 1 and that the full `Cd·A` path is exercised — not just the Stokes
      asymptote, which is what STOP 40 currently covers. Add alongside STOP 40 in
      `test_el_kernel_forces`.

## PE semi-implicit integrator (C++)

- [ ] **`ELApplyParticleForces=No` ⇒ PE velocity unchanged.** Largely guaranteed by
      construction: the driver arms PE only in the pre-advance pass, and
      `clearELHydroStates()` runs at the end of every PE step, so an un-armed body takes
      the explicit branch with `getForce() = 0`. A full-application integration test
      (one fluid step, assert PE velocities unchanged) would still be worth having as a
      guard against that construction being broken later.

## Transfer / force integration on the synthetic Q2/P1 mesh

- [ ] **Constant velocity + constant ∇p ⇒ expected `drag` and `pressure = −V·∇p`** from
      kernel-sampled fields. Extend `test_el_transfer`: seed `val_p` with a linear field,
      call `EL_INTEGRATE_PARTICLE` then `EL_COMPUTE_PARTICLE_FORCES`, compare.
- [ ] **Linear velocity field ⇒ expected vorticity/lift inputs away from walls**
      (Saffman–Mei branch), with wall-adjacent normalization preserved.
- [ ] **Spread-conservation retargeted to `feedback_force`:** assert
      `Σ_a b_a = −Σ_i feedback_i`, and that the deposited field excludes pressure and
      gravity. The Phase-1 conservation test still passes as-is; this adds the
      feedback-specific assertion.
- [ ] **MPI serial / 2-rank / 8-rank variants**, registered like `el-transfer-*`.

## Runtime acceptance

- [ ] **One-sphere settling in quiescent fluid** via `q2p1_el_pipeflow`: set
      `Properties%Gravity` (currently 0 in the `.dat`) and check terminal velocity
      against the force-law value using the **net submerged weight**
      `(ρ_p − ρ_f)·V·g`. See the gravity/buoyancy ownership note in
      `libs/pe/doc/technical-notes/euler-lagrange-solver-conventions.md` — under the
      active E-L solver the driver supplies both gravity and buoyancy, so a test that
      omits the Archimedes term will disagree with the force law.
- [ ] **Regression guard:** re-run the `el-*` transfer tests, the
      `pe-el-semi-implicit-drag` test, one existing DNS/FBM smoke target, and a
      `q2p1_el_frozen_trace` trajectory check, to confirm the restored `α = 0.35`
      relaxation keeps it within tolerance.

## Related

- [`el_methods_inventory.md`](el_methods_inventory.md) — the implemented closures and their literature origins.
- [`el_validation_report.md`](el_validation_report.md) — the validation campaign and its verdicts.
- [`euler_lagrange_development.md`](euler_lagrange_development.md) — development notes and staged roadmap.
