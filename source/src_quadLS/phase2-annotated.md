# Euler-Lagrange Phase 2: One-Way Particle Physics — Annotated Plan

> **About this document.** This is the Phase 2 plan with review annotations folded
> in. Annotations are marked `> **Review note:**` and are grounded in the current
> code (`source/src_el/*`, `source/src_particles/dem_query.f90`, and the PE
> collision systems under `libs/pe/pe/core/collisionsystem/`). Where an annotation
> *changes* the original intent, the original bullet is struck through and the
> corrected version follows.

> **Architecture note (PE-side changes — read first).** All PE-side Phase 2 work
> goes in the **existing** E-L collision-system specialization, not in
> `HardContactAndFluid.h`. The pattern is already in place and already active:
> - `pe/core/collisionsystem/HardContactEulerLagrange.h` — full
>   `CollisionSystem< C<CD,FD,BG,response::HardContactEulerLagrange> >`
>   specialization (~4400 lines), with its own `initializeVelocityCorrections`
>   (line 4296) and velocity-prediction loop (lines 1916-1926).
> - `pe/core/response/HardContactEulerLagrange.h` — the matching response solver.
> - Selected at compile time by `pe/config/Collisions.h:92`
>   (`#define pe_CONSTRAINT_SOLVER pe::response::HardContactEulerLagrange`);
>   `HardContactAndFluid` is commented out at line 97.
>
> So the "less invasive, isolated E-L extension" you want is the established
> design: extend `HardContactEulerLagrange.h` only; `HardContactAndFluid.h` stays
> untouched. Any earlier annotation in this doc that referenced
> `HardContactAndFluid.h` line numbers should be read against
> `HardContactEulerLagrange.h` instead — and note the two make **opposite**
> gravity/buoyancy decisions (see the gravity section).

## Summary

Implement Phase 2 only: replace the provisional Stokes/Schiller-Naumann particle
forcing with the planned one-way E-L force model, while keeping particle feedback
to the fluid diagnostic-only until Phase 3. The new implementation adds Di Felice
drag, pressure-gradient force, configurable lift, and a real semi-implicit
particle drag update through the PE interface.

> **Review note (de-risking).** The transfer layer already samples everything
> Phase 2 needs. `EL_SAMPLE_SIZE = 18` (`el_quadrature.f90:120-127`) and
> `EL_INTEGRATE_PARTICLE` already accumulates, kernel-weighted and
> conservation-reduced:
>
> | slot      | quantity            |
> |-----------|---------------------|
> | `1`       | normalization       |
> | `2:4`     | velocity            |
> | `5:13`    | ∇u (3 rows)         |
> | `14`      | pressure            |
> | `15:17`   | ∇p                  |
> | `18`      | epsilon_f           |
>
> Phase 1 consumes only `sample(1)` + `sample(2:4)` today
> (`el_transfer.f90:104`). So Phase 2 is mostly **additive force math on data
> already in hand** — no quadrature/halo changes. This removes the riskiest class
> of change.

## Key Changes

### Force modeling in `source/src_el`

- Add `EL_COMPUTE_PARTICLE_FORCES(...)` returning `drag`, `pressure`, `lift`,
  `particle_total`, `feedback_force`, `drag_B`, and sampled `u_f`.
  > **Review note.** Keep `gravity` in the returned set (renamed `grav_buoy`): for
  > the active `HardContactEulerLagrange` solver, gravity+buoyancy are the driver's
  > responsibility (see gravity section). Decomposition:
  > - `grav_buoy   = (ρ_p − ρ_f)·V_i·g`            → net submerged weight.
  > - `particle_total = drag + pressure + lift + grav_buoy`  → the force mapped to
  >   PE (the EL solver applies it verbatim via `getForce()`).
  > - `feedback_force = drag + lift`               → the Newton's-3rd-law reaction
  >   spread to the fluid (Phase 3). Pressure, gravity, and buoyancy are **not**
  >   transmitted to the fluid.
- Implement `ELDragModel = difelice | stokes | schiller_naumann`; default
  `difelice` for `q2p1_el_pipeflow`.
  > **Review note.** `el_config.f90` currently whitelists only
  > `{stokes, schiller_naumann}` (`el_config.f90:30-37`). Add `difelice` there and
  > to `EL_PRINT_CONFIG`. The existing `EL_DRAG_FORCE` (`el_forces.f90`) already
  > computes the Stokes/SN path; extend rather than replace it so the legacy models
  > stay available.
- Use the current convention that `Properties%Viscosity(1)` is **kinematic**
  viscosity and compute `mu_f = rho_f * nu_f` inside E-L force routines.
  > **Review note — confirmed.** `EL_DRAG_FORCE` already does
  > `dynamic_viscosity = density*viscosity` (`el_forces.f90:23`); the momentum
  > diffusion coefficient uses `Viscosity(1)` bare (`QuadSc_corrections.f90:52`).
  > Convention holds.
- Decode existing transfer samples as: normalization, velocity, velocity-gradient
  rows, pressure, pressure gradient, and epsilon.
  > **Review note — already implemented.** See the slot table above. Add named
  > parameters/accessors for these offsets so the force code does not hard-code
  > `sample(15:17)` etc.

### Force closures

- Di Felice drag exactly as specified, with `Re_p >= 1e-12`, `epsilon_f` clipped
  by existing `ELEpsFMin`, and `B_drag = |F_drag| / max(|u_rel|, tiny)` for
  semi-implicit use.
  > **Review note.** `B_drag` is collinear-valid: Stokes/Di Felice drag is along
  > the slip, so `B·u_rel` reproduces `F_drag`. Assert this in tests (cheap, guards
  > sign/normalization of the exact quantity PE consumes).
- Pressure-gradient force `F_gradp = -V_i grad(p)` when `ELPressureForce = Yes`.
  > **Review note.** This is the *hydrodynamic* pressure force (resolved dynamic
  > pressure gradient), distinct from buoyancy. Keep it. In a quiescent settling
  > case `grad(p) ≈ 0`, so it contributes nothing there — the settling driver is
  > the driver-supplied `grav_buoy` term (see gravity section). Consistent.
- ~~Particle gravity force `m_i g` from `Properties%Gravity`; do not add separate
  Archimedes `-rho_f V_i g`.~~ → **Driver supplies net submerged weight
  `(ρ_p − ρ_f)·V·g` (gravity + buoyancy) in the mapped force.**
  > **Review note — CORRECTED, and this supersedes an earlier draft.** An earlier
  > version of this note read buoyancy from `HardContactAndFluid.h` and concluded
  > "remove gravity from the Fortran side." **That was the wrong collision
  > system.** The *active* solver is `HardContactEulerLagrange`, which makes the
  > opposite decision and states it explicitly
  > (`HardContactEulerLagrange.h:1901-1903`):
  > ```cpp
  > // Euler-Lagrange body-force semantics such as gravity and buoyancy are
  > // handled by the outer driver, not by this PE collision system.
  > ```
  > Its velocity prediction applies **only** the mapped external force —
  > `dv = invMass*dt*getForce()` (`:4300`) — with **no** buoyancy term
  > (`v_[j] = body->getLinearVel()`, `:1920`). Contrast `HardContactAndFluid.h:2005`,
  > which *does* apply `(ρ_p − ρ_f)·V·g` inside PE; that path is inactive.
  >
  > Consequences for the E-L path:
  > - The Fortran driver **must** include gravity **and** buoyancy in the force it
  >   maps to PE, as the net submerged weight `F_grav+buoy = (ρ_p − ρ_f)·V_i·g`.
  >   (So the original "no Archimedes" instruction is wrong for a submerged
  >   particle — buoyancy `−ρ_f·V_i·g` must be included.)
  > - `g` is the configured gravity; `ρ_f` is the carrier density already passed to
  >   the force routines. This keeps **all** body-force semantics on the driver
  >   side, exactly as the EL solver intends, and needs **no** PE change for
  >   gravity.
  > - `Properties%Gravity` (or the EL config's own gravity vector) is the source of
  >   `g`; the **fluid** momentum equation still carries no gravity, by our
  >   convention.
  >
  > **DECIDED — (A) driver-side gravity+buoyancy.** Compute
  > `grav_buoy = (ρ_p − ρ_f)·V_i·g` in the Fortran force routine and add it to
  > `particle_total`. This matches the active EL solver's stated contract
  > (`HardContactEulerLagrange.h:1901-1903`) and requires **no** PE change for
  > gravity — `HardContactEulerLagrange::initializeVelocityCorrections` applies the
  > mapped force verbatim. The rejected alternative (porting `HardContactAndFluid`'s
  > PE-side buoyancy term into the EL specialization) is **not** pursued: it
  > contradicts the file's design intent and widens the PE surface.
  > So `particle_total = drag + pressure + lift + grav_buoy`, while
  > `feedback_force = drag + lift` is unchanged (gravity/buoyancy are never spread
  > to the fluid).
- Lift controls: `ELLiftModel = none | saffman_mei | saffman_mei_wall`, default
  `none`; implement the isolated coefficient/wall-factor functions with literature
  comments, validation focused on switchability and deterministic values.
- Add `ELMagnus = No` config now: parsed and printed, left off and unimplemented;
  enabling it stops with a clear unsupported message.
  > **Review note.** Put the unsupported-stop in `EL_VALIDATE_CONFIG`, not deep in
  > the force loop, so a misconfigured run fails immediately at startup.

### PE-side semi-implicit E-L update path

> **Review note — scope flag.** This is the largest and highest-risk item; it
> crosses into the `libs/pe` submodule, but stays inside the **already-isolated**
> `HardContactEulerLagrange` specialization — `HardContactAndFluid.h` is not
> touched. `tParticleData` (`dem_query.f90:12-24`) carries only
> `position, velocity, force, torque, density, radius, systemIdx` — there is **no
> slot for `B_drag`, `u_f`, or `F_other`.** Land this as its own gated milestone
> (gate = the analytic semi-implicit test) *before* the force-model work, so the
> integrator change is validated in isolation.

- Extend the Fortran/C particle force bridge with E-L hydrodynamic state fields
  (`B_drag`, `u_f(3)`, `F_other(3)`) sufficient for PE to apply
  `U_{n+1} = (U_n + dt_sub/m·(B·u_f + F_other)) / (1 + dt_sub·B/m)` before PE
  contact/lubrication substeps.
  > **Review note.** Touches: Fortran `tParticleData`, the matching C struct,
  > `setForcesMapped`/`c2f_interface.cpp`, **and** the integrator in the EL
  > specialization. `F_other` here = pressure + lift + grav_buoy (everything
  > except drag); drag enters implicitly via `B`/`u_f`. Wire it into
  > `HardContactEulerLagrange::initializeVelocityCorrections`
  > (`HardContactEulerLagrange.h:4296-4312`), which is currently the fully explicit
  > `dv = invMass*dt*getForce()` (`:4300`). Replace that with the semi-implicit
  > form so `dv` becomes `dt/m·(B·u_f + F_other)/(1 + dt·B/m)` minus the current
  > velocity contribution as appropriate. Because the EL solver puts **all**
  > body forces (incl. grav_buoy) on the driver, there is no separate PE buoyancy
  > term to compose with — `F_other` already contains it. Keep the change confined
  > to this one routine; the surrounding velocity-prediction loop
  > (`:1916-1926`) stays as-is.
- Disable the legacy `applyFluidForces(fullStepSize, alpha=0.35)` smoothing for
  the E-L pipeflow path once semi-implicit mode is active.
- Keep torque from hydrodynamics zero; PE contacts/lubrication remain
  authoritative for torque.
  > **Review note — matches current code.** `el_transfer.f90:127` already sets
  > `owned_particles(i)%torque = 0.0d0`.
- Preserve `ELApplyParticleForces = No` as a diagnostic mode that computes fields
  and forces without modifying PE state.
  > **Review note.** Already wired as `el_apply_particle_forces`
  > (`el_transfer.f90:124-129`). Keep that guard around the new bridge call too.

> **Review note — feedback ≠ total plumbing.** Today the *same* `force` is both
> applied to PE and spread to the fluid via
> `EL_BROADCAST_FROM_OWNERS → EL_DEPOSIT_PARTICLE` (`el_transfer.f90:131-141`).
> In Phase 2 these diverge. Cleanest mapping without growing the `(4,n)` result
> arrays: apply `particle_total` to PE locally in the owner loop, and set
> `owned_result(2:4) = feedback_force` so **only feedback is broadcast/deposited**.
> Make this explicit, or the Phase 1 spread-conservation test silently starts
> asserting the wrong vector.

### Phase 3 out of scope

- `el_field_data%force_rhs` continues to be assembled/exported for
  diagnostics/tests only.
- No insertion of `force_rhs` into the Navier-Stokes RHS.
- No semi-implicit drag term in the fluid matrix.

## Test Plan

### Force unit tests (extend `test_el_kernel_forces` or add `test_el_phase2_forces`)

- Di Felice reduces to Stokes within 1% for `Re_p < 1e-3`, `epsilon_f = 1`.
  > **Review note.** Also assert the limit form directly at `Re_p → 0, ε_f = 1`
  > so the test cannot pass trivially because both sides are tiny.
- Di Felice force increases monotonically as `epsilon_f` decreases at fixed slip.
- Pressure force equals `-V_i gradp` componentwise.
- ~~Gravity force equals particle mass times configured gravity.~~
  > **Review note — REVISED.** Assert the driver `grav_buoy` term equals the **net
  > submerged weight** `(ρ_p − ρ_f)·V_i·g` componentwise (not `m_i·g`), and that
  > `particle_total` includes it while `feedback_force` does not. This is the
  > driver-side quantity the EL solver consumes verbatim.
- Feedback force contains only drag plus lift, never pressure or gravity.
- Lift model `none` gives exactly zero; Saffman-Mei paths return finite
  deterministic values for a fixed velocity gradient.
- **(added)** `B_drag` consistency: `B·u_rel` reproduces `F_drag`
  magnitude/direction at the sampled slip.

### Transfer/force integration tests on the synthetic Q2/P1 mesh (Phase-1 style)

- Constant velocity + constant pressure gradient produces the expected force from
  kernel-sampled fields.
- Linear velocity field produces expected vorticity/lift inputs away from walls.
- Wall-adjacent particles retain normalized `epsilon_f`, pressure-gradient, and
  velocity sampling behaviour.
- MPI serial/2/8-rank variants remain registered like `el-transfer-*`.
- **(added)** Spread-conservation now targets `feedback_force`: assert
  `Σ_a b_a = -Σ_i feedback_i` and that the deposited field excludes
  pressure/gravity.

### PE interface coverage

- Analytic semi-implicit check: one particle, several `dt`, against the
  closed-form velocity. **This is the gate for the PE-interface milestone.**
- Stiff drag stability for `tau_p << dt` with no velocity overshoot.
- `ELApplyParticleForces = No` leaves PE velocity unchanged.
- **(added)** Grav/buoy sanity: with zero hydrodynamic force, mapping
  `force = grav_buoy = (ρ_p − ρ_f)·V·g` produces a PE velocity increment of
  `(ρ_p − ρ_f)/ρ_p · g · dt` through `dv = invMass·dt·getForce()`
  (`HardContactEulerLagrange.h:4300`), confirming the driver-side gravity path is
  correct and not double-counted (PE adds no buoyancy of its own).

### Runtime acceptance

- Lightweight settling-style regression: one sphere in quiescent fluid via the
  E-L pipeflow executable.
  > **Review note.** Terminal velocity balances **net submerged weight**
  > `(ρ_p − ρ_f)·V·g` (driver-supplied) against drag — **not** `m·g`. Check against
  > the force-law value using `(ρ_p − ρ_f)`. Confirm `ρ_f` used in `grav_buoy`
  > matches the carrier density, and that fluid gravity is off.
- Re-run existing E-L transfer tests and one existing DNS/FBM smoke/build target
  to confirm no FBM regression.

## Assumptions (revised)

- This phase implements Phase 2 with focused tests modeled after Phase 1; it does
  not start Phase 3 fluid feedback.
- The semi-implicit implementation modifies the **EL** collision specialization
  (`HardContactEulerLagrange::initializeVelocityCorrections`,
  `HardContactEulerLagrange.h:4296`) rather than approximating it through the
  mapped-force-only path. `HardContactAndFluid.h` is not modified.
- **Gravity and buoyancy are applied by the driver**, not PE. The active
  `HardContactEulerLagrange` solver explicitly delegates body-force semantics to
  the outer driver (`HardContactEulerLagrange.h:1901-1903`) and applies only the
  mapped `getForce()`. So the Fortran E-L force maps
  `drag + pressure + lift + (ρ_p − ρ_f)·V·g` to PE. The fluid momentum equation
  carries no gravity.
- Saffman-Mei formulas coded with citations in comments; full Segré-Silberberg
  validation remains Phase 5.

## Suggested sequencing

1. **Config + force closures** on already-sampled data (Di Felice, pressure, lift;
   `EL_COMPUTE_PARTICLE_FORCES`) with the force/integration unit tests.
2. **Feedback/total split** through deposit (`owned_result(2:4) = feedback_force`),
   update the spread-conservation test.
3. **PE semi-implicit interface** as its own gated milestone (bridge fields +
   integrator change in `HardContactEulerLagrange.h`), gated by the analytic +
   stiff-stability tests.
4. **Settling acceptance** case last, once buoyancy (PE) and drag (coupled) are
   both in place.

## Phase 2 test status (tracking)

Legend: **[done]** implemented and registered · **[prepared]** specified here,
not yet coded.

### Force closures (Fortran)
- **[done]** Di Felice → Stokes limit within 1% (`test_el_kernel_forces`, STOP 40).
- **[done]** Di Felice monotonic in ε_f at fixed slip (STOP 47).
- **[done]** Pressure force `= −V·∇p` componentwise (STOP 41).
- **[done]** `grav_buoy = (ρ_p−ρ_f)·V·g` componentwise (STOP 42).
- **[done]** Lift `none` ⇒ exactly zero (STOP 43).
- **[done]** `particle_total = drag+pressure+lift+grav_buoy` (STOP 44).
- **[done]** `feedback_force = drag+lift` only (STOP 45).
- **[done]** `B_drag·slip = drag` collinearity (STOP 46).
- **[done]** Saffman-Mei / -wall return finite deterministic values (STOP 48/49).
- **[prepared]** Direct Di Felice limit *form* at ε_f=1, moderate Re (verify the
  voidage factor is exactly 1 and the full `Cd·A` path — not just the Stokes
  asymptote — is exercised). Add alongside STOP 40.

### PE semi-implicit integrator (C++)
- **[done]** Defining implicit relation solved exactly; closed-form match over a
  range of `dt` (`pe-el-semi-implicit-drag`).
- **[done]** Stiff stability `τ_p ≪ dt`: no overshoot, monotone contraction,
  `B→∞ ⇒ u_f`.
- **[done]** Pure-forcing limit `B=0`.
- **[done]** Convergence to the exact exponential as `dt→0`.
  > These cover the math via the extracted pure helper
  > `pe/core/collisionsystem/ELSemiImplicitDrag.h`, which the integrator now
  > calls — so the unit test exercises the exact code path without standing up
  > MPI/world/collision machinery.
- **[prepared]** `ELApplyParticleForces=No` ⇒ PE velocity unchanged. Largely
  guaranteed by construction now: the driver arms PE only in the pre-advance
  pass, and `clearELHydroStates()` runs at the end of every PE step, so an
  un-armed body takes the explicit branch with `getForce()=0`. A full-app
  integration test (one fluid step, assert PE velocities unchanged) would
  confirm the wiring end-to-end; needs the staged pipeflow harness.

### Transfer/force integration on the synthetic Q2/P1 mesh (Fortran)
- **[prepared]** Constant velocity + constant ∇p ⇒ expected `drag` and
  `pressure = −V·∇p` from kernel-sampled fields (extend `test_el_transfer`:
  seed `val_p` with a linear field, call `EL_INTEGRATE_PARTICLE` then
  `EL_COMPUTE_PARTICLE_FORCES`, compare).
- **[prepared]** Linear velocity field ⇒ expected vorticity/lift inputs away from
  walls (Saffman-Mei branch), wall-adjacent normalization preserved.
- **[prepared]** Spread-conservation retargeted to `feedback_force`:
  `Σ_a b_a = −Σ_i feedback_i`, and assert the deposited field excludes
  pressure/gravity. (Phase-1 conservation test still passes as-is; this adds the
  feedback-specific assertion.)
- **[prepared]** MPI serial/2/8-rank variants, registered like `el-transfer-*`.

### Runtime acceptance
- **[prepared]** One-sphere settling in quiescent fluid via `q2p1_el_pipeflow`:
  set `Properties%Gravity` (currently 0 in the `.dat`), check terminal velocity
  against the force-law value using **net submerged weight** `(ρ_p−ρ_f)·V·g`.
- **[prepared]** Regression guard: re-run `el-*` transfer tests + the
  `pe-el-semi-implicit-drag` test + one existing DNS/FBM smoke target, and a
  `q2p1_el_frozen_trace` trajectory check to confirm the restored `α=0.35`
  relaxation keeps it within tolerance.
