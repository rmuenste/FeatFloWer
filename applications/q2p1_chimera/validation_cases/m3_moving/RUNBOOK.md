# M3 — Moving bodies with the Chimera component (Phase 6)

Design: `chimera-integration-design.md` §10 row 6 ("Moving"); Phase-6 notes
in §3. Builds on the Phase-5 periodic Hasimoto case
(`validation_cases/m2_hasimoto/RUNBOOK.md`). First executed 2026-09-04 on
the single-core sandbox (9 oversubscribed ranks).

## What moves and how

- `SimPar@ChimeraMotion = prescribed | free`. Bodies carry X_k, U_k, Ω_k
  (body table columns after H: `ux uy uz wx wy wz rho_s`).
- **No submesh ALE.** Each submesh is solved in the translating frame of
  its body: the mesh stays at its initial fit, the frame change
  `u' = u − U_k` is exact and adds only the fictitious force
  `−ρ dU_k/dt` (a uniform body force, the Phase-5 path); rotation enters
  through the inner Dirichlet data `Ω × r`. Consequences: locator and
  Stokes-frozen factorisation are reused, the Robin sample points are
  shifted by `X_k − X0_k` and wrapped, the fringe/penalty data go back
  to the lab frame (`u' + U_k`; rigid velocity in the hole).
- Per coupling update: `X^{n+1} = X^n + Δt U^n`, `a^n = (U^n − U^{n−1})/Δt`,
  markers re-classified (strong) or the lumped penalty re-tabulated on
  all levels (weak), donors and lab sample points rebuilt; after the
  solve free spheres integrate
  `U^{n+1} = U^n + Δt (F + F_ext)/(m_s + c ρ_f V + Δt κ 6πμR)` (virtual-mass
  factor `c = ChimeraAddedMass`, default 0.5; implicit Stokes-drag
  linearisation `κ = ChimeraDragImplicit`, default 1) and
  `Ω^{n+1} = Ω^n + Δt T/(I + Δt κ 8πμR³)`. Explicit staggered coupling,
  first order in time (H13 in-step outer iterations and the PE hand-off
  H14 are open).
- Lines: `ChimeraBody<k>: time X Y Z Ux Uy Uz Wx Wy Wz` every step.

## Cases (all: 6³ periodic box, L2 = D/h 4, level-2 atmosphere, μ = ρ = 1,
dt 0.01 to t = 4.01, `ChimeraSubStokes`)

1. **Galilean check** (`_data/q2p1_param_galilean*.dat`,
   `_data/hasimoto_moving.dat`): the Hasimoto sphere translates at the
   prescribed velocity `−U_sup` of the static case with the body force
   scaled to f = 1 (`q2p1_param_galilean_static.dat` is that reference:
   F_s 0.97260993, U_sup 0.17523049, Re ≈ 0.06). Galilean invariance
   demands the same steady force and a vanishing lab-frame cell average;
   the sphere crosses 8.4 background cells in the run, so the force
   history measures the **continuity across grid crossings**.
2. **Free sedimentation** (`_data/q2p1_param_sediment*.dat`,
   `_data/hasimoto_free.dat`): a sphere of density 2 released at rest
   under `ChimeraBodyGravity = (0,0,−0.515664)` (so that
   `(ρ_s − ρ_f) V g = 1e-2`) with the fluid counter-force
   `ConstantForcing = +1e-2/V_fluid` (total momentum conserved: the
   sedimenting-suspension setup of the FF-EL closures). Steady state: the
   sphere settles at U_s, the fluid drifts up, and the relative
   superficial velocity `U_rel = ⟨u⟩_cell − U_s` must satisfy
   `F_ext = 6πμR K U_rel` with the K of the static case (1.80–1.82 at
   this resolution).

Post-processing: `python3 validation_cases/m3_moving/moving_stats.py
prot.txt --tmin 2 --h 0.083333 --R 0.1666665 --fref <F_ref>`.

## Staging

As in the M2 runbook (periodic box, axis partition, 9 ranks); the
Phase-6 decks and body tables are vendored in `_data/`. The static anchor
`q2p1_param_hasimoto.dat` must stay bit-identical (checked: `P6_H_S_L2`).

## Results (2026-09-04)

| Case | Variant | ⟨F_z⟩ (t ≥ 2) | p2p/⟨F⟩ | ⟨u_z⟩_cell (lab) | U_rel | K(F_ref) | wall |
|---|---|---|---|---|---|---|---|
| G0_S static f = 1 | S | 0.97260993 | — | 0.17523049 | 0.17523 | 1.8020 meas / 1.8165 bal | 4.4 min |
| G_S prescribed −U_sup | S | 0.98117 | 1.03 % | 9.42e-4 (0.5 % of U_sup) | 0.17617 | 1.757 (F_ref = static F_s) | 5.8 min |
| G_W prescribed −U_sup | W (γ 1e5) | 0.98078 | 0.63 % | 1.13e-3 | 0.17636 | 1.756 | 4.6 min |
| SED_S free (c_am 0.5, plain rotation) | S | — NaN at t = 0.35: Ω grows ×1.7 per step (explicit update vs stiff viscous torque, I/(8πμR³) = 0.0037 < Δt/2); translation only oscillated (decaying) | | | | | |
| SED4_S free, `ChimeraDragImplicit 1` | S | 9.99363e-3 | 0.015 % | −2.08e-4 | 1.882e-3 | 1.691 (F_ref = F_ext) | 4.7 min |
| SED4_W free, `ChimeraDragImplicit 1` | W | 9.99505e-3 | 0.018 % | −1.59e-4 | 1.866e-3 | 1.706 | 4.7 min |

Momentum-balance reference for the prescribed case: f V_fluid = 0.98061.
For the sedimentation cases the Hasimoto-convention force is
F_H = F_s + f_c V_solid = 1.0191e-2 (mean-gradient buoyancy of the fluid
counter-force f_c = 1e-2/V_fluid), giving K_meas = 1.724 (S) / 1.740 (W).

## Reading

- **Galilean invariance holds.** The lab-frame cell average of the
  translating case is 0.5 % of U_sup (exactly 0 for a Galilean-invariant
  discretisation); the mean force over 8.4 cell crossings is within
  +0.06 % (S) / +0.02 % (W) of the momentum-balance value f V_fluid, i.e.
  the moving sphere transmits the momentum that the static coupling
  leaks (static F_s is 0.8 % low). The relative superficial velocity is
  0.5–0.6 % above the static one, so K is 2.5 % below the static K —
  the moving-frame discretisation error at D/h = 4.
- **Force continuity across grid crossings** (the M3 metric): 1.0 %
  peak-to-peak for the strong variant (hole/fringe sets change as the
  sphere crosses cells), 0.63 % for the weak variant — Chimera-W is the
  smoother coupling, as the paper claims, but the strong variant is far
  from the FBM staircase jumps.
- **Free motion needs the implicit drag/torque linearisation.** With
  `ChimeraDragImplicit = 1` (Stokes slope, exact implicit Euler in
  creeping flow) both variants are stable at Δt = 0.01 although the
  viscous relaxation times of the sphere (0.007 translational, 0.0037
  rotational) are below Δt; the virtual-mass factor 0.5 was sufficient
  for the translational coupling.
- **Sedimentation reaches force balance to 0.06 %** but the sphere and
  the fluid drift down together (U_s from −1.81e-3 at t = 0.5 to
  −2.09e-3 at t = 4, ⟨u_f⟩ ≈ −2e-4): the periodic cell conserves total
  momentum only if the coupling does, and the static 0.8 % momentum leak
  of the overlap (surface traction below f V_fluid) is integrated in
  time — the leak rate 8e-5 per unit time matches the observed drift.
  The relative velocity therefore reads 7 % high (K 1.69–1.74 vs the
  static 1.80–1.82). Conclusion: the motion machinery (frame change,
  reclassification, integrator) is validated by the Galilean case; the
  dynamics accuracy item is the conservative (variationally consistent)
  force evaluation, not Phase 6.

## Open

- H13 (in-step outer iterations, time-accurate Chimera-S) and H14 (PE
  hand-off); the internal Newton–Euler integrator is the Phase-6
  dynamics.
- ten Cate sedimentation in a walled box (needs an all-`Wall` box
  background and a finer ladder; the periodic sedimentation case above
  is the Phase-6 dynamics gate).
