# Torque-path review: FF FBM torque calculation → pe torque application

Reviewed 2026-09-03 as D6 (non-spherical particle) preparation, on owner
request, before starting the D/h=16 resolution rung. Scope: the per-particle
hydrodynamic force/torque pipeline — FEM integration in FeatFloWer, MPI
reduction, handoff through the C interface, and application inside the pe
rigid-body update (serial-PE mode, `HardContactAndFluid` solver). The
viscometer instrument torques (`VISC_*`, bob-wall `BndrForce` sum) are outside
this scope — they carry their own eight-fold five-digit metrology record in
the datasheet.

Motivation: every campaign validation so far used spheres, for which the
rotational half of the pipeline is nearly unobservable (isotropic inertia, no
orientation dependence in the α field, zero mean torque in most gates). D6
makes rotation first-class physics; this review is the audit of that path
while it is still cheap to fix.

## Path audited

1. `source/src_quadLS/QuadSc_force_serial.f90` —
   `ForcesLocalParticlesSerial_Standard` / `_KVEL`: element-wise volume-form
   traction `t = -p n + ν(∇u+∇uᵀ)·n` with `n = -∇α_h`, torque
   `∫ (x - x_c) × t`, `x_c` the live particle center; `Prop@ForceScale`
   factors applied per force/torque component; one `COMM_SUMMN(forceArray, 6N)`.
2. `source/src_particles/dem_query.f90` — `setForcesMapped` →
   `setLocalParticle2` (struct-coupled force+torque+identity path).
3. `libs/pe/src/interface/object_queries.cpp` — `setPartStruct`:
   `body->setForce(f); body->setTorque(t)` (overwrite semantics);
   `synchronizeForces()` → `body->applyFluidForces(stepsize)`.
4. `libs/pe/pe/core/rigidbody/RigidBody.h` — `applyFluidForces`:
   `v += dt·F/m`, `w += dt·(R·Iinv_body·Rᵀ)·T`, optional under-relaxation of
   BOTH force and torque (campaign default: fully explicit), then
   `resetForce()`.
5. `libs/pe/pe/core/collisionsystem/HardContactAndFluid.h` —
   `resolveContacts` velocity seeding (gravity/buoyancy + gyroscopic Euler
   term), contact impulses via world-frame `getInvInertia()`,
   `integratePositions` (exponential-map quaternion update, `R_` rebuilt).

## Verdict: sound for spheres; three defects block trustworthy non-spherical runs

(Two found by the source read-through; a third — the volume/mass factor — was
caught by the new unit test's mass check during the fix, see D-3.)

### Sound (verified, not assumed)

- Torque moment arm is taken about the current particle center each step, and
  the traction integrand is the same one the drag gates certified (D1
  Hasimoto five-digit lattice corrections, D3 Beetstra, D2 ten Cate).
- No double application of the fluid wrench: `setForce`/`setTorque` overwrite,
  `applyFluidForces` consumes into (v,w) and resets the accumulators, so the
  solver's own `initializeVelocityCorrections` (which re-reads
  `getForce`/`getTorque`) sees zeros from the fluid side.
- Frame handling is correct everywhere torque meets inertia:
  `getInvInertia()` returns `R·Iinv_body·Rᵀ` (RigidBody.h:915-918) and is used
  in `applyFluidForces`, in contact impulses, and in the velocity seeding.
- The seeding includes the gyroscopic term
  `w += dt·I⁻¹((I·w)×w)` — required for spheroid tumbling (Jeffery), present
  for all body types.
- Orientation integration: `dq = Quat(φ, |φ|)`, `q ← dq·q`, `R` rebuilt —
  standard first-order exponential map, adequate at campaign dt.
- Under-relaxation (`applyFluidForces` relaxation parameter) treats force and
  torque identically; campaign runs use the fully explicit default.
- `ForceScale` deck factors scale force and torque components independently
  (symmetry decks depend on this).

### Defect D-1 (pe): ellipsoid inertia tensor wrong by 25% on the a-b axis

`pe/core/rigidbody/EllipsoidBase.h`, `calcInertia()`:

```
I_[0] = 0.2  m (b² + c²)   // correct
I_[4] = 0.2  m (a² + c²)   // correct
I_[8] = 0.25 m (b² + a²)   // WRONG — solid ellipsoid factor is 1/5, not 1/4
```

The commented-out block directly above has the correct `0.2` on all three
axes; the live line looks like an experiment that was never reverted.
Consequence: I_zz overestimated 25% for every ellipsoid; rotational response
about the c-axis correspondingly too sluggish. Invisible to spheres (own
class, correct `2/5 m r²`). Must be fixed (with a unit test pinning all three
principal values) before any D6 run that lets an ellipsoid rotate.

### Defect D-2 (pe): ellipsoid gets full gravity — no buoyancy reduction

`pe/core/collisionsystem/HardContactAndFluid.h`, `resolveContacts` velocity
seeding (~line 2031-2110): sphere, capsule, cylinder, and box each seed

```
v += dt · g · vol·(ρ_body − ρ_liquid)/m      (buoyancy-scaled gravity)
```

but ellipsoid (and non-fixed triangle mesh) falls through to the generic
`else` branch, which seeds `v += dt·g` — full gravity, no `(ρ−ρ_f)` scaling.
At the campaign density ratio 1.1 that is an ~11× overestimate of the net
driving force for a settling ellipsoid: a D6.1 Oberbeck gate run this way
fails at leading order before hydrodynamics enters. Fix: give the ellipsoid
the same buoyancy branch (or better, collapse the four copy-pasted per-type
branches into one generic `body->getVolume()`-based path so new shapes cannot
miss it).

### Defect D-3 (pe): ellipsoid volume and mass low by 25%

Found while pinning D-1 with the new unit test (its mass check failed after
the inertia fix): `EllipsoidBase::getVolume()`, `calcVolume()`, and
`calcMass()` all returned `pi*a*b*c` — the correct `(4/3)*pi*a*b*c` lines sat
commented out beside them, same fingerprint as D-1. Consequence: mass, weight,
buoyancy volume, and (through the mass factor) the whole inertia tensor were
25% low. Introduced in the same creep-bench bring-up lineage (c19ed8f); the
only `createEllipsoid` caller anywhere is commented-out code in
`pe/interface/setup_creep.h`, so nothing live ever depended on the wrong
values. Residual left in place: `EllipsoidBase::calcDensity` still carries a
sphere-signature copy-paste `(0.75 m / pi r^3)`; it has no callers.

### Notes (not defects, staging constraints)

- N-1: the FF torque moment arm `x − x_c` has **no periodic minimum-image
  convention**. A body crossing a periodic boundary gets a corrupted torque
  (force is unaffected). All campaign cases to date had fixed particles or
  walled domains. D6.1 in the d11-style periodic cell must therefore use a
  fixed body at the cell center (D1/D3 discipline: body force on the fluid,
  measure the wrench) or fix the arm before allowing free motion.
- N-2: `synchronizeForces()` (object_queries.cpp:254-255) contains a no-op
  `Vec3 tau = body->getTorque(); body->setTorque(tau);` in the
  non-sphere branch — leftover scaffolding, remove on the next pe touch.
- N-3: `RigidBody` has a linear DOF mask (`linearDofMask_`) but no angular
  counterpart; a "translation-locked but free to rotate" D6 configuration is
  available via `setFixed`/mass tricks only. Fine for planned D6.1/D6.2.
- N-4: the sphere branch of the velocity seeding carries an unconditional
  per-step stdout diagnostic gated on `cfdRank == 1`; it did not fire in any
  d52 log (0 hits), but it is a spam hazard for future configs.

## Consequences for the campaign

- The **D/h=16 resolution rung (spheres) is unaffected** by both defects and
  proceeds on the frozen certified build.
- FIXED 2026-09-03, same day: D-1 + D-2 + D-3 repaired in pe commit 6971b13
  with new test `pe_ellipsoid_inertia_test` (volume, mass, all three principal
  inertia values, I*Iinv identity, sphere degeneracy against the Sphere
  class); full pe serial suite 15/15; e4_l3 bitwise twin PASS (job 142801,
  particle streams byte-identical — datasheet row pe_ellfix_twin). Pin bump
  pending owner push.
- D6.1 design should use a fixed spheroid + body-force-driven flow (N-1),
  which also sidesteps D-2 entirely; D-2 still must be fixed for any free
  ellipsoid later (D6.2+).
