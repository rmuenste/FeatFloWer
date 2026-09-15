# D6.2 — Jeffery orbit in a planar Couette box (CASE_SPEC)

Drafted 2026-09-08; CLOSED 2026-09-15 (rows d62_v0b_spin, d62_v1b_orbit, d62_v2_clearance). Second non-spherical family; the first live test of FREE
ROTATION end-to-end: FBM torque -> pe angular update (gyroscopic term, the
D-1-corrected inertia) -> orientation integration. D6.1 exercised only the
torque NULL on a fixed body; here the torque drives the physics.

Owner decisions taken:
- Purpose-built PLANAR Couette box - NOT the legacy `Couette_20x4x10` (that
  is the old cylindrical viscometer: curved shear + radially graded layers)
  and not the V2 annulus. Jeffery's solution is for linear shear; curvature
  would contaminate the gate.
- Shear Reynolds number <= 0.05 for the entire first pass; a finite-Re
  ladder only after the low-Re orbit is certified.
- Wall clearance is a FIRST-CLASS systematic: the owner has measured a real
  wall influence on this problem before ("I needed to put more space between
  the wall and the particle"). The clearance ladder (V2) is therefore a
  required gate, not an option, and the default box errs generous.
- Runs on Fritz (workload split).

Reference: `literature/jeffrey_1922.pdf` - G.B. Jeffery, Proc. R. Soc. A 102
(1922) 161-179, title-page verified (README note: filename keeps the
provider's misspelling). Finite-Re literature (Ding & Aidun 2000) NOT in the
folder - not needed while Re <= 0.05; acquire when the Re ladder starts.

---

## 1. Physics and gates

A neutrally buoyant prolate spheroid (r_e = a/b = 2), translation-locked but
FREE TO ROTATE, at the center of a planar Couette box: walls at z = +-H/2
moving at u = +-U x^, periodic in x (streamwise) and y (spanwise). Linear
shear gammadot = 2U/H, vorticity +gammadot along +y (corrected 2026-09-10;
the draft said -y).

Jeffery (1922), axis started in the shear plane (x-z) -> pure tumbling
orbit:

```
period      T = 2 pi (r_e + 1/r_e) / gammadot          [T*gammadot = 15.70796 at r_e = 2]
waveform    |dphi/dt| = gammadot (cos^2 phi + r_e^2 sin^2 phi) / (r_e^2 + 1)
            -> at r_e = 2: |dphi/dt| ranges gammadot/5 (axis || flow, slow)
               to 4 gammadot/5 (axis || gradient, fast) - a 4:1 modulation
```

Convention (pinned 2026-09-10 after the review found the printed formula
swapped slow/fast for the analyzer's angle): phi = atan2(axis_z, axis_x),
measured FROM THE FLOW AXIS in the x-z plane, read directly from the
per-step DNS_PART_AXIS record. With u = gammadot z x^ the vorticity is
+gammadot y^, so the physical spin is omega_y = +gammadot/2 for a sphere
and phi DECREASES (dphi/dt = -omega_y). Gates quote |dphi/dt|; the sign
is checked separately (analyzer line "rotation sense"). The formula
above in the gradient-axis convention reads gammadot (r_e^2 cos^2 +
sin^2)/(r_e^2+1) - same physics, different origin of phi.

| gate | quantity | band | rationale |
|---|---|---|---|
| G-sphere (V0, PRIMARY control) | sphere spin in shear: omega_y / gammadot = +1/2, i.e. (dphi/dt)/gammadot = -1/2 in the analyzer's angle, steady | +-1% | exact, shape-free; certifies the torque->rotation chain before any orbit |
| G-orient (V1, added 2026-09-10) | \|dphi/dt\| binned by phi vs Jeffery: rate at flow alignment / rate at the gradient axis = 1/r_e^2, rms residual over the bins | +-10% placement; rms < 3% of gammadot | extrema alone accept a phase-shifted orbit (review finding 1) |
| G-period (V1, PRIMARY) | measured T*gammadot vs 15.70796 | +-3% at default clearance, tightened by V2 | the integral test of the whole rotational path |
| G-waveform (V1) | dphi/dt(phi) vs Jeffery, incl. the 4:1 slow/fast ratio | ratio +-5% | shape test - catches compensating errors the period integral hides |
| G-wall (V2, REQUIRED) | period shift between the two clearances | measured, monotone toward Jeffery with clearance | the owner-observed wall effect, quantified not assumed |
| G-plane | axis y-component stays ~0 through >= 1 orbit | |axis_y| < 0.02 | in-plane orbit is neutrally stable at Stokes; drift flags integration error |

## 2. Geometry, units, numbers

Campaign box units (d11/d61 conventions): rho = 1, nu = 1, mu = 1.

- Body: a = 0.5, b = c = 0.25 (r_e = 2, length 2a = 1). Center at the box
  center, axis initially +x (slow phase of the orbit - gentlest start).
  particleDensity_ = 10.0 (AMENDED 2026-09-09, see "Rotational coupling
  stability" below; the original neutrally-buoyant 1.0 is unusable at this
  Re/dt). With translation locked and gravity off, density enters ONLY
  through the rotational inertia, i.e. the rotational relaxation time
  tau_rot ~ rho_p b^2 / (15 nu) = 0.042 t.u. at rho_p = 10 - against the
  orbit period 78.5 t.u. (tau_rot * gammadot = 0.008) the body still follows
  the zero-inertia Jeffery solution; the +-3% period gate is unaffected.
- Shear: gammadot = 0.2 -> Re_a = gammadot a^2 / nu = 0.05 (the owner cap;
  definition recorded). Wall speed U = gammadot H / 2.
- Default box (V1): L_x x L_y x H = 8 x 6 x 8 (16 body lengths of x-period,
  12 spanwise, wall-to-center clearance 4 = 8 semi-major axes). U = 0.8.
- Clearance rung (V2): H = 4 (clearance halved, U = 0.4 to keep gammadot),
  same L_x, L_y - one parameter moves.
- Duration: transient ~5 t.u. + >= 1.25 periods; T = 78.54 at gammadot=0.2
  -> run to t = 105 (dt = 0.01 -> 10500 steps).

Resolution (thin axis 2b = 0.5): direct-writer coarse box 19 x 15 x 19
(h ~ 0.42), production at level 4: h = 0.0526 -> 2b/h = 9.5 (certified
class), 2.77M elements. Level-3 smoke first (2b/h = 4.75, minutes).

**Rotational coupling stability (found 2026-09-09, V0/V1 first attempt).**
The FBM-PE coupling is explicit: one torque exchange per step, omega
advanced by dt * I^-1 * T. The fluid torque relaxes omega toward the local
half-vorticity on the viscous time scale tau_rot = I / C_T, for a sphere
rho_p a^2 / (15 nu). The explicit update is stable only for
g = dt / tau_rot * u < 2 (u ~ 1.5-2.5, the short-time unsteady torque
enhancement a / sqrt(pi nu dt)). At rho_p = 1, nu = 1, a = 0.25, dt = 0.01:
tau_rot = 0.0042, g ~ 3.6 -> torque sign-alternates and grows ~2.5x per step
(measured: 0.058, -0.154, ... -2.5 at t=0.1, +23.7 at t=0.15, -94 at
t=0.2), NaN by t = 0.24 (sphere) / 0.33 (spheroid). Nondimensionally
g = 15 u (gammadot dt) / (rho_r Re_b): at the Re cap the only levers are
steps per period N and rho_r, g = 18850 u / (N rho_r). Choice: rho_r = 10
at dt = 0.01 (N = 7850): g ~ 0.36 (sphere), smaller for the spheroid,
monotone convergence with tau_rot = 0.042 t.u. (converged before t = 0.3).
The alternative dt = 0.002 costs 5x the wall time for the same effect.
This diagnosis stands on its own evidence (sign-alternating torque growth
at the predicted gain, cured by the predicted inertia change). It is NOT
the ten Cate "dt floor": that one was refuted as PE/CFD step
desynchronization (row dt_stability_refuted) - the draft's analogy is
withdrawn.

## 3. Fixture requirements (the build work)

1. **Mesh**: `gen_d62_shearbox.py` in the d52-writer style (closed-form node
   placement, uniform cells - no grading). z walls tagged as moving walls;
   x/y faces periodic. Periodic-face encoding to be copied from the QBOX9
   fixture (the d11/d31 PERIODIC_COMM machinery); G0 verifies pair counts
   the way d31_smoke did.
2. **Partition**: axis_uniform grid (periodic case - METIS invalid), e.g.
   5 x 4 x 7 = 140 subdomains -> 141 ranks on 2 fritz nodes. STAGING RISK
   (memory `axis-uniform-tolerance-defect`): legacy partpy `-4` is the
   correct tool; pe_partpy's span-relative tolerance mis-partitions. G0
   includes partition verification before any physics.
3. **Moving-wall BC**: constant-velocity wall profile (+-U x^) - an
   Inflow-type entry in the app's GetVeloBCVal, mirroring the viscometer's
   rotating-wall profile. Deck-parametrized U.
4. **pe setup - rotation-only body**: the D6.1 xyz path calls
   setFixed(true), which kills rotation. New json key
   `particleMotion_`: "fixed" (default, = D6.1) | "rotationOnly"
   (setLinearDofMask(0,0,0), NOT fixed, angular DOFs free). Sphere path and
   D6.1 decks byte-identical with the key absent; e4_l3 twin gates the pin.
5. **Analysis**: `tools/d62_jeffery_analysis.py` - phi(t) from
   DNS_PART_AXIS (unwrapped atan2), period from successive zero crossings,
   dphi/dt(phi) waveform against Jeffery, omega_y/gammadot for V0.
6. Lubrication OFF (sphere-only model; D6.1 guard already refuses it for
   ellipsoids). PE serial mode; json stepsize_ = deck dt = 0.01.

## 4. Run matrix (all Fritz)

| run | body | box | level | duration | est. cost |
|---|---|---|---|---|---|
| G0 | ellipsoid, rotationOnly | 8x6x8 | L3 | 50 steps | minutes - creation, periodic pairs, axis record, translation lock |
| V0 | r_e=1 spheroid (semiAxes 0.25^3, so DNS_PART_AXIS fires), rotationOnly | 8x6x8 | L4 | to t=10 (1 rad of spin at omega = gammadot/2; tau_rot = 0.04) | one segment - the omega = -gammadot/2 control |
| V1 | r_e=2 | 8x6x8 | L4 | t=120 (3 pi-crossings) | DONE 2026-09-13 (row d62_v1b_orbit): PASS all gates - period +0.30%, modulation +0.8%, orientation-resolved waveform 0.48% rms; three 24h segments at ~21 s/step with orientation continuation |
| (V0/V1 first attempt, rho_r = 1) | | | | NaN at t = 0.24 / 0.33 | rotational coupling instability - see section 2; kept as evidence (fritz jobs 4199314/15) |
| V2 | r_e=2 | 8x6x4 | L4 | t=120 | DONE 2026-09-15 (row d62_v2_clearance): PASS - period +0.81% at H=4 vs +0.30% at H=8, wall effect +0.51 pp shrinking with clearance, clearance-free extrapolation +0.1..+0.25%; D6.2 CLOSED |
| V3 (optional) | r_e=2, axis +y start | 8x6x8 | L4 | shorter | log-rolling state - only after V1/V2 gate |

Cost of the required set: roughly 3-6 node-days - between D6.1 and one
viscometer segment.

## 5. What this family deliberately defers

Finite-Re orbit physics (period growth, orbit-plane stability - needs Ding
& Aidun in the literature folder first), orbit-constant drift over many
periods, log-rolling vs tumbling selection, and any suspension of rotating
bodies (D6.3 scope). The clearance ladder measures the wall systematic; it
does not attempt an analytic wall correction.

## 6. Relation to the rest of the campaign

Independent of the D/h=16 rung (v25f) - different domain, observable, and
nodes; runs in parallel. Builds directly on D6.1's setup keys, the
DNS_PART_AXIS record, and the twin-gate discipline. The pe torque-path
review (docs/md_docs/dns_torque_path_review.md) is the code-level
foundation: what D6.1 verified as a null, D6.2 promotes to the observable.
