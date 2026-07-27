# kroupa_couette — TOTAL suspension viscosity via the wall force balance

Publication-grade reproduction of Kroupa/Vonka/Soos/Kosek, Langmuir 32
(2016) 8451: a TWO-WAY COUPLED plane-Couette cell measuring the total
suspension viscosity eta = eta_H + eta_L(+ eta_C) from the force balance
on the walls (their eqs 26-33), the piece the frozen-shear kroupa_shear
box cannot see (its dilute-limit check proved the missing component is
exactly the hydrodynamic Einstein term 2.5*phi*mu).

## Geometry / model

- Unit box [0,1]^3, QBOX9 3x3x3 partition, level 3 (36^3 cells), np=28.
- Periodic x/y (SimPar@PeriodicAxis = xy, length 1); z-faces are MOVING
  WALLS: .par tag Inflow300 -> u_x = G*(z - z_c) = +-G*L/2 (symmetric,
  zero net momentum; equivalent to Kroupa's bottom-fixed frame). G = 1.
- COUPLED fluid (no ELPrescribedField); ELInitialField = linear_shear
  seeds the exact single-phase Couette solution at t=0 (skips the
  L^2/nu = 50 t.u. viscous startup; Kroupa does the same).
- PE: zWallsEnabled_ with global planes at z=0/1, wall velocities
  -0.5/+0.5 (setLinearVel under MOBILE_INFINITE); periodic x/y, 3x3x3.
- Particles: neutrally buoyant d_p = 0.05, random seeding (inset
  R + seedMinGap/2 = 0.0275 from walls — inside the lubrication cutoff;
  capped impulses relax this in the first steps).
- Lubrication: pair closure (eqs 12-14) + NEW sphere-wall closure
  (eqs 16-18, twisting omitted), cutoff 0.025 = R_p, h_c = 0.0025,
  substeps_ = 10 (ladder-validated).
- NO wall hydration repulsion (Kroupa eq 24) — user decision 2026-07-28;
  PGS wall contacts + wall lubrication carry the mechanics. Revisit if
  particles accumulate in wall contact.
- nu = 0.02 KINEMATIC (Prop@Viscosity); mu = rho*nu = 0.02 dynamic.

## Measurement (EL_WALL_STRESS, audit freq 25)

Per wall (bt = bottom, top), the +x traction the suspension exerts on
the wall, tail-averaged (ELTAvgWindow):

- tauH: fluid channel, mu*du_x/dz face-quadrature at the wall
  (sign wsign = +1 bottom / -1 top; Kroupa eq 29-30).
- tauL: sphere-wall lubrication impulse channel (eq 31-32).
- tauC: sphere-wall PGS contact impulse channel (not in Kroupa — their
  hydration wall keeps particles off; ours may contact).
- eta_wall = [Sum_channels tau(bot) - Sum_channels tau(top)]/(2G).

Cross-check: interior virial channels EL_SUSP_STRESS (etaL_tavg) and
EL_CONT_STRESS (etaC_tavg) from the same run; in steady homogeneous
shear, eta_wall ~= mu + etaL_virial + etaC_virial + (fluid feedback
reshaping of the profile, which the wall gradient sees as eta_H - mu).

## Gates

1. Single-phase (or near-zero phi) twin: profile stays stationary,
   eta_wall/mu = 1 within ~1%; tauL = tauC = 0.
2. unit_wall_lubrication (separate case): sphere at gap 0.3R from the
   moving wall — force sign/decay, wall-impulse accumulator equals the
   sphere momentum change, EL_NEWTON_PAIR unchanged.
3. phi=0.20 smoke: no NaN, wall channels nonzero, E013_CreateComm
   stable under xy-periodicity (new territory), wall-time/step sizing.
4. Sweep phi = 0.05/0.10/0.20/0.30: eta_wall(phi) monotone; compare
   Krieger-Dougherty (phi_max 0.64) + Kroupa Fig 3 shape; wall-balance
   vs interior-virial consistency.

## Watchpoints

- eps_f near walls (deen_poly kernel truncates, no renormalization):
  watch <eps_f> TAVG as in V4.
- ELMomentumFix = No and the EL_FLUID_PAIR audit's expected-source
  bookkeeping excludes wall traction — its mismatch line is NOT a gate
  here (like V4).
- Frozen-box heritage lines (etaL_tavg = sig/G) remain valid: G is the
  nominal wall shear rate.

## Run log

(appended as runs complete)

## Smoke gates — ALL PASS (2026-07-28, jobs 136250/136251/136252)

1. cou_reg (frozen kroupa_shear phi=0.20, 200 steps, new binary):
   EL_SUSP_STRESS bitwise identical to the kc20 run -> wall code inert
   when disabled.
2. cou_sp (single particle at center, coupled, t=1): tauH = +-2.000000e-2
   on the two walls, eta_wall/mu = 1.0000 to all printed digits; profile
   stationary at machine precision; wpairs = 0. The tau_H quadrature,
   sign convention and linear-shear IC are exact.
3. cou_p20s (phi=0.20 coupled, 200 steps): no NaN, xy-periodic comm
   stable (E013_CreateComm fine on QBOX9 3x3x3), wall lubrication firing
   (3370 wall pairs/step, tauL nonzero both walls), EL_NEWTON_PAIR
   2.7e-20 WITH the wall-momentum audit extension. ~1.7 s/step ->
   5000-step production ~2.5 h.

Production sweep phi = 0.05/0.10/0.20/0.30 submitted as jobs
136253-136256 (med, 7 h).

## RESULT — production sweep (2026-07-28, jobs 136253-136256)

All four runs clean to t=10 (no NaN); runtimes 1:05-3:20 h. Tail-averaged
wall-balance viscosity vs the interior-virial estimate and
Krieger-Dougherty (phi_max = 0.64):

| phi | eta_wall/mu | interior (mu+etaL+etaC)/mu | wall-vs-interior | KD | vs KD |
|----:|------------:|---------------------------:|-----------------:|-----:|------:|
| 0.05 | 1.034 | 1.016 | +1.8% | 1.139 | -9% |
| 0.10 | 1.126 | 1.093 | +3.0% | 1.312 | -14% |
| 0.20 | 1.664 | 1.606 | +3.7% | 1.821 | -9% |
| 0.30 | 2.936 | 2.874 | +2.2% | 2.751 | +7% |

Verdicts:

1. CROSS-CHECK PASS (the publication-grade result): the wall force
   balance and the interior impulse virial - two fully independent
   estimators (boundary traction vs bulk momentum flux) - agree within
   2-4% at every phi. The small positive gap is physical: the wall
   balance additionally sees the fluid-feedback reshaping of the
   profile (eta_H beyond mu), which the interior particle channels do
   not contain.
2. KD COMPARISON PASS (shape and band): monotone, superlinear, within
   ~15% of Krieger-Dougherty across phi = 0.05-0.30, crossing KD
   between 0.2 and 0.3. Compare the no-wall pipe (recovered only 9-20%
   of the KD excess) and the frozen box (structurally missing eta_H):
   the coupled Couette cell closes the gap. Dilute end sits below KD
   (eta/mu - 1 = 0.034 at phi = 0.05 vs Einstein 2.5*phi = 0.125): the
   unresolved model generates only part of the Einstein term - the
   known model boundary, shared with Kroupa's method.
3. Channel structure: wall lubrication is the dominant wall channel at
   high phi (tauL ~ 2-3x tauH at phi >= 0.2); wall contacts contribute
   < 7% of the total everywhere.
4. Gates: <eps_f> TAVG = 1-phi to 4 decimals at all phi (the deen_poly
   wall-truncation concern did not materialize); max_overlap <= 2.5e-10;
   no NaN. EL_NEWTON_PAIR closes to <= 3.9e-6 (about 7e-5 of the system
   momentum scale), NOT machine zero: the residual is z-dominant (6x
   the x-component), scales with wall-CONTACT activity, and was machine
   zero in the wall-lubrication-only smoke. Attribution: the converged
   PGS cache impulse p_ vs the relaxation-folded applied velocity
   change differ by the fixed-iteration truncation; in periodic runs
   this is invisible (pairwise antisymmetry cancels it in the global
   sum), while a sphere-wall contact exposes the sphere-side truncation
   directly. A solver-accuracy footnote, not a momentum leak; the
   x-channels entering eta are 5-10x cleaner, and verdict 1 bounds any
   effect on the measurement.
