# s07_pipe_viscosity - viscosity convention of the fluid operator (developed pipe flow)

Purpose: show that the fluid's viscous term carries the dynamic viscosity. The solver's
diffusion matrix takes `Properties%Viscosity` without a density factor against mass and
convection matrices that carry the density, i.e. it expects mu. Euler-Lagrange decks give the
kinematic viscosity (`Prop@Viscosity`, or `Prop@DynVisc`, which the parser divides by the
density), because the drag closures form mu = rho nu themselves. For rho_f /= 1 the fluid's
viscous term was therefore too small by the factor rho_f. Key (2026-11-15):
`SimPar@ELFluidViscosity = kinematic | dynamic`; default `dynamic` since 5a882280 (owner decision
2026-11-15), `kinematic` = campaign behaviour.

Branch feature/el-fluid-viscosity (worktree FeatFloWer-visc), commits eb41fc84 (key, default
still kinematic), ad9a93a1 (review fixes), 5a882280 (default dynamic); binary of the recorded runs build-el-e (ad9a93a1, dirty 0), pe 43e8fb5.
Run directories: /data/warehouse17/rmuenste/runs/visc07.

## Set-up

Column mesh of case s05 (pipe R 13.85 mm, length 55.4 mm, O-grid, level 3, 27 z slabs, 28 ranks;
mesh not committed, see the s05 RUNBOOK), uniform inlet U = 0.002 m/s (boundary type Inflow975,
constant), open outlet, rho 10, `Prop@DynVisc = 1.5e-3` (nu 1.5e-4), Re = rho U D / mu = 0.37.
Single-phase flow: 30 resting spheres from the floor of the s05 bed only because the application
needs particles; `ELApplyParticleForces = No`, `ELApplyFluidFeedback = No`,
`ELVolumeAveraged = none`. dt 2 ms, 400 steps (0.8 s; viscous time R^2/nu = 1.3 s with mu,
12.8 s with the campaign coefficient). Decks `q2p1_param_kinematic.dat` /
`q2p1_param_dynamic.dat` (one line apart; since the default flip the kinematic deck sets the key
explicitly and the dynamic deck would not need it), `example.json`, `particles.xyz`.

Reading: `tools/poiseuille_dp.py <rundir> <frame> --uin 0.002 --mu 1.5e-3 --radius 0.01385`:
layer-mean pressure against z, least-squares slope over the second half of the pipe
(27.7 ... 49.7 mm), against Hagen-Poiseuille dp/dz = -8 mu U / R^2 = -0.12512 Pa/m; maximum axial
velocity over U (2 when developed).

## Result (jobs 160654 kinematic, 160655 dynamic; 4 min 45 s each)

| Convention | frame | dp/dz, Pa/m | / Hagen-Poiseuille | max u_z / U |
|---|---|---|---|---|
| kinematic (default) | 0.4 s | -1.6053e-2 | 0.1283 | 1.509 |
| kinematic (default) | 0.8 s | -1.3795e-2 | 0.1103 | 1.760 |
| dynamic | 0.4 s | -1.25132e-1 | 1.0001 | 1.9997 |
| dynamic | 0.8 s | -1.25113e-1 | 1.0000 | 2.0001 |

- `dynamic`: developed after 0.4 s, pressure gradient 1.0000 of Hagen-Poiseuille (-0.003 %).
- Default: the profile is still developing at 0.8 s (effective kinematic viscosity nu/rho); the
  gradient falls towards 0.100 = 1/rho of the reference. The defect is the factor rho.
- Identical numbers on the first build of the key (eb41fc84, jobs 160634 / 160635).

## What else moves with the key

Start-up case s06 (rho 1000) with `ELFluidViscosity = dynamic`, all four formulations (jobs
160656-160659): time constants 8.8440 / 4.4931 / 7.0569 / 7.2214 ms, unchanged in all digits
shown against the default; the particle viscous force of Bfull / A grows from -6.1e-9 N to
-6.0e-6 N (factor 980) and stays negligible. The suspension start-up does not depend on the
fluid's viscous term. Not rerun with the key: the Ergun column s05 (drag-dominated bed; the wall
shear of the free column above the bed changes), the column settling cases.

## Supported settings and guards (negative tests)

Refused at start-up with exit 1: invalid value (job 160661), `MatrixRenewal` S1 (160662: the
stress matrices are built from nu at initialisation), `Prop@PowerLawExp` 0.8 (160663); also
non-Newtonian / viscoelastic, D0. A property reload at run time (`Reload_Prop` of the process
control) would restore nu: the run then stops at the next Euler-Lagrange pass (guard in
`EL_CLOSURE_VISCOSITY`; not exercised by a run). The deck convention is unchanged (nu or
`Prop@DynVisc`); `fluidViscosity_` of the particle json is an independent input (mu, as before).

## Certification of ad9a93a1

- Default path: ten certified EL configurations bitwise against runs/coupling_p4p5/ergun_gates2,
  strict comparison, jobs 160644-160653 (runs/coupling_p4p5/visc_gates2); first round on
  eb41fc84: 160624-160633, bitwise. Bfull at the default key bitwise against the certified s06
  run (160660).
- With the key: the four s06 runs bitwise between eb41fc84 and ad9a93a1.
- DNS twins (build-dns2, certified DNS recipe): e4_l3 job 160664, DKT job 160665 bitwise against
  the certified heads (md5 ca98e7af225b / 89abce88b675 and 82411f64f661 / cfb8b42d0969); first
  round 160642 / 160643 bitwise.
- Unit tests: test_el_momentum_model (viscosity convention: three checks; force split under the
  dynamic convention), test_el_kernel_forces, test_el_va_flux, test_el_transfer pass.
- Reviews of eb41fc84 (GPT-6 Astra high; an independent Claude agent): no finding on the default
  path; every viscosity consumer the EL decks reach confirmed consistent under the key. Fixed in
  ad9a93a1: property reload reverts the conversion unnoticed (both) - guarded stop; S matrices
  of MatrixRenewal S1 keep nu (both) - refused, with D0; PowerLawExp /= 1 not refused (GPT) -
  refused; conversion state public and written from outside (agent) - private behind a setter;
  tests did not reach the particle forces under the key (both) - added; echo glued to a
  progress line, stop messages from every rank (agent). Left and stated: the conversion routine
  itself has no unit test (covered by this case); the wall-stress diagnostic
  (`el_diagnostics.f90`) always forms rho nu, so under the default it reports rho times the
  stress the operator carries; `init_q2p1_ext` in the application file is not used by the
  application (marked); with a thermodynamic density model the parser's division and the
  conversion could use different densities (not an EL configuration).

## Default flipped to dynamic (5a882280, 2026-11-15)

Owner decision. Binary build-el-g (5a882280, dirty 0). The echo line `EL fluid viscosity = ...`
is now printed in every run.

- Pipe case: deck without the key 1.0000 of Hagen-Poiseuille (job 160695), deck with
  `ELFluidViscosity = kinematic` 0.1103 at 0.8 s (job 160694); both identical to the table above.
- Certified EL configurations at the new default (jobs 160680-160689,
  runs/coupling_p4p5/visc_gates3): the six with rho_f = 1 (tier-2 twin, semi-implicit twin, P4
  frozen explicit / semi, P5 with and without feedback) bitwise; the four with rho_f = 10 (s01
  two-way level 2, VB_L2 default, T9a, T10a4) not bitwise, as intended.
- The same four with `SimPar@ELFluidViscosity = kinematic` added to the deck (jobs
  160690-160693): bitwise. A certified case with rho_f /= 1 is reproduced with this one key.
- DNS twins (build-dns3): e4_l3 job 160699, DKT job 160700 bitwise against the certified heads.
- Consequence for decks: a deck with MatrixRenewal S /= 0, D0, `Prop@PowerLawExp` /= 1 or a
  non-Newtonian fluid is refused at start-up unless it sets `kinematic`.
- Ergun column reruns with the new default: s05 RUNBOOK.

## Open

- Size of the change in the four rho_f = 10 certified cases (s01, s04 VB) at the new default:
  not evaluated; their datasheet rows belong to the `kinematic` convention.
- Reruns of rho_f /= 1 cases with the key where the fluid's viscous term matters (wall-bounded
  columns).
