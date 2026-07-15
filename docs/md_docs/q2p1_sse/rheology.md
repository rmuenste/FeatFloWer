# q2p1_sse Rheology

Primary files:

- `source/src_util/ReadExtrud3DParameters.f90`
- `source/src_quadLS/QuadSc_user.f90`
- `applications/q2p1_sse/app_init.f90`

## State Contract

Material and rheology setup is parsed into `myMultiMat`. The viscosity functions
then select the active material by index and use its `tRheology` description.

## Main Viscosity Entry

`AlphaViscosityMatModel(NormShearSquare, iMat, Temperature)` is the main
SSE-relevant viscosity entry found in `QuadSc_user.f90`.

It supports temperature shift factors when `Temperature` is present, including
the configured `AtFunc` variants, then evaluates the configured rheology model.

Observed model branches include:

- Carreau/Paderborn-style branch
- normalized power-law branch
- Polyflow Carreau branch
- Ellis branch
- Hogen power-law branch
- Bingham branch

## Initialization Use

`General_init_ext` prints flow-curve tables for every material around
`myProcess%T0` and a range of shear rates. This is an early validation point:
if material parsing or rheology changes are wrong, this table is often the first
visible symptom.

## Automatic Timestep Use

When `mySetup%bAutomaticTimeStepControl` is enabled, initialization estimates a
characteristic viscosity and timestep:

- For `SSE`, `TSE`, and `XSE`, the characteristic size and velocity are derived
  from barrel dimensions and screw speed.
- For `DIE`, the characteristic size and velocity are derived from extrusion
  gap size and extrusion speed.
- `AdjustTimeStepping` updates `TSTEP`, `DTGMV`, and `TIMEMX`.

This makes rheology part of startup control flow, not only the nonlinear
assembly/solver path.

