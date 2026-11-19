# Preparation of cross-code case 11 (minimum fluidisation) in q2p1_el_pipeflow (2026-11-19)

Branch `feature/el-min-fluidization` (worktree `FeatFloWer-fluid11`, base c76ca927 = master + periodic
kernel + periodic coarse pressure operator), pe pin 43e8fb54 (unchanged). Case spec (peer session,
read-only): `liggghts-dem/LIGGGHTS-PUBLIC/md_docs/case_specs/11_min_fluidization/SPEC.md`. Case directory:
`applications/q2p1_el_pipeflow/validation_cases/s11_min_fluidization/` (RUNBOOK, decks, json, tools). Runs:
`/data/warehouse17/rmuenste/runs/fluid11/<run>`.

Marks: **[V]** verified here (code read, unit test or run), **[I]** inferred, not verified.

## 1. Drag closures: table, Gidaspow, Wen-Yu, Ergun

`SimPar@ELDragModel` now has a table in `EL_CONFIG` (`el_drag_names`, ids `EL_DRAG_*`, column
`el_drag_standard_form`). The validation accepts exactly the names of the table; the momentum-model
refusal of the campaign closure `difelice` reads the column; `EL_FORCES:EL_DRAG_CLOSURE` looks the id up
and dispatches to one function per closure family. A further closure is one table entry, one function and
one CASE line. The existing closures were moved verbatim into `EL_DRAG_B_DIFELICE` and
`EL_DRAG_B_SINGLE_SPHERE` (default path bitwise, section 5) [V].

New closures [V, code and unit test]:

| name | per-particle coefficient drag_B (F_d = drag_B (u_f - v)) |
|---|---|
| `ergun` | V_p [150 (1 - eps) mu / (eps d^2) + 1.75 rho_f \|u - v\| / d] |
| `wenyu` | V_p 3/4 (C_d \|u - v\|) rho_f eps eps^-2.65 / d, C_d = 24/Re (1 + 0.15 Re^0.687) for Re < 1000, 0.44 above, Re = rho_f eps \|u - v\| d / mu |
| `gidaspow` | `ergun` for eps <= 0.8, `wenyu` above |

u - v is the interstitial slip (the kernel-sampled fluid velocity, interstitial under the continuity stage
that every momentum model implies). Zero slip: finite without a special branch (Ergun 25 pi mu d (1 - eps)/eps;
Wen-Yu: C_d |u - v| is formed as 24 mu/(rho eps d)(1 + 0.15 Re^0.687), limit 3 pi mu d eps^-2.65, reached
continuously; at slip 1e-12 m/s the excess is exactly 0.15 Re^0.687 = 1e-8).

**Convention.** beta is the interphase coefficient per unit mixture volume of the two-fluid model A
(Gidaspow): the particle phase receives beta (u - v) besides its share -(1 - eps) grad p. With
n = (1 - eps)/V_p, drag_B = V_p beta / (1 - eps). This is the form of every closure of the table
(`difelice_std` is the same form: 0.5 rho C_d A |u - v| eps^(2 - chi) = V_p beta/(1 - eps) with Di Felice's
beta) [V, code]. The momentum models consume it as follows (`EL_CONFIG` model table, `EL_MM_DRAG_SCALE`) [V, code]:

| model | particle force | drag form |
|---|---|---|
| A, Bfull | F_d + pressure force -V_p grad p (+ viscous force) | V_p beta/(1 - eps) |
| B | F_d / eps, no pressure force | V_p beta/(eps (1 - eps)) |
| none | F_d, pressure force per the campaign keys | V_p beta/(1 - eps) |

Fixed-bed identity (unit test): with particles at rest and superficial velocity U the model-A drag per
volume is n F_d = beta U/eps; with the pressure-gradient force the total is n F_d/eps = the Ergun gradient
150 mu (1 - eps)^2 U/(eps^3 d^2) + 1.75 rho (1 - eps) U^2/(eps^3 d), to 1e-12 at eps 0.40 / U 0.10 and eps
0.36 / U 0.25 (gidaspow) [V]. The spec's reference script `umf_reference.py` uses the same relation
(`beta * ur / eps`).

**CFDEM and MFiX [V, sources read under /data/warehouse17/rmuenste/code/liggghts-dem, read-only]:**

- CFDEM `GidaspowDrag.C`: `betaP` is beta/(1 - eps) ("per m3 PARTICLE"); `dragCoefficient = Vs*betaP`,
  divided by voidfraction only for `modelType "B"`; dense branch for voidfraction <= `switchingVoidfraction`
  (default 0.8), Wen-Yu above with the Schiller-Naumann C_d and no cap at Re 1000; Re = d eps |Ur|/nu.
- MFiX `model/drag_gs.f` `DRAG_GIDASPOW`: `lDgA` is beta/(1 - eps); Ergun for `EPg <= 0.8`, Wen-Yu above
  with C_d = 0.44 for Re > 1000; Re = d |u - v| rho eps / mu. `model/des/drag_gp_des.f`: `F_gp = DgA*PVOL`
  (model A), `/EP_G` for model B. MFiX returns zero drag at exactly zero slip (`RE == 0`); here the limit.

Both codes therefore use exactly this per-particle form for model A and the 1/eps form for model B. The
switch is at eps <= 0.8 (Ergun at 0.8 itself) as in both codes; the brief asked for Ergun below 0.8 and
Wen-Yu from 0.8 on; the difference is the single value eps = 0.8. The branches do not join at 0.8:
Ergun/Wen-Yu = 1.1533 at zero slip and 1.1427 at slip 0.5 m/s for the case-11 particle [V, unit test].
Above Re 1000 the codes differ (CFDEM keeps the Schiller-Naumann C_d); irrelevant here (Re about 9 at U_mf).

U_mf from weight gradient = drag (bisection as `umf_reference.py`): gidaspow at eps 0.40 0.210383 m/s
(spec table, column Ergun: 0.210), wenyu 0.169723 m/s (column Wen-Yu: 0.170) [V, unit test].

## 2. Inlet velocity table

`SimPar@ELInletVelocityTable = "<file>"` (quotes optional; the value is taken as a path, not read
list-directed, because a `/` would end a list-directed read). Rows "t U" (seconds, superficial m/s), blank
lines and lines starting with `#` or `!` skipped, blank or comma separated, text after the two numbers
ignored. U(t) linear between rows, constant beyond the ends. Refused with a message and `STOP 1`: missing
or unreadable file, a row without two numbers, a non-finite value, times not strictly increasing, no rows,
and the table together with any of `ELInletVelocityStart / End / RampTime` (key registry of the parser).
Module `EL_INLET_TABLE` (table code, no MPI); `EL_CONFIG:EL_VALIDATE_CONFIG` loads it. Every rank parses
the deck itself, so every rank reads the same file and either holds the same rows or stops with the same
message (agreed form instead of read-and-broadcast; a broadcast would add a collective to the deck parser,
which is also entered on a deck reload) [V, code; unit test of the refusals]. The echo (`EL_PRINT_CONFIG`)
names the file and prints every row as `EL_INLET_TABLE t U`. Without the key `EL_INLET_VELOCITY` evaluates
the ramp expression unchanged [V, unit test 41 and gate E2w5]. The table counts as a configured inlet, so
the superficial-to-interstitial conversion of the inlet value under the continuity stage applies to it
as to the ramp.

Tool: `s11_min_fluidization/tools/make_inlet_table.py` writes the protocol of the spec (U0 0.30, steps of
0.02 m/s every 0.5 s down to 0.04 m/s after 0.5 s of preparation, 20 ms linear transitions, end 7.5 s) and
the step list; its rows equal the `csv` output of the spec's `inlet_table.py` [V, diff].

## 3. Bed audit line

`SimPar@ELBedAuditFreq = N` (default 0: no line, no collective, nothing else changes). Module
`EL_BED_AUDIT`, called by the workers after the fluid solve of every N-th step. Line

    EL_BED_AUDIT t= step= p_in_out= p_in p_out area_in_out= A_in A_out dp_hydro_dyn= dp dp_hydro dp_dyn
                 n= N weight_z_floor_z= W_z F_floor z50_95_99_max= z50 z95 z99 zmax eps_bed= e window= z_lo z_hi
                 gradz_win= G

- Plane pressures: the solver pressure is P1 discontinuous, per element p0 at the vertex mean c and the
  gradient; p(x) = p0 + grad p . (x - c), the form the EL sampler uses. Integrated over every element face
  in the bottom (zmin) and top (zmax) plane of the EL domain box with 2 x 2 Gauss points on the bilinear
  face (exact for this integrand), summed over the ranks, divided by the summed area. No vertex
  projection [V, unit test on the unit-cube mesh: affine field exact, discontinuous field = face-centroid sum].
- Hydrostatic part: `dp_hydro = -rho_f g_z (zmax - zmin)` when the fluid carries gravity
  (`ELFluidGravity = Yes`), else 0. Every momentum model implies `ELFluidGravity = No` (dynamic pressure,
  buoyancy on the particle side), so in case 11 the pressure contains no hydrostatic part and
  dp_dyn = dp [V, code of EL_RESOLVE_MOMENTUM_MODEL; confirmed in the probes: dp_hydro = 0 and dp_dyn
  tracks the particle force]. The 0.98 Pa of the gas column over 100 mm is therefore not in the number.
- Bed weight: sum of (rho_p - rho_f) V_p g_z over the particles (the EL_FORCE_BUDGET grav_buoy z).
- Floor contact: -(contact impulse of the pe bottom z wall in the last macro step)_z / dt.
- Bed heights from a histogram of the particle centres (4000 bins over the column: 0.025 mm at 100 mm),
  linear inside a bin.
- Porosity of the bed window [zmin + 2 d, z99 + r - 2 d] from the positions, exact sphere-slab overlap
  (the spec's definition: bottom and top two diameters excluded) [V, unit test of the slab volume].
- Pressure gradient of the bed window: -<dp/dz> over [z_lo, z_hi], the element-wise constant gradient of
  the P1 pressure weighted with the element volume inside the window (exact for elements bounded by two
  z planes, as here; an estimate otherwise) [V, unit test]. It is the gradient of the same bed core as
  eps_bed and the gate value of F2 in the post tool (commit 87474a41, added after the first fixed-bed
  probe showed that dp over the whole bed carries the looser end layers: 0.918 against 0.962 in the core).
- The forces on the particles (drag, pressure force, buoyant weight) are already in `EL_FORCE_BUDGET`
  (every step, transfer pass with the fluid state before the solve of the step); the post tool joins the
  two lines by step.

Post tool `s11_min_fluidization/tools/fl11_post.py`: per protocol step mean and std of dp_dyn over the last
0.3 s and over 0.2-0.35 s / 0.35-0.5 s (flagged when they differ by more than 2 % of W/A), force-based dp,
bed height, window porosity, floor share; gates F1-F5 of the spec (F2 with the measured dp/dz =
dp_dyn / H_eq, H_eq = N V_p/(A (1 - eps_bed)), against the closure at eps_bed).

## 4. Viscosity convention

`SimPar@ELFluidViscosity` defaults to `dynamic` on this branch [V, el_config.f90]. The deck keeps the EL
convention: `Prop@DynVisc = 1.8d-5,1d0` (the parser divides by the density: nu = 1.8e-5 m2/s for
rho_f 1.0); the closures receive nu and form mu = rho nu, the fluid operator receives mu = rho nu =
1.8e-5 Pa s. No `ELFluidViscosity` key is needed (default). For rho_f = 1.0 the two conventions are
numerically identical. The pe json `fluidViscosity_` holds mu = 1.8e-5 (independent input).

## 5. Verification

- Unit tests (Slurm 163806 at 7f05a069, 163820 at 87474a41; build-fluid11): `el-drag-closures-serial`,
  `el-inlet-table-serial`, `el-bed-audit-serial` and the nine existing EL / pco tests: 12/12 passed.
- Default-path gate (RUNBOOK of the case, section 5): base binary of c76ca927 and the branch binaries
  built the same way on longhouse1a, the same decks without any new key: D_s6 (s10 periodic, model A
  weighted, difelice_std), G_tier2 (stokes), VB40 (difelice), E2w5 (Inflow975 ramp, cylinder): all EL
  streams, solver lines bitwise (jobs 163807, 163822).
- Probes at density ratio 860 (RUNBOOK sections 6-7): stable at dt 2.5e-4 s, 0.55-0.64 s per step for bed R
  and about 1.9 s for bed A on 37 ranks.

Further change found necessary by the probes (commit ab23bb6e): the application closed every box with
`SimPar@PeriodicAxis` (`bNoOutflow` forced to true), so a laterally periodic column with an outlet was
refused by the continuity stage. An explicit `SimPar@NoOutflow = No` is now kept; without the key or with
Yes nothing changes (every committed periodic deck sets Yes; gate D_s6 bitwise).

## 6. Limits and open points

- The three feature commits were built and tested together at 7f05a069; the intermediate commits
  1f736fc9 and 978ed6bf were not built on their own [I: they compile, the split follows file boundaries].
- The audit pressure is the area mean over the boundary planes; the spec's CFDEM functionObject averages the
  patch values as well. The P1 pressure next to the inlet carries the discretisation error of the first
  element layer [I].
- `floor_z` reads the pe wall-contact impulse of the bottom z wall; it is zero when the particles rest on
  each other only and is meaningful only with `zWallsEnabled_` in the json.
- The audit evaluates the window porosity with the mean particle radius (monodisperse here).
