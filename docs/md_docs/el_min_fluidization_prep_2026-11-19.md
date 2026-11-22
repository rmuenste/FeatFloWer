# Preparation of cross-code case 11 (minimum fluidisation) in q2p1_el_pipeflow (2026-11-19)

Branch `feature/el-min-fluidization` (worktree `FeatFloWer-fluid11`, base c76ca927 = master + periodic
kernel + periodic coarse pressure operator), pe pin 43e8fb54 (unchanged). Case spec (peer session,
read-only): `liggghts-dem/LIGGGHTS-PUBLIC/md_docs/case_specs/11_min_fluidization/SPEC.md`. Case directory:
`applications/q2p1_el_pipeflow/validation_cases/s11_min_fluidization/` (RUNBOOK, decks, json, tools). Runs:
`/data/warehouse17/rmuenste/runs/fluid11/<run>`.

Marks: **[V]** verified here (code read, unit test or run), **[I]** inferred, not verified.

## 1. Drag closures: table, Gidaspow, Wen-Yu, Ergun

`SimPar@ELDragModel` now has a table in `EL_CONFIG` (`el_drag_names`, ids `EL_DRAG_*`, column
`el_drag_standard_form`). The validation accepts exactly the names of the table and resolves the id
once (`el_drag_closure_id`, `EL_RESOLVE_DRAG_CLOSURE`; review r1 finding 14: no per-particle name lookup);
the momentum-model refusal of the campaign closure `difelice` reads the column; `EL_FORCES:EL_DRAG_CLOSURE`
dispatches on the id to one function per closure family, with an explicit CASE for every registered
closure (Stokes and Schiller-Naumann included) and a fatal `CASE DEFAULT` (MPI_Abort 915), so a closure
added to the table without its CASE stops the run instead of falling back to the single-sphere law (unit
test 62 of `test_el_drag_closures` evaluates every table entry). A further closure is one table entry,
one function and one CASE line. The existing closures were moved verbatim into `EL_DRAG_B_DIFELICE` and
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
Module `EL_INLET_TABLE`; `EL_CONFIG:EL_VALIDATE_CONFIG` loads it. Review r1 (finding 7) replaced the
agreed form (every rank reading the file twice, list-directed) by: rank 0 of MPI_COMM_WORLD reads ONE
snapshot (one open, every line into memory, then parsed), each row must start with two tokens that are
decimal numbers (strict token check, so `0.5 /` is refused instead of being read list-directed with the
velocity of the row before), the end of the file is distinguished from a read error, a line that does not
fit the 1023-character buffer is refused; status, message and rows are broadcast
(`EL_INLET_TABLE_LOAD_SHARED`, only on the table path: no collective without the key), and every rank
stops with the same message on a refusal. An empty path (`SimPar@ELInletVelocityTable = ""`) is refused.
Parsing, distribution and evaluation are separate routines [V, unit tests 36-46 of
`test_el_inlet_table`]. The deck reload `Reload_SimPar` goes through the generic parser, which knows no EL
key: the table (like every EL setting) is a start-up setting; the echo of the table says so in the
protocol (finding 8; no atomic replacement implemented). The echo (`EL_PRINT_CONFIG`)
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
                 gradz_win= G floor_bal_z= F eps_kern_slip_win= e s win_margin= m gradz_win_dyn= Gd

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
- Porosity of the bed window [zmin + m, z99 + r - m] from the positions, exact sphere-slab overlap
  [V, unit test of the slab volume]. Margin m = MAX(2 d, kernel support radius of the mean particle,
  `EL_KERNEL_SUPPORT_RADIUS` as in the transfer), printed as `win_margin=` (review r1 finding 4, amended
  spec: the window excludes two diameters and the mapping support; 2.5 d for the case-11 kernel). Before
  round 1 the margin was 2 d (logs without `win_margin=`).
- Pressure gradient of the bed window: -<dp/dz> over [z_lo, z_hi], the element-wise constant gradient of
  the P1 pressure weighted with the element volume inside the window, the part inside taken as the
  height fraction of the element: EXACT ONLY for elements bounded by two z planes (the case-11 mesh); for
  other elements, even a tapered hex with horizontal top and bottom faces, an estimate [V, unit test].
  `gradz_win` is the total gradient; `gradz_win_dyn = gradz_win - dp_hydro/(zmax - zmin)` (appended in
  round 1) is the dynamic one, the counterpart of dp_dyn (equal when the fluid carries no gravity, as in
  every momentum model). It is the gradient of the same bed core as
  eps_bed and the gate value of F2 in the post tool (commit 87474a41, added after the first fixed-bed
  probe showed that dp over the whole bed carries the looser end layers: 0.918 against 0.962 in the core).
- Geometry and unavailable values (review r1 finding 12): a positive `ELBedAuditFreq` needs
  `ELDomainType = box` (refused at start-up on every rank otherwise). Unavailable quantities are printed
  as NaN instead of 0, -1 or -HUGE: a plane without element faces (p_in / p_out, dp), no particle
  (heights, window, eps_bed, gradz_win), an empty window, no particle centre in the window
  (eps_kern_slip_win), and `floor_bal_z` outside its assumptions (below).
- `floor_bal_z` (finding 6) is formed only when the CFD side armed pe with the forces
  (`ELApplyParticleForces = Yes`; the expected impulse is then exactly the impulse of the armed drag,
  pressure, lift and buoyant weight); with `No` nothing of the expected impulse acted and the line prints
  NaN (before: an upward "floor force" equal to the computed weight). It is the force of ALL walls, and it
  assumes that pe adds no force of its own (no pe gravity). The momentum balance is kept only with
  `ELBedAuditFreq > 0` (default path stores nothing). Unit tests 80-89 of `test_el_bed_audit` (margin,
  balance incl. the force-off contact-free fixture, kernel-sampled window sums with the bin-centre
  selection at both ends).
- The forces on the particles (drag, pressure force, buoyant weight) are already in `EL_FORCE_BUDGET`
  (every step, transfer pass with the fluid state before the solve of the step); the post tool joins the
  two lines by step.

Post tool `s11_min_fluidization/tools/fl11_post.py`: per protocol step mean and std of dp_dyn over the last
0.3 s and over 0.2-0.35 s / 0.35-0.5 s (flagged when they differ by more than 2 % of W/A), force-based dp,
bed height, window porosity, floor share; gates F1-F5 of the spec. (Wording of 2026-11-19, superseded:
F2 was then dp_dyn / H_eq, H_eq = N V_p/(A (1 - eps_bed)), against the closure at eps_bed. Since the
amendment of 2026-11-20 the tool gates F2 on the dynamic window gradient against the closure at the
mapped gas fraction, and only for a window outside the mapping support; sections 7.2 and RUNBOOK 10-12.)

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

## 7. Velocity spikes of the first protocol run and the HRZ weights of the implicit drag

**Symptom [V].** The first bed-R protocol run (orchestrator, `runs/fluid11/R_proto`, job 163872, binary
87474a41) stopped at step 11 315 (t 2.83 s) with a NaN: the velocity multigrid diverged (rate 14.6, then
5e51) while the pressure solves stayed at 2 cycles. The logged maximum velocity norm was 0.7-0.8 m/s most
of the time but had intermittent episodes from step 2618 on (up to 15 m/s at step 3550, 13-20 m/s in the
last 30 steps). Integral quantities did not follow: dp, superficial flux (EL_VA_FLUXMOM eps u_z / U),
continuity residual and bed heights were normal through every episode; the velocity nonlinear loop went
from 1 to 3 iterations during the episodes. A local defect, not a global balance error.

**Location [V].** Diagnostic line `EL_UMAX_AUDIT` (`SimPar@ELVelocityAuditFreq`, module
EL_VELOCITY_AUDIT, commit d3d23e4b) in a rerun of the first 4 800 steps (`R_diag`, job 163992, the EL
streams bitwise with R_proto: same trajectory). Every maximum above 1.2 m/s sits at a Q2 vertex dof at
z = 30.0 or 32.5 mm, i.e. at the bed surface (z99 about 29-30 mm), where the transport porosity is
0.77-1.0 (a few particles); 89.3 % of them (1 674 of 1 875; corrected in review r1, first written
as 95 %) at a dof whose logged implicit drag coefficient B_i is NEGATIVE (counts in the RUNBOOK;
the logged field is the post-advance diagnostic value of the owner rank, not the assembled
coefficient the solver used, see 7.1). The velocity there grows by 1.5-3.6 % per step, which is
thstep |B_i| divided by the lumped (eps-weighted) fluid mass of a vertex dof.

**Cause.** The semi-implicit drag adds thstep B_i to the diagonal of the fluid matrix and +B u_f to the
right-hand side; B_i = sum_p B_p int K_p phi_i / norm_p with the Q2 basis functions phi_i. A Q2 vertex
function is negative over part of its elements, so a node reached only by those parts of a particle
kernel gets B_i < 0: anti-damping. In a liquid the lumped fluid mass is orders larger than thstep |B_i|
and the term is harmless (all earlier cases); in the gas of case 11 (rho_f 1, rho_p/rho_f 860) it is a
few percent of the mass per step, and a negative B_i of growing size ends in a non-positive diagonal.

Hypotheses of the orchestrator checked with the same data:
- (a) inlet: the largest inlet value over the 4 800 steps is 0.76 m/s, the smallest inlet transport
  porosity 0.396; no maximum on the inlet plane (z <= 1 mm: 0 of 4 800 audit lines). Excluded.
- (b) low porosity: smallest element porosity over the run 0.359 (floor 0.2, no clipping); the spikes
  sit where eps is 0.77-1.0. Excluded as the cause.
- (c) advective CFL: before an episode the velocity is 0.7-0.8 m/s and the nonlinear loop converges in
  one iteration; the growth is exponential at the rate thstep |B_i| / m_i, not a CFL breakdown (the CFL
  grows only as a consequence). Not the cause.
- (d) outlet: no maximum near the outlet (z >= 95 mm: 0). Excluded.
- (e) stiffness of the semi-implicit drag: closest, but the issue is the sign of the nodal coefficient,
  not its size.

**Remedy (commit 77ffd070).** `SimPar@ELDragImplicitWeights = rowsum | hrz` (default rowsum: every
existing log bitwise). Under `hrz` the implicit coefficient B and its explicit split partner -B u_f are
deposited with the diagonal-scaling (Hinton-Rock-Zienkiewicz) weights of each element,
w_i = int K phi_i^2 * (sum_j int K phi_j) / (sum_j int K phi_j^2) / norm >= 0, which keep the element
total of the row-sum weights (so sum_i B_i = B_p, and the split -B (u_i - u_f) is exact for a uniform
field); the drag impulse the particle receives and every other reaction keep the row-sum weights. This is
the standard remedy for higher-order elements whose row-sum lumping is not positive; it changes only how
the linearisation correction of the drag is distributed inside an element. Unit test
(`test_el_transfer`): a small kernel gives a row-sum minimum of -4.5 % of B, the HRZ weights are
non-negative with the same total, the split partner is deposited with them. The case deck sets
`SimPar@ELDragImplicitWeights = hrz`.

Verification (job 163999, `runs/fluid11/R_hrz`): the first 1.3 s of the protocol with `hrz` ran without
a spike: largest velocity norm 1.05 m/s over all 5 200 steps (start-up collapse), one nonlinear
iteration per step, plateau 0.97 W/A with 2.3 % of the weight on the floor. Default-path gate of the
binary bitwise (D_s6, G_tier2, E2w5). Not verified: the full protocol (job 164004 running at the time of
writing), and whether `hrz` changes the fixed-bed gradient of the probes (the row-sum and HRZ weights
distribute only the implicit correction B (u_i - u_f), which vanishes in a converged steady state with
u_i = u_f; in a fixed bed u_i differs from the sampled u_f at element scale, so a small change is
possible) [I].

**General property (2026-11-20).** The negative row-sum weight is not specific to case 11: every
semi-implicit Euler-Lagrange run deposits B_i = sum_p B_p int K_p phi_i / norm_p with the Q2 basis, and a
node reached only by the negative lobes of phi_i gets B_i < 0 (unit test: -4.5 % of B_p for a small
kernel). In a liquid the anti-damping thstep |B_i| is negligible against the lumped fluid mass, so the
runs are stable, but the nodal weights are formally wrong there too (an anti-damping coefficient where
the physical drag damps). `hrz` keeps the element totals and changes only the distribution of the
implicit correction B (u_i - u_f) inside an element. Whether it should become the default is an owner
decision. Its effect on a liquid gate [V, job 164099, `runs/fluid11/hrzliquid`]: the D_s6 deck (s10
periodic settling bed, rho_f 1000, mu 0.05, model A weighted, 18 335 particles) run for 200 steps (0.1 s)
with rowsum and with hrz: the vertical drag per step (`EL_FORCE_BUDGET`, steps 1-200) agrees to 1.8e-4
relative (largest at step 193; mean difference 4.9e-5; it is fixed by the bed weight), the mean
particle velocity is 0.08 % lower with hrz at step 20, 0.24 % at 60, 0.33 % at 100 and 0.42 % at 200
(slip alike), i.e. a small change that has not levelled off after 0.1 s; the long-run effect on the
s10 numbers (core 1.047) is not measured.

### 7.1 Momentum consistency of the hrz split (review r1, finding 2, 2026-11-21)

**Defect [V, unit test, job 164663].** In the version above the explicit partner of the implicit
drag was -B u_f, deposited with the HRZ weights h_i, while the particle samples
u_f = sum_i w_i u_i^n with the row-sum weights w_i. With the theta scheme of `QuadSc_def`
(matrix diagonal theta dt B_i, old-velocity defect -(1-theta) dt B_i u_i^n, B_i = B h_i) the
drag terms of the fluid momentum balance of one step were

    dt [ -w_i F_d + B h_i u_f - B h_i (theta u_i^{n+1} + (1-theta) u_i^n) ],

and in a field that does not change over the step (u^{n+1} = u^n) their sum is
-F_d + B (u_f - sum_i h_i u_i): a net fluid force without a particle counterpart whenever the
field is not uniform, because sum_i h_i u_i != sum_i w_i u_i. Unit test
`test_el_transfer` (steady non-uniform Q2 field, small kernel): residual -1.85 % of F_z.

**Remedy.** Under `hrz` the partner is formed by the fluid from the OLD NODAL velocity with the
same diagonal, `EL_FIELDS:EL_APPLY_DRAG_SPLIT_PARTNER` (called in `QuadSc_main` next to the
explicit source, before the boundary values of the new step are imposed): + dt B_i u_i^n on
every dof, B_i the captured coefficient `drag_B_source` (distributed form, like the explicit
source). The deposit no longer carries a partner (`EL_DEPOSIT_PARTICLE` argument `hrz` instead
of `split`; the transfer result rows are 5 again). The drag matrix acts in the MOMENTUM SOLVE,
i.e. on the momentum predictor u~ (the velocity before the pressure projection and the velocity
correction that follow it; review r2, P2); the discrete momentum equation per dof i is

    (M_i + theta dt (A + B_i)) u~_i = (M_i - (1-theta) dt A) u_i^n + dt (-w_i F_d + f_i^other)
                                      + theta dt B_i u_i^n + (pressure gradient, body forces),

i.e. the drag terms of the momentum solve are dt [ -w_i F_d - theta B_i (u~_i - u_i^n) ] with
B_i = sum_p B_p h_ip >= 0 (A the remaining operator, F_d the drag impulse per dt that pe applies
over its sub-cycled semi-implicit steps with the frozen u_f, f^other the non-drag reactions; case 11
runs BE, theta = 1: -dt B_i (u~_i - u_i^n)). The projection then changes u~ to u^{n+1} without a
drag term.

What is and is not conserved:
- u~ = u^n (a field the momentum solve leaves unchanged), the configuration the unit tests exercise:
  they evaluate the deposit, the partner and the diagonal for a GIVEN steady field (no solve, no
  projection). There the split correction vanishes node by node (not only in the sum) and the fluid
  receives exactly -F_d dt, the mirror of the particle drag impulse [V, unit tests
  `test_el_transfer` on the regular mesh with a small and a wide kernel and on a distorted mesh,
  `test_el_periodic_mpi` on 4 and 8 ranks with the periodic seam]. This is the whole conservation
  claim; it is not a statement about the steps of a run.
- In a run the momentum solve receives in addition the impulse -theta dt sum_i B_i (u~_i - u_i^n),
  which the particle does not see (pe integrates with the frozen u_f of the start of the step). It
  is a non-negative nodal damping of the predictor increment (B_i >= 0), the stabilising part of the
  linearisation; it is NOT momentum exchange and is not mirrored. It is not measured: the velocity
  audit and the output fields carry the projected u^{n+1}, and two equal successive PROJECTED
  velocities do not show that u~ = u^n within the step. A stepwise statement would need the term
  evaluated at the predictor, before the projection (open [I]).
- The default `rowsum` has the corresponding term B (u_f - sum_i w_i u~_i^theta), which vanishes in
  the sum for a steady field (u_f is sampled with the same w_i) but not node by node, and whose nodal
  coefficient can be negative.
- `rowsum` (default) is unchanged bit for bit [V, gates D_s6, G_tier2, E2w5, VB40 of the round-1
  binary 627caf47, job 164668, `s11_min_fluidization/runs_versions.txt`].

**Consequence observed in the smoke runs [V, jobs 164669-164673, `runs/fluid11/R_smoke_r1*`].** The
former partner -B h_i u_f together with the diagonal acted as a relaxation B h_i (u_f - u_i^{n+1}) of
every nodal velocity towards the kernel average u_f: besides the net force above it damped the nodal
velocity modes that the kernel average does not see (sum_i w_i u_i = 0), which feel no drag in the
Euler-Lagrange coupling. The consistent split removes that damping (only the increment is damped).
Case 11, bed R, constant U 0.30 from the rsa bed, same deck and json (friction 0.15 per material):

| run | binary | max \|u\| steps 1-400 | max \|u\| steps 1000-2000 | dp* at 400 / 1000-2000 |
|---|---|---|---|---|
| R_smoke_r1(_2000), hrz consistent | 627caf47 | 2.21 m/s (step 390, bed top) | 1.57 m/s | 1.364 / 0.96-1.02 |
| R_smoke_r1_oldhrz, former hrz | c8541cfc | 1.00 m/s | (R_hrz: 0.87) | 1.366 / (R_hrz 0.96-1.00) |
| R_smoke_r1_rowsum | 627caf47 | 1.09 m/s | - | 1.364 / - |

No exponential growth, no negative coefficient at the maxima (dragB_fsrc_solver > 0), the integral
quantities agree to 0.2 %; but the largest nodal gas speed is about 1.6-2.2 times that of the former
hrz and of rowsum. Whether these sub-kernel velocity modes need their own damping (e.g. a drag form
that samples the implicit part with the same nodal weights for the particle, so that both sides see
B h_i (u_i - v)) is open [I].

Changed diagnostics under `hrz`: `EL_FLUID_SOURCE` / `EL_FEEDBACK_CONSERVATION` report the
deposited explicit reaction only (the partner is no longer deposited), and `EL_FLUID_PAIR`
compares the fluid momentum change with dt times the deposited source; the partner
+ theta-part of the diagonal belong to the implicit side there, as for `rowsum`.

**Velocity audit (finding 10).** `EL_UMAX_AUDIT` appends `dragB_fsrc_solver= Bs fs`: the drag
coefficient and the z component of the explicit source captured in the PRE-advance pass (the
values the fluid solve of the step used), summed over all ranks that share the dof (one
`E013Sum2` on copies, entered by every worker on every audit step). The older field
`dragB_fsrc=` keeps its meaning (post-advance diagnostic pass, owner rank only) for the
comparison with older logs; the counts in the RUNBOOK of negative coefficients at the spike
locations (1 674 / 1 875 = 89.3 % of the logged maxima, not 95 % as written above) refer to that
older field.

### 7.2 Fix round 2 (review r2, 2026-11-22)

Implementation of the external review r2 (findings 1-10 of the brief; record files RUNBOOK,
example json and datasheet are kept by the campaign record, not here). Fortran commits c0d71b65,
95e5923a, 729f77c9, 3e801eab, d27a74fb; binary `q2p1_el_pipeflow.d27a74fb-pece2a7f7` (build job
164949); post tool 1fbb28e4, movie job 2cfeb048. Every claim was checked in the code first.

**Post tool `fl11_post.py` (findings 1-3) [V, job 164960].**
- F2 is gated only when every contributing audit window excludes max(2 d, mapping support) at both
  ends. The support is the kernel width factor of the start-up echo times d (EL_CONFIG:
  EL_KERNEL_SUPPORT_RADIUS = 2 x factor x r for every kernel); the margin is measured on every audit
  line of the fixed-step windows: lower end minus the inlet plane (`Global CFD domain extents`),
  z99 + r minus the upper end, and `win_margin` (2 d for binaries before 627caf47), within the print
  precision of the line. Otherwise the ratios are printed as a measurement and the verdict is
  `unavailable (window inside the mapping support)`.
- F2 uses `gradz_win_dyn`; a log without it uses `gradz_win` only when `dp_hydro` = 0 on every audit
  line used (dp_hydro = -rho_f g_z (zmax - zmin) is zero exactly when the fluid carries no gravity or
  rho_f g_z = 0, i.e. no hydrostatic gradient), else unavailable. All three section-11 runs: zero.
- Every ratio needs a positive finite denominator (no inlet plane, empty bed: nan, `unavailable`, no
  exception); F1-F3 and F5 are evaluated independently; a step is complete when the audit samples its
  window densely (interval at most a tenth of the window, no gap wider than one interval); a steps
  file is compared with the echoed protocol (missing and foreign steps named, foreign ones excluded)
  and cannot move the F5 target, which is the lowest plateau of the ECHO; F5 runs over the constant
  part of that plateau after an audit-coverage check, the 20 ms transition into it is reported as a
  labelled stricter number. `--self-test`: 12 synthetic logs (reference, missing inlet plane, empty
  bed, incomplete and foreign steps file, short window, sparse audit, margin and gradient rules).

| run (`post_r2.txt`) | F1 | F3 ratio | F2 worst \|ratio_k - 1\| | window margin | F2 verdict | F5 (plateau / with transition) |
|---|---|---|---|---|---|---|
| R_r1_proto | 0.978732 FAIL | 1.0474 PASS | 3.49 % | 2.5 d | PASS | 23.81 / 23.81 mm PASS |
| R_r1_oldsplit | 0.982573 PASS | 1.0490 PASS | 2.38 % (measurement) | 2.0 d < 2.5 d | unavailable | 24.00 / 24.00 mm PASS |
| R_proto_movie | 0.981421 PASS | 1.0420 PASS | 1.62 % (measurement) | 2.0 d < 2.5 d | unavailable | 24.03 / 24.03 mm PASS |

The two historical 2 d windows lie inside the 2.5 d support: their F2 numbers are measurements, not
gate results. The F2 comparison between the two splits (3.49 % vs 2.38 %) also mixes two windows.

**`ELDragImplicitWeights = hrz` and the non-Newtonian defect (finding 4) [V, job 164952].** Verified:
with ELMomentumModel = none, ELFluidViscosity = kinematic, semi-implicit feedback, FlowType =
non-Newtonian and MatrixRenewal S = 0 every former check passed, the partner + dt B_i u_i^n went to
the right-hand side, and the hand-built defect of `Matdef` (idef = -1: mass, convection, STRESS) has
no -theta dt B_i u_i, so the partner stays in the residual of an unchanged field (the assembled-matrix
defects of the Newtonian path and of non-Newtonian S >= 1 carry B on the diagonal). That combination
is now refused at start-up (`EL_CONFIG:EL_DRAG_HRZ_DEFECT_SUPPORTED`, message naming FlowType,
MatrixRenewal and ELDragImplicitWeights) and the hrz branch of the momentum assembly stops on it (a
FlowType or MatrixRenewal changed later by Reload_SimPar, whose generic parser reads both). Old binary:
accepted, 3 steps; new: refused before the first step. Follow-up [I]: the missing defect term itself
(no test case); the same branch also lacks the theta dt B_i u_i of the default `rowsum` semi-implicit
drag, which this round does not refuse (pre-existing, outside the brief).

**Velocity audit (finding 5) [V, job 164952].** The transfer captures `drag_B_source` and
`fluid_feedback_source` in every configuration; the solve applies B only under semi-implicit fluid
feedback, the source only with fluid feedback, neither for a prescribed field. `EL_AUDIT_SOLVE_APPLIES`
mirrors these guards and the audit prints NaN for a term the solve did not apply (collectives
unchanged). Old binary, explicit drag with feedback (G_tier2): Bs 1.81e-3 printed; new: NaN, fs
unchanged. Feedback off (VB40): old Bs 2.0e-6, fs -6.0e-7; new NaN NaN. E013Sum2 sums for the
communication schemes 3 and 4 only: ELVelocityAuditFreq > 0 with another scheme is refused at
start-up (the pattern of the other E013Sum2 consumers; iCommSwitch has no deck key, 3 is the default,
so the case is not reachable from a deck).

**Inlet definition at a reload (finding 6) [V, job 164952; unit test checks 48-51].** Start-up only
stays the contract: Reload_SimPar runs the generic parser, which reads no SimPar@EL* key. The reload
now reads the inlet keys of the reload file the way the application parser reads them and reports a
difference (changed table path, table removed, table <-> ramp, ramp values changed or removed, both
together, unreadable value) as a WARNING on the protocol and the terminal; the run keeps its
definition. A ramp-only deck gets its start-up notice (`EL inlet velocity ramp: start-up setting ...`,
not an `EL_` line). E2w5 with a reload deck that swaps the ramp for a table: old binary silent, new
`WARNING: EL inlet definition of the reload file _data/reload_table.dat ignored: ramp -> table ...`.

**Inlet reader (finding 7) [V, failing before: job 164948, STOP 47].** Tabs become blanks before the
blank-line and comment tests; tab-only and tab-indented comment lines are skipped, data rows
tokenise as before.

**Movie render job (finding 8) [V, jobs 164957-164959].** Fresh frame directory per render, each
worker waited for, a failed worker stops the encode, the movie is replaced only on success.

**Design note (finding 9).** Section 7.1 now writes the drag terms with the momentum predictor u~ and
restricts the conservation claim to the given-field configuration of the unit tests.

**P3.** `el_bed_audit.f90`: the window-gradient weighting is exact for constant horizontal
cross-section, stated at the integration routine; one closure module `fl11_closures.py` for
`fl11_post.py` and `core_deficit.py` (identical values to both former copies, job 164960); the
unit-test job prints the checkout HEAD and the tree hashes of `source/` and of the application.

**Gates [V, job 164951; `runs/fluid11/gate/pc_compare_d27a74fb-pece2a7f7.txt`].** (a) Default path: D_s6,
G_tier2, E2w5, VB40 against the `__base` runs of job 163807 and the first 400 steps of case 11 with
the restitution keys absent (`R_fix2_off400`, inputs of R_smoke_r1) against R_smoke_r1 (627caf47):
all common EL_ streams bitwise (strict), XxX 40/60/10/40/800 bitwise, case-11 particle logs 27/27
identical; every `EL_` line as text bitwise against the previous binary 20867701-pece2a7f7 (against
the older references one extra line, the start-up `EL_CONTACT_PARAMS` echo of pe ce2a7f7). (b) Case 11
with restitution on (`R_fix2_id400`, deck and json of R_id400_new, 400 steps) against R_id400_new:
every `EL_` line bitwise (13 077, among them EL_RESTITUTION_STATS 400, EL_BED_AUDIT 40, EL_UMAX_AUDIT
40), XxX 800/800, particle logs 27/27 identical. The record run (Fortran 627caf47) therefore stands
for the final head on the case-11 path.

### 7.3 Fix round 3 (review r3, 2026-11-22)

Commits 46d06ce6 (Fortran, finding 8), 70306954 and a45e9303 (post tool, findings 1-7), b9e8f410
(movie tool and render job, findings 9-10); binary `q2p1_el_pipeflow.b9e8f410-pece2a7f7` (build job
164964, unit tests 164965: 15/15). Every reviewer counterexample was reproduced first.

**Post tool `fl11_post.py` [V, jobs 164968, 164972; `runs/fluid11/unit/post_r3/old_vs_new.txt`].**
- The official gates F1-F5 use the ECHOED protocol only, each plateau once. A `--steps` file is
  validated in all three columns and for multiplicity; a deviating file gives a separately labelled
  steps-file summary (NOT gated). Only a log without a table echo is evaluated from a steps file, and
  every verdict then carries `[protocol from the steps file, not from the echo]`. On the unchanged
  R_r1_proto log the tool of 1fbb28e4 printed F1 0.980867 PASS for a duplicated 0.18 row and 0.980269
  PASS for an omitted 0.30 row and accepted a 0.18 row with a shifted start silently; the new tool
  prints the recorded 0.978732 FAIL in all three and names the deviation.
- F2: each contributing mapped sample must be finite with 0 < eps < 1 (the -1 sentinel of old
  binaries and an unclipped 1 are refused), and every fixed step needs a positive closure and a finite
  ratio before the maximum (zero-velocity plateau included). The old tool printed `worst 0.0000 PASS`
  for a mapped porosity 1.0; now unavailable with the reason.
- Inlet history with the solver rule (EL_INLET_TABLE_VALUE): the last value is held up to the
  evaluation horizon min(MaxSimTime, start + MaxNumStep dt) (deck values when the log's time line
  agrees within its 4-digit print precision, the log otherwise); without a horizon F5 is withheld for
  that hold. A table without its redundant last knot: old F5 tested 0.10 m/s (PASS), new 0.05 m/s (FAIL
  with the 60 mm centre of the fixture).
- A missing floor_bal_z is nan for any number of audit lines (old: TypeError on a one-line log); the
  EL_NEWTON_PAIR reconstruction needs the matching line and a positive finite dt.
- F2 margin: the inlet elevation and a finite lower, upper and reported-or-legacy margin on every
  contributing line are required; missing geometry is its own reason (old: `met`, PASS).
- 'The run reached t' (an audit line or the step count of the last time line of the log) is checked
  separately from the sampling (distinct samples, first / last within one cadence of the window ends,
  no gap wider than one cadence, print-precision tolerance of the time). A uniform 30 ms audit now
  qualifies (old: no window), 40 ms gaps in a 30 ms audit do not, a dropped last audit line of a run
  that reached its end gives a consistent complete result. F5 qualifies the particle count over the
  whole run (start to horizon) separately from the height coverage of the lowest plateau; a log with
  only the final rows gives `unavailable (count over the run not established)` (old: PASS).
- Absolute z: bed height z99 + r - zmin, containment against zmin + Lz/2 from the echoed extents
  (unavailable without them; `--lz` only cross-checked). A column translated by -20 mm: old 40 mm
  PASS, new 60 mm above the floor, FAIL.
- `--self-test`: 30 synthetic logs, every counterexample of reviews r2 and r3 among them;
  `old_vs_new.txt` runs each through the tool of 1fbb28e4 as well.

`post_r3.txt` of R_r1_proto, R_r1_oldsplit and R_proto_movie: every gate number identical to
`post_r2.txt` (F1 0.978732 FAIL / 0.982573 / 0.981421, F3 1.0474 / 1.0490 / 1.0420, F2 3.49 % PASS /
unavailable / unavailable, F5 23.81 / 24.00 / 24.03 mm PASS with the count covered over 0-7.5 s).

**HRZ refusal (finding 8) [V, job 164967].** `EL_CONFIG:EL_FLUID_SOLVE_ACTIVE` is the one fluid-solve
predicate (transport, velocity-audit terms, start-up refusal). The prescribed-field deck (fb4 +
`ELPrescribedField = linear_shear`): refused by d27a74fb, runs 3 steps with the new binary; the active
non-Newtonian S0 solve is still refused, the message no longer suggests rowsum. The post-reload
backstop aborts MPI_COMM_WORLD (EL_DRAG_HRZ_RUNTIME_CHECK): a Newtonian hrz run whose Reload_SimPar
switches FlowType to non-Newtonian ends after step 1 on every rank with MPI_ABORT, exit 149 (d27a74fb:
ERROR STOP on the workers only; in this run the launcher also ended the job, exit 1, but the new path
does not rely on that).

**Movie (findings 9-10).** `fl11_movie.py` takes its step table and audit columns from
`fl11_post.prepare` / `step_table`; the render job keeps the check log and the temporary movie in its
job directory and exits with the status of the final rename [V, job 164969].

**Gates [V, job 164966; `runs/fluid11/gate/pc_compare_b9e8f410-pece2a7f7.txt`].** (a) Default path: D_s6, G_tier2, E2w5,
VB40 against `__base` and `R_fix3_off400` (restitution keys absent) against R_smoke_r1: all common EL_
streams bitwise (strict), XxX 40/60/10/40/800, particle logs 27/27; every `EL_` line bitwise against
the round-2 binary d27a74fb-pece2a7f7 (against the older references one extra line, the
`EL_CONTACT_PARAMS` start-up echo of pe ce2a7f7). (b) `R_fix3_id400` (restitution on, 400 steps)
against R_id400_new: every `EL_` line bitwise (13 077; EL_RESTITUTION_STATS 400, EL_BED_AUDIT 40,
EL_UMAX_AUDIT 40), XxX 800/800, particle logs 27/27.

### 7.4 Structure of the case-11 post tool (review r4, 2026-11-22)

The fourth review showed verdict reversals on equivalent inputs (a redundant table knot moved F1 from
0.978732 FAIL to PASS) and on degraded inputs (one dropped audit row: PASS from the remaining subset).
The cause was structural: holds, coverage, statistics, branch membership and extrema came from
different views of the input. `fl11_post.py` (commit d70850f8; `fl11_closures.py` unchanged) is now built
from five pieces, each defined once; `fl11_movie.py` (616696cc) uses the same evaluation.

**A. Record reader.** One tokenizer reads every `EL_` record the tool uses (EL_BED_AUDIT,
EL_FORCE_BUDGET, EL_NEWTON_PAIR as `key= values` records with the required and optional keys and their
value counts; EL_INLET_TABLE as exactly two numbers). A record is complete or it is nothing; no field is
taken by position. A malformed record after the last complete record, or a last line without its line
end, is an unfinished trailing record: ignored, noted on stderr, so the report is unchanged. A malformed
interior record is diagnosed in the report and leaves a hole. The same pass reads the start-up echo
(kernel width factor, drag model, fluid gravity, bed audit frequency, the row count of the table echo,
domain extents), the time lines and the finish marker. The deck is read with the application's syntax
(last `@`, last `=`, first list-directed token of the value, so a trailing comment is not part of it);
a missing key and an unreadable value are told apart.

**B. Protocol.** The solver's rule U(t) (EL_INLET_TABLE_VALUE: linear between the rows, the first and
the last value held beyond the ends) is applied to the echoed table (complete against its echoed row
count) and reduced to maximal constant-velocity holds: adjacent equal-valued segments, the holds beyond
both ends included, are merged before the 0.1 s threshold; the terminal hold ends with the run (open
when the run has not finished). A step ends where its hold ends and starts where the previous hold ends
(the start of the transition into it). Redundant knots, subdivided holds and a dropped redundant last
knot give byte-identical output. Without an echo a valid steps file (no duplicate or overlapping rows)
is the protocol, every verdict labelled; with an echo a steps file is only validated against it (three
columns, multiplicity) and gives a labelled summary when it differs.

**C. Clock and audit sequence.** The expected audit rows are the multiples of the echoed audit frequency
up to the last completed step; their times are t0 + k dt where a constant time step is established (one
printed dt on every time line; the exact deck TimeStep when it agrees within the 4-digit precision of the
time line, else the audit rows). An observed row is validated against its step (multiple of the
frequency, time within the print precision of t0 + k dt, particle count and buoyant weight changing
together); identical copies are collapsed, conflicting copies rejected with the reason, a missing or
rejected expected row is a coverage hole. The run ends with the last completed step (the solver finishes
the step that crosses MaxSimTime; a Reload_SimPar that moves MaxSimTime is followed through the time
lines and the observed steps); without the finish marker the run has not finished.

**D. Per-step evaluation.** The evaluation window (the last 0.3 s of the step) must lie inside the hold
and be covered by validated rows (at least ten expected rows); statistics, sub-windows and extrema use
the same rows, selected by their step time. Buoyant weight, inlet area and particle count must be
constant in the window: they normalise that step (no run-wide normalisation). A step is `measured`, or
`unmeasurable` with the reason (short hold, coverage hole, non-finite dp_dyn, non-constant weight / area
/ count, open terminal hold).

**E. Verdicts.** Every gate is PASS, FAIL or unavailable with the reason. F1-F3 are unavailable when a
protocol step is unmeasurable (it could belong to the fluidised or the fixed set); the subset value is
printed as a measurement. F2 and F3 use the closure of the drag-model echo; a differing `--closure` is a
labelled comparison. F2 needs a finite positive mapping support, finite margins on every contributing
line (tolerance from the print precision of both coordinates, undecidable otherwise), valid mapped samples
and a positive closure with a finite ratio on every fixed step. F5 count (every expected row validated
over the whole finished run, one count) and F5 height (all tied lowest holds covered, no centre above
zmin + Lz/2; a centre at mid-height is not above) are independent of F1-F3 and of each other; an observed
exceedance or count change is a FAIL even with holes. Physical command-line arguments must be finite and
positive. The report ends with one `verdicts:` line.

Evidence: [V, job 164976; `runs/fluid11/unit/post_r4/`] `--self-test` 52/52 (every counterexample
of reviews r2-r4; `old_vs_new4.txt` runs the r4 cases and two real-log variants through the tool of
9984944d as well: the redundant knot (3.80, 0.18) gave F1 0.980674 PASS and the dropped 0.8500 s row
0.980269 PASS there, now 0.978732 FAIL and unavailable). `--self-test-real R_r1_proto`, 300 seeded
variants: 100 equivalence transforms (knot in a hold or a ramp, subdivided hold, last knot removed,
unfinished trailing record) with byte-identical output, 200 degrading transforms (rows deleted,
duplicated, conflicting, a field NaN, log truncated, audit thinned, count or weight perturbed, finish
marker removed) without a flipped verdict, 170 of them with at least one verdict unavailable.
`post_r4.txt` of R_r1_proto, R_r1_oldsplit, R_proto_movie and R_r2_rest: every gate number and verdict
of `post_r3.txt`. `fl11_movie.py --stats` on R_id400_new: both round-3 caches rebuilt, a fresh cache reused.

## 8. Floor force of the bed audit: the pe wall-contact impulse is the unrelaxed PGS impulse

**Finding [V].** In the fixed branch of the protocol run (R_proto_hrz, job 164004) dp* + floor/W with the
pe wall-contact impulse (`getElWallContactImpulse`, audit field floor_z) was 1.016 at 0.16 m/s rising to
1.089 at 0.04 m/s, while a static bed must give 1. The ratio floor_pe / (W - F) is 10/9 = 1.111 in every
fixed step (1.1111-1.1113) and also in the fluidised steps (1.108-1.111): a constant factor, not a
geometric effect. The EL_NEWTON_PAIR audit shows the same: the particle momentum is constant (dp ~ 0)
while the wall impulse exceeds the impulse of the applied forces by 11.1 % (its mismatch column).

**Cause [V, pe source read].** `HardContactEulerLagrange` folds every hard-contact velocity correction
into the bodies with the uniform under-relaxation `relaxationParam_ = 0.9` (v += 0.9 dv; the comments
there explain why both sides of a cross-rank contact must get the same fold). The converged contact
impulses p_[i] therefore over-state the impulse actually applied to the bodies by 1/0.9. The wall
channel (`wallContactImpulse_`, accumulated from p_[i] after the PGS loop) and the contact virial
(`contactVirial_`, same loop, EL_CONT_STRESS) both report p, not 0.9 p. The dynamics are correct (the
bed is static and its momentum is conserved); only the two diagnostic channels are scaled. An earlier
suspicion (shadow copies at the periodic planes counted twice) does not fit: the factor is the same in
fluidised and fixed states and equals 1/0.9 to four digits; contacts of remote bodies with global
bodies are rejected by the contact filter.

**Needed in pe (not changed here, pin 43e8fb54 kept):** accumulate relaxationParam_ * p_[i] (the applied
impulse) in wallContactImpulse_ and contactVirial_, or document both channels as unrelaxed. Affected
FeatFloWer diagnostics: EL_NEWTON_PAIR (its mismatch carries 0.11 x the wall impulse whenever particles
rest on a z wall), EL_WALL_STRESS tauC, EL_CONT_STRESS (contact stress, Kroupa campaign) - all scaled by
1/0.9 where contacts carry load [I for the last two: same accumulation].

**FeatFloWer remedy (commit c8541cfc).** The bed audit line gets `floor_bal_z = (dp - expected)_z / dt`,
the wall force from the momentum balance of the particles over the PE advance (EL_NEWTON_PAIR_END keeps
dp and the applied-force impulse every step). It is the floor force when the floor is the only wall the
particles touch (top wall far above the bed in case 11). Over the existing protocol log (same two
numbers from the EL_NEWTON_PAIR lines): dp* + floor/W = 1.0003-1.0016 in the fixed branch, 0.993-1.002
in the fluidised steps; the fluidised floor shares become 1.2-2.0 % (were 1.3-2.3 %). The post tool
uses floor_bal (from the audit line, or from EL_NEWTON_PAIR for older logs).

Verification [V]: `runs/fluid11/R_static` (job 164100, binary c8541cfc): the settled bed of R_proto_hrz
at t 7.5 s (positions from the particle logs; json `seedAllowContact_ true`, a first attempt in job
164099 was refused by the seed-gap check) at a constant 0.04 m/s for 400 steps: dp* + floor_bal/W =
1.0002-1.0011, floor_pe / floor_bal = 1.111 throughout, dp* 0.2013-0.2014.

## 9. The 5 % deficit of the core gradient against the closure

Measured (R_proto_hrz, fixed branch 0.04-0.14 m/s): -<dp/dz> over the bed window [2 d, H - 2 d] is
0.946-0.954 of Gidaspow (= Ergun) at the exact window porosity 0.369. Decomposition at four output frames
(0.04, 0.06, 0.10, 0.14 m/s; `s11_min_fluidization/tools/core_deficit.py`: the code's kernel sample replicated
from the VTK fields - deen_poly kernel, 27-point Gauss quadrature per element, periodic images, P0
porosity exact, velocity trilinear from the vertices; replicated total drag within 0.2 % of
EL_FORCE_BUDGET at all four frames) [V]:

| factor | 0.04 | 0.06 | 0.10 | 0.14 |
|---|---|---|---|---|
| (a) kernel-sampled porosity at the window particles: mean 0.3722-0.3731, +0.0035 above the slab value -> drag | 0.9855 | 0.9857 | 0.9862 | 0.9866 |
| (b) nonlinearity: mean of drag(eps_i) against drag(mean eps) | 1.0020 | 1.0020 | 1.0020 | 1.0019 |
| (c) sampled slip u_z 2.3-2.4 % below U/eps (particles sample the slower flow next to them) | 0.9779 | 0.9775 | 0.9770 | 0.9765 |
| (d+e) rest: model-A identity n<F_d>/eps against the measured window gradient | 0.9883 | 0.9843 | 0.9813 | 0.9801 |
| product = measured / closure | 0.9544 | 0.9504 | 0.9474 | 0.9461 |

Of the rest, the porosity the fluid momentum balance carries (deposited element porosity 0.3698-0.3705,
+0.001 above the slab value) accounts for 0.3 %; the remainder (0.9-1.7 %, growing with U) is not split
further (pressure force at the particle positions against the window mean, viscous term of model A,
P1 gradient measure) [I].

**Where it sits [V].** Restricted to the interior of the window (more than one kernel radius, 2.5 d,
from both window ends), the sampled porosity is 0.0004 BELOW the exact value (the self-voidage sense, a
+0.3 % drag effect) and the sampled slip is 0.996 of U/eps. A gradient from the vertex pressure between
the vertex planes 5 and 20 mm is 1.014-1.024 of the closure at the exact porosity of that slab, between
2.5 and 22.5 mm 0.956-0.966. Factors (a) and (c) come from the particles within one kernel radius of the
window ends: their kernels reach into the loose floor layer (and are cut by the floor) and into the bed
top, which the 2 d exclusion of the spec does not remove for a 2.5 d kernel. The self-voidage bias of
this code (eps_eff about 0.975 for a single particle) does not apply in a packed bed: there the
particle's own volume is part of a field that is already dense, and the interior bias is -0.0004.

**Like-for-like reference.** The closure evaluated at the porosity the drag actually uses: the mean
kernel-sampled porosity of the window particles, now in the audit line (eps_kern_slip_win, commit
c8541cfc) and in the post tool as ratio_k. It removes factor (a) only; (c) is a property of the
resolved flow (flow through the looser regions), (b) is +0.2 %. In R_static (0.04 m/s, settled bed) the
code's own window mean of the sampled porosity is 0.37205 (the replication above gave 0.37220) and the
sampled slip 0.1063 m/s (replication 0.1060): ratio_w 0.954, ratio_k 0.991.
