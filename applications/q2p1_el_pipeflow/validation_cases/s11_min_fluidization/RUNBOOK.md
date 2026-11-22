# s11_min_fluidization - minimum fluidisation of a gas-solid bed (cross-code case 11)

Record of 2026-11-19 (preparation and cost probes; the protocol has NOT been run). Spec and reference
(peer session, read-only): `liggghts-dem/LIGGGHTS-PUBLIC/md_docs/case_specs/11_min_fluidization/SPEC.md`,
`reference/umf_reference.py`, `reference/inlet_table.py`. Origin: MFiX PIC07 / NETL experiment.
Branch `feature/el-min-fluidization` (worktree `FeatFloWer-fluid11`), pe pin 43e8fb54. Design note of the
code: `docs/md_docs/el_min_fluidization_prep_2026-11-19.md`. Runs: `/data/warehouse17/rmuenste/runs/fluid11`.

Marks: **[V]** verified here (code read, unit test, run output), **[I]** inferred or computed, not verified
by a run, **[P]** taken from the peer's spec.

## 1. Purpose

Gas flows up through a bed of spheres. Fluidised: the time-averaged pressure drop equals the buoyant bed
weight per area (closure-free, gate F1, 1 +- 0.02). Fixed bed: dp/dz follows the drag law at the bed
porosity (F2, +-5 %); U_mf against the closure at the code's own porosity (F3, +-5 %); against the
measured 0.182 m/s reported only (F4); particle count and bed containment (F5) [P].

## 2. Set-up

| Item | Value |
|---|---|
| Particles | d 0.87 mm, rho_p 860 kg/m3 (pe `benchRadius_` 4.35e-4) [P] |
| Gas | rho_f 1.0 kg/m3, mu 1.8e-5 Pa s; deck `Prop@Density = 1d0,1d0`, `Prop@DynVisc = 1.8d-5,1d0` (nu = mu/rho for the closures; the fluid operator gets mu under the default `ELFluidViscosity = dynamic`; identical numbers for rho_f = 1); json `fluidViscosity_` 1.8e-5 = mu [V] |
| g | 9.81 m/s2 on the particles (buoyant weight (rho_p - rho_f) V g); the fluid carries no gravity (implied by the momentum model): pressure without hydrostatic part [V] |
| Bed R (first) | 15 x 15 mm, laterally periodic, N 9 789 (static height 25 mm at eps 0.40), column 100 mm |
| Bed A | 25 x 25 mm, N 27 190, same column |
| Mesh | coarse cubic cells of 5 mm (bed R: 3 x 3 x 20) or 4.17 mm (bed A: 6 x 6 x 24), `MaxMeshLevel 2`: 2.5 mm = 2.9 d (bed R, as the peer's CFDEM mesh 6 x 6 x 40) / 2.08 mm = 2.4 d (bed A) |
| Boundaries | x, y `Periodic` (`SimPar@PeriodicAxis = xy`, `PeriodicLength = L`); floor `Inflow975` (uniform superficial inflow, `SimPar@ELInletVelocityTable`); top `Outflow` (`SimPar@NoOutflow = No`) |
| Particles in pe | periodic x, y (`periodicX_/Y_`), z walls at the floor and at the column top (`zWallsEnabled_`; the fluid outlet is not open for particles, the top wall is 65-75 mm above the bed) [V, json]; restitution 0.9, friction 0.3 [P]; 10 pe sub-steps per fluid step (as s05, s10) |
| Partition | 3 x 3 x 4 axis-uniform (pe needs >= 3 processes per periodic axis, fluid and pe own the same slabs): 36 workers + master = 37 ranks |
| Coupling | `ELMomentumModel = A`, `ELMMPressureOperator = weighted` (as s10; implies continuity form C, semi-implicit drag, pressure force, no fluid gravity), `ELVACompatibility = none` (open outlet), `ELVASourceTimeScheme = BDF2`, kernel `deen_poly` 2.5 d (radius 2.175 mm) |
| Drag | `ELDragModel = gidaspow` (common closure of the spec; model-A form, see the design note) |
| Pressure | coarse level 1 (`MGMinLev = MGMedLev = 1`, needed by the weighted operator), `MGCrsSolverType = 2` (UMFPACK, safe for open boxes; types 3/4 are refused on weakly open components), periodic coarse operator `gathered` by default for this app |
| Protocol | `make_inlet_table.py`: 0.30 m/s for 0.5 s preparation + 0.5 s step 0, then 13 steps of -0.02 m/s every 0.5 s down to 0.04 m/s, 20 ms linear transitions; end 7.5 s [V, identical to the spec's `inlet_table.py`] |

Choice of the mesh and kernel [I]: the kernel radius 2.5 d is the one of s10 (d 1 mm, cells 3.33 mm,
ratio radius / cell 0.75) and s05; here radius / cell is 0.87 (bed R). Cells of 2-3 d keep the
kernel supported (a 1 x 1 x 40 mesh with cells four times the kernel left particles without fluid support in
s10) and match the CFDEM mesh of the peer (6 x 6 x 40). Bed A cannot keep 5 mm coarse cells with a 3 x 3
partition of 25 mm, hence 4.17 mm.

Momentum model [I]: model A with the weighted pressure matrix is the configuration certified on the s05
column (G6) and on the s10 periodic bed; the plateau F1 is the momentum bookkeeping that model A was built
for. No reason against it was found. Bfull and B are available with the same deck (one key).

## 3. Keys of this case (new on this branch)

- `SimPar@ELDragModel = gidaspow` (also `wenyu`, `ergun`).
- `SimPar@ELInletVelocityTable = "_data/inlet_table.txt"` (rows "t U"; excludes the three ramp keys).
- `SimPar@ELBedAuditFreq = 10` (EL_BED_AUDIT line every 10 steps; default 0 = off).

## 4. Staging and submission

    T=/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-fluid11/applications/q2p1_el_pipeflow/validation_cases/s11_min_fluidization/tools   # absolute: the next line changes directory
    cd /data/warehouse17/rmuenste/runs/fluid11
    # protocol run of bed R (rsa bed = the peer's CFDEM bed, seed 1)
    $T/stage_s11.sh R_proto                                # S11_BED=A for bed A
    $T/fl11_submit.sh R_proto /data/warehouse17/rmuenste/runs/fluid11/bin/q2p1_el_pipeflow.<commit> 24:00:00 long
    # constant-velocity probe
    S11_INLET=const:0.30 $T/stage_s11.sh R_P30 SimPar@MaxNumStep=400 SimPar@ELParticleLogFreq=100
    # evaluation
    python3 $T/fl11_post.py R_proto --steps R_proto/_data/steps.txt

`stage_s11.sh` writes the coarse mesh (`make_s11_mesh.py`), partitions it (`PyPartitioner.py 36
axis_uniform x3-y3-z4 FL36P fl_box/file.prj`, copy of the pe_partpy branch
fix/axis-uniform-boundary-precision in `runs/batch10/tools/pe_partpy_axisfix`), writes the bed
(`make_bed_s11.py`: `rsa` = the peer's `make_bed.py` with PERIODIC=1, phi 0.25 over 60 mm, md5-identical
for seed 1 [V]; `lattice` = jittered simple cubic, solid fraction 0.43, top 34.6 mm, a dense start for
fixed-bed probes), the inlet table and step list, the deck (`q2p1_param_s11.dat` with the mesh name,
part count and period filled in, further keys from the command line) and the json (`example_s11.json`,
`stepsize_` and `timesteps_` following the deck). `fl11_submit.sh` writes `job.sbatch` and
`versions.txt` (commit, binary md5, bed md5, deck keys, table md5) and submits. Jobs point to binary
copies in `runs/fluid11/bin`.

## 5. Gate of the default path

Default path (none of the new keys) against the base binary of c76ca927, both built in
`build-fluid11` on longhouse1a by the recipe of `runs/pcoarse/tools/build_pc.sh` (gcc 14.3, OpenMPI 4.1.6,
USE_PE, FBM acceleration; copies in `runs/fluid11/bin`), the same decks, compared with
`runs/pcoarse/tools/pc_compare.sh` (EL streams strict, solver lines, particle logs) [V]:

| Deck (clone of) | what it covers | job | base vs 7f05a069 | base vs 87474a41 (final) |
|---|---|---|---|---|
| D_s6 = `pcoarse/round6/D_s6_gath_chk` (s10 periodic, 20 steps, 28 ranks) | difelice_std, model A weighted, wrap kernel, periodic coarse operator | 163807 / 163822 | bitwise | bitwise |
| G_tier2 = `round6/E_gates/G_tier2_twin__unset` (60 steps) | stokes | 163807 / 163822 | bitwise | bitwise |
| VB40 = `round6/E_gates/G_VB_L2_default_kin__unset`, 40 of 200 steps | difelice (campaign form) | 163807 | bitwise | not rerun (the commits after 7f05a069 touch the parser for periodic decks and the audit only) |
| E2w5 = `runs/ergun05/E2p_ramp0.4_A_w2`, 5 steps | Inflow975 ramp, cylinder wall, model A weighted | 163807 / 163822 | bitwise | bitwise |

Base and 7f05a069 ran on summergames, the final binary on nx-05. D_s6 of the base binary is also bitwise
with `pcoarse/round6/D_s6_gath_chk` (their c76ca927 binary). The fixed-bed probe of this case gives the
same EL_FORCE_BUDGET stream with ab23bb6e (R_P10_dt25) and 87474a41 (R_P10w_dt25): the audit does not
touch the solution [V]. Unit tests (Slurm 163806 at 7f05a069, 163820 at 87474a41): 12/12 passed
(`el-drag-closures`, `el-inlet-table`, `el-bed-audit` and the nine existing EL / pco tests) [V].
Schiller-Naumann has no gate deck here (unit test `el-kernel-forces` only).

## 6. Probes at density ratio 860 (bed R)

All with the binary of the branch (ab23bb6e; R_P10w and A_P10 with 87474a41), 37 ranks on nx-05, dt
2.5e-4 s (fluid) with 10 pe sub-steps, constant inlet velocity, starting from the generated bed (not
settled). Run directories `runs/fluid11/<run>` (versions.txt, staged.txt). The first submission of
R_P30/R_P10 (jobs 163808/163809, binary 7f05a069) stopped at the start: the parser closed every periodic
box (`bNoOutflow` forced), fixed in ab23bb6e (keeps an explicit `SimPar@NoOutflow = No`).

| Run | job | bed, U | steps | stable | pressure MG | continuity after the corrector | max overlap |
|---|---|---|---|---|---|---|---|
| R_P30_dt25 | 163812 | rsa (phi 0.25 to 60 mm), 0.30 m/s | 400 (0.1 s) | yes, no NaN, no warning, no clipping | 2.01 cycles/solve, none at the cap of 20, DefFinal/DefInit 7e-5 | res1_L1 6e-8 (predictor 1e-5) | 1.12 % d (0.083 s, collapse) |
| R_P30L_dt25 | 163817 | same, 0.30 m/s | 2000 (0.5 s) | yes | 2.17 cycles/solve, none at the cap, 1.7e-4 | res1_L1 median 1.8e-8, max 1.1e-7; superficial flux eps u_z / U = 1.0002 (0.989 ... 1.006) after 0.2 s | 1.12 % d at 0.083 s (impact of the collapsing loose bed), <= 0.56 % d after 0.2 s |
| R_P10_dt25 / R_P10w_dt25 | 163813 / 163822 | lattice (phi 0.43 to 34.6 mm), 0.10 m/s | 400 (0.1 s) | yes | 2.05 cycles/solve, none at the cap | res1_L1 7.7e-10 at the end | 0.41 % d |
| A_P10_dt25 | 163821 | bed A, lattice (phi 0.48 to 31 mm), 0.10 m/s | 200 (0.05 s) | yes | 2.10 cycles/solve | res1_L1 2.6e-8 | 0.82 % d |

No alternating or growing quantity was seen in the continuity audit (predictor and corrector values of a
step agree, flux within 0.6 % of U) [V]. Density ratio 860 needs no smaller step than 2.5e-4 s in these
probes; 1e-4 and 5e-5 were not needed and not run. A larger step (5e-4) was not tried.

**Pressure drop and hydrostatics [V].** `dp_hydro = 0` in every line (the fluid carries no gravity under
model A), and dp_dyn equals the force-based value (drag + pressure force summed over the particles,
EL_FORCE_BUDGET) / A to 0.2 % in the median over the audit steps after 0.1 s (R_P30L_dt25: 90 % of the
steps within 1.0 %, largest 1.4 % at step 1770; R_proto_hrz: median 0.2 %, 90 % within 0.7 %, largest
1.9 %, mean signed difference 0.00 %): the plane pressures close the momentum balance of the column in the
mean, single steps of the slugging bed deviate by up to 2 %. (An earlier version of this paragraph said
"0.1-0.3 % at every audit step"; corrected after the review of 2026-11-21.) The 0.98 Pa of the gas column are not contained.

**Fluidised, 0.30 m/s (R_P30L) [V].** The loose bed collapses (z99 59 -> 26 mm in 0.1 s, window porosity
0.75 -> 0.42), dp overshoots to 1.24 W/A at 0.1 s, the bed then lifts and bubbles: z50 oscillates between
12.8 and 15.9 mm. Means (W/A = 126.41 Pa, W = 2.844e-2 N):

| window | dp_dyn A / W | std | (drag + pressure force) / W | floor contact / W | z99 | window porosity |
|---|---|---|---|---|---|---|
| 0.20-0.50 s | 0.985 | 0.049 | 0.989 | 0.019 | 28.7 mm | 0.489 |
| 0.20-0.35 s | 0.989 | 0.040 | 0.992 | 0.019 | 28.4 mm | 0.483 |
| 0.35-0.50 s | 0.980 | 0.056 | 0.984 | 0.019 | 29.1 mm | 0.496 |

dp + floor contact = 1.004 W: the 1.5 % below the weight is carried by the floor (particles touching the
inlet plane in the bubbling bed). Not a plateau measurement yet: 0.3 s after a start-up transient, 5 %
fluctuation, period of the oscillation not resolved by the window [I].

**Fixed bed, 0.10 m/s (R_P10w) [V].** Settled by 0.075 s: window porosity 0.4269 (z 1.74 ... 24.7 mm),
z99 26.0 mm, floor contact 0.77 W, dp_dyn 0.3148 W/A (39.8 Pa). Against Gidaspow (= Ergun at this
porosity) at the measured porosity:

| reading | measured | closure | ratio |
|---|---|---|---|
| window gradient -<dp/dz> vs closure(eps 0.4269) | 1591 Pa/m | 1654 Pa/m | 0.962 |
| dp / H_eq (uniform bed at 0.4269, 26.16 mm) | 1524 Pa/m | 1661 Pa/m | 0.918 |
| dp / H (H = z99 + r = 26.4 mm, eps 0.4325) | 1508 Pa/m | 1560 Pa/m | 0.967 |

The bed core is 3.8 % below the closure (inside the F2 tolerance of 5 %); the whole-bed readings differ by
the looser end layers (floor and top). Wen-Yu at the same porosity would give 1853 Pa/m. Gidaspow U_mf at
this porosity: 0.255 m/s [V, fl11_post closure functions].

**Bed A (A_P10, 0.05 s, not settled) [V].** Runs; the bed is still collapsing at the end (window porosity
0.45, floor 1.1 W): no pressure statement.

## 7. Cost and projection

Wall time per step (Statistics.txt, overall time / steps, start-up included; nx-05, 37 ranks) [V]:

| Run | regime | steps | overall | per step | fluid solver share |
|---|---|---|---|---|---|
| R_P30_dt25 | loose bed falling | 399 | 164 s | 0.41 s | 14 % |
| R_P30L_dt25 | bubbling bed at 0.30 m/s | 1999 | 1102 s | 0.55 s | 11 % |
| R_P10w_dt25 | collapse + fixed bed | 399 | 255 s | 0.64 s | 9 % |
| A_P10_dt25 | bed A, collapsing | 199 | 388 s | 1.95 s | 7 % |

The rest of a step (86-93 %) is the particle work (EL transfer passes, pe with 10 sub-steps,
communication); `Statistics.txt` does not split it further. The pressure solve of the fluid is 6-9 %.

Projection of the protocol (7.5 s = 30 000 steps at 2.5e-4 s) [I]:

- Bed R: 0.55-0.64 s/step -> 4.6-5.3 h. One job (med partition 8 h, or long) fits with margin.
- Bed A: 1.9-2.0 s/step measured in the collapse phase (ratio A/R about 3.0 for 2.78 times the particles
  and 4.8 times the cells) -> 16-17 h, about 14-20 h with the spread between the regimes of bed R. Fits the
  24 h limit with about 30 % margin; a slow node (s05 saw one node 1.5 times slower) would not fit.
- A restart re-inserts the particles from particles.xyz (pe `resume_`/checkpointer not used by this
  application; project memory), so the protocol cannot be chained over jobs: a run that hits the limit is
  lost. If bed A has to be safe: more ranks along z (3 x 3 x 5 = 46 ranks puts the bed on two pe slabs;
  coarse 6 x 6 x 24 does not divide by 5, needs nz 25 or 30), a larger step (5e-4, untested), or a
  particle-state restart (not available).

Recommendation: bed R as the common bed (5 h per protocol run, room for a closure control run and a step
check in a day); bed A only after bed R works in every code, as one single-shot run of about 16 h.

## 7a. Velocity spikes of the first protocol run and the remedy (2026-11-19, evening)

First protocol run of bed R (orchestrator: `R_proto`, job 163872, binary 87474a41, deck as committed then):
NaN at step 11 315 (t 2.83 s, U 0.22 m/s, bed fluidised), the velocity solve diverged. The maximum
velocity norm had spikes of 2-15 m/s from step 2618 on, while dp, flux and continuity stayed normal [V].

Diagnosis (design note section 7): rerun of the first 4 800 steps with the new line
`EL_UMAX_AUDIT` every step (`R_diag`, job 163992, binary d3d23e4b; EL_FORCE_BUDGET bitwise with R_proto
over all 4 801 lines, i.e. the same trajectory). Max 15.45 m/s at step 3549 (the run's 15 at ~3550).
Of the 1 875 audit lines with |u| > 1.2 m/s, 1 674 (89 %) have the maximum at a dof with a NEGATIVE
deposited implicit drag coefficient B_i, all at z = 30.0-37.5 mm (bed surface, transport porosity
0.77-1.0); none on the inlet plane (largest inlet w 0.76 m/s, smallest inlet porosity 0.39), none near
the outlet, smallest element porosity 0.355 (no clipping). The rest sit at neighbouring freeboard dofs
(B = 0) of the same episode. Cause: Q2 row-sum weights of the implicit drag go negative at nodes reached
only by the negative lobes of a basis function; in the gas the anti-damping thstep |B_i| is a few
percent of the lumped fluid mass of a vertex dof per step [V; growth rate consistent with the logged
B_i]. Orchestrator hypotheses (a) inlet, (b) low porosity, (d) outlet excluded by the same data; (c) CFL
is a consequence, not the cause; (e) is the closest (sign, not size, of the semi-implicit coefficient).

Remedy: `SimPar@ELDragImplicitWeights = hrz` (commit 77ffd070; HRZ lumping of the implicit coefficient
and its split partner, non-negative, same element total; default rowsum bitwise: D_s6, G_tier2, E2w5
bitwise against the base binary, job 163999). The case deck now sets it.

Verification (`R_hrz`, job 163999, binary 77ffd070, the protocol deck + `hrz`, velocity audit every 10
steps, 5 200 steps = 1.3 s: 0.30 m/s for 1 s, then 0.28 m/s) [V]: completed, no NaN, no warning. Largest
velocity norm over ALL 5 200 steps 1.048 m/s (step 388, collapse of the start bed; afterwards 0.8-0.9 m/s,
the interstitial 0.3/0.38 inside the bed); no audit line above 1.2 m/s, no maximum at a dof with B < 0
(none exists now), the maxima sit inside the bed. Velocity nonlinear loop 1 iteration in every step,
largest velocity multigrid rate 2e-13. Smallest element porosity 0.356, inlet w at most 0.76 m/s.
dp_dyn A/W 0.973 +- 0.068 (0.7-1.0 s) and 0.970 +- 0.067 (1.0-1.3 s) with 2.3-2.4 % of the weight on the
floor; (drag + pressure force)/W 0.979 / 0.975. Cost 0.56 s per step (2 933 s), as before.
The trajectory differs from R_proto from step 2 on (different distribution of the implicit drag).

Full protocol resubmitted: `runs/fluid11/R_proto_hrz`, job 164004 (partition long, 10 h limit, binary
77ffd070, deck = R_proto + `ELDragImplicitWeights = hrz` + `ELVelocityAuditFreq = 10`; expected about
4.7 h). What to look at: `python3 tools/umax_summary.py R_proto_hrz/run_slurm.log 1.5` (runs/fluid11/tools),
`grep 'maximum norm of velocity'`, and the evaluation `fl11_post.py R_proto_hrz --steps
R_proto_hrz/_data/steps.txt`.

## 8. Results of the protocol (bed R, 2026-11-19/20)

**Status after the review of 2026-11-21 (section 10): the runs of this section used a contact friction
of 0.6 (not the specified 0.3) and the first form of the hrz split. They are kept as the record of that
configuration; the rerun with the corrected deck and binary is in section 11. The gate table below is
the evaluation by the rules in force on 2026-11-20; the amended rules are applied in section 10.**

Runs (`runs/fluid11`, details in runs_versions.txt) [V]:

- `R_proto` (job 163872, binary 87474a41, deck before `hrz`): FAILED at step 11 315 (t 2.83 s, U 0.22 m/s):
  velocity spikes from step 2618 (negative implicit drag weights at the bed surface, section 7a), velocity
  solve diverged. No gate evaluated.
- `R_proto_hrz` (job 164004, binary 77ffd070, deck = R_proto + `ELDragImplicitWeights = hrz` +
  velocity audit): COMPLETED in 6 h 10 min (0.74 s per step on this node; the probes gave 0.55-0.64),
  largest velocity norm 1.05 m/s. Evaluation `fl11_post.py R_proto_hrz --steps R_proto_hrz/_data/steps.txt`
  (post tool of 2026-11-20: F1 and F3 by the agreed rules, floor force from the particle momentum balance).

Per step (0.2-0.5 s of each step; dp* = dp_dyn A / W):

| U m/s | 0.30 | 0.28 | 0.26 | 0.24 | 0.22 | 0.20 | 0.18 | 0.16 | 0.14 | 0.12 | 0.10 | 0.08 | 0.06 | 0.04 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| dp* | 0.973 | 0.983 | 0.975 | 0.979 | 0.984 | 0.987 | 0.990 | 0.870 | 0.752 | 0.636 | 0.523 | 0.413 | 0.305 | 0.201 |
| floor / W (balance) | 0.020 | 0.017 | 0.019 | 0.016 | 0.016 | 0.012 | 0.012 | 0.132 | 0.250 | 0.366 | 0.478 | 0.588 | 0.695 | 0.799 |
| dp* + floor / W | 0.993 | 1.000 | 0.994 | 0.995 | 0.999 | 1.000 | 1.002 | 1.001 | 1.002 | 1.002 | 1.001 | 1.001 | 1.001 | 1.000 |
| core porosity | 0.505 | 0.480 | 0.459 | 0.434 | 0.411 | 0.387 | 0.370 | 0.370 | 0.370 | 0.369 | 0.369 | 0.369 | 0.369 | 0.369 |
| core gradient / closure | | | | | | | | | 0.946 | 0.947 | 0.947 | 0.949 | 0.951 | 0.954 |

The floor force of the first evaluation (pe wall-contact impulse) read 10/9 of the applied force
(dp* + floor 1.016-1.089 in the fixed branch); it is the unrelaxed PGS impulse of pe, design note
section 8. The table uses the balance value.

| Gate | rule | FeatFloWer (R_proto_hrz) | verdict |
|---|---|---|---|
| F1 plateau | mean dp* over the steps with dp* > 0.95 (0.30-0.18 m/s, 7 steps) | 0.981 (dp* + floor 0.998) | PASS (1 +- 0.02) |
| F2 fixed branch | core gradient / Gidaspow at the core porosity 0.369, 6 steps 0.14-0.04 | 0.946-0.954 (worst -5.4 %); interior 5-20 mm 1.014-1.024; whole bed (H = z99 + r) 0.964-0.974 | FAIL (+-5 %), split in design note section 9 |
| F3 U_mf | Ergun form dp* = a U + b U^2 through the origin over the fixed steps, at dp* = 1 / at the plateau | 0.1813 / 0.1783 m/s against Gidaspow 0.1642 at 0.369: 1.104 / 1.086 (linear fit, old rule: 0.1820, 1.109) | FAIL (+-5 %) |
| F4 | U_mf / 0.182 | 0.996 / 0.980 | reported |
| F5 | particle count, highest centre at 0.04 m/s | 9789 throughout; 24.0 mm < 50 mm | PASS |

Other codes on the same bed and protocol [P, peer messages of 2026-11-20, not recomputed here].
**Historical: the CFDEM row is the top-driven run that the peer withdrew on 2026-11-21 (gas entered as
jets) and the MFiX row is an interim run; the reading below that FeatFloWer packs denser no longer holds.
Current peer values are in section 10.**

| Code | plateau | fluidised steps | packed porosity | gradient / closure | U_mf / closure |
|---|---|---|---|---|---|
| CFDEM (bed C) | 0.942 (0.932-0.952) | 0.22-0.30 | 0.390 | 0.869 (whole bed) | 1.07 / 1.13 |
| MFiX-DEM (interim) | 0.959, 0.976 | 0.28-0.30 only | 0.392 | 0.754 (core) | 1.33 (unexplained) |
| FeatFloWer (this record) | 0.981 | 0.18-0.30 | 0.369 | 0.946-0.954 (core), 0.964-0.974 (whole bed) | 1.104 / 1.086 (Ergun form), 1.109 (linear) |

Readings [I]: the defluidised FeatFloWer bed packs denser (core 0.369) than the spec table assumes (0.40)
and than the CFDEM / MFiX beds (0.39), so its closure U_mf is lower (0.164 m/s); that the measured U_mf
(0.181) equals the NETL value is a coincidence of packing and closure, not a validation. The core
deficit against the closure sits at the window ends (kernel 2.5 d against a 2 d exclusion); the window
interior is 1-2 % above the closure. Bed A not run.

## 9. Open points

1. Protocol range: answered by the run - the defluidised bed packs at 0.369, U_mf about 0.18 m/s, seven
   fluidised steps (0.18-0.30 m/s) qualify for F1. (The probe-based worry of 2026-11-19 about 0.30 being
   too close to U_mf was based on the looser collapsed beds of the probes.)
2. (Settled 2026-11-20, section 10: F2 is taken at the mapped gas fraction, and only from a window that
   excludes the mapping support.) F2 definition and the window: with a 2.5 d kernel the 2 d exclusion of the spec leaves the kernel
   of the window-edge particles reaching the loose end layers; core 0.946-0.954, interior 1.014-1.024
   (design note section 9). A window of kernel radius + 1 d, or the like-for-like closure at the
   kernel-sampled porosity (ratio_k, audit field eps_kern_slip_win, from 2026-11-20 runs on), was to be
   agreed with the peer (wording of 2026-11-20, kept as history).
3. pe: the wall-contact impulse and the contact virial report the unrelaxed PGS impulse (1/0.9 of the
   applied one); fix needed in pe (design note section 8); the bed audit uses the momentum balance.
4. `ELDragImplicitWeights = hrz` as default: owner decision; effect on the s10 liquid gate D_s6 over
   200 steps measured (design note section 7, RUNBOOK of this case section 7a).
5. (Settled 2026-11-20, section 10: the F3 reference is the closure at the whole-bed porosity; this
   run then gives 1.042, PASS.) F3 FAIL under the first rule (+8.6 / +10.4 %): dp* is the pressure drop of the whole bed, whose mean porosity over
   H = z99 + r is 0.376 against the core 0.369 (looser floor and top layers). Gidaspow U_mf at 0.376 is
   0.174 m/s, against which the measured 0.181 is +4 %, close to the whole-bed gradient ratio
   0.964-0.974. Which porosity the F3 reference takes was a definition to agree with the peer [I]
   (wording of 2026-11-20, kept as history).
6. Bed A cost (section 7) and the missing particle-state restart; bed A not run.
7. Not run: a step-size study, models B / Bfull, closures wenyu / ergun as controls.

## Particle movie of the protocol run

`/data/warehouse17/rmuenste/runs/fluid11/movie/fluid11_R_proto_movie_particles.mp4` (52 s, 2560 x 1440,
30 fps, 74 MB; not tracked). Run `R_proto_movie` (job 164399, binary c8541cfc, 5 h 59 min): the deck of
`R_proto_hrz` with `SimPar@ELParticleLogFreq = 20` (one particle frame every 5 ms) and
`SimPar@OutputFreq = 0.0249d0`; its per-step table is identical to `R_proto_hrz` to the printed
digits [V]. One image per logged particle frame, no interpolation (1 501 frames); largest displacement
of a particle between two frames 1.37 d after t = 0.2 s, 3.2 d during the collapse of the loose start
bed (`movie/check_R_proto_movie.log`) [V]. Tool `tools/fl11_movie.py` (`--check`, `--stats`,
`--frames`), render job `tools/fl11_movie_render.sbatch` (final render job 164578, 5 min 45 s, Slurm; the first render, job
164533, carried the CFDEM markers of the peer's void top-driven run and was replaced).
Content: x-z and y-z cuts and a top-down slab (1 d thick, discs = largest section of the sphere in
the slab, colour = vertical particle velocity, +-0.10 m/s), solid fraction in 1 mm bins, traces of
U, dp*, dp* + floor / W and the bed top with a time cursor; MFiX (run F_none) and CFDEM (run fluid_C_gid_inlet2,
uniform inlet) step means as peer-reported markers [P]. What it shows: in this 17 d wide periodic box the fluidised bed slugs
(a plug rises, particles rain through a void band below it, about ten cycles per second), it does
not bubble; the regime label follows the F1 rule, so the 0.18 m/s step reads "fluidised" with the
particles at rest.

## 10. External review of the branch (GPT-6.2 Astra, high; 2026-11-21) and what it changed

Verdict CHANGES REQUIRED, 15 findings; brief and report in the session scratchpad
(`review_brief_fluid11.md`, `codex_review_fluid11_r1.md`). Checked here against code and logs [V]:

1. **Contact friction was 0.6, not 0.3.** pe adds the friction coefficients of the two materials of a
   contact (`libs/pe/src/core/Materials.cpp`: `csfTable_(mat,mat) = csf + csf`, mixed pairs
   `csf_i + csf_j`), and the EL setup passes the json value to every material. `staticFriction_ = 0.3`
   therefore gave 0.6 for particle-particle and particle-floor contacts in every run of sections 5-8.
   The json now carries 0.15 per material (commit d2743338). The convention is general in pe and applies
   to every earlier case that set a friction value (not re-evaluated here; owner decision).
2. **The restitution 0.9 of the spec is not used.** The hard-contact EL solver (`HardContactEulerLagrange`)
   is inelastic; it never reads the material restitution. Deviation from the spec, stated to the peer.
3. **`hrz` split not momentum-consistent** in a steady non-uniform field: the fluid received
   `B (u_f - sum_i h_i u_i)` per particle beyond the drag reaction. Fixed in the fix round (design note
   section 7); global balance of the recorded run: dp* + floor = 1.000-1.002 in the fixed bed.
4. **Gates by the amended spec** (post tool of 2026-11-21): F3 against the whole-bed porosity 0.3761:
   Gidaspow U_mf 0.1740 m/s, ratio **1.042, PASS** (at the plateau 1.025; at the core porosity 1.104 /
   1.086 as before). F2 against the closure at the mapped gas fraction: ratio 0.984-0.993 on
   `R_proto_movie` (the log of `R_proto_hrz` has no mapped value), **not gated**: the audit window of
   these runs excludes 2 d, not the 2.5 d kernel support the amended spec asks for, so the amended F2
   is unavailable for them (corrected 2026-11-22 after the second review, section 12; this item first
   called it a PASS). The fix round widens the window.
5. **Overlap limit of the spec (1 % of d) exceeded during preparation only**: `R_proto_hrz` 1.30 % of d
   at 0.0825 s (four log lines above 1 %, 0.0775-0.09 s, landing of the loose start bed); after 0.2 s
   the largest is 0.67 % of d, in the evaluation windows (after 0.7 s) 0.61 %.
6. Post tool: step starts of the echo-only step list (were the transition ends), disjoint sub-windows,
   F5 over the whole plateau of the lowest velocity and only for a complete protocol. The recorded F5
   stands (largest centre height in that plateau 24.03 mm).
7. Fortran findings (inlet-table reader and reload, `floor_bal_z` outside its assumptions, velocity audit
   coefficient, bed-audit geometry, closure dispatch, undefined logical under `PeriodicAxis`): fix round,
   listed with commits in the design note.

Peer values of 2026-11-21 [P, not recomputed]: CFDEM with a uniform inlet packs to 0.3695, plateau 0.955,
U_mf 0.191 m/s (1.11 of the closure at whole-bed porosity) at 40 cell rows and 0.976 / 0.362 / 0.168
(1.03) at 80 rows (not mesh-converged); MFiX (run F_none) 0.995 / 0.370 / 0.173 (1.02). The first CFDEM
and DPMFoam runs (top-driven) are void.

## 11. Reruns after the review: friction 0.3 per contact (2026-11-21/22)

Same deck, bed, mesh, inlet table and partition as the protocol run of section 8 (dense output as in
`R_proto_movie`), json friction 0.15 per material = 0.3 per contact, contacts still inelastic (e = 0).
Evaluation `fl11_post.py` (rules of the amended spec), output in `post.txt` of each run directory [V];
re-evaluated with the post tool of fix round 2 (section 12) in `post_r2.txt`, which is the source of
the F2 column. The three runs have the same physical deck; the output settings differ (`R_r1_proto` and
`R_proto_movie`: particle log every 20 steps, fields every 0.0249 s; `R_r1_oldsplit`: 400 steps, 0.5 s).

| run (job) | binary | drag split | F1 plateau | dp* + floor | U_mf (dp* = 1) | whole-bed porosity | F3 ratio | F2 worst | packed core porosity | F5 |
|---|---|---|---|---|---|---|---|---|---|---|
| `R_proto_movie` (164399), friction 0.6 | c8541cfc | first `hrz` | 0.981 PASS | 0.998 | 0.181 | 0.3761 | 1.042 PASS | 1.6 %, not gated (window inside the mapping support) | 0.370 | PASS |
| `R_r1_proto` (164674) | 627caf47 | consistent `hrz` | **0.9787 FAIL** | 0.9935 | 0.1770 | 0.3726 | 1.047 PASS | 3.5 % PASS | 0.366 | PASS |
| `R_r1_oldsplit` (164675) | c8541cfc | first `hrz` | 0.9826 PASS | 0.9975 | 0.1764 | 0.3720 | 1.049 PASS | 2.4 %, not gated (window inside the mapping support) | 0.365 | PASS |

- **Friction** (0.6 -> 0.3, same split): the fixed bed packs denser (core 0.370 -> 0.365, bed top
  24.07 -> 23.9 mm), U_mf falls from 0.181 to 0.176 m/s, the fluidised bed stands 0.1-0.5 mm higher.
- **F2 of the two runs with the first split**: their binaries exclude 2 d (1.74 mm) at the window ends,
  the mapping support is 2.5 d (2.175 mm), so under the amended spec F2 is unavailable for them; the
  ratios are measurements. `R_r1_proto` excludes the support (binary 627caf47) and is gated. The 3.5 %
  against 2.4 % therefore also compares two windows, not only two drag splits.
- **Drag split** (same friction): the plateau means agree to 0.004 in dp* (single steps differ by up to
  0.0064: 0.9769 against 0.9833 at 0.28 m/s) and U_mf to 0.0006 m/s, but the
  plateau of the consistent split is 0.9787, 0.001 outside the F1 gate (recorded as FAIL), and the
  global balance dp* + floor force is 0.9935 on the plateau against 0.9975 (fixed bed: 1.0001 against
  1.0005-1.0013, so the consistent split closes the fixed-bed balance better).
- **Peak gas speed** (`EL_UMAX_AUDIT`, mean of the audited maximum per velocity step): consistent split
  1.22 m/s at U = 0.30, 0.89 at 0.18, 0.23 at 0.04; first split 0.80 / 0.52 / 0.12. In the packed bed at
  0.04 m/s the interstitial speed is 0.109 m/s: the first split shows it, the consistent split shows
  twice that at single (Q2 vertex) dofs. The consistent split therefore carries a node-to-node
  oscillation of the gas velocity inside the bed. Second review (section 12): the saved fields confirm
  it (vertex values of w in the packed bed from -0.028 to 0.222 m/s, six of 288 vertices negative,
  against 0.101 to 0.115 with the first split); the corrected diagonal acts on the velocity increment
  and gives no restoring force for a stationary nodal pattern, while the first split damped the pattern
  at the price of the steady force error of review round 1. The run is not a validation of local gas
  velocities, and the oscillation can bias the gates (F1 moved by 0.004). Remedy (damping term or
  sampling) is an owner decision, not done here.
- Overlap: the 1 % limit of the spec is exceeded during the preparation, 1.17 % of d at 0.07 s (collapse
  of the start bed); at most 0.75 % of d after 0.2 s (`R_r1_oldsplit` 1.19 % / 0.81 %).
- Cost 6 h 23 min each on 37 ranks.

Restitution (e = 0.9) needs the pe change `feature/el-restitution` (pe pull request #59, seven external
review rounds on 2026-11-21/22); the run with it will be recorded in section 13 (job 164928, not finished at this revision). `R_r1_rest` (164701, first
version of the pass, no rest velocity), `R_r2_rest` (164872, review round 2, budget 24) and `R_r3_rest`
(164910, review round 3) are development runs and not part of the record; `R_r4_rest` (164919) was
cancelled (its binary doubled the normal increment in the friction limit of the pass).

## 12. Second external review of the branch (GPT-6.2 Astra, high; 2026-11-22) and fix round 2

Reviewed head c430b4d2, verdict CHANGES REQUIRED. Report and brief are in the session record; the
findings and what was done [V unless marked]:

1. **Post tool awarded an amended F2 PASS to runs whose window lies inside the mapping support** (P0):
   F2 is now gated only if every audit line of the fixed-step windows excludes max(2 d, mapping support)
   at both ends; otherwise the ratios are printed as a measurement and the verdict is unavailable.
   Section 11 corrected accordingly (1fbb28e4).
2. F2 uses the dynamic window gradient (`gradz_win_dyn`); an older log falls back to the total gradient
   only if its hydrostatic part is zero on every audit line used (true for all runs of this case).
3. Post tool: unavailable audit states no longer stop the evaluation, F1-F3 and F5 are evaluated
   independently, a step counts only with a densely sampled window, a `--steps` file is checked against
   the echoed protocol, F5 covers the constant part of the lowest plateau (the 20 ms transition is
   printed as a labelled stricter number). `fl11_post.py --self-test`: 12 synthetic logs.
4. `ELDragImplicitWeights = hrz` is refused at start-up where the momentum defect is hand-built without
   the drag diagonal (non-Newtonian flow with `MatrixRenewal` S = 0); case 11 takes the Newtonian
   assembled-matrix branch (729f77c9). Follow-ups, not done: the missing defect term itself; the default
   `rowsum` semi-implicit drag reaches the same branch and is still accepted (as before this branch).
5. Velocity audit: a term the solve did not apply (explicit drag, feedback off, prescribed field) is
   printed as NaN; the audit is refused for a communication scheme other than 3 or 4 (3e801eab).
6. A reload reports a changed inlet definition (start-up-only stays the contract); a ramp-only deck gets
   its start-up notice (95e5923a). Tab-only and tab-indented comment lines of the inlet table (c0d71b65).
7. Render job: fresh frame directory per render, each worker waited for, no encode after a failed
   worker, the old movie replaced only on success (2cfeb048).
8. Record: this section, the corrections in sections 9 and 11, datasheet rows (superseded gate rules
   labelled), design note 7.1 (the drag diagonal acts on the momentum predictor; the conservation
   statement is limited to the given-field unit tests) and 7.2.

Evidence of fix round 2 (binary d27a74fb-pece2a7f7, build job 164949): EL unit tests 15/15 (job 164950,
`unit/unit_r2-164950.out`), failing-before runs of items 4-6 (job 164952, `fb_r2/`) and of the tab case
(job 164948); gates (job 164951, `gate/pc_compare_d27a74fb-pece2a7f7.txt`): D_s6, G_tier2, E2w5, VB40 and
400 steps of case 11 without restitution bitwise against the references (all common streams, particle
logs 27/27), and 400 steps of case 11 WITH restitution bitwise against `R_id400_new` (13 077 `EL_` lines
including the bed, velocity and restitution audits, 800 solver lines, 27/27 particle logs). The last gate
is what will let the record run of section 13 (Fortran sources 627caf47) stand for this head.

The reviewer's judgement of the gas-velocity oscillation is in section 11. Still explicit after both
reviews: F1 of `R_r1_proto` fails by 0.001, the overlap limit is exceeded during the preparation, and
the friction-only runs are inelastic.
