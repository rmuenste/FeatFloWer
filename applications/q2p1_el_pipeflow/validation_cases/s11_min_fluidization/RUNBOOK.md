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

    T=applications/q2p1_el_pipeflow/validation_cases/s11_min_fluidization/tools
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
EL_FORCE_BUDGET) / A to 0.1-0.3 % at every audit step: the plane pressures close the momentum balance of
the column. The 0.98 Pa of the gas column are not contained.

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

## 8. Results table (to be filled by the protocol runs)

| Run | bed | dt | F1 plateau | fixed-branch slope (Pa s/m) | eps_bed | U_mf fit | U_mf closure at eps_bed | F2 worst | F3 | F4 (/0.182) | F5 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| (not run) | R | | | | | | | | | | |
| (not run) | A | | | | | | | | | | |

## 9. Open points

1. Packing and the protocol range [I, from the probes]: the pe beds settle at a core porosity of 0.427
   (collapsed lattice) and 0.418 (collapsed loose bed, transient), not 0.40 as the spec's table assumes;
   Gidaspow U_mf at 0.427 is 0.255 m/s. With U0 = 0.30 m/s only the first steps are above U_mf and
   possibly none above 1.2 U_mf (the F1 selection). The bed of the protocol is the defluidised one and may
   pack differently. To be raised with the peer before the protocol run: a higher U0 (e.g. 0.40 m/s) or
   F1 over the steps above the measured U_mf.
2. F1 plateau: the 0.3 s window after start-up gave 0.985 +- 0.049 with 1.9 % of the weight on the floor;
   whether the long-time plateau of the protocol reaches 1.00 +- 0.02 is open (bubbling, 0.5 s per step).
3. F2 definition: gradient of the bed core (window, gate value of fl11_post) against dp over the whole bed
   (0.962 vs 0.918 / 0.967 in the fixed-bed probe). The spec measures dp over the column and the porosity
   of the core; the definition should be agreed with the peer.
4. Bed A cost (section 7) and the missing particle-state restart.
5. Not run: the protocol itself (by instruction), a step-size study, bed A beyond 0.05 s, models B / Bfull,
   closures wenyu / ergun as controls.
6. The pe json keeps `boxWallsEnabled_ false` and z walls at the column ends; the top wall at 100 mm was
   never approached (zmax 59.5 mm at the start, 30 mm when bubbling).
