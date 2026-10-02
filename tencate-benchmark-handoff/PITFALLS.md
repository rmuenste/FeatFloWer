# PITFALLS — every difficulty the FeatFloWer campaign hit on this benchmark

Format: **symptom → cause → what to check in a new code.** Sources: datasheet
rows (`reference/datasheet_rows_tencate.csv`), the practitioner's guide
(`campaign_docs/dns_practitioners_guide.md`, §n), the campaign owner's memory
notes (translated here; FeatFloWer key names kept in parentheses so the
original can be recognised), and the paper. Ordered roughly from "silently
wrong physics" to "operational".

## A. Input and coupling consistency

**P1 — Viscosity convention (dynamic vs kinematic) is invisible at ρ = 1 and
nothing checks it.**
Symptom: a settling sphere lands at a plausible-looking but wrong terminal
velocity; in the pre-campaign E4 run the deck had μ = 0.053 instead of
0.058 Pa·s and the peak came out +8.7 % fast (row `d0_visc`). In the EL
(point-particle) path of the same code, putting μ into a slot read as ν gave a
960×-overdamped drag at ρ_f = 960 — "a suspiciously small, perfectly
force-balanced velocity" (memory note `ff-viscosity-convention`).
Cause: FeatFloWer's parameter parser stores `Prop@Viscosity` as-is and derives
the same slot from `Prop@DynVisc`/density; a source comment
(`QuadSc_main.f90:364`) documents the slot as kinematic for the EL closures,
while the resolved DNS/FBM decks and force routine use it as dynamic
(README §9 item 2). Two conventions in one code, no runtime check that the
fluid solver's and the rigid-body solver's fluid properties agree.
Check: (1) print the effective ν and μ at startup and compare with Table I;
(2) assert that the fluid solver and the particle/contact solver read the same
ρ_f, μ_f (row `d0_visc`: FeatFloWer's deck and PE json carry both copies);
(3) at rest, the hydrodynamic force must equal (ρ_p − ρ_f) V g — a check that
catches density but *not* viscosity errors; (4) check u_peak/u_∞ against
Table II for E4 first (fast, 16 min at L3): a μ error of −8.6 % shows as +3.6 %
in u_t by Schiller–Naumann (row `d0_visc`).

**P2 — Fluid and rigid-body integrators silently stepping with different dt
produced a fake "stability floor" and a fake order-of-accuracy study.**
Symptom: reducing dt from 0.5 to 0.25 ms made E4 blow up (growing sawtooth in
particle velocity from t ≈ 0.03 s, |v| to O(1) m/s, run "finishes
successfully") at both L2 and L3; it was recorded as an added-mass instability
of the partitioned coupling (rows `dt_stability`, `dt_stability_l2`, with a
Causin–Gerbeau–Nobile 2005 attribution) and a dt rule "dt has a LOWER bound"
was written. All dt ≠ 1 ms results — a temporal-order fit, a dt ladder — were
also affected.
Cause: the serial rigid-body engine advanced with its own config step
(json `stepsize_` = 1 ms) and ignored the CFD dt passed to it (row
`pe_stepsize_mismatch`: "stepSimulationSerial: fullStepSize =
config.getStepsize()"); at CFD dt = 0.25 ms the body moved 4× too far per CFD
step. Peaks are quasi-static and looked plausible (E4 L3 transient +37 % at
t = 0.02 s decaying to +1.8 % by t = 0.4 s).
Resolution: with both steps equal, dt = 0.25 ms is fully stable (row
`dt_stability_refuted`: 0 sawtooth warnings, +1.94 %); "the dt window is
accuracy-bounded only" (guide §3). A fatal first-step guard now aborts on
mismatch (row `stepsize_guard`).
Check: a single source of truth for dt shared by both solvers, or a startup
assertion that they are equal; a watchdog on step-to-step velocity
sign-alternation (FeatFloWer `DNS_SAWTOOTH_WARNING`, fires during the growth
phase, row `sawtooth_watchdog`). General lesson: when a partitioned coupling
gets *worse* as dt decreases, verify that both halves actually use that dt
before invoking added-mass theory.

**P3 — A default contact/constraint solver that silently discards fluid
forces.**
Symptom: force and velocity exactly 0 for all time, no error message (row
`d0_solver`).
Cause: the rigid-body library's default solver was an Euler–Lagrange variant in
which hydrodynamic forces enter only through a velocity-correction path; forces
written through the FBM interface were a dead end. The DNS builds need a
different solver selected at configure time (memory note
`dns-pe-solver-selection`; guide §1).
Check: a 10-step smoke run whose particle *moves*; a zero-force/zero-velocity
trajectory is a configuration signature, not a physics result. Make the
solver choice visible in the build/run record.

**P4 — Angular velocity zeroed every step by a leftover debug line.**
Symptom: `DNS_PART_STATE` angular velocities exactly 0 for all time under
non-zero torque; a DKT (two-sphere) "frictional stall" that was really the
inability to roll (row `hcaf_angvel_reset`, guide §6).
Cause: `w = Vec3(0,0,0)` in the integrator (an old debug artefact).
Impact on this benchmark: second order — a single sphere on the box axis has
spin ≈ 0 anyway; the certified ten Cate results stand (row `hcaf_angvel_reset`
"certified results stand"), and the fixed binary reproduces the E4 L3 log
byte-for-byte (row `d22_g0_twins`).
Check: verify ω ≠ 0 under torque in any rotation-coupled test before trusting
contact-mechanics conclusions.

**P5 — Stale per-body velocity corrections for fixed bodies (latent).**
Row `pe_stale_dv_fixed_bodies`: correction vectors resized but never cleared;
fixed bodies skip the initialisation, so after a body-index reshuffle a wall can
inherit a mobile body's correction and impart ghost impulses. Campaign impact
none (serial mode, bodies never reorder); relevant to any code with
index-recycled body storage. Check: momentum bookkeeping after
creating/destroying bodies next to fixed walls.

**P6 — Gravity applied twice (fluid body force + particle weight).**
Memory note (D6.4 G0 findings, item 4): FeatFloWer applies the deck gravity to
the *fluid* and the rigid-body side seeds the net weight; setting both gives
double buoyancy. The certified decks have `Prop@Gravity = 0` and gravity only
in the PE json. Check: force at rest = (ρ_p − ρ_f) V g exactly; if it is off
by a factor near 2 or by ρ_f V g, gravity is counted twice or buoyancy is
missing.

**P7 — Symmetry-domain force reconstruction factors.**
Memory note (`ff-deck-staging-pitfalls` item 2): the quarter-box decks carry
`Prop@ForceScale = 0,0,4,0,0,0` (z-force × 4, lateral forces and torques
zeroed). Used in a full domain this is a ×4 loop gain on the vertical force and
produces an exponential coupling divergence (~35 per time unit) that is
*insensitive to dt and density* — "that insensitivity is the fingerprint
distinguishing it from CFL or added-mass". It was the sole cause of a series of
failed viscometer runs that had been blamed on dt and ρ_p/ρ_f = 1.0 (same note,
2026-08-25 update). Check: any force scaling/symmetry factor must be tied to
the mesh file, not to the deck template; if a coupling diverges exponentially
regardless of dt, look for a gain error before a stability theory.

## B. Resolution, time step and what is a converged result

**P8 — Grid-crossing force noise is first order in h (guide §4).**
Symptom: force jitter on a moving sphere; rms of detrended F_z / mean F_z =
3.84 / 1.79 / 0.88 % at D/h = 11.4 / 23.9 / 49.1 (constant-velocity protocol,
row `d12_prescribed_ladder`; free-fall plateaus 4.1 / 2.4 / 0.9 %, row
`d12_plateau_noise`); velocity jitter 6e-6 … 5e-5 m/s.
Cause: the moving interface re-classifies indicator DOFs as it crosses cells —
intrinsic to immersed/fictitious-boundary methods.
Check: tolerances tighter than ~1 % at D/h ≈ 24 measure noise, not physics;
compare peaks from plateau averages, not single samples; expect halving per
refinement level.

**P9 — The resolved force follows lubrication only down to gap ≈ 2 cells; the
final approach to the wall is not a DNS result (guide §5).**
Measured against Brenner's exact sphere-to-plane solution at constant velocity
(Re 0.78, row `d21_prescribed`; figure `reference/figures/d21_brenner_crossover.png`):
the FBM force is within ~10 % of the lubrication divergence down to gap
≈ 2–3 h, −20 % at ≈ 1.1–1.7 h, and the departure curves of D/h = 24 and 49
collapse in gap/h. In the free-fall E1 approach the force departs > 20 % from
Brenner at gap = 1.20 h (L3) / 1.32 h (L4) (row `d21_free_approach`). Below
~2 h the dynamics are the *contact model* (guide §5–6), and the resting gap is
a numerical equilibrium (CASES.md, last section: 0.46 h / 0.23 h / 0 at
L2/L3/L4). The paper had the same problem below one lattice spacing and added
an explicit lubrication force F_w = −6πμ r_p u_⊥ (r_p/h − r_p/Δ_0), Δ_0 = 1
grid spacing [p.4017 eq. (10)], which improved the velocity decay but made the
sphere approach the bottom for "an unrealistically long time" [p.4023].
Check: report the approach phase, do not gate on it; state the gap below which
the contact/sub-grid model takes over; compare codes above gap ≈ 2 h only.

**P10 — Adding a full analytic lubrication force on top of a partly resolved
film double-counts (rows `d22_g2_brenner`, `d22_g2b_deficit`).**
FBM alone: −14.8 % / −25.9 % vs Brenner in the 1h–2h / sub-1h bands; FBM + the
full resistance set: +74.6 % / +67.4 %; FBM + a *deficit* form (analytic minus
its value at the activation gap, i.e. only the unresolved part): +7.2 % /
+24.0 %. With the deficit form the E1/E2 approach reproduces the paper's Fig. 13
behaviour and lands in finite time (row `d22_g3_tencate`). Check: any sub-grid
lubrication model must subtract what the mesh already resolves and switch on at
a mesh-tied gap (here 2 h_min), not at a fixed physical distance.

**P11 — The FBM sphere behaves as if narrowed by ≈ 0.14 h: a_eff = a − 0.14 h
(guide §2; rows `d11_rh_collapse`, `d11_aeff_sign_erratum`, `d11_coarse_probes`).**
Symptom: drag systematically low (−6 … −2 % at D/h 6–24 in a periodic-array
benchmark, −9 … −12 % at D/h 3–4), converging first order in h.
Cause: the discrete no-slip constraint is under-enforced between velocity
nodes; flow penetrates ~0.14 cells into the nominal solid. (A sign erratum
was recorded: v2 of the guide said +0.14 h; the measured deficits force the
minus sign.) The paper's LBM has the opposite bias — its sphere appears
*larger* (hydrodynamic radius r_h/r_0 = 1.12 at r_0 = 4 lu, Table II; Sec. III.C
p.4017) — and the paper calibrates r_h with Hasimoto's periodic-array drag
before every run; without calibration "the velocity ratio is underpredicted
some 20 %" [p.4022]. Check: measure your own effective radius with a
Hasimoto/Stokes drag probe at the production D/h; expect the uniform
L3→L4 shift of −2.1 … −2.4 pp across E1–E4 (row `d13_matrix`) to be an
interface-representation effect, not a flow-regime effect.

**P12 — The good dt = 1 ms, D/h = 24 numbers ride on error cancellation
(guide §3; rows `e4_ladder_l4`, `e4_l3_dt_ladder_sync`).**
Spatial error at L3 is positive (+1.5 … +1.9 pp for E4), the temporal error at
1 ms negative (−1.1 … −1.5 pp); the sum +0.8 % looks converged but neither
term is small. Halving dt moves the peak *away* from the experiment at L3
(+0.81 → +1.40 → +1.94 %) and toward it at L4 (−1.56 → −0.74 %). The temporal
order could not be pinned from the peak metric (row `e4_l3_dt0p7`: "ORDER
PINNING INCONCLUSIVE"). Check: never certify on one (h, dt) pair; run at least
one dt-halving and one refinement and fit err = S(h) + T(dt)
(`tools/tencate_error_decomposition.py` is a 100-line template).

**P13 — Resolution bookkeeping bugs.** (a) A parallel reduction took the
min-over-ranks of each rank's local *maximum* cell size, so the reported h_min
was mode-dependent by 4× (row `d0_hmin`; fixed, now
1.31425429e-3 m / D/h 11.413 in both modes). (b) DOFs inside the particle are
double counted on subdomain-interface planes when rank-local counts are summed
(row `d0_smoke_twin`: 1108 serial vs 1190 parallel, 6.9 %; memory note D6.4
item 5). Check: compute h and D/h from the mesh file, not from a reduction you
have not tested in both serial and parallel; use a multiplicity-aware DOF count.
"D/h counts Q2 elements; nodal spacing is h/2 — halve literature values quoted
for LBM/IBM node counts when comparing" (guide §2).

**P14 — Internal fluid inertia at ρ_p/ρ_f ≈ 1.15.** The paper states that for
an immersed-boundary-type method "the inertia of the internal fluid may affect
the particle motion in a nonphysical manner" and subtracts the rate of change of
internal-fluid momentum from the total force [p.4015, p.4017 eq. (8)], which let
them simulate at density ratio 1.15. Not traced in FeatFloWer for this
hand-off. Check: know how your method treats the fictitious fluid inside the
body; test with a density-ratio sweep (the campaign's DKT probes at 1.10 and
1.02 were clean once P2 was fixed — guide §3).

**P15 — Spheres below D/h ≈ 1.5 are hydrodynamically transparent (row
`d11_ctrl`).** Irrelevant at the resolutions above, but a smoke-test grade of
D/h ≈ 11 already gives +3.4 … +4.4 % (guide §2 "D/h ≈ 11 is smoke-test grade").

## C. Reference discipline

**P16 — Digitised curves vs printed tables (row `tc_ref_audit`, guide §7).**
The digitised E1/E2 velocity minima are +3.4 / +3.9 % faster than the paper's
printed u_max/u_∞; E2 at L3 read −3.1 % (FAIL) against the curve and +0.6 %
(PASS) against the table (row `e2_l3`). Gate peaks on printed values, use the
curves for shape. Pin the reference's own definitions (here: Re on the Abraham
u_∞; St = Re ρ_p/(9 ρ_f)) from its text before believing a discrepancy (guide §7
"third convention pin of the campaign").

**P17 — Confined vs unbounded velocity.** The box slows the sphere to ≈ 95 %
of u_∞ (paper p.4020). Comparing with u_∞ hides an 8.7 % error as 3.8 %
(`tc-ref_README.md`).

**P18 — Time of peak is not an observable at low Re.** Broad plateau; the
digitised E4 curve has 16 samples; the FeatFloWer–PIV RMS drops from 16.4 % to
4.8 % of the peak with a −40 ms time-origin shift (row `e4_ladder_l3`); row
`e3_l3` records a +0.050 s lag. Compare curve shapes with an explicit
time-origin alignment and say so.

**P19 — E1 (Re 1.5) cannot be matched to the experiment within 2 % by any
simulation on record.** ten Cate's LBM −5.6 %, FeatFloWer L4 −5.2 %, Richardson
0.886–0.899 (rows `sim_v_sim`, `err_decomposition`); the paper: "Increased
resolution does not result in an improvement" [p.4021]. Set the E1 reference
band to include the paper's S-series and do not spend resolution on it (guide
§2).

**P20 — Quarter-box symmetry assumption untested against the full box.**
Plan §11 item 5 promised a quarter-vs-full twin row; none exists (README §9).
A full-box run at matched D/h would settle whether the symmetry planes bias the
peak — the sphere is on the axis and the flow is nominally axisymmetric, but
the box is square. If your code runs the full box, you provide that twin.

## D. Operational traps (library-agnostic reading of FeatFloWer-specific notes)

**P21 — Inline comments on numeric input lines** (memory note item 1; guide
§1): a list-directed Fortran read took "… D/h=8 (CASE_SPEC 4)" and set the
mesh level to 8 → out of memory. Keep numeric lines bare; validate parsed
values by echoing them.

**P22 — Restart semantics.** (a) A fluid-field restart re-created the particles
from the input files — initial position, orientation, zero velocity — so every
chained segment re-inserted the sphere at rest; the torque trace looked
continuous across the seam, so "seam continuity is NOT evidence of state
continuity — check positions" (memory note, 2026-09-11; fixed 2026-09-15 by a
particle checkpoint written with every fluid dump). (b) Legacy dumps lacked one
field and poisoned the restart (NaN at the first residual). (c) Restart dumps
written only at output events: a wall-clock timeout wrote nothing (memory note
2026-09-04). (d) Hard-linked dump slots shared between cloned run directories
were overwritten in place through the shared inode (2026-09-01). (e) A chain
driver with a segment cap stopped one segment short when the end time was
raised (2026-09-23). Check: restart = fluid *and* particle state; copy, never
hard-link, restart files; dump on a wall-clock schedule as well as on steps;
restart tests that compare positions, not just forces.

**P23 — Staging by template.** Partition count must match the partition
directory (`SubMeshNumber`); boilerplate particle files must exist even when
superseded; closed boxes need the "no outflow" switch or the pressure
multigrid diverges at step 2; an application default made the box periodic
when a key was left at its sentinel value (segfault); run directories staged
under a build tree lost 60 GB of VTK (stage tool docstring). Memory items
3–4 and D6.4 G0 items 1–2. Check: stage by cloning a *running* certified case
and diff against it (guide §1: "Stage cases with `tools/stage_tencate_case.py`
… rather than hand-editing decks"; the tool applies count-checked
substitutions and refuses to stage into build trees).

**P24 — Output volume.** VTK on a once-refined output mesh is ~8× larger;
gate runs run with VTK off and particle logs on (stage tool defaults
`OutputFreq = 1d6`, `OutputLevel MAX-1`). L4 runs need ~60 GB RAM on 32 ranks
(`job.sbatch`), L3 settling boxes ≈ 3.7 GB/rank (memory note D6.4 item 3).

**P25 — Legacy ad-hoc wall-force code.** FeatFloWer carried a pre-campaign
`ENABLE_LUBRICATION` wall-sliding force from an earlier suspension study; the
owner ruled it dead code ("start clean", memory note
`enable-lubrication-legacy-kroupa`), and the ad-hoc near-wall thresholds
(5·h_min, a hard-coded 0.07) were replaced by the measured 2 h rule (guide §5).
Lesson: near-wall corrections must be calibrated against an exact solution
(Brenner) at your resolution, not inherited.

**P26 — Treat contact parameters as physics, not stabilisers (guide §6).**
Everything below ~2 cells of gap *is* the contact model; for liquid-immersed
contacts prefer near-zero friction or a lubricated contact. For this benchmark
(no rebound, St ≤ 4.13) a hard contact with the bottom plane sufficed.
