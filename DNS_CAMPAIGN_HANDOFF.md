# DNS campaign handoff — state as of 2026-09-10

Orientation document for anyone (human or review agent) joining the DNS
(fictitious-boundary, FBM) validation campaign cold. It says what exists,
what is settled, what is in flight, and the rules. It does not repeat the
numbers: every verdict lives in the datasheet, and this file points at rows
by their `case` id. Previous edition (2026-08-03) is in git history.

Branch `feature/dns-validation` (~250 commits since fork; the owner pushes,
nobody else). Companion: `EL_CAMPAIGN_HANDOFF.md` (Euler-Lagrange sibling,
untracked, 2026-08-13).

## 0. Read these, in this order

| # | document | what it is | freshness |
|---|---|---|---|
| 1 | `docs/md_docs/dns_validation_datasheet.md` (+ `.csv`) | **the authority.** 155 verdict rows: suite, case id, quantity, expected + source, measured, rel_error, tolerance, verdict. Row prose carries the reasoning; the last ~20 rows are the current state. The two files share the row set (case/quantity keys) but not the ordering or every word of historical prose; the CSV is what tools parse, the MD is what people read | current (rows `d61_review_corrections`, `d62_conventions_pinned`, 2026-09-10) |
| 2 | `applications/q2p1_dns_drag/validation_cases/*/` | per-family specs: `d11_hasimoto/RUNBOOK.md`, `d22_lubrication/DESIGN_SPEC.md`, `d51_viscometer/CASE_SPEC.md` + `CASE_SPEC_V2.md`, `d61_oberbeck/CASE_SPEC.md`, `d62_jeffery/CASE_SPEC.md`. Gates, geometry, run matrices, dated amendments | current for D5/D6 |
| 3 | `docs/md_docs/dns_torque_path_review.md` | FF traction/torque integration and pe torque application reviewed line by line; ellipsoid defects D-1..D-4 with fix records; notes N-1..N-4 | 2026-09-05 |
| 4 | `docs/md_docs/dns_practitioners_guide.md` v2.3 | the distilled rules: D/h, dt, force noise floor, near-contact crossover, contact parameters, reference discipline, array drag law, job costs | v2.3, 2026-09-10 (§12 Fritz + segmented restarts, §13 ellipsoids) |
| 5 | `dns-validation-campaign-plan.md` | the design document, stages and literature canon (959 lines) | 2026-08-03 design; **stage numbering differs from the datasheet, see §1** |
| 6 | `docs/artifact_snapshots/dns_campaign_artifact_live.html` | the dashboard page as last published (one-page status with figures) | 2026-09-02 |
| 7 | `docs/md_docs/dns_figures/` | committed analysis figures per family (`setup_scenes_units.md` documents the setup scenes) | rolling |
| 8 | `tools/` | `d61_oberbeck_analysis.py`, `d62_jeffery_analysis.py` (with `--plot`), `stage_tencate_case.py`, `el_slurm_submit*` | current |
| 9 | `build_guide.md`, `CLAUDE.md`, `docs/md_docs/featflower_partitioner_usage_guide.md` | building, PE modes, partitioning | stable |

Not for this review: the other root handoffs (`EL_CAMPAIGN_HANDOFF.md`,
`dashboard-handoff.md`, `benchmark-regression-handoff.md`, `mpi_prf_*`,
`stage0-acceptance-battery-plan.md`, `inspection-artifacts-handoff.md`) belong
to other threads; `PR_DESCRIPTION*.md` describe branch history up to August.

Dashboard artifacts (owner-published, results only):
DNS https://claude.ai/code/artifact/d2db32bb-9f15-4f07-9098-07d63032f4fc,
EL https://claude.ai/code/artifact/ff5bca87-80f1-4d22-b088-c89bf412ae01.

## 1. Stage crosswalk (plan vs datasheet)

The plan (Aug 3) and the datasheet `suite` column diverged when the
non-spherical work was pulled forward. Use the datasheet names.

| plan stage | datasheet suite | content |
|---|---|---|
| D0 | `d0_infra`, `d0_solver`, `infrastructure` | twins, pin bumps, guards, cluster bring-up |
| D1 | `d1_metrology` | Hasimoto ladder, periodic-coupling fix, ten Cate matrix, grid-crossing noise |
| D2 | `d2_pairs` | sphere-wall/sphere-sphere approach, lubrication (D2.2), DKT (D2.3) |
| D3 | `d3_collective` | random-array drag law (D3.1), hindered settling (D3.2) |
| D5 (EL cross-validation) | not started | see §4 |
| D6.1 Kroupa box (plan) | `d5_rheology` (d51/d52 viscometer) | suspension viscosity in a numerical viscometer |
| D4 non-spherical (plan) | `d6_nonspherical` (d61 Oberbeck, d62 Jeffery) | single ellipsoid, fixed then free rotation |
| D6.2 ATC (plan) | not started | capstone |

## 2. What is settled (headline per family, row ids)

- **D1.1b periodic coupling (the August headline).** Q2 vertex-share table
  had no periodic pairs, so periodic faces were solved traction-free.
  Fixed in `parentcomm.f90` (`d11b_periodic_comm`: 702 = 27x26 pairs on the
  torus). Hasimoto ladder now converges to the Stokes drag with O(h)
  interface widening (`d11_rh_collapse`); the phi-anomaly closed
  (`d11_phi_small_postfix`, -1.5%). Rule: periodic runs need axis-uniform
  Cartesian partitions, METIS is invalid.
- **D1.2 force noise floor is first order in h** (`d12_prescribed_ladder`:
  3.8/1.8/0.9% at D/h 11/24/49).
- **D1.3 ten Cate** matrix complete; pe/CFD stepsize mismatch found and
  fatal-guarded; no dt stability floor (the earlier floor was the desync).
- **D2.2 lubrication (Kroupa 2016 model ported into HardContactAndFluid,
  pe PRs #25/#26/#28/#29).** Deficit formulation verified against Brenner
  (`d22_g2b_deficit`, +7% band vs +6.7% predicted); ten Cate landing finite
  with lubrication on (`d22_g3_tencate`). Lubrication is sphere-only and
  the setup refuses it for ellipsoids.
- **D2.3 DKT** closed (`d23_result`): both contact models complete
  drafting-kissing-tumbling; the earlier frictional "stall" was the
  `hcaf_angvel_reset` defect (omega zeroed in integratePositions), removed.
- **D3.1 random-array drag law** (CHERD manuscript basis): FBM reproduces
  Beetstra at N=54 to +0.4% (`d31_n54_verdict`); all nine L5 refinement
  checks measured (`d31_l5_re9_p010`); matched-Re transport made
  multiplicative (`d31_matched_re_operation`); alpha-gradient functional
  carries a level-independent +2% overshoot (`fbm_functional_closure`).
- **D3.2 hindered settling**: the 6d-column exponent is a confinement law
  driven by collective return flow (`d32_wide_attribution`).
- **D5.1 numerical viscometer (d52, Fritz L4 rung).** Instrument: annulus
  torque T(0) at L4 = 83.77503 vs analytic, -0.001% (`d52_v24f_baseline_hr`;
  L3 was +0.54%). phi ladder at L3: 0.05 Einstein-consistent (+0.6%,
  `d52_v21_einstein`), 0.10 Batchelor (-0.8%, `d52_v22_phi10`), 0.20
  Krieger-Dougherty (-0.9%, `d52_v23_phi20`); sub-grid lubrication adds
  +0.76% / +2.71% at 0.10 / 0.20 (`d52_v22L_settled`, `d52_v23L_settled`).
- **D6.1 Oberbeck (fixed prolate spheroid, r_e=2) CLOSED**
  (`d61_v5_halfsize`): drag-ratio Y/X = 1.14148 vs 1.14532 (-0.34%) once
  the finite-cell lattice term is removed by the half-size rung; absolutes
  -1.1/-1.5% raw; torque nulls 1e-8. Full ladder in rows `d61_v0_anchor` ..
  `d61_v3b_transverse`. Four pe ellipsoid defects found by the review and
  by gate G0 (inertia 0.25->0.2, missing buoyancy seeding, volume 4/3,
  containsPoint z-term) - all fixed, unit-tested, twin-gated
  (`pe_ellfix_twin`, `pe_isell_twin`).
- **Restart-dump machinery repaired** (v24f seg-2 NaN x3): `BackUpFreq`
  counts OUTPUT FRAMES not steps; legacy dump wrote scratch as coordinates
  and no MaterialDistribution; `ProcCtrl` Dump_Out now uses the complete
  writer; coordinate sanity guards at both readers (the executed restart
  path is `init_sol_same_level`). Details: row `d52_v24f_baseline_hr`,
  memory `ff-deck-staging-pitfalls`.

## 3. In flight (2026-09-10)

| run | where | job | state | what closes it |
|---|---|---|---|---|
| d52 v25f, phi=0.05 at L4 (D/h=16 rung) | Fritz, 6 nodes | 4203909 (segment 5, t 240->250) | running, ~t=244 | eta_L4 = T(phi)/83.77503 over the t>=230 plateau vs composite Einstein; rung verdict row; decides whether phi=0.20 needs an L4 rerun |
| d62 V0b, r_e=1 spin control | Fritz | 4200094 | DONE, PASS (`d62_v0b_spin`, omega/gammadot = -0.50495, +0.99%) | - |
| d62 V1b, r_e=2 Jeffery orbit, t->120 | Fritz, 2 nodes | 4200095 (seg 1) + self-chaining `chain_v1b.sbatch` (4209002 queued afterany) | seg 1 at t~26, ~21 s/step, 3-4 segments | `tools/d62_jeffery_analysis.py <concatenated logs> --gammadot 0.2 --re 2 --tmin 0.5 --plot`; gates: T*gammadot = 15.708 +-3% from pi-crossings, 4:1 waveform +-5%, axis_y < 0.02; preview at t=25 already on the Jeffery curve (figure NOT committed yet, owner's call) |

D6.2 finding this week (row `d62_v01_rho1_unstable`, spec §2): at
rho_p = 1 the explicit torque exchange is unstable (rotational relaxation
time 0.004 t.u. < dt = 0.01, gain ~3.6, NaN by t=0.3). Cure: rho_p = 10
at unchanged dt (tau_rot = 0.04 t.u. against a 78.5 t.u. period, the
zero-inertia Jeffery limit is intact). Stands on its own evidence; it is
NOT the refuted ten Cate "dt floor" (that was PE/CFD desync).
Reviewed 2026-09-10 (`DNS_CAMPAIGN_REVIEW.md`, untracked): five findings,
all addressed the same day - see rows `d61_review_corrections` and
`d62_conventions_pinned`.

## 4. Open list (priority order)

1. D6.2: V1b verdict, then V2 clearance ladder (H=4 box, new mesh via
   `d62_mesh_v1/gen_d62_shearbox.py --h 4` + axis-uniform partition;
   the owner's earlier observation is that walls measurably slow the orbit,
   so V2 is REQUIRED before D6.2 closes); optional V3 log-rolling start.
2. v25f verdict; if eta_L4 - eta_L3 is large, the phi=0.20 rung at L4.
3. Website: D6 non-spherical results are HELD OFF the site until the owner
   declares the picture clear (owner directive 2026-09-04). Site = clean
   results only, ledger consumes datasheet rows verbatim.
4. pe housekeeping for the next push: `rotationOnly` commit is local on
   `feature/ellipsoid-rigidbody-fixes` (ahead 1); banner "volume fraction"
   assumes a unit cube (cosmetic); upstream defect list items 3 and 5
   (frozen-field test guard, checkpoint deadlock) still open.
5. D6.1 optional: half-size 45-degree variant (quantitative off-diagonal),
   r_e = 3.
6. Longer: D5 (DNS-EL twins, DNS-informed closures), the plan's D6.2 ATC
   capstone, Ding & Aidun 2000 for a finite-Re Jeffery ladder, the
   `VISC_TORQUE_DNA` -> `VOL` log-tag rename (routine already renamed
   2026-08-27, tag kept until D5.1 closes), FF PR to master.

## 5. Rules that bind (each has a memory file; violations cost jobs)

- **Pushing.** Claude never pushes FF or pe. Owner pushes; then the Fritz
  clone is updated with `git pull && git submodule update`. Commits carry
  `Co-Authored-By` + `Claude-Session` trailers.
- **Pin bumps.** Every `libs/pe` pin bump is gated by an e4_l3 bitwise twin
  (instrument logs `particle_force.log`, `bulk_flow.log` byte-identical)
  before any production run. Pin lineage: f9b7115 (PR #28) -> 6971b13
  (ellipsoid D-1/2/3) -> 7853aeb (D6.1 setup keys) -> 1b0dda2 (D-4
  predicates) -> 6dc2261 (isTypeEllipsoid, current FF pin).
- **Builds.** DNS builds need `-Dpe_CONSTRAINT_SOLVER=pe::response::HardContactAndFluid`
  (the default silently zeroes FBM dynamics), USE_PE + USE_PE_SERIAL_MODE +
  ENABLE_FBM_ACCELERATION via the 3-step cmake, SED_BENCH=OFF. Frozen
  local trees: `build-dns-pe-serial/-parallel/-lubaddon` read-only.
- **Decks.** Clone a certified running rundir, never a template. Inline
  comments on numeric lines poison the parser; `SimPar@TimeStep` must equal
  json `stepsize_`; `SubMeshNumber` = partition count; full-domain
  `ForceScale` = 1,1,1,1,1,1; `Prop@Viscosity` is KINEMATIC; `start/`
  boilerplate must exist; `BackUpFreq` counts output frames (set 1 with
  `OutputFreq` = dump cadence); copy restart dumps with `cp -a`, never
  hardlink; a segmented trajectory never restarts across machines.
- **Fritz** (NHR@FAU): `SSH_AUTH_SOCK=/run/user/3086/ssh-agent.socket ssh fritz`
  (owner unlocks the agent per boot). Workspace
  `/anvme/workspace/k115ce12-featflower`: `FF/FeatFloWer` clone (at the
  pushed tip plus explicitly listed patches), `mesh_repo`, `rundirs/`,
  `build-dns-gcc14` (gcc 14.2 + OpenMPI 4.1.8, certified by
  `fritz_build_twin`). Whole 72-core nodes are billed; 24 h wall limit
  everywhere, so long runs are segmented through the dump slots.
- **Mesh dirs and rundirs stay untracked**; figures for analysis are
  committed under `docs/md_docs/dns_figures/`, presentation renders and
  `.pvsm` scenes are not.
- **Datasheet discipline.** One row per verdict, appended to both `.md` and
  `.csv`, verdict vocabulary PASS / FAIL / RECORDED / RESOLVED (plus
  qualified forms such as FAIL-diagnosed); a superseded row is never
  deleted, the superseding row names it.
