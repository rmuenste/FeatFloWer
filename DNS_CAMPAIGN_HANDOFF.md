# DNS campaign handoff — state as of 2026-08-03

Companion to `EL_CAMPAIGN_HANDOFF.md`. Branch `feature/dns-validation`
(~60 commits, owner pushes). Records: `docs/md_docs/dns_validation_datasheet.md`
(+.csv, authoritative), `docs/md_docs/dns_practitioners_guide.md`,
`dns-validation-campaign-plan.md` (design, stages D0–D6). Artifacts:
DNS https://claude.ai/code/artifact/d2db32bb-9f15-4f07-9098-07d63032f4fc (🌊),
EL https://claude.ai/code/artifact/ff5bca87-80f1-4d22-b088-c89bf412ae01 (🫧).

## Done

- **D0**: twin gates, runtime keys replace SED_BENCH, HardContactAndFluid
  solver selection (memory: HCEL zeroes FBM), h_min fix, viscosity fix,
  sawtooth watchdog (validated both directions), PE/CFD stepsize fatal
  guard (validated), staging tool `tools/stage_tencate_case.py`.
- **D1.3 complete**: E1–E4 × L2/L3/L4 at dt=1ms vs printed Table II
  ratios (digitized E1/E2 curves are +3.4/+3.9% wrong — tc-ref/README).
  Spatial error geometry-dominated; finest configs match ten Cate's own
  LBM to +0.03…+0.9% (E1's −5% is sim-vs-experiment, in the 2002 paper
  itself). **pe_stepsize_mismatch** found+fixed (serial PE steps with
  json `stepsize_`, ignores CFD dt; guard in `fbm_main.f90`): all old
  dt≠1ms rows superseded; synced ladder L3 +0.81/+1.41/+1.94%
  (1/0.5/0.25ms), **no stability floor** (earlier floor = 4× desync
  artifact), sub-linear temporal order ~0.5–1.
- **D1.1 Hasimoto + D1.1b root cause (the headline)**: drag ladder
  converges (order ~1.7) to **−8.9% asymptote** — wrong answer. Code
  trace + u_x plane test (ratio 0.102, job 137506) prove: **Q2 comm has
  NO periodic pairs**; `myVERTLINK`/`VerticeCommunicationScheme`
  (`source/src_mpi/parentcomm.f90:496-586`) built with non-periodic
  octree search, gates all E013 exchanges + parallel pressure → periodic
  faces solved **traction-free**. Face-level (mg_mpi) matching works.
  Also: 2× periodic face double-count (`parentcomm.f90:787-844`), MUMPS
  label-swap (`Create_GlobalNumbering`). **Affects every "periodic" run
  in BOTH campaigns** (EL momentum box, RZ settling, kroupa xy) —
  re-interpretation pending. Partitioning rule (owner): periodic runs
  need axis-uniform Cartesian partitions; METIS invalid.
- **D2.3 DKT wired**: `applications/q2p1_dkt` + json-driven
  `setupDraftKissTumbSerial` (libs/pe branch `dkt-setup`, pushed;
  pointer bumped). Mesh `dkt_mesh_v1/` (owner, untracked; committed
  README+REPRODUCE). D/h=8: drafting ✓ kissing ✓ (t_kiss≈18.1);
  offset run: tilt onset then **stall at 4.2°** in exact contact —
  candidates: contact friction, D/h=8 wake, time. dt-floor probes at
  ratios 1.10/1.02: clean (no floor).

## Open / next

1. **Fix the periodic Q2 coupling** (D1.1b): periodic search in
   myVERTLINK construction; acceptance test = existing 4-run Hasimoto
   ladder → 1.832. Then re-measure r_h curve (much of it was
   confinement shortfall), revisit φ-anomaly (likely closes), revisit EL
   periodic-run interpretations (incl. RZ FAIL attribution).
2. DKT D/h=16 rung (mesh ready, needs checkpoint cadence set) — tumble
   stall discriminator; frictionless-contact probe if stall persists.
3. D1.2 grid-crossing noise; E2/E3 dt points; practitioner guide v1.
4. Housekeeping: ~66G buried EL rundir viz deleted, logs zstd'd;
   `.gitignore` guard; disk shared/volatile — keep viz off for gates.

## Traps (memories exist for all)

Deck/json must agree: viscosity+density (d0_visc), **TimeStep==stepsize_**
(fatal-guarded), solver=HardContactAndFluid compile flag, no inline
comments on numeric deck lines, `_mesh/<MeshFolder>` + `_data/MG.dat` +
`start/` staging trio, ForceScale 4× only for quarter-box, D/h counts
elements (nodal = 2×; halve LBM/IBM literature values).
