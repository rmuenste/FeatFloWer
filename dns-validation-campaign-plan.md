# DNS (FBM) validation campaign — design document

Drafted 2026-07-31. Companion to the completed EL campaign
(`docs/md_docs/el_validation_report.md`); reuses its discipline, tooling and
report/datasheet/RUNBOOK format. Target end state: the resolved-particle DNS
method (Fictitious Boundary Method + PE rigid-body backend) is understood,
documented and validated to the same depth as the EL method — so both can be
used selectively, and DNS can serve as in-house ground truth for EL closures
("DNS-informed EL").

The deeper motivation (R. Münster, 2026-08-01): pre-campaign DNS studies
with this codebase — the Kroupa-box effective-viscosity study and the
Archimedes-tube-crystallizer slug — *did* produce correct physics, but
only through painful trial-and-error on resolution, dt, collision
handling and meshing. The campaign's yardstick is that such cases become
set-up-able in an informed, structured way from written guidelines
(D1–D3), first-time-right at predictable cost. That claim is tested
literally: both pain-cases return as the D6 capstone ("Endgegner")
stage.

---

## 0. What we are validating (method map, from code survey 2026-07-31)

The DNS method in FeatFloWer is the Turek/Wan FBM:

- **Representation**: per-Q2-DOF 0/1 indicator `FictKNPR`, rebuilt from
  scratch every step by geometric point-in-body queries
  (`source/src_quadLS/QuadSc_boundary.f90:60` → `fbm_updateFBMGeom`,
  `source/src_fbm/fbm_main.f90`). PE backend classifier supports **spheres,
  capsules, boxes, ellipsoids, cylinders, triangle meshes**
  (`libs/pe/src/interface/object_queries.cpp`). Optional HashGrid
  acceleration (`SimPar@UseHashGridAccel`).
- **No-slip**: internal Dirichlet by row filtering — solution values set to
  the rigid-body field v = V + ω×(x−x_c), defect zeroed, matrix rows
  unit-diagonalized (`Boundary_QuadScalar_Val/Def/Mat`,
  `QuadSc_boundary.f90:542-728`). No penalty, no forcing term.
- **Hydrodynamic force/torque**: volume integral of (−pI + μ(∇u+∇uᵀ))·(−∇α)
  over the interface band (`ForcesLocalParticles`,
  `source/src_quadLS/QuadSc_force_extension.f90:24`; KVEL-accelerated serial
  variant in `QuadSc_force_serial.f90`).
- **Rigid-body update + collisions**: PE `HardContactEulerLagrange` (PGS,
  velocity-level contacts) with optional Kroupa lubrication; serial-PE mode
  (per-rank full-domain PE world, forces summed over the CFD comm) or
  parallel-PE mode (rank 0 excluded, Cartesian decomposition, shadow
  copies). FullC0ntact is a second backend (not in scope for this campaign
  except where an app forces it).
- **Time stepping**: fixed dt from `SimPar@TimeStep`. `TimeStepCtrl` exists
  but every DNS call site is commented out; `bGapAdaptiveTimeStep` is a
  declared, unused stub. `ComputeCFL` prints fluid + particle CFL but
  nothing acts on it.
- **Mesh**: fixed uniform multigrid hierarchy (1:8 refinement,
  `SimPar@MinMeshLevel/MaxMeshLevel`). No dynamic refinement around
  particles; r-adaptivity machinery (`UmbrellaSmoother` with
  distance-to-body monitor) exists but is not wired into the DNS loop.

Known validation state going in:

- ten Cate E4 sedimentation: serial −0.132921 / parallel −0.132877 m/s
  (0.033% apart), +3.8% vs the 0.128 published value (guides 02/07). The
  committed regression baseline (`tools/featflower_test/testcases/`) is a
  **self-consistency** check against an under-resolved L2 run whose source
  file (`E4_lvl2.txt`) is not in the repo. The "L3 within ~2% of PIV" claim
  has **no in-repo evidence** (missing `tc-ref/` PIV CSVs, missing
  guide_02_part2).
- DFG cylinder (fixed body): C_D 5.5795 / C_L 0.010619 baselines, green.
- Nothing else: no multi-particle physics benchmark, no fluidization physics
  gate (guide 06 is compile-only), no DKT app (PE has `setupDraftKissTumb`
  but no deck), no drag-correlation comparison in `q2p1_dns_drag`, no
  non-spherical validation, no resolution/dt guidance anywhere.

---

## 1. Campaign principles (inherited from the EL campaign, binding)

1. **Probe-first**: one short smoke gated before any sweep. Every case has a
   RUNBOOK.md (design → gates → job numbers → verdicts); every quantitative
   claim becomes a datasheet row (expected / source / measured / tolerance /
   verdict). FAILs recorded honestly with attribution, never tuned around.
2. **Independent estimators beat tighter tolerances** (the EL Couette
   lesson): wherever possible each stage carries two routes to the same
   number (e.g. force-integral vs momentum-balance; wall force vs bed
   weight).
3. **Bitwise twins + A/B worktrees** settle "pre-existing vs ours".
4. **Master/worker asymmetry discipline**: rank 0 has no PE in parallel
   mode; diagnostics must be collective-safe in rank-asymmetric paths.
5. Mesh fixtures need `git add -f` (`.gitignore` excludes `*.tri`/`*.par`);
   per-case pinned `example.json` (never the shared
   `libs/pe/pe/interface/example.json` scratch file).
6. The user pushes; agent commits only. FF commits end with the standard
   Co-Authored-By line. Never push `libs/pe`.

Campaign artifacts (mirror the EL set):

- `docs/md_docs/dns_validation_report.md` — the science record.
- `docs/md_docs/dns_validation_datasheet.md` + `.csv` — same columns as EL
  (`suite | case | quantity | expected | expected_source | measured |
  rel_error | tolerance | verdict`).
- Per-case RUNBOOKs under `applications/<app>/validation_cases/<case>/`.
- **`docs/md_docs/dns_practitioners_guide.md`** — the usability deliverable:
  resolution guidelines, dt guidelines, mode selection (serial vs parallel
  PE), lubrication-threshold rules, staging/submission recipes. This is a
  first-class output of Stages 1–2, not an afterthought.
- Dashboard extension (`docs/dashboard/`) once the datasheet exists.

---

## 2. Stage overview

| Stage | Name | Depends on | Headline deliverable |
|---|---|---|---|
| D0 | Infrastructure, recertification, usability debt | — | Reproducible ten Cate truth chain; unified diagnostics; generalized staging tools |
| D1 | Single-sphere metrology (resolution + dt guidelines) | D0 | D/h and dt guidelines with error bars; practitioner's guide v1 |
| D2 | Pair interactions (approach, lubrication ladder, DKT) | D1 | DKT app + validated pair dynamics; resolved-vs-modeled gap crossover map |
| D3 | Collective: fixed random arrays → hindered settling → fluidization | D1 (D2 for contact-heavy runs) | DNS drag law vs Beetstra/Tenneti; wall-bounded hindered settling; fluidization pressure-drop gate |
| D4 | Non-spherical single particle | D0 (parts of D1) | Jeffery orbit + spheroid settling validated; non-sphere diagnostics fixed |
| D5 | DNS↔EL cross-validation, DNS-informed EL | D1–D3 | Closure-comparison memo; DNS-measured hindrance correction for EL |
| D6 | Capstone ("Endgegner") re-runs: Kroupa box, ATC slug | D1–D3 (D2 for contact-heavy regimes) | Both prior pain-cases reproduced with setup taken straight from the practitioner's guide |

D2, D3, D4 are parallelizable after D1; D5 consumes their outputs; D6 is
the closing acid test of the guides themselves. Each
stage below lists prerequisites, cases, gates, deliverables, and an exit
gate that must be green (or FAIL-attributed) before dependent stages start.

---

## 3. Stage D0 — Infrastructure, recertification, usability debt

Goal: make DNS runs stageable, observable and reproducible the way EL runs
are, and put the existing single-sphere claim on defensible footing.

### D0.1 Tooling generalization
- Generalize `tools/el_stage_rundir.sh` / `el_slurm_submit.sh` into
  app-agnostic `tools/ff_stage_rundir.sh` / `ff_slurm_submit.sh` (or add an
  `--app` flag): DNS cases stage a different fileset (`q2p1_param.dat`,
  pinned `example.json`, `cube.json`/`sampleRigidBody.xml` where needed,
  `_adc`/mesh partition tree) and variable rank counts (np = px·py·pz + 1
  in parallel-PE mode, not always 28).
- Adopt the maintained `featflower-partition` (METIS 5) for graph
  partitions; document and wrap the legacy `tools/partpy` `-4` Cartesian
  path (required whenever parallel PE is used, since PE's `MPI_Cart_create`
  decomposition must align spatially with the CFD partition). Decide and
  record: migrate `-4` into `featflower_partitioner` vs. pin the legacy
  script with a usage note.

### D0.2 Diagnostics unification (grep-able vocabulary, both PE modes)
- Today the trajectory channel differs by mode: `SED_BENCH_VEL` (serial
  only, `QuadSc_force_serial.f90`) vs PE's `Position:`/`Velocity:` prints
  (parallel). Introduce a unified, collective-safe tag set printed
  identically in both modes, e.g.:
  - `DNS_PART_STATE t= id= x y z vx vy vz |omega|` (per particle, or capped
    subset for large N)
  - `DNS_PART_FORCE t= id= Fx Fy Fz Tx Ty Tz` (post-`sync_forces`)
  - `DNS_RESOLUTION t= dofs_per_particle_global h_min D_min/h` — fix the
    parallel Dofs-per-Particle print (currently rank-1-local by design;
    make it a proper reduction on a rank-symmetric path)
  - `DNS_CFL t= cfl_fluid cfl_particle gap_min/h` (extend `ComputeCFL`)
- These become the parse targets for all later gates (same pattern as
  `EL_MEAN_SLIP` etc.).

### D0.3 Retire per-case compile-time hacks
- `-DSED_BENCH` (`QuadSc_force_extension.f90:424`: zeroes Fx/Fy, ×4 Fz) must
  become runtime config — `Prop@ForceScale` already exists and expresses
  exactly this. One campaign binary must be able to run every case.
  Bitwise/tolerance twin against the flag-built binary on the 10-step smoke
  before retiring the flag.
- Pin per-case `example.json` in each case dir (disease #7 from the EL
  campaign); stop staging the shared scratch copy.

### D0.4 ten Cate truth chain (recertification)
- ~~Digitize the ten Cate PIV trajectories~~ **DONE 2026-07-31**: digitized
  E1–E4 curves provided by the user and committed as `tc-ref/`
  (`case_E*_h.csv` = t vs gap h/d; `ref_E*.dat` = t vs v_z [m/s]; format
  and sanity checks in `tc-ref/README.md`). This repairs the dangling
  references in
  `pipemesh_v1/handoff_euler_lagrange_drag_validation_ten_cate.md`.
  Sharpened gate: the confined-PIV E4 peak is −0.1230 m/s, so the
  documented serial DNS −0.1329 is **+8.1% vs PIV** (not +3.8% vs the
  unbounded u_∞ = 0.128) — the honest baseline the D1 resolution matrix
  must explain or close.
- Commit the mesh fixtures this campaign owns (`git add -f`):
  `benchSym/` quarter-box + `benchSym/mesh12/NEWFAC` (1×1×12 Cartesian),
  and the full-box `ten_cate_mesh_v1/` (already mirrored in the EL v1b
  case). Record explicitly that quarter-box and full-box are different
  discretizations of the same experiment.
- Re-run E4 serial and parallel with the unified diagnostics; regression
  baseline moves from "match an undocumented L2 run" to "match the pinned
  L3 trajectory + PIV within stated tolerance".

### D0.5 Reproducibility floor
- Re-green the two existing featflower_test cases on the campaign branch.
- One bitwise twin: 10-step sedimentation smoke, campaign binary vs current
  branch head, byte-identical `SED_BENCH_VEL`/`DNS_PART_STATE` lines (after
  D0.3 this becomes the walls-off-style regression anchor).

### D0.6 DNS constraint-solver decision (architecture gate, agreed 2026-07-31)
- **DNS and EL do not share a CollisionSystem.** `HardContactEulerLagrange`
  is EL-specific at its core (hydro forces enter as `dv_`/`dw_` velocity
  corrections inside the PGS sweep; the FBM path instead injects resolved
  forces via `setForcesMapped` through the standard force→integration
  pathway). DNS gets its own solver — a dedicated `HardContactDNS`,
  starting from `HardContactSemiImplicitTimesteppingSolvers`,
  `HardContactLubricated` or `HardContactAndFluid`; the choice is recorded
  as a short design note in the report.
  *Empirical confirmation (2026-07-31, campaign day one)*: a fresh build
  tree inheriting pe master's default (`HardContactEulerLagrange`) produced
  an all-zero sedimentation trajectory — force and velocity exactly 0 for
  10 steps, no error raised — because `setForcesMapped` forces are a dead
  end in that solver. Interim campaign selection until `HardContactDNS`
  exists: `HardContactAndFluid` via
  `cmake -DCMAKE_CXX_FLAGS="-Dpe_CONSTRAINT_SOLVER=pe::response::HardContactAndFluid"`
  (the define is `#ifndef`-guarded; no submodule edit needed). The
  zero-trajectory signature is now a documented failure mode.
- **Share components, not the class**: (a) the Kroupa lubrication closure
  functions become a solver-agnostic pure-math component both solvers call
  (with different policy — see D2.2); (b) the accumulator/audit machinery
  (contact virial, per-wall impulses, Newton-pair bookkeeping) becomes
  free-standing, exposed through the `pe/interface/el_optional_api.h`
  detector pattern so the Fortran query chain works against whichever
  solver implements it. The detector refactor at the end of the EL branch
  is exactly what makes this cheap.
- **Equivalence twin before any campaign numbers**: E4 sedimentation
  (10-step smoke + full trajectory) under `HardContactEulerLagrange` vs
  the DNS solver, same inputs. Bitwise identity is not expected across
  solvers — the gate is trajectory agreement within the D1 noise floor,
  with the wall-contact phase compared and RECORDED.
- Operational consequence: one build = one `pe_CONSTRAINT_SOLVER`, so DNS
  and EL campaign binaries come from separate build trees; every D5
  cross-comparison row records which binary/build produced each side.

**Exit gate D0**: staged-from-scratch E4 run (both PE modes) reproduces the
guide 02/07 numbers through the new tooling and unified diagnostics;
datasheet rows exist for serial/parallel cross-check (≤0.1%) and vs-PIV
(tolerance set after D1 resolution study — provisionally RECORDED); DNS
solver decision recorded and the cross-solver equivalence twin green.

Cost: small — mostly code/tooling; runs are the 36 s smoke and two ~1 h E4
trajectories.

---

## 4. Stage D1 — Single-sphere metrology: resolution & dt guidelines

Goal: answer "how many cells per diameter, what dt, and what error do I get
for it" — quantitatively, per Re regime. This stage produces the guidelines
the user asked for and calibrates every later stage's cost model.

### D1.1 Fixed-sphere resolution ladder (no moving boundary — cleanest signal)
- Case: periodic simple-cubic sphere array at low φ, Stokes/low-Re flow,
  `q2p1_dns_drag` machinery (serial-PE only today — acceptable; note as
  constraint). Reference: Hasimoto (1959) dilute-array expansion; Zick &
  Homsy (1982) / Sangani & Acrivos (1982) tables for general φ.
- Sweep mesh level (D/h ladder, measured by the fixed `DNS_RESOLUTION`
  diagnostic) at fixed φ ∈ {0.05, 0.2}; observable: drag vs reference.
- Deliverable: error-vs-D/h convergence curve for the α-integral force on a
  *stationary* interface; the asymptotic order of the FBM force (expected
  ~1st order in h for staircase indicator methods) measured, not assumed.

### D1.2 Grid-crossing noise (translating sphere, frozen dynamics)
- Case: sphere dragged at constant velocity across the mesh
  (`SimPar@skipFBMDynamics = Yes`, prescribed motion), uniform flow frame.
- Observable: force time-series oscillation amplitude vs D/h as the
  indicator re-classifies DOFs. This characterizes the intrinsic FBM force
  noise floor — needed to set honest tolerances for every moving-particle
  gate later (and explains trajectory jitter practitioners see).

### D1.3 ten Cate E1–E4 resolution × dt matrix (the core study)
- Full-box mesh (`ten_cate_mesh_v1`), all four regimes (Re 1.5 / 4.1 /
  11.6 / 31.9), levels L2/L3/L4, dt ladder ×{2, 1, 1/2, 1/4} around the
  guide value at the middle level.
- Observables per run: full v(t) trajectory vs PIV (peak velocity, time of
  peak, approach-to-wall deceleration), plus `DNS_CFL` statistics.
- Gates: L-convergence monotone toward PIV; quantify Δu_t(D/h, dt, Re).
  Bottom-wall approach phase analyzed separately (it is lubrication/contact
  dominated — belongs to D2's crossover map, recorded not gated here).
- Cost estimate before launch (probe-first): one L4 probe run timed before
  committing to the full matrix; expected O(20–30) jobs, minutes (L2) to
  ~day (L4 E1, viscous time scale) each.

### D1.4 Free-fall cross-check at moderate Re (optional, decision point)
- Uhlmann & Dušek (2014) sphere-settling regimes give published DNS with
  documented resolution requirements (D/h up to ~24) — one matched case
  (e.g. Ga ≈ 144, vertical regime) would anchor the guidelines against an
  independent DNS. Include only if the D1.3 cost model says an adequate
  mesh fits the cluster (decide at stage midpoint).

**Deliverables D1**: `dns_practitioners_guide.md` v1 with (a) recommended
D/h per target accuracy and Re band, (b) dt rule (particle CFL bound +
grid-crossing rate), (c) measured FBM force convergence order and noise
floor, (d) cost-per-(D/h, level) table for sizing future runs. Datasheet
rows for every matrix cell (PASS against convergence-model prediction, or
RECORDED).

**Exit gate D1**: E4 at recommended resolution within stated tolerance of
PIV with an error bar that the D1.1/D1.2 analysis explains; guidelines
published.

---

## 5. Stage D2 — Pair interactions: approach, lubrication ladder, DKT

Goal: validate two-body hydrodynamics and map where resolved hydrodynamics
hands over to the modeled (PE lubrication / hard-contact) layer — the DNS
analogue of the EL campaign's lubrication arc, and the direct bridge to it.

### D2.1 Sphere–wall approach (resolved lubrication ladder)
- Case: sphere settling onto a plane wall (the tail of every ten Cate run,
  isolated): prescribed approach and free approach variants.
- Reference: Brenner (1961) / Cooley & O'Neill (1969) drag divergence
  F ∝ 1/gap for gap → 0.
- Sweep D/h; observable: measured drag amplification vs gap/h. Deliverable:
  the **crossover map** — for given D/h, below which gap/h the resolved
  force falls off the analytic curve and the PE lubrication correction
  (threshold currently ad-hoc: `5·h_min` computed but the setter is
  commented out in `q2p1_bench_sedimentation/app_init.f90:121`; SRR
  hard-codes 0.07) must take over. Output: a defensible
  lubrication-threshold rule as a function of h — closing a loop the EL
  campaign left open (ε_c sensitivity).
- Re-enable/gate the `set_lubrication_threshold_from_mesh` path per this
  rule.
- **Legacy code note (do not reuse):** the `#ifdef ENABLE_LUBRICATION`
  blocks (`sliding_wall_force`, `total_lubrication`) in
  `QuadSc_force_serial.f90`/`QuadSc_main.f90` are leftovers of a
  pre-campaign Kroupa-box effective-viscosity study (rough preliminary
  wall-sliding forces feeding the total-stress → effective-viscosity
  path; results matched Krieger–Dougherty, but the implementation is
  explicitly not to be revived). Any D2 near-contact treatment starts
  clean from the Brenner/Cooley–O'Neill gates. HCAF itself has no
  lubrication (API stubs return 0/no-op).

### D2.2 Sphere–sphere approach
- Same protocol, two spheres, head-on and shear (tangential) motion, vs
  the two-sphere lubrication asymptotics (Jeffrey & Onishi resistance
  functions; same Kroupa closure family the PE layer implements —
  consistency check DNS ↔ PE model in the resolved regime, divergence
  expected only below the crossover gap).
- **Shared-closure extraction (code work, implements D0.6)**: extract the
  Kroupa closure functions into the solver-agnostic component both
  `HardContactEulerLagrange` and the DNS solver call. DNS policy differs
  from EL policy in two ways: the correction activates only below the
  mesh-resolvable gap identified by the D2.1 crossover map (threshold from
  h, replacing the ad-hoc `5·h_min`), and it applies only the *unresolved
  remainder* — never double-counting the portion of the lubrication force
  the FBM integral already resolves. EL policy is unchanged; shared
  coefficients guarantee the D5 comparisons hold the physics contract
  fixed.

### D2.3 Drafting–kissing–tumbling (the qualitative flagship)
- Prerequisite build-out: PE has `setupDraftKissTumb`
  (`commf2c_dkt_`, see `docs/md_docs/pe_initialization.md:64`) but **no
  application, deck, or mesh**. Deliverable: `q2p1_dkt` case (likely as a
  `q2p1_fc_ext`/`q2p1_bench_sedimentation` variant rather than a new app —
  decide during implementation) with committed mesh + pinned config.
- References: Fortes, Joseph & Lundgren (1987) experiments; 3D sphere DKT
  computations (Glowinski et al. 2001; Apte, Martin & Patankar 2009 —
  select the primary quantitative reference during case design and pull
  the PDF into the repo root for the dashboard).
- Gates (robust to chaos — DKT is symmetry-breaking): phase sequence
  observed (drafting → kissing → tumbling), kissing onset time and
  trailing-sphere approach velocity within literature band; post-tumbling
  divergence RECORDED qualitatively, not gated numerically.
- Two-mode discipline: run serial-PE and parallel-PE twins; agreement
  through kissing (pre-contact chaos amplification) is the internal gate.

**Exit gate D2**: crossover map published in the practitioner's guide;
DKT case green on phase-sequence + kissing-time gates; lubrication
threshold rule adopted (or divergence FAIL-attributed).

Cost: moderate — approach ladders are short runs; DKT trajectories are
E4-scale ×2 spheres, O(10) jobs.

---

## 6. Stage D3 — Collective behavior: arrays, hindered settling, fluidization

Goal: multi-particle validation at increasing coupling strength, with the
EL campaign's RZ failure informing the experimental design (wall-bounded,
not tri-periodic, for settling).

### D3.1 Static random arrays — DNS drag law (highest value/cost ratio)
- Case: frozen random sphere arrays in the periodic unit cell
  (`q2p1_dns_drag` + a seeding tool for random non-overlapping packs; the
  current setup enforces 1-diameter gaps — relax to literature-standard
  minimum gaps, noting the resolution requirement between close spheres
  from D2).
- Sweep φ ∈ {0.1, 0.2, 0.3, 0.4}, Re ∈ {≈0.2, 10, 50}, ≥3 random seeds
  each (report seed spread as the error bar).
- Reference: Beetstra, van der Hoef & Kuipers (2007); Tenneti, Garg &
  Subramaniam (2011); Ergun/Wen–Yu as engineering anchors.
- Gate: mean drag within the spread of the two DNS correlations (they
  differ by up to ~15% — the gate is the band, not one curve).
- **This is the primary DNS-informed-EL input** (see D5): the same F(φ, Re)
  table that validates DNS calibrates/replaces Di Felice hindrance in EL.

### D3.2 Wall-bounded hindered settling
- Case: N ≈ 20–50 spheres settling in a closed/walled column (sized by the
  D1 cost model; resolution per D1 guidelines is the feasibility limit —
  state N honestly rather than under-resolving).
- Reference: Richardson–Zaki n(Re) band; qualitative comparison with the
  EL v2_rz_settling result — DNS in the walled column is expected to
  recover hindrance where tri-periodic EL showed cluster-induced
  enhancement. Even a two-point φ sweep discriminates.
- Gates: ⟨v_settle⟩/u_t decreasing in φ (mechanism), exponent within RZ
  band (quantitative, tolerance set generously and honestly given N and
  box effects); volume/mass conservation and force-symmetry audits at
  machine precision (build the DNS analogue of the Newton-pair audit:
  Σ particle hydrodynamic force vs fluid momentum + wall stress budget —
  an independent-estimator pair as per principle 2).
- Explicit decision point at design time: if the D1 cost model says
  N ≥ 20 at adequate D/h exceeds the cluster, descope to N = 8–16 and
  reframe the quantitative RZ row as RECORDED (mechanism gate stands).

### D3.3 Fluidization (stretch tier — go/no-go after D3.1)
- `q2p1_bench_fluidization` exists but has never had a physics gate.
- Minimal defensible physics gates that don't require full bed statistics:
  (a) fixed-bed pressure drop vs Ergun at the D3.1-validated resolutions;
  (b) Δp ≈ bed weight/area at incipient fluidization; (c) minimum
  fluidization velocity vs Ergun-derived u_mf within a stated band.
- Serial-PE SRR variant (hard-coded 0.07 lubrication threshold) is
  reviewed against the D2 threshold rule first.
- Go/no-go: only if D3.1 shows the required per-particle resolution is
  affordable at bed scale; otherwise RECORDED as future work with the cost
  analysis as the artifact.

**Exit gate D3**: D3.1 datasheet rows complete (this is the mandatory
core); D3.2 mechanism gate green; D3.3 verdict (results or attributed
descope) recorded.

Cost: D3.1 is the big-ticket item — O(36) runs (4φ × 3Re × 3 seeds) but
steady-state, so short walltimes at fixed arrays; D3.2 is a few long jobs.

---

## 7. Stage D4 — Non-spherical particles, single body first

Goal: certify the already-present non-spherical FBM path and give it
diagnostics, then validate against the two cleanest analytic references.

### D4.1 Non-sphere plumbing audit (prerequisite, code work)
- The classifier and force integral are shape-agnostic, but: `ComputeCFL`
  reads `theParticles(1)%radius` (sphere-centric), several PE query helpers
  are sphere-only, lubrication filters non-spheres out, and there is no
  orientation output. Deliverables:
  - equivalent-radius / bounding-radius handling in CFL + diagnostics,
  - `DNS_PART_ORIENT t= id= q0 q1 q2 q3` (quaternion trace) in the unified
    diagnostic set,
  - a classification unit smoke per shape (ellipsoid, capsule, box):
    indicator volume vs analytic body volume at the D1 resolutions
    (volume-conservation gate, machine-precision style).

### D4.2 Jeffery orbit (prolate spheroid in Couette shear)
- Reference: Jeffery (1922) — orbit period T = 2π(r + 1/r)/G for aspect
  ratio r; analytic, parameter-free. The moving-wall Couette machinery from
  the EL campaign (`Inflow300` BC in `QuadSc_user.f90`, kroupa_couette mesh
  retagging recipe) is directly reusable for the DNS box.
- Torque-free rotation: gates on orbit period (few %, tolerance informed by
  D1 noise floor at the chosen a/h), orientation trajectory shape, and
  drift over ≥2 periods (RECORDED — FBM staircase torque noise will drift
  the orbit; measuring that drift rate is itself a guideline output).
- Sweep aspect ratio r ∈ {2, 4} and one resolution ladder point.

### D4.3 Settling spheroid (Stokes and inertial)
- Stokes tier: Oberbeck (1876) / Happel & Brenner analytic drag for prolate
  and oblate spheroids, broadside and edgewise — terminal velocity gate at
  low Re in a large box (confinement correction applied or box-size ladder
  RECORDED).
- Inertial tier: orientation stability — at finite Re a settling spheroid
  turns broadside-on; gate is the mechanism (correct attractor from a
  tilted release), with the transient RECORDED against literature DNS
  (reference selection during case design; pull the PDF).
- `q2p1_creep` (ellipsoid.obj, `setupCreep`) is the starting app; decide
  during design whether to run OBJ triangle mesh vs analytic
  `Ellipsoid` PE body — **run both once**: agreement between the analytic
  and triangulated classifier for the same spheroid is itself a
  cheap, valuable gate on the mesh-particle path (which q2p1_span/drill
  production cases rely on).

**Exit gate D4**: Jeffery period + Oberbeck terminal velocity green;
orientation-attractor mechanism green; non-sphere diagnostics merged;
practitioner's guide gains a non-spherical section (resolution measured per
minor axis, not diameter).

Cost: small-moderate — single-particle runs at D1-informed resolutions.

---

## 8. Stage D5 — DNS ↔ EL cross-validation and DNS-informed EL

Goal: the strategic payoff — use each method where it is strong, and use
DNS to close EL's documented model boundaries.

### D5.1 Same-physics twins (method-comparison memo)
- ten Cate E1–E4: DNS (this campaign, full-box + quarter-box) vs EL v1b
  (existing) vs PIV — one figure per case, with the EL self-voidage bias
  (ε_eff ≈ 0.975, +6% co-flow) annotated. States explicitly which method
  to use at which D/h and cost.
- Sphere–wall approach: DNS resolved curve (D2.1) vs EL's behavior (EL has
  no particle-wall lubrication in pipe geometry — documented boundary);
  quantifies what EL misses near walls.

### D5.2 DNS-informed EL closures (the deliverable that changes EL)
- Drag/hindrance: fit F(φ, Re) from D3.1 and compare against EL's Di Felice
  ε^−χ hindrance over the same grid. If deviation exceeds the EL
  campaign's documented bias band, implement a configurable closure in
  `source/src_el/el_forces.f90` (new `ELDragClosure = difelice |
  dns_table` key, default unchanged) and re-run one EL validation case
  (RZ-style or kroupa_couette φ point) with it — one datasheet row showing
  the DNS-informed closure's effect.
- Self-voidage: DNS gives the true undisturbed fluid velocity at a particle
  (no self-induced field problem). Use a matched dilute settling case to
  measure the EL kernel's self-voidage bias directly against DNS truth,
  replacing the analytic ε_eff ≈ 0.975 estimate with a measured curve vs
  δ/d_p (kernel width factor) — potentially a corrected sampling rule.
- Lubrication threshold: D2.1's crossover rule feeds back into the EL/PE
  `lubricationCutoff_`/ε_c choice (the deferred ε_c sensitivity item from
  the EL campaign, now anchored by resolved data).

### D5.3 Selection guide
- Final chapter of the practitioner's guide: decision table (N particles,
  d_p/L, Re, φ, wall proximity, shape) → DNS / EL / EL-with-DNS-closures,
  with cost estimates per row from the campaign's measured job costs.

**Exit gate D5**: report + datasheet complete; memo and selection guide
published; any EL closure changes validated and gated.

---

## 9. Stage D6 — Capstone ("Endgegner") cases

Both cases were solved before the campaign, by brute trial-and-error.
They return at the end as the acid test of the campaign's central
promise: a practitioner armed with the D1–D3 guidelines sets them up in
an informed, structured way — no parameter fishing. Success is measured
as much by the *setup process* (every choice traceable to a guideline
table) as by the physics gates; every place the guide proves
insufficient is logged and fed back into it. That feedback is the real
deliverable.

### D6.1 Kroupa box — effective viscosity of a sheared suspension
- Prior art: pre-campaign in-house DNS agreed with Krieger–Dougherty
  predictions and the Kroupa & Šoóš simulations (the Kroupa length scale
  differed, but for effective viscosity the Re regime and solid fraction
  dominate). Reference: Kroupa & Šoóš 2016 [in repo], §13.
- Observable: total stress → effective viscosity μ_eff(φ); gate:
  Krieger–Dougherty band + consistency with the prior in-house result.
- Rule: fresh implementation only — the legacy
  `ENABLE_LUBRICATION`/`sliding_wall_force` path (D2.1 note) stays
  retired.

### D6.2 Archimedes tube crystallizer (ATC) — suspension regime vs RPM
- Prior in-house study: a single liquid slug with a **rotation boundary
  condition** standing in for the full ATC helix (full-helix DNS is out
  of reach); the rotating BC "drives" the fixed-mesh slug through the
  ATC; ~5% particle volume fraction.
- Observable: particle-suspension regime as a function of rotation
  speed — at low RPM particles settle out; above a critical RPM they
  suspend evenly (the engineering target).
- Gate: reproduce the regime boundary (settling ↔ even suspension) of
  the prior study; mechanism-level observables (vertical concentration
  profile vs RPM) rather than point values.
- Prerequisites beyond D1–D3: rotation-BC recertification (audit the
  existing rotating-boundary support before staging); D3.2 diagnostics
  (concentration profiles) reused.

**Exit gate D6 (= campaign exit)**: both cases green (or
FAIL-attributed) with setup performed strictly from the practitioner's
guide; guide-defect log folded back into the guide.

---

## 10. Case ledger (initial — IDs are datasheet keys)

| ID | Case | App | Observable | Reference | Gate type | Stage |
|---|---|---|---|---|---|---|
| D0-SMOKE | 10-step sedimentation smoke, bitwise twin | bench_sedimentation | diagnostic lines | branch head | byte_identity | D0 |
| D0-E4-XMODE | E4 serial vs parallel PE | bench_sedimentation | u_t | internal | ≤0.1% | D0 |
| D0-E4-PIV | E4 vs digitized PIV | bench_sedimentation | v(t) | ten Cate 2002 PIV | provisional (tol from D1) | D0 |
| D1-ARR-CONV | SC array drag, D/h ladder | dns_drag | F/F_Stokes | Hasimoto; Zick–Homsy | convergence order | D1 |
| D1-NOISE | translating-sphere force noise | fc_ext variant | ΔF amplitude vs D/h | self (characterization) | RECORDED | D1 |
| D1-TC-MATRIX | E1–E4 × L2/L3/L4 × dt ladder | bench_sedimentation | v(t), u_t | ten Cate PIV | per-cell tol model | D1 |
| D1-UD14 | moderate-Re free fall (optional) | bench_sedimentation | regime + u_t | Uhlmann–Dušek 2014 | decision point | D1 |
| D2-WALL | sphere–wall approach ladder | fc_ext variant | F(gap) | Brenner/Cooley–O'Neill | crossover map | D2 |
| D2-PAIR | sphere–sphere approach | fc_ext variant | F(gap), normal+shear | lubrication asymptotics | crossover map | D2 |
| D2-DKT | drafting–kissing–tumbling | new dkt case | phase seq., t_kiss | Fortes 1987; Apte 2009 | mechanism + band | D2 |
| D3-RAND | random-array drag F(φ,Re), 3 seeds | dns_drag | ⟨F⟩ ± spread | Beetstra 2007; Tenneti 2011 | within correlation band | D3 |
| D3-HS | wall-bounded hindered settling | new case | ⟨v⟩/u_t vs φ | Richardson–Zaki band | mechanism + RECORDED | D3 |
| D3-FLU | fixed-bed Δp; u_mf (go/no-go) | bench_fluidization | Δp, u_mf | Ergun; bed weight | band | D3 |
| D4-VOL | shape indicator volume smokes | geom/creep | V_indicator/V_exact | analytic | precision | D4 |
| D4-JEFF | Jeffery orbit, r = 2, 4 | new couette case | orbit period | Jeffery 1922 | few-% + drift RECORDED | D4 |
| D4-OBER | spheroid Stokes settling | creep | u_t (2 orientations) | Oberbeck/Happel–Brenner | few-% | D4 |
| D4-ORIENT | inertial orientation attractor | creep/bench | attractor | literature DNS (TBD) | mechanism | D4 |
| D5-TWIN | DNS vs EL vs PIV, E1–E4 | both | v(t) overlays | campaign data | memo | D5 |
| D5-CLOS | Di Felice vs DNS F(φ,Re) | el_pipeflow | closure deviation | D3-RAND | band + EL re-run | D5 |
| D5-SV | EL self-voidage vs DNS truth | both | bias vs δ/d_p | D-campaign data | measured curve | D5 |
| D6-KROUPA | sheared-suspension μ_eff(φ) (Kroupa box) | new box case | μ_eff(φ) | Krieger–Dougherty; Kroupa–Šoóš 2016; prior in-house | band + setup-process audit | D6 |
| D6-ATC | ATC slug: suspension regime vs RPM | new rotating-BC case | regime boundary, concentration profile | prior in-house study | mechanism + regime boundary | D6 |

Full citations for every reference in this table, with acquisition status
(dashboard convention: cited literature PDFs live at repo root), are in
§13.

---

## 11. Risks and standing decision points

1. **Cost blow-up** is the main risk: resolved DNS scales as (D/h)³ per
   particle per level. Mitigation: D1's cost table is a *gate* for D3.2/D3.3
   scoping; every stage sizes N and resolution from measured cost, not
   hope. Descopes are recorded, not hidden.
2. **`q2p1_dns_drag` is serial-PE-only** (parallel path aborts by design).
   Fine for fixed arrays (no PE dynamics needed), but rank-count-limited
   via the CFD side; if L4 arrays are needed, either extend the setup to
   parallel PE or accept the resolution ceiling — decision at D3 design.
3. **Grid-crossing noise may dominate tolerances** at affordable D/h; if
   D1.2 shows a large noise floor, gates shift from instantaneous values to
   time-averaged/filtered observables (state the filter in the RUNBOOK).
4. **DKT chaos**: post-contact trajectories are not gateable; the design
   pre-commits to pre-contact observables so a FAIL cannot be argued away.
5. **Two ten Cate meshes** (quarter benchSym vs full box) are different
   discretizations; every cross-comparison states which is used. The
   quarter-box symmetry assumption itself gets one full-box vs quarter-box
   twin row in D1.
6. **Parallel-PE Cartesian partitioning** depends on the unmaintained
   METIS-4 `tools/partpy` path — D0.1 must produce a supported recipe
   before any parallel-PE campaign runs.
7. **FullC0ntact backend is out of scope** except D4's OBJ-mesh twin;
   recorded so the campaign's claims are clearly scoped to the PE backend.
8. **Existing DNS regression baselines were produced under
   `HardContactEulerLagrange` builds** (including the committed
   featflower_test sedimentation baseline). The D0.6 solver split
   therefore needs its own equivalence twin, and baselines are re-pinned
   under the DNS solver once that twin is green — never silently reused
   across solvers.

## 12. First tactical moves when execution starts

1. `git checkout -b feature/dns-validation-phase1` (branch off master per
   EL-campaign convention; PR targets `master`).
2. D0.2 diagnostics + D0.3 SED_BENCH retirement (code, small, testable via
   the 10-step smoke + bitwise twin).
3. D0.4 PIV digitization + mesh fixture commits (`git add -f`).
4. D0.6 solver decision note + `HardContactDNS` skeleton + cross-solver
   equivalence twin (pe work happens on a new pe branch off
   `feature-el-extensions`; never pushed by the agent).
5. D0.1 tooling generalization; stage E4 from scratch through it (exit
   gate D0).
6. Open the D1 RUNBOOKs and run the first probe (one L3 E4 + one L4 timing
   probe) before any sweep.

---

## 13. Literature references

Full citations, grouped by campaign role. Status tags: **[in repo]** = PDF
available locally — DNS-campaign papers live in `literature/` (untracked;
index in `literature/README.md`), EL-era papers remain at repo root until
the dashboard's path check is updated; **[wanted]** = please provide the
PDF (the EL campaign showed having the originals on hand pays off during
case design); **[optional]** = useful context, not gate-defining;
*(verify)* marks bibliographic details recalled from memory that should be
checked against the actual paper. First batch of 12 PDFs received
2026-07-31; the three *(verify)*-tagged entries among them were checked
against the title pages and confirmed.

### Method lineage (FeatFloWer FBM + PE)
- S. Turek, D. Wan, L.S. Rivkind, "The fictitious boundary method for the
  implicit treatment of Dirichlet boundary conditions with applications to
  incompressible flow simulations", in *Challenges in Scientific
  Computing*, LNCSE 35, Springer (2003). **[optional]** *(verify)*
- D. Wan, S. Turek, "Direct numerical simulation of particulate flow via
  multigrid FEM techniques and the fictitious boundary method", *Int. J.
  Numer. Meth. Fluids* 51 (2006) 531–566. **[in repo]**
  (`literature/wan_turek_2006.pdf`, citation confirmed) — the method paper
  behind `ForcesLocalParticles`' volume-integral force; also the origin of
  the Glowinski-style collision model examined there.
- D. Wan, S. Turek, "Fictitious boundary and moving mesh methods for the
  numerical simulation of rigid particulate flows", *J. Comput. Phys.* 222
  (2007) 28–56. **[optional]** *(verify)*
- R. Münster, O. Mierka, S. Turek, "Finite element-fictitious boundary
  methods (FEM-FBM) for 3D particulate flow", *Int. J. Numer. Meth.
  Fluids* 69 (2012) 294–313. **[in repo]**
  (`literature/munster_mierka_turek_2012.pdf`, citation confirmed) — the
  3D FBM lineage of this codebase; covers static AND adaptively aligned
  meshes (grid deformation), relevant to the dormant r-adaptivity
  machinery noted in §0.
- K. Iglberger, U. Rüde, "Massively parallel rigid body dynamics
  simulations", *Comput. Sci. Res. Dev.* 23 (2009) 159–167. **[optional]**
  — PE library origin.
- T. Preclik, U. Rüde, "Ultrascale simulations of non-smooth granular
  dynamics", *Comput. Part. Mech.* 2 (2015) 173–196. **[optional]** — the
  non-smooth hard-contact/PGS formulation family behind
  `HardContactEulerLagrange` and the planned `HardContactDNS`.
- M. Kroupa, M. Vonka, M. Soos, J. Kosek, "Utilizing the discrete element
  method for the modeling of viscosity of concentrated suspensions",
  *Langmuir* 32 (2016) 8451–8460. **[in repo]** (`kroupa_soos_2016.pdf`) —
  lubrication closures shared by EL and (per D0.6/D2.2) the DNS solver.

### D0/D1 — single sphere, resolution & dt
- A. ten Cate, C.H. Nieuwstad, J.J. Derksen, H.E.A. Van den Akker,
  "Particle imaging velocimetry experiments and lattice-Boltzmann
  simulations on a single sphere settling under gravity", *Phys. Fluids*
  14 (2002) 4012–4025. **[in repo]** (`ten_cate_piv.pdf`) — E1–E4
  definitions and the PIV trajectories to digitize in D0.4.
- M. Schäfer, S. Turek, "Benchmark computations of laminar flow around a
  cylinder", in *Flow Simulation with High-Performance Computers II*,
  Notes Numer. Fluid Mech. 48, Vieweg (1996) 547–566. **[optional]** — DFG
  C_D/C_L reference already regression-pinned.
- H. Hasimoto, "On the periodic fundamental solutions of the Stokes
  equations and their application to viscous flow past a cubic array of
  spheres", *J. Fluid Mech.* 5 (1959) 317–328. **[in repo]**
  (`literature/hasimoto_1959.pdf`) — D1-ARR-CONV dilute-array drag
  expansion.
- A.A. Zick, G.M. Homsy, "Stokes flow through periodic arrays of spheres",
  *J. Fluid Mech.* 115 (1982) 13–26. **[in repo]**
  (`literature/zick_homsy_1982.pdf`) — D1-ARR-CONV tabulated drag at
  general φ.
- A.S. Sangani, A. Acrivos, "Slow flow through a periodic array of
  spheres", *Int. J. Multiphase Flow* 8 (1982) 343–360. **[optional]** —
  cross-check on Zick–Homsy.
- M. Uhlmann, J. Dušek, "The motion of a single heavy sphere in ambient
  fluid: a benchmark for interface-resolved particulate flow simulations
  with significant relative velocities", *Int. J. Multiphase Flow* 59
  (2014) 221–243. **[in repo]** (`literature/uhlmann_dusek_2014.pdf`) —
  includes documented resolution requirements (D/h up to ~24), directly
  comparable to our guideline study; D1.4 go/no-go now unblocked on the
  reference side.
- N. Mordant, J.-F. Pinton, "Velocity measurement of a settling sphere",
  *Eur. Phys. J. B* 18 (2000) 343–352. **[optional]** — experimental v(t)
  time series, alternative free-fall anchor.

### D2 — pair interactions
- H. Brenner, "The slow motion of a sphere through a viscous fluid towards
  a plane surface", *Chem. Eng. Sci.* 16 (1961) 242–251. **[in repo]**
  (`literature/brenner_1961.pdf`) — D2-WALL analytic drag divergence.
- M.D.A. Cooley, M.E. O'Neill, "On the slow motion generated in a viscous
  fluid by the approach of a sphere to a plane wall or stationary sphere",
  *Mathematika* 16 (1969) 37–49. **[in repo]**
  (`literature/cooley_oneill_1969.pdf`) — near-contact asymptotics
  complementing Brenner.
- D.J. Jeffrey, Y. Onishi, "Calculation of the resistance and mobility
  functions for two unequal rigid spheres in low-Reynolds-number flow",
  *J. Fluid Mech.* 139 (1984) 261–290. **[in repo]**
  (`literature/jeffrey_onishi_1984.pdf`) — D2-PAIR normal + tangential
  two-sphere resistance functions.
- S. Kim, S.J. Karrila, *Microhydrodynamics: Principles and Selected
  Applications*, Butterworth-Heinemann (1991). **[optional]** — closed-form
  compilation of the above.
- A.F. Fortes, D.D. Joseph, T.S. Lundgren, "Nonlinear mechanics of
  fluidization of beds of spherical particles", *J. Fluid Mech.* 177
  (1987) 467–483. **[in repo]**
  (`literature/fortes_joseph_lundgren_1987.pdf`) — the DKT experiment.
- R. Glowinski, T.-W. Pan, T.I. Hesla, D.D. Joseph, J. Périaux, "A
  fictitious domain approach to the direct numerical simulation of
  incompressible viscous flow past moving rigid bodies: application to
  particulate flow", *J. Comput. Phys.* 169 (2001) 363–426. **[in repo]**
  (`literature/glowinski_pan_hesla_joseph_periaux_2001.pdf`) — canonical
  DKT computation.
- S.V. Apte, M. Martin, N.A. Patankar, "A numerical method for fully
  resolved simulation (FRS) of rigid particle–flow interactions in complex
  flows", *J. Comput. Phys.* 228 (2009) 2712–2738. **[in repo]**
  (`literature/apte_martin_patankar_2009.pdf`, citation confirmed) — 3D
  sphere DKT with quantitative trajectories; candidate primary D2-DKT
  reference.
- W.-P. Breugem, "A second-order accurate immersed boundary method for
  fully resolved simulations of particle-laden flows", *J. Comput. Phys.*
  231 (2012) 4469–4498. **[in repo]** (`literature/breugem_2012.pdf`) —
  alternative modern DKT reference with published trajectories; pick Apte
  or Breugem as primary during D2-DKT design.

### D3 — collective behavior
- R. Beetstra, M.A. van der Hoef, J.A.M. Kuipers, "Drag force of
  intermediate Reynolds number flow past mono- and bidisperse arrays of
  spheres", *AIChE J.* 53 (2007) 489–501. **[in repo]** (`literature/beestra_2007.pdf`)
  — D3-RAND primary correlation.
- S. Tenneti, R. Garg, S. Subramaniam, "Drag law for monodisperse
  gas–solid systems using particle-resolved direct numerical simulation of
  flow past fixed assemblies of spheres", *Int. J. Multiphase Flow* 37
  (2011) 1072–1092. **[in repo]** (`literature/tenneti_2011.pdf`) —
  D3-RAND second correlation (the gate is the band between the two).
- M.A. van der Hoef, R. Beetstra, J.A.M. Kuipers, "Lattice-Boltzmann
  simulations of low-Reynolds-number flow past mono- and bidisperse arrays
  of spheres: results for the permeability and drag force", *J. Fluid
  Mech.* 528 (2005) 233–254. **[optional]** — low-Re limit of D3-RAND.
- J.F. Richardson, W.N. Zaki, "Sedimentation and fluidisation: Part I",
  *Trans. Inst. Chem. Eng.* 32 (1954) 35–53. **[in repo]**
  (`literature/richardson_1954.pdf`) — D3-HS exponent band (also closes
  the loop on the EL campaign's RZ reframe).
- N.-Q. Nguyen, A.J.C. Ladd, "Sedimentation of hard-sphere suspensions at
  low Reynolds number", *J. Fluid Mech.* 525 (2005) 73–104. **[optional]**
  — DNS hindered-settling reference at modest N.
- X. Yin, D.L. Koch, "Hindered settling velocity and microstructure in
  suspensions of solid spheres with moderate Reynolds numbers", *Phys.
  Fluids* 19 (2007) 093302. **[optional]** — finite-Re hindered settling
  DNS, closer to our achievable regime than pure-Stokes references.
- S. Ergun, "Fluid flow through packed columns", *Chem. Eng. Prog.* 48
  (1952) 89–94. **[in repo]** (`literature/ergun_1952.pdf`) — D3-FLU
  fixed-bed Δp and u_mf anchor.
- C.Y. Wen, Y.H. Yu, "Mechanics of fluidization", *Chem. Eng. Prog. Symp.
  Ser.* 62 (1966) 100–111. **[optional]**
- M.A. van der Hoef, M. van Sint Annaland, N.G. Deen, J.A.M. Kuipers,
  "Numerical simulation of dense gas–solid fluidized beds: a multiscale
  modeling strategy", *Annu. Rev. Fluid Mech.* 40 (2008) 47–70.
  **[optional]** — framing for the DNS→EL closure ladder (§8/D5).

### D4 — non-spherical particles
- G.B. Jeffery, "The motion of ellipsoidal particles immersed in a viscous
  fluid", *Proc. R. Soc. Lond. A* 102 (1922) 161–179. **[in repo]**
  (`literature/jeffrey_1922.pdf`, title page verified) — D4-JEFF orbit
  period, parameter-free.
- A. Oberbeck, "Über stationäre Flüssigkeitsbewegungen mit
  Berücksichtigung der inneren Reibung", *J. Reine Angew. Math.* 81 (1876)
  62–80. **[optional]** — original spheroid Stokes drag; the practical
  formulas come from Happel & Brenner.
- J. Happel, H. Brenner, *Low Reynolds Number Hydrodynamics*, Noordhoff
  (1973 ed.). **[in repo]** (`literature/happel_brenner_1973.pdf`, full
  book) — D4-OBER drag formulas for prolate/oblate spheroids, both
  orientations.
- M.N. Ardekani, P. Costa, W.-P. Breugem, L. Brandt, "Numerical study of
  the sedimentation of spheroidal particles", *Int. J. Multiphase Flow* 87
  (2016) 16–34. **[in repo]** (`literature/ardekani_2016.pdf`, arXiv
  1602.05769 preprint, title page verified) — candidate primary for
  D4-ORIENT (inertial orientation attractor, resolved DNS); also treats
  spheroid pair DKT (bridges D2.3/D4).
- W.W. Willmarth, N.E. Hawk, R.L. Harvey, "Steady and unsteady motions and
  wakes of freely falling disks", *Phys. Fluids* 7 (1964) 197–208.
  **[optional]** — if a disk/oblate case is added.
- K.O.L.F. Jayaweera, B.J. Mason, "The behaviour of freely falling
  cylinders and cones in a viscous fluid", *J. Fluid Mech.* 22 (1965)
  709–720. **[optional]** *(verify)*

### D5 — EL side (already in the EL campaign's canon)
- R. Di Felice, "The voidage function for fluid-particle interaction
  systems", *Int. J. Multiphase Flow* 20 (1994) 153–159. **[wanted]** —
  the closure D5-CLOS tests against DNS truth.
- I.M. Krieger, T.J. Dougherty, "A mechanism for non-Newtonian flow in
  suspensions of rigid spheres", *Trans. Soc. Rheol.* 3 (1959) 137–152.
  **[optional]** — already the EL campaign's KD anchor.

**Priority acquisition list — remaining** (second batch received
2026-07-31: Beetstra 2007, Tenneti 2011, Richardson–Zaki 1954,
Ergun 1952, Jeffery 1922 [verified], Happel–Brenner full book,
Ardekani 2016 [verified; arXiv preprint — also covers spheroid DKT,
bridging D2.3/D4]; filenames in `literature/README.md`): only
**Di Felice 1994** is still outstanding.
