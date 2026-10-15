# EL volume-averaged continuity — stage record (branch feature/el-volume-averaged)

Worktree `FeatFloWer-va`, base feature/dns-validation fa6c242d + merged
feature/el-rz-postfix; pe feature/el-box-walls-on-master c02cf28 (= master
1f60291 + boxWallsEnabled_). Build `build-va/` (certified EL recipe, gcc 14.3).

## Why

The CFDEM survey (liggghts-public-34, 2026-10-14) and our column logs agree:
the EL fluid solves single-phase div u = 0, so the volume-averaged return flow
through a settling cloud does not exist (sampled u_f at the particles ≈ 0
where −φu_p/ε would be 0.008 / 0.023 at N = 40 / 120), and the particle gets
ρ_f buoyancy only. Each costs one power of ε in the Di Felice hindrance
exponent (ceiling χ−1 ≈ 2.1 at Re ≈ 80 vs measured 1.72 vs DNS 4.57).

## What was implemented (commit 34ead1f6 + guard db857540)

`SimPar@ELVolumeAveraged = continuity` imposes div u = g,
g = −(∂ε/∂t + u·∇ε)/ε, through the pressure Poisson right-hand side (G/Δt on
the constant P1 entry of each element, after `Matdef_General_LinScalar(...,1)`
in `Transport_q2p1_UxyzP_fluid_core`; derivation in the reconnaissance notes:
C δp = −(1/k)Bᵀũ + G/k). ∇ε = −∇α_p from a new cell-averaged kernel-gradient
field deposited with α_p (`EL_KERNEL_DERIVATIVE`); element velocity = vertex
mean; the volume mean of g is removed (constant nullspace, closed/periodic
boxes). Momentum operator unchanged. New diagnostic `EL_VA_SOURCE`
(∫g, ḡ, L1 of the two terms). `history_valid` seeds ε_old at the first
finalisation (the initial ε_old = 1 gave a spurious −α/Δt kick at step 1).

## Gates (tier-2 momentum_conservation, QBOX9, 28 ranks, 60 steps)

| run | job | result |
|---|---|---|
| key unset, before the guard | 149705 | BITWISE vs branch reference (job 149704) on all six EL streams |
| key unset, after the guard | 149712 | BITWISE |
| key on, before the guard | 149706 | step-1 ∫g = 6.8e-2 (the kick), particle dp −4.62e-9 vs −4.78e-9 off |
| key on, after the guard | 149709 | step-1 ∫g = 3.3e-8; drift 1.8e-12; Newton pair 6e-30; dp −4.7185e-9 vs −4.7835e-9 |

The remaining 1.4% in dp (1.1% in ⟨u_p⟩) with the source on is the particle
sampling its own displacement flow at the kernel scale, a second
self-induced effect to be removed by the undisturbed-field stage; the
collective displacement flow of a cloud is the physics the stage is for.

## In flight

Column (v2w decks + key): u_t job 149710, N=40 seed 1 job 149711 (level 2,
109 ranks). Expected: sampled u_f at the particles moves toward −φu_p/ε and
U/u_t at N=40 drops below the 0.928 of the single-phase run.

## Fritz offload (2026-10-14)

EL worktree on NHR@FAU Fritz, separate from the DNS clone:
`/anvme/workspace/k115ce12-featflower/FF/FeatFloWer-el` on
feature/el-volume-averaged (1fd573e2, pe c02cf28; fetched with
`--recurse-submodules=no` because the merged history references an unpushed
pe commit). Build `build-el-fritz/` (Unix Makefiles — no ninja on Fritz;
gcc 14.2 + OpenMPI 4.1.8 from `~/ff-env-gcc14.sh`, parallel PE, serial mode
OFF, ENABLE_FBM_ACCELERATION=ON, `Q2P1_MESH_DIR` = workspace mesh_repo);
ldd prefixes as the guide §8 requires. Mesh fixture copied by rsync to
`FeatFloWer-el/dkt_mesh_v1/coarseB` (+ DKT_B_108), case links recreated.
Helpers: `tools/fritz_el_submit_np.sh` (`--export=NONE`, env sourced in the
job), `tools/stage_fritz_rungs.sh` (variants cont / std / stdcont / all ×
u_t + 12 clouds = 52 rundirs under build-el-fritz/applications, 639 MB).
Twin job 4331912 (tier-2 momentum, 28 ranks, singlenode) pending.
Workspace expires 2026-10-28 (10 extensions left).

### Fritz twin PASS and ladder submission (2026-10-14)

Tier-2 momentum twin on Fritz (job 4331912, 28 ranks) vs the branch
reference log: EL_MEAN_SLIP, EL_CONTACT_STATS, EL_FEEDBACK_CONSERVATION,
EL_NEWTON_PAIR, EL_VOLUME_CONSERVATION bitwise on every non-roundoff field;
EL_MOMENTUM_ELEMINT worst relative difference 3.4e-7 (a printed 1e-8
quantity, last digit); roundoff-level fields differ by ≤ 1e-16 absolute.
Same signature as the DNS `fritz_build_twin` row (gcc 14.3 local vs 14.2
Fritz, different OpenMPI). Certified for independent runs; no trajectory is
ever restarted across machines.

Submitted on Fritz (owner go-ahead): 52 jobs 4332057–4332108, multinode,
109 ranks / 2 nodes / 04:00:00 each — variants cont, std, stdcont, all ×
(u_t + the 12 DNS clouds). Ids in
`build-el-fritz/applications/v2wF_ids.txt` on Fritz. Each variant is
analysed against its own u_t with tools/d32_ladder_analysis.py after
tools/el_particle_log_merge.py.

### Local continuity column pair RESULT (2026-10-14, jobs 149710 / 149711, med, 109 ranks, 1 h 26 / 1 h 34)

`ELVolumeAveraged = continuity` on the VA binary (build-va, ee1ae990), level
2, δ = 2.5 d, campaign Di Felice form, no pressure force; u_t + the N = 40
seed-1 cloud (byte-identical to the campaign and pressure-force twins).
Clean: no NaN, volume error ≤ 1.3e-15, Newton pair ≤ 1.6e-15; EL_VA_SOURCE
ḡ removed to 2.8e-7 every step, ∫g before removal 2.4e-4 at N = 40.

| run | u_t or U | U/u_t | ε seen | u_f at particles (z) |
|---|---|---|---|---|
| vac_ut | 0.35349 (campaign 0.35734, −1.1 %) | – | 0.9994 | 2.1e-4 |
| vac_n40_s1 | 0.3146 ± 0.0016 | **0.890** | 0.976 | **7.9e-3** |
| campaign twin (v2w_n40_s1) | 0.3342 | 0.928 | 0.976 | 5.6e-4 |
| pressure-force twin (v2wpf_n40_s1) | 0.3202 | 0.903 | 0.976 | 5.4e-4 |

Reading:
- The return flow now exists at the particles: the kernel-sampled fluid
  velocity at N = 40 rises from 5.6e-4 to 7.9e-3, which is exactly
  −φ_seen u_p/ε_seen = 0.024 × 0.317 / 0.976 = 7.8e-3. The volume-averaged
  continuity delivers the return flow the sampled voidage implies.
- Single-cloud exponent estimate ln(U/u_t)/ln(1−φ_cloud), φ_cloud 0.047:
  campaign 1.55, pressure force 2.11, continuity **2.41**; increments +0.56
  and +0.86 against the bookkeeping's +0.65 (pressure) and +1.0 (return
  flow). Same ordering and size; the six-cloud exponents come from Fritz
  (v2wF_cont predicted ≈ 3.4 before the sampled-input offset of ≈ 0.6;
  v2wF_all runs the standard drag form and is predicted ≈ 2.1 — the
  earlier ≈ 4.0 written here was the campaign-form value, corrected
  2026-10-15).
- The lone particle settles 1.1 % slower, as in the tier-2 gate: with
  div u = g it samples its own displacement flow. This is the self-induced
  velocity that the undisturbed-field stage (U1) is to remove; it is
  bounded by the 0.7 % the pressure force alone produces, and both are
  small against the collective effect.

Row `v2w_cont_n40_probe` RECORDED. The pressure-force rung itself is
recorded on `feature/el-rz-postfix` (RUNBOOK section 11, commit 8df5a228,
n = 2.36 over six clouds, +0.64 vs predicted +0.65).

### Fritz u_t preview (2026-10-14 17:30, four u_t runs + two N = 20 clouds finished, 46 to go)

| variant | u_t | vs campaign 0.35734 | note |
|---|---|---|---|
| cont | 0.35349 | −1.1 % | identical to the local vac_ut to all printed digits (cross-machine) |
| std | 0.36964 | +3.4 % | = Dallavalle ε = 1 fixed point 0.3686 within 0.3 %: the standard form's ε_eff^(−1.1) nearly cancels the lone-particle self-voidage bias that the campaign form's ε_eff^(−3.1) carries (ε_eff ≈ 0.975 → ratio ε^−2 = 1.05 in drag → +3.4 % in u_t at Re ≈ 80) |
| stdcont | 0.36565 | +2.3 % | = std × cont (0.3696 × 0.989) |
| all | 0.36295 | +1.6 % | = stdcont × pressure force (0.3656 × 0.993) |

The three single-particle effects are independent and multiplicative at the
0.1 % level. Final-time (not windowed) N = 20 seed-1 previews: cont
U/u_t 0.894 (campaign seed mean 0.965, DNS 0.892), std 0.952. Six-cloud
exponents per variant follow when the set completes.

### Local standard-form column pair RESULT (2026-10-14, jobs 149717 / 149718, med, 109 ranks, 1 h 26 / 1 h 35)

`ELDragModel = difelice_std` on the VA binary, continuity off, no pressure
force; u_t + the N = 40 seed-1 cloud. Clean: no NaN, volume error ≤ 8.5e-16,
Newton pair ≤ 2.5e-15.

| run | u_t or U | U/u_t | u_f at particles (z) |
|---|---|---|---|
| vas_ut | 0.36964 (campaign 0.35734, +3.4 %; Fritz std_ut identical to all digits) | – | 6.7e-6 |
| vas_n40_s1 | 0.3612 ± 0.0016 | **0.977** | 5.3e-4 |
| campaign twin | 0.3342 | 0.928 | 5.6e-4 |

Reading: the standard form alone weakens the hindrance, single-cloud
exponent 0.47 vs 1.55 campaign (−1.08; bookkeeping −1.3). This is the
expected sign: the ε² that the campaign form carries in the drag stands in
for the two volume-averaged effects (pressure-gradient share and return
flow) that the fluid equations omit. Removing it without adding them
leaves almost no collective hindrance. The standard form is therefore only
meaningful together with `continuity` and `ELPressureForce` (Fritz
variants stdcont and all). Row `v2w_std_n40_probe` RECORDED.

Local N = 40 seed-1 summary (same cloud, own u_t, windowed 15–25):
campaign 0.928 → pressure force 0.903 → continuity 0.890 → standard form
alone 0.977; DNS seed mean 0.826.

### Fritz 4 × 13 factorial RESULT (completed 2026-10-15 00:55, jobs 4332057–4332108, all COMPLETED, 3 h 05 – 3 h 46 each, no timeouts)

Four variants on the VA binary (build-el-fritz, twin-certified), level 2,
δ = 2.5 d, each with its own u_t and the twelve byte-identical DNS clouds;
analysed locally (tools/el_particle_log_merge.py + d32_ladder_analysis.py,
window 15–25; Fritz python has no numpy). Keys: cont = `ELVolumeAveraged =
continuity`; std = `ELDragModel = difelice_std`; stdcont = both; all = both
+ `ELFluidGravity = No` + `ELPressureForce = Yes`. Clean: no NaN, volume
error ≤ 1.7e-15, Newton pair ≤ 8e-15 everywhere except the bottom-wall
arrivals noted below.

| variant | u_t | U/u_t seed means N = 20 / 40 / 80 / 120 | n (rms) | Δn vs campaign: measured / bookkeeping |
|---|---|---|---|---|
| campaign (local, 2026-10-03) | 0.3573 | 0.965 / 0.935 / 0.870 / 0.805 | 1.72 (0.005) | – |
| pressure force (local, 6 clouds, RUNBOOK §11) | 0.3547 | – / 0.911 / – / 0.746 | 2.36 (0.007) | +0.64 / +0.65 |
| cont | 0.3535 | 0.947 / 0.899 / 0.808 / 0.722 | **2.64** (0.006) | +0.92 / +1.00 |
| std | 0.3696 | 0.983 / 0.978 / 0.953 / 0.925 | 0.61 (0.005) | −1.11 / −1.30 |
| stdcont | 0.3657 | 0.968 / 0.941 / 0.885 / 0.827 | 1.52 (0.004) | −0.20 / −0.30 |
| all | 0.3630 | 0.956 / 0.916 / 0.838 / 0.764 | 2.16 (0.005) | +0.44 / +0.35 |
| DNS (same clouds) | 0.4061 | 0.892 / 0.826 / 0.689 / 0.588 | 4.57 | |

Seed spans ≤ 0.015 at every rung; every variant is a clean power law
(rms ln-residual ≤ 0.007).

Reading:
- All five measured increments agree with the finite-Re bookkeeping
  (docs/md_docs/el_drag_convention_note.md, Re ≈ 80, χ 3.10, a 0.46) to
  within 0.2 in n, in sign and size: pressure force +0.64 (+0.65), return
  flow +0.92 (+1.00), standard form −1.11 (−1.30), standard + return flow
  −0.20 (−0.30), standard + both +0.44 (+0.35). The bookkeeping is
  validated as the accounting of what the fluid equations and the drag
  convention each contribute to the hindrance exponent.
- The bookkeeping's absolute level is 2.31 for the campaign set-up against
  1.72 measured: a constant offset of ≈ 0.6 that no variant moves. It is
  the sampled-input bias established by the kernel-width and level-3 rungs
  (diluted voidage ε seen 0.927 vs 0.888 at N = 120, self-disturbed slip),
  and it is the target of the undisturbed-field stage.
- The physically complete campaign-form set — continuity + pressure force
  with the campaign ε^(−χ) drag — was not part of the factorial. The
  bookkeeping predicts 1.72 + 1.65 = **3.4** for it (3.96 absolute); with
  the sampled-input offset removed, ≈ 4.0 against DNS 4.57. It is the
  single most informative remaining rung (13 jobs); staged as `contpf`,
  not submitted (owner decision).
- Standard-form variants: `all` (the CFDEM-style combination) reaches 2.16,
  short of the campaign form with the same physics because the ε²
  difference between the conventions is worth 1.3 in n at this Re. Which
  convention is right is a question for the drag-correlation literature
  (Di Felice 1994 defines F₀ at the superficial flux), not for this ladder;
  the ladder says the equations and the convention must be chosen
  together.

Bottom-wall arrivals (audit note, not a FAIL): the fastest particle of
some clouds reaches the bottom wall at t ≈ 23.2–24.0 (campaign runs
included; same clouds, same window). When it bounces the audit stays
clean; when it comes to rest (z = 0.50, v = 0) `EL_NEWTON_PAIR` shows a
fixed z mismatch of 4.1e-5 per audited step, and the impact step itself
can show up to 3.5e-2 (std_n80_s2, a 0.31 bounce impulse). Six runs:
all_n20_s1, stdcont_n20_s1, stdcont_n40_s3, std_n40_s2, std_n40_s3,
std_n80_s2. The resting-contact wall impulse is not reported by
`getElWallContactImpulse`, so the audit undercounts it; the particle
dynamics are unaffected (it rests). Effect on the windowed U ≤ 0.2 %,
below the seed span. The bounce height (0 to +0.19) varies with the
penetration at detection, i.e. the bounce is pe's position correction,
not restitution. To be fixed in the wall-impulse accounting before any
case that relies on particles resting on a wall.

Rows `v2w_cont_exponent`, `v2w_std_exponent`, `v2w_stdcont_exponent`,
`v2w_all_exponent` RECORDED.

### Peer corrections and the pressure-force feedback defect (liggghts-public-1e, 2026-10-15; verified in this tree)

1. **The particle pressure force was never fed back to the fluid.**
   `source/src_el/el_forces.f90`: `particle_total = drag + pressure + lift
   + grav_buoy` but `feedback_force = drag + lift`. With `ELPressureForce
   = Yes` the particle receives −V∇p and the fluid keeps its full −∇p over
   the particle volume: an O(φ) momentum imbalance, ≈ 5 % in the column.
   The Newton-pair audit could not see it because its expected impulse is
   built from the particle-side sum (`el_transfer.f90`,
   `el_pair_expected_impulse = dt*(drag_impl + lift + pressure +
   grav_buoy)`), not from what was deposited. Consequence: the
   pressure-force rung (RUNBOOK §11, n = 2.36) and the `all` variant
   (2.16) ran an inconsistent momentum set; their agreement with the
   bookkeeping increments is not evidence that the set was right.
   Fix (this commit): new key `SimPar@ELPressureFeedback = Yes` adds the
   pressure force to the deposited feedback, which makes the
   cell-averaged fluid pressure force −ε∇p, the Model-A form (the
   alternative, ε-weighting ∇p in the operator, is the same thing and
   must not be combined with it). Default No keeps every certified log
   bitwise. Gates: key-unset twin bitwise vs the branch reference; key-on
   tier-2 with the audit; then the pressure-force rung and `contpf` are to
   be re-run with the key on, and the §11 / `all` rows re-labelled
   "inconsistent set" in the datasheet.
2. **CFDEM's Di Felice uses χ(Re) by default** (DiFeliceDrag.C:103
   `setSwitches(17,true)`; χ = 2 only if `voidageFunctionDiFelice` is set
   false). The earlier survey note that CFDEM defaults to χ = 2 is
   withdrawn; the ε² convention point is unchanged.
3. **VTK output drops the EL fields at `OutputLevel = MAX+1`**
   (`source/OutputProfiles.f90`: a field is written only if `SIZE(field)
   == NoOfElem`, and at MAX+1 NoOfElem is the refined output mesh's
   count). Use `OutputLevel = MAX` when ELAlphaP / ELEpsilonF / ELDepsFDt
   are wanted.
4. Peer's absolute bookkeeping (Re ≈ 80): camp 2.35, pf 3.01, cont 3.35,
   std 1.02, stdcont 2.02, all 2.68 — the same increments as ours
   (+0.66, +1.00, −1.33, −0.33, +0.33). Model-B references for the
   unbounded exponent: DPMFoam Wen-Yu with the volume source 3.09 at
   Ar 6338, Garside & Al-Dibouni 3.07; Richardson–Zaki with the wall term
   at d/D = 1/6 gives 4.77, consistent with DNS 4.57 being a confinement
   law. A DPMFoam twin of v2w on our clouds is running on their side.
5. CFDEM gate for the continuity stage: with the flow solved, Σ ε U_z V
   matches N V_p |w| to 0.7 % (their voidfractionTest, divided void
   fraction; Gauss aborts, trilinear puts 1e10 into boundary cells).
   Specs: LIGGGHTS-PUBLIC md_docs/case_specs/ (03_column,
   04_voidfraction).

**Gate for `ELPressureFeedback` (job 149907, 2026-10-15, tier-2 box, 28 ranks, three
chained runs of 60 steps):** key unset → `run_g0_twin_va5.log` bitwise against
`reference/run_g0_twin_va_2026-10-14.log` on all six streams. Pressure force
on without / with feedback: `EL_FEEDBACK_CONSERVATION` rhs z 9.0587e-6 →
9.1420e-6 (+0.92 %, the deposited pressure force), residual 1e-20 both;
Newton pair unchanged (particle side); no NaN. `EL_MOMENTUM_ELEMINT` cannot
see the change: in a closed box ∫u dV = ∮ x (u·n) dS = 0 by incompressibility,
so the "fluid" entry is a discretisation residual whatever the forcing. The
column re-run is the physical test. Staged locally on the VA binary, NOT
submitted: `build-va/applications/q2p1_el_rundir_vapf_{ut,n40_s1..3,n120_s1..3}`
(pressure force + feedback, campaign drag, continuity off), the like-for-like
re-run of RUNBOOK §11; `contpf` on Fritz needs the key added and a rebuilt
binary (after push).

**Peer reply on the inconsistent set and the DPMFoam twin (liggghts-public-1e,
2026-10-15, after the factorial was sent; arithmetic checked here):**

- Their pf column assumes the consistent pressure gradient (particle feels
  −V∇p with |∇p| = φΔρg, which Model A with fluid −ε∇p and our unweighted
  −∇p + `ELPressureFeedback` both give). For a steady homogeneous
  suspension the inconsistent set (feedback = drag only) gives
  |∇p| = nF_d, so F_d = VΔρg/(1+φ) instead of VΔρg(1−φ): a φ² difference,
  effective pressure share k_eff = ln(1+φ)/(−ln(1−φ)) = 0.98/0.96/0.92/0.89
  at φ = 0.02/0.04/0.08/0.12 (checked). At Re ≈ 80 a unit of k is worth
  ≈ 0.65 in n, so the inconsistent set lowered the fitted n by ≈ 0.06
  (peer's corrected figure; their first message said 0.02–0.04). Prediction
  for the staged re-runs: pf and `all` rise by ≈ 0.06, i.e. 2.36 → ≈ 2.42
  and 2.16 → ≈ 2.22. The gap to DNS is untouched by the fix; the fix
  matters in a bed (φ ≈ 0.6: Δp comes out at W/(1+φ) ≈ 0.63 W).
- DPMFoam R1, homogeneous column, Wen–Yu, 12×12×4 cells: u_t 0.3967,
  U/u_t 0.925/0.860/0.735/0.625 at φ 0.02/0.04/0.08/0.12, n = 3.69
  (rms 0.002) against its own closure's 3.09. DPMFoam also sits 0.6 off its
  bookkeeping, steeper rather than shallower. Clean window t = 8–20 only;
  past t ≈ 25 at φ ≥ 0.08 the column layers (their reading: one-step
  face-flux spikes from the return flow).
- DPMFoam R2, our D3.2 clouds, 6d column, no-slip walls, 2×2×8 cells:
  u_t 0.399 off-axis (the axis probe sits on a 4-cell corner and oscillates);
  seed-mean U/u_t 0.917/0.850/0.717/0.552 at N = 20/40/80/120 (DNS
  0.892/0.826/0.689/0.588; ours 0.965/0.935/0.870/0.805). The cloud
  stretches upward (N = 120 window φ_cloud 0.076 vs 0.117 initial), so the
  exponent depends on the φ definition: 6.29 with the window φ_cloud of
  `d32_ladder_analysis.py`, 4.48 against the initial φ (DNS 4.57). For the
  same-cloud comparison use U/u_t at equal N with the window φ_cloud next
  to it. Open item for us: state which φ our n values use (the ladder tool
  uses φ_cloud = N V_p/(A H_cloud) in the window) and report the initial-φ
  fit alongside.
- Ergun (for E2/G6): CFDEM fixed-bed slopes come out at 0.67–0.76 of its
  own closures, plateaus within 3 %, U_mf 40–45 % high; mechanism open.
  Treat the CFDEM Ergun numbers as "reference as run".
- Details: LIGGGHTS-PUBLIC md_docs/case_specs/{03_column,05_ergun}/SPEC.md
  §12; their RUNBOOK runs 61–88.
