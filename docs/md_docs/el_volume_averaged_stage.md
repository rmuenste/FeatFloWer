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
  (v2wF_cont, and v2wF_all for both together, predicted ≈ 3.4 and ≈ 4.0
  before the sampled-input offset of ≈ 0.6).
- The lone particle settles 1.1 % slower, as in the tier-2 gate: with
  div u = g it samples its own displacement flow. This is the self-induced
  velocity that the undisturbed-field stage (U1) is to remove; it is
  bounded by the 0.7 % the pressure force alone produces, and both are
  small against the collective effect.

Row `v2w_cont_n40_probe` RECORDED. The pressure-force rung itself is
recorded on `feature/el-rz-postfix` (RUNBOOK section 11, commit 8df5a228,
n = 2.36 over six clouds, +0.64 vs predicted +0.65).
