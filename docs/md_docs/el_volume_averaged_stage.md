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
