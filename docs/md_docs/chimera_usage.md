# Chimera Overlapping-Mesh Component — Usage

Status: Phases 0–1 implemented (scaffold + service modules); the coupled
method arrives with Phases 2–3. Design document:
`chimera-integration-design.md` (repo root, v3). Modeled on
`fbm_acceleration_usage.md`, including its verification protocol.

## What it is

The Chimera method (arXiv:2506.22831) augments the background Q2/P1
fictitious-domain solve with small body-fitted "atmosphere" submeshes
around each particle. Submeshes solve their own Navier–Stokes problems
(rigid-body Dirichlet inside, Robin coupling outside) and provide
surface-stress force/torque integration — substantially more accurate
and smoother than the FBM volume-integral force at the same background
resolution. Variants: Chimera-S (strong hole/fringe Dirichlet
constraints; accuracy champion for stationary configurations) and
Chimera-W (weak interior-penalty coupling; smooth forces for moving
particles).

## Activation model — no build flag

There is **no CMake option and no preprocessor flag**. The component is
always compiled; the runtime key `SimPar@ChimeraEnable` (default `No`)
is the only switch. Consequences:

- One binary serves both modes; a deck without Chimera keys behaves
  bit-identically to the pre-Chimera code (this is regression-gated:
  `q2p1_fc_ext_cylinder` baseline at tolerance 0, and no Chimera echo
  lines appear in `prot.txt` unless enabled).
- `ChimeraEnable = Yes` requires an application that calls
  `Chimera_Initialize` (the dedicated `q2p1_chimera` app, Phase 3). Any
  other application aborts at the first time step with a clear message.

## Parameters

See `parameter_reference.md`, section "Chimera Overlapping-Mesh
Component", for the full key table (`ChimeraVariant`, `ChimeraOuterBC`,
`ChimeraParticleFile`, `ChimeraSubmeshFile`, `ChimeraSubmeshLev`,
`ChimeraRobinAlpha`, `ChimeraGammaMax`, `ChimeraOuterIters`,
`ChimeraSubNL`, `ChimeraWriteVTK`). Invalid values abort at parse time
(`CHIMERA_VALIDATE_CONFIG`).

## Verification protocol (the fbm_acceleration_usage.md pattern)

1. **Off-regression** (any phase): run the default deck of
   `q2p1_fc_ext_cylinder` on a current build — results and protocol
   output must match the committed baseline exactly. This is the
   machine-checkable form of "the standard operational mode is not
   disturbed".
2. **Unit tests**: `ctest -R chi-` in a `BUILD_TESTING=ON` build runs
   the Phase-1 service tests (`chi-geometry-serial`, `chi-eval-serial`,
   `chi-locator-serial`, `chi-sparse-direct-serial`).
3. **Physics acceptance** (Phase 3+): `q2p1_chimera` on the steady DFG
   flow-around-cylinder configuration; `ChimeraForce:` C_D/C_L pinned
   against the DFG reference band, with the FBM baseline as cross-check
   on identical background meshes.

## Component layout

`source/src_chimera/` (see its README.md for the layer rules and file
inventory); hooks into existing code are limited to one guarded call in
`Transport_q2p1_UxyzP_fluid_core`, the parser cases, and CMake wiring —
everything else is behind the `CHIMERA_API` facade.
