# D2.2 sub-grid lubrication — design specification (owner-approved direction)

Status: direction approved by owner 2026-08-20 ("this is how it should be
done"); implementation NOT started. Prerequisite: the pe merge-repair branch
`fix/lubrication-merge-resolution` must land upstream first (see the
lubrication-merge audit, datasheet-adjacent memory, and
`~/nobackup/code/FF-EL/pe-lubfix-worktree/PR_DESCRIPTION.md`).

## 1. Philosophy (owner-set)

Lubrication is contact-scale particle physics and lives in pe as a
**switchable C++ add-on component usable by any collision pipeline** —
unit-testable in isolation, configured from json + mesh metadata.
No physics on the Fortran side; FeatFloWer supplies only grid metadata
(h_min) through one plumbing call. FF's legacy `ENABLE_LUBRICATION`
ifdef stays dead.

## 2. Disposition of pe PR #23 (three layers)

| layer | content | verdict |
|---|---|---|
| 1 | `pe/core/lubrication/LubricationModel.h` + 4 unit tests: pure closure, Kroupa 2016 eqs 12–15 pair + 16–19 wall, Vinogradova slip, epsCritical saturation, semi-implicit integrator, modelLegacy regression path | **KEEP** — this is the component |
| 2 | detection/AABB machinery: `MaxContacts` lubrication-contact branches (sphere–sphere, sphere–plane), `aabbPadding()` inflation | **KEEP, PROMOTE** the gate from compile-time `PE_LUBRICATION_CONTACTS` to the runtime config switch so any pipeline can use it |
| 3 | `HardContactLubricated` isolated pipeline + default-solver flip in `Collisions.h` | **RETIRE** once the stage is extracted; re-point its tests at the runtime switch (also resolves the fix-branch test-skip issue) |

Wholesale revert rejected: it would discard tested layers 1–2 to remove
layer 3.

## 3. The refactor PR (pe-side, one PR)

1. **Extract the application stage** from `HardContactLubricated.h:2221-2472`
   (config snapshot → loop tagged contacts → apply impulses) into a shared
   `pe/core/lubrication/LubricationStage.h`, invoked after
   `synchronizeVelocities()`, gated on `lubrication::isEnabled()`.
2. **Hook `HardContactAndFluid`** with that call. Switch off ⇒ provably
   identical behavior (Gate 0 below). Switch on ⇒ Kroupa forces in the
   campaign solver — which also dissolves the "json switch silently does
   nothing" footgun.
3. **Runtime-gate the detection branches** in `MaxContacts.h` (replace the
   compile-time macro), both sphere–sphere and sphere–plane.
4. **AABB padding = the EFFECTIVE cutoff**: min(cutoffFactor·aRef,
   clampFactor·meshDx) when enabled, 0 when disabled — never inflate more
   than the physics can act on. (Also removes the legacy unconditional
   `lubricationThreshold = 1e-2` absolute-units padding at HEAD,
   `SphereBase.h:222` — confirmed present and un-overridden in all DNS
   setups; practically negligible in d=1 cases, > radius in SI cases.)
5. **Loud refusal** if `lubricationEnabled_: true` reaches a solver/build
   that cannot honor it.
6. **Retire `HardContactLubricated`**; keep `Collisions.h` default solver
   neutral (`HardContactEulerLagrange`).
7. **One serial-mode + CFD-coupled test** (the coverage gap: zero such
   tests exist). The stage behaves like contact response, which serial PE
   already handles deterministically.

## 4. FF-side wire (plumbing only)

One interface call at setup passing h_min (FF already computes it — the
`DNS_RESOLUTION` print) → `lubrication::setMeshDx()`. This activates the
already-implemented-but-dead `lubricationMeshClampFactor_` clamp
(`setMeshDx` currently called only from a unit test). With clamp factor
c = 2 the activation gap is **gap < 2h automatically per resolution** —
the D2.1 crossover rule (resolved FBM tracks Brenner to ~10% down to
2–3 cells, −20% at ~1.5 cells; departure collapses in gap/h) wired in,
instead of a hand-computed cutoff factor per case.

## 5. Setup procedure (the "robust by construction" requirement)

- Activation: mesh clamp, default c = 2 (deck-overridable).
- Physical parameters (viscosity, slip length) from json; must mirror the
  deck fluid properties — add a startup consistency check in the spirit of
  the pe_stepsize_mismatch guard (nothing checks fluid-property agreement
  today, `d0_visc`).
- Model selection: `kroupa2016` (full) | `legacy` (ten Cate Eq.-10-class
  regression) | later `kroupaDeficit` (see §7).
- Defaults: lubrication OFF everywhere; enabling is an explicit json act.

## 6. Validation ladder

| gate | what | pass criterion |
|---|---|---|
| G0 | lubrication-OFF twins with refactored pe, vs certified logs — AMENDED 2026-08-20 after the hcaf_angvel_reset fix (a71c34d) joined the branch: (a) **e4_l3** stays the BITWISE instrument (settling, spin ≈ 0, ω-fix inert); (b) **DKT 20-step** is now a physics-DIFF case — the ω restoration changes rotating trajectories by construction, so its twin is EXPECTED to differ; its deviation is the ω-restoration signal and feeds the D2.3 frictional rerun (datasheet hcaf_angvel_reset) | e4_l3 bitwise; DKT diff attributed to rotation (report, don't gate); SI ten Cate family remains the AABB-sensitive one |
| G1 | unit tests: existing 4 re-pointed at the runtime switch + the new serial/CFD-coupled test | all pass |
| G2 | **Brenner wall-approach rerun** (D2.1 fixture, constant-V protocol, Re 0.78) with lubrication ON | sub-2h force follows Brenner's exact solution where FBM alone was −20%; ALSO measure the 2–3h overlap band overshoot = the double-counting number that decides §7 |
| G3 | **lubricated ten Cate**: reproduce Fig. 13 (pairs S16/S18 at Re 1.5, S17/S19 at Re 4.1) — first `legacy` model (like-for-like with the paper's Eq. 10), then full Kroupa | early velocity-decay improvement reproduced; Kroupa's epsCritical saturation avoids the paper's documented time-to-contact overprediction ("the abrupt stop remains" pathology) |

## 7. Physics follow-up (option, not blocker)

Full-vs-deficit closure: in resolved FBM the fluid already carries most of
the lubrication force down to ~2h, so the full model double-counts in the
overlap. Deferred by design: G2 measures the actual double-count at
D/h = 8 before deciding. If warranted, add `modelKroupaDeficit` —
Kroupa resistances minus their value at the activation gap (reduces to
ten Cate Eq. 10 in the leading normal term, keeps the saturation fix).
D2.1's inverse-collapse-curve idea remains the more ambitious variant.

## 8. Performance & reproducibility (recorded 2026-08-20 discussion)

- Pair-building rides the existing HashGrids broad phase (AABB inflation →
  candidate pairs → exact gap in fine detection): O(N)-class, no naive
  distance loops. The switch must control inflation + detection + stage
  coherently (it does, once §3.3 lands).
- AABB-size → hierarchy-level → enumeration-order effects: irrelevant for
  the monodisperse low-N d=1 benchmark class (single size class, sparse
  contacts, single-contact DKT provably order-free). Production caveats:
  polydispersity populates multiple grid levels; dense packings are
  statically indeterminate so the sequential solver's force distribution
  is order-dependent — production gates there must be statistical
  (seed ensembles), never single-trajectory.

## 9. Sequencing

1. Owner pushes `fix/lubrication-merge-resolution`, opens upstream PR (fix
   the broken merge first — refactoring wants a compiling base).
2. pe refactor PR (§3) on top; review + unit tests.
3. FF wire (§4) + setup guard (§5) — small FF-side commit.
4. Validation ladder G0→G3; G2's overlap measurement feeds the §7 decision.
5. Campaign records: datasheet rows per gate; practitioner's guide §5/§6
   update (the gap ≲ 2h regime stops being "contact model only").

Out of scope: EL-solver (HardContactEulerLagrange) changes; FF-side
lubrication physics of any kind; sphere–sphere lubrication validation
beyond the unit tests (no benchmark fixture yet — D2.2's original
sphere–sphere approach case can follow the same ladder later).
