# Runbook: d0_smoke_twin — SED_BENCH retirement cross-binary equivalence (campaign D0-SMOKE)

Goal: prove that retiring the compile-time `-DSED_BENCH` hack in favor of
runtime `Prop@ForceScale` + `SimPar@PrintParticleState` changes nothing in
the physics, byte for byte.

## Configuration

- Case: ten Cate E4 quarter-box sedimentation, 10-step smoke
  (`SimPar@MaxNumStep = 10`, dt = 0.001, level 2), serial-PE mode, np = 32
  (31 METIS parts), mesh `_mesh/NEWFAC` staged from the
  `build-hcaf-pe-serial` reference app dir.
- Old binary: `build-hcaf-pe-serial/applications/q2p1_bench_sedimentation`
  (built 2026-07-31 08:14 from `master` @ 5de38823, `SED_BENCH=ON`,
  `pe_CONSTRAINT_SOLVER = HardContactAndFluid` via in-tree Collisions.h
  state).
- New binary: `build-dns-pe-serial/...` (campaign branch @ d0b9f0f6,
  `SED_BENCH` retired,
  `-Dpe_CONSTRAINT_SOLVER=pe::response::HardContactAndFluid` via
  `CMAKE_CXX_FLAGS`).
- Decks identical except: new deck adds `SimPar@PrintParticleState = Yes`
  (the key the old binary does not know). Both carry
  `Prop@ForceScale = 0d0,0d0,4d0,0d0,0d0,0d0` (ignored by the old binary,
  applied by the new one — the hack and the runtime path express the same
  arithmetic).
- Rundirs (evidence, do not delete): `q2p1_dns_rundir_twin_old/`,
  `q2p1_dns_rundir_twin_new/` at repo root; logs `run_twin_old.log`,
  `run_twin_new.log`.

## Gates and results (2026-07-31)

| Gate | Result | Verdict |
|---|---|---|
| Old binary re-run vs head of 1300-step reference log (`run_sed_hcaf_np32.log`) | 10/10 `SED_BENCH_VEL` lines byte-identical | PASS (determinism) |
| New binary vs old binary | 10/10 `SED_BENCH_VEL` lines byte-identical | PASS (equivalence) |
| Unified tags emitted | `DNS_PART_STATE`/`DNS_PART_FORCE` present, trajectory matches legacy lines at ES16.8 | PASS |

## Finding recorded along the way (D0-SOLVER)

The first new-binary run used a tree inheriting pe master's default
`pe_CONSTRAINT_SOLVER = HardContactEulerLagrange` and produced an
**all-zero trajectory** (force and velocity exactly 0 for all 10 steps,
no error raised): FBM forces written via `setForcesMapped` are a dead end
in the EL solver, which consumes hydro input only as `dv_` velocity
corrections. This empirically confirms campaign decision D0.6 (DNS and EL
must not share a CollisionSystem) and defines the failure signature:
a silently motionless particle means *wrong solver*, not a physics bug.
Interim DNS solver until `HardContactDNS` exists: `HardContactAndFluid`,
selected per build tree with
`cmake -DCMAKE_CXX_FLAGS="-Dpe_CONSTRAINT_SOLVER=pe::response::HardContactAndFluid"`.

## Known gaps (documented, not tuned around)

- `SED_BENCH_FORCE` x/y columns may print `-0.000000E+00` under the
  runtime path (0d0 × negative partial force) where the old hack printed
  hard zeros; velocity/position lines are unaffected. Not observed in the
  10-step window; noted as a permissible formatting-only difference.
- The `id=` column of `DNS_PART_STATE` is the global list index in
  serial-PE mode but the PE systemID low byte in parallel-PE mode;
  identical for single-particle cases. Unify when a multi-particle
  parallel case first needs it.
