# Runbook: unit_migration_hydro (mid-step ownership migration Newton-pair audit)

**Purpose:** fast unit-style reproducer for the PE momentum leak caused by
particle ownership migration (fixed in the same libs/pe commit as the
cross-rank contact Newton-pair fix). A single sphere (r = 0.005) is seeded at
z = 0.333533, just above the z = 1/3 PE partition plane of the shared QBOX9
3x3x3 periodic box, and settles under Prop@Gravity z through the plane
(crossing around step 7). Fluid feedback is OFF, difelice drag,
`ELMomentumAuditFreq = 1`, 20 steps of dt = 1e-3 with 2 PE substeps.
No contacts are involved at all.

**np = 28** (27 workers + 1 master). Run exactly like
`unit_straddling_contact` (see its RUNBOOK), staging this case instead.

## What it catches (broken PE, before the fix)

1. On migration the old owner's body is demoted to a shadow copy but kept its
   armed Euler-Lagrange hydro state; `clearELHydroStates()` clears OWNED
   bodies only, so the stale armed shadow kept computing the full semi-implicit
   hydro dv every substep and shipping it to the new owner as a velocity
   correction: the particle received ~1.9x gravity/drag forever after the
   crossing (audit: constant mismatch_z = -4.6e-6 per step from step 8 on).
2. At the crossing step itself, the promoted body was unarmed for the second
   substep (the CFD driver arms only once per PE step, on the pre-step owner),
   giving a +2.6e-7 blip at step 7.

## Acceptance

`EL_NEWTON_PAIR mismatch` at machine zero (|mismatch| < 1e-15, typically
<= 3e-20) for EVERY step, including the crossing step, and the per-step
velocity increment stays ~ -9.8e-3 (semi-implicit g*dt) across the crossing.

## Reference results

- Broken PE: step 7 mismatch_z = +2.57e-7, steps 8..20 mismatch_z = -4.6e-6
  each; particle Delta-v jumps to -1.86e-2/step (1.9x) after crossing.
- Fixed PE: all 20 steps |mismatch| <= 3e-20; Delta-v = -9.80e-3/step
  throughout.
