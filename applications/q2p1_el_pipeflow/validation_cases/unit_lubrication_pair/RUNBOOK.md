# unit_lubrication_pair — cross-rank lubrication pair visibility + fold-once

Unit test for the pairwise lubrication closure (pe commits 916600a +
24295f5): two spheres (d_p = 0.05) straddle the z = 1/3 PE rank boundary
(3x3x3 decomposition) at surface gap 0.0075 = 0.3*R — inside the
lubrication cutoff (0.025) but NOT in contact, and with NEITHER sphere
geometrically overlapping the neighbor subdomain (center distance to the
plane 0.02875 > R). Without the shadow-copy margin the pair is invisible
to both ranks; with the margin (= cutoff) both ranks see it. Frozen
linear shear G = 1 drives relative sliding through the drag closure, so
the sliding lubrication term produces a nonzero virial.

## Gates
1. EL_SUSP_STRESS pairs > 0 every audit step  -> shadow-copy margin works.
2. EL_SUSP_STRESS sig nonzero                 -> virial accumulates.
3. EL_NEWTON_PAIR at machine zero             -> fold-once antisymmetry
   exact across the rank boundary.
4. No contact fires (max_overlap = 0, gap stays positive).

## Run log
(appended as runs complete)

### 2026-07-26 — job 135994 (np=28, 8 s): ALL GATES PASS
- pairs = 20 every audit step (1 pair x 10 substeps, both boundary ranks
  counting) — without the shadow-copy margin this is 0, since neither
  center is within R of the boundary plane.
- sig_xz ~ 7e-7 nonzero, dissipative sign, sliding-dominated structure.
- EL_NEWTON_PAIR worst 1.8e-22 — fold-once antisymmetry exact.
- max_overlap 0.0 throughout (no contact; pure lubrication pair).
