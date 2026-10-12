# FBM query testing, Layer 1: conclusion

Response to `fbm_acceleration_review_2026-10-11.md`. Layer 1 tests the **classification
query** — `pointInsideParticlesAccelerated()` against `pointInsideParticles()` on the same
frozen state — serially and under MPI. Harness: `tools/fbm_query_stress/`.

The layer exists because of one structural fact: every downstream consumer (the KVEL
candidate sets, the alpha field, the forces) reads the *output* of this query, so a
Standard-versus-KVEL force comparison cannot detect a classification miss. Neither could
the campaign's own twins: e4_l3 has one particle and DKT two, below pe's
`gridActivationThreshold = 32`, so the grid is never active in either.

## What was found and fixed

| finding | status |
|---|---|
| P1 stale spatial cells | fixed: `HashGrids::refreshSpatialIndex()`, wired into FeatFloWer after `step_simulation()` |
| P2 overlap ownership | fixed: minimum system id in both query paths |
| `*inpr` written differently by the two paths | fixed: both write the solid flag |
| the Fortran verifier compared only the top byte of the system id | fixed: compares all eight |

**P1 is a broad-phase defect, not a shape problem.** A body is filed by its AABB *minimum*
(`HashGrids.h:909`) and a point query visits 27 cells (`:223`). A body whose minimum sits
just below a cell boundary therefore has its far surface nearly two cells from where it is
filed, and `findContacts()` files bodies *before* contact resolution integrates positions.
An interior point near the advancing surface comes back as fluid. The failure is governed by
where the minimum sits relative to the cell lattice, **not by how far the body moved**: it
reproduces after a displacement of 0.035 % of a body diameter. No CFL-style or
displacement-bound argument can rule it out; the fix has to be geometric.

Sphere, box, capsule, ellipsoid and cylinder all lose the same four points at the same
coordinates once their AABB fills the cell span. The cylinder appeared immune in the first
version of the test only because its AABB was 0.990 of the span against 1.000 for the others
— a measurement artefact, now removed by placing bodies by their real AABB minimum.

**Why re-filing, and what it costs.** `HashGrid::update()` recomputes a body's hash and
returns immediately when it is unchanged, so a refresh is one hash and one compare per body
plus cell moves only for bodies that actually crossed. For 900 bodies that is ~900 hashes
per step. The alternative — a validity guard with an exhaustive fallback — degenerates:
with many moving bodies *some* body crosses a cell every step, so the fallback is permanent,
at ~100k DOF x 900 bodies of containment tests. Because the refresh leaves the index current,
the next `findContacts()` finds nearly every hash unchanged, so the work is moved rather
than duplicated. The refresh is the existing update phase factored out and called by
`findContacts()` itself, so there is one code path.

**Ownership.** A point inside several bodies received "the first one found", which depended
on traversal order — and the two paths traverse differently, so they disagreed on 21 of 42
sampled points inside two overlapping spheres. Both now take the minimum system id, which is
order-free. Neither path can stop at the first hit any more.

## Parallel

The same harness under MPI at 2 and 4 ranks, each stage run twice (control without the
refresh, verdict with it):

| stage | control | refreshed |
|---|---|---|
| M1 bodies away from the cut | ok | ok |
| M2 body on the cut, neighbour holds a shadow | ok | ok |
| M3 body crossing, queried before the next contact search | 9 missed (2 ranks), 42 (4 ranks) | **0** |
| M4 after migration | ok | ok |
| M5 five shapes, rotated, on the cut | ok | ok |

So shadow copies are indexed correctly and migration is handled; the only parallel failure
was the same staleness, amplified because a body crosses several cells in one step. The
single refresh call closes it at both decompositions with no per-rank special casing, and
FeatFloWer already calls it in the shared Fortran path, so parallel PE inherits the fix.

## Lessons that constrain Layers 2 and 3

1. **Assert the preconditions of the test.** The grid-active check caught three cases where
   the harness passed while testing nothing: background bodies spread across a cut left both
   ranks under the 32-body threshold, and a lattice wider than a 4-rank slab did it again.
   The threshold is **per rank**. A Layer-2 gate must assert that it is exercising the code
   it claims to.
2. **Uniform randomisation is not a generator here.** Stages A, C and twelve random
   configurations pass with zero mismatches over ~20000 queries; the defect needs the
   cell-boundary coincidence. A random generator must be biased toward boundaries.
3. **Compare identities, not booleans.** Ownership differences are common rather than exotic
   (half the sampled overlap points), and the existing verifier could not see them.
4. **No rank's silence is agreement.** Verdicts reduce over all ranks; this is the discipline
   the review asked for in P4 and it is already implemented in the MPI harness.
5. **Headers and library must come from the same pe source.** A test TU that disagrees with
   the library links cleanly and is undefined behaviour.
6. **The Fortran verifier is dead code**: it compiles only under `VERIFY_HASHGRID`, which no
   normal build defines. Layer 2's first task is making the gate exist before making it
   strict.

## Not addressed here

The P2 *representation* question — a DOF inside several bodies still stores one owner.
A sketch for an exception list (canonical owner plus a short CSR side table of multiply-owned
DOFs) is recorded in the session notes; it keeps the hot path unchanged and would let a
Layer-2 gate compare owner *sets* rather than winners.

Roadmap item 6 (seed-and-flood) and the conservative AABB enumeration the review argues for
are untouched, as is the campaign A/B (`UseHashGridAccel` on/off on a viscometer dump), which
is deliberately last.
