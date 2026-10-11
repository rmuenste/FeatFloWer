# Layer-1 stress test for the FBM classification query

`pointInsideParticlesAccelerated()` (HashGrid) against `pointInsideParticles()`
(exhaustive) on the **same frozen world state**. Everything downstream — the KVEL
candidate sets, the alpha field, the forces — consumes the result of this query, so a
force comparison cannot detect a classification miss. Only this comparison can.

Every stage asserts pe's grid is ACTIVE (`gridActivationThreshold = 32` bodies); a
stage that silently runs without the grid tests nothing, which is why the single- and
two-particle campaign cases never exercised this path.

```bash
# against any FF build tree (libs/pe inside it)
bash tools/fbm_query_stress/fbm_query_stress.build.sh build-twin-fbmaccel /tmp/fbm_stress
(cd /tmp && /tmp/fbm_stress)

# with meshes and DistanceMaps: needs a pe built with CGAL
WITH_CGAL=1 bash tools/fbm_query_stress/fbm_query_stress.build.sh build-petest-cgal /tmp/fbm_stress_cgal
```

A CGAL-enabled pe is built standalone with
`-DPE_USE_CGAL=ON -DCGAL_DIR=/sfw/cgal/6.1.2/gcc13.2.0/lib64/cmake/CGAL`. The test TU
**must** carry pe's own CGAL flags — `-DPE_USE_CGAL -DCGAL_DISABLE_GMP=1 -isystem
<cgal>/include` — because `PE_USE_CGAL` changes the `TriangleMesh` layout (a TU that
disagrees with the library links but is undefined behaviour), and `CGAL_DISABLE_GMP`
is what keeps CGAL's CORE headers compiling under gcc 13.

## Stages

| stage | what it builds | assertion |
|---|---|---|
| A | 12 monodisperse spheres, settled and after integration | must pass |
| B | polydisperse spheres, AABB minima just below cell boundaries, 0.035 % motion | must pass |
| C | sphere, box, capsule, cylinder, ellipsoid, union; rotated; before/after motion | must pass |
| D | triangle meshes, one with `enableDistanceMapAcceleration`, plus a fixed mesh | must pass |
| D2 | the same box geometry with and without a DistanceMap | characterises the band |
| E | overlapping bodies | characterises ownership |
| F | randomised shapes, sizes, rotations, velocities over 3 seeds x 4 steps | must pass |
| G | the adversarial cell-boundary construction, per shape | must pass |

Stages E and D2 report rather than fail: ownership between overlapping bodies is
undefined today, and a DistanceMap is an approximation whose disagreement band is the
thing worth measuring.
