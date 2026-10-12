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

# testing a MODIFIED pe: the headers must come from the same source as the library
WITH_CGAL=1 PE_SRC=../pe-worktrees/<branch> \
  bash tools/fbm_query_stress/fbm_query_stress.build.sh build-petest-fixes /tmp/fbm_stress
```

`PE_SRC` is not a convenience. A test TU that includes different pe headers than the
library was built from links cleanly and is undefined behaviour; the first run of this
test against a modified pe only failed loudly because the change added a new method.

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


## The parallel harness

`fbm_query_stress_mpi.cpp` asks the same question under MPI, where it gains three axes
no serial stage can reach: shadow copies (the exhaustive path walks local bodies *and*
`ghostBodies()`, the accelerated path returns grid order), the `synchronize()` seam, and
migration across a subdomain cut.

```bash
mpicxx ... -Itools/fbm_query_stress tools/fbm_query_stress/fbm_query_stress_mpi.cpp \
   build-petest-mpi/lib/libpe.a ...          # pe built with PE_USE_MPI=ON
mpirun -np 2 ./fbm_query_stress_mpi          # also valid at 4
```

| stage | what it builds |
|---|---|
| M1 | bodies far from the cut — the serial question, in parallel |
| M2 | a body centred on the cut, so the neighbour holds a shadow copy |
| M3 | a body crossing in one step, queried before the next contact search |
| M4 | after migration completes: ownership has changed sides |
| M5 | all five shapes, rotated, straddling the cut |

Every stage runs **twice**: a control pass without the index refresh, which documents the
defect and does not decide the exit code, and a pass with `refreshSpatialIndex()` after
motion, which does. Verdicts are reduced with `MPI_Allreduce`, so no rank's silence counts
as agreement, and the grid-active check is reduced with `MPI_MIN`, so a stage where *any*
rank lacks a grid is reported as having tested nothing.

That last check has caught the harness three times: background bodies straddling the cut
left both ranks under pe's 32-body threshold, and later a body lattice wider than a 4-rank
slab did the same. The bodies are now a compact 4 x 4 x 3 block placed inside each rank's
own slab. **The threshold is per rank, not global.**
