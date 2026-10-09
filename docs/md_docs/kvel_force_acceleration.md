# KVEL Force Computation Acceleration

## Overview

The FBM (Fictitious Boundary Method) force computation in serial PE mode
(`ForcesLocalParticlesSerial`) originally iterates over **all** mesh elements
on every MPI rank for every particle each timestep to find the boundary elements
that contribute to the hydrodynamic force integral. For a particle that occupies
only a small fraction of the domain this is highly wasteful.

The KVEL acceleration replaces the brute-force element loop with a **candidate
set**: only elements that are topologically adjacent to DOFs known to be inside
the particle are tested. This reduces the work from O(NEL × N_particles) to
O(N_boundary_elements × N_particles).

The KVEL code is compiled with `ENABLE_FBM_ACCELERATION`, which is **ON by
default** whenever `USE_PE=ON` (see `fbm_acceleration_usage.md`). It is used
only in **serial PE mode** (`USE_PE_SERIAL_MODE=ON`): the parallel PE force path
(`ForcesLocalParticles` / `ForcesRemoteParticles`) integrates over all elements
and does not read the cache, so the cache is not built there.

---

## Key Data Structures

### `FictKNPR` — Binary Solid/Fluid Flag

```fortran
INTEGER, ALLOCATABLE :: FictKNPR(:)   ! (var_QuadScalar, QuadSc_var.f90)
```

Stores a **binary** occupancy flag for every DOF on the local mesh:

| Value | Meaning |
|-------|---------|
| `0`   | DOF is in fluid domain |
| `1`   | DOF is inside *some* particle |

**Important:** In the PE code path (`fbm_getFictKnprFC2`), `FictKNPR(i)` is
always 0 or 1 — it does **not** encode *which* particle the DOF belongs to.
This is in contrast to the non-PE Wangen path (`GetFictKnpr_Wangen`) where
`FictKNPR(i)` is set to the particle loop index `IP`.

### `FictKNPR_uint64` — Full PE System ID

```fortran
type(tUint64), allocatable :: FictKNPR_uint64(:)  ! (var_QuadScalar)
```

where `tUint64` is defined in `source/src_util/types.f90`:

```fortran
type tUint64
  integer(c_short), dimension(8) :: bytes
end type tUint64
```

This stores the **64-bit PE system identifier** of the particle that owns
each DOF, encoded as 8 × `c_short` (16-bit signed integers). For DOFs in
the fluid domain all bytes are set to `-1`. This representation exists because
Fortran has no native unsigned 64-bit integer type, whereas the PE C++ library
uses `uint64_t` particle IDs.

### `ParticleVertexCache` — Cached DOF Index Lists

```fortran
TYPE tVertexCache
  INTEGER :: nVertices
  INTEGER, ALLOCATABLE :: dofIndices(:)
  INTEGER :: particleID
  TYPE(tUint64) :: longId          ! PE system id of particle IP at build time
END TYPE tVertexCache

TYPE(tVertexCache), ALLOCATABLE :: ParticleVertexCache(:)
```

Allocated as `ParticleVertexCache(1:N_particles)`. For each particle it holds
the list of local DOF indices (corner vertices, edge midpoints, face midpoints
and element centres, in ascending DOF order) that are inside that particle.
This allows the force routine to look up "which DOFs belong to particle IP"
in O(1) rather than scanning the full DOF array.

### `FictKNPR_IP` — Dense Local Particle Index per DOF

```fortran
INTEGER, ALLOCATABLE :: FictKNPR_IP(:)   ! (var_QuadScalar)
LOGICAL :: bKVEL_IndexValid
```

Built together with the cache: `FictKNPR_IP(i) = IP` if `FictKNPR_uint64(i)`
holds the system id of particle `IP` in the `getAllParticles()` ordering,
`0` otherwise (fluid DOF, or a body that is not in the particle list, e.g. a
box). `bKVEL_IndexValid` is `.TRUE.` only after a successful build in the
current timestep. `FictKNPR` and `FictKNPR_uint64` are not modified; all
other consumers (boundary conditions, output, `dem_query`) keep using them.

---

## The `longIdMatch` Function

Defined in `source/src_particles/dem_query.f90`:

```fortran
logical function longIdMatch(idx, longFictId)
  use var_QuadScalar, ONLY : FictKNPR_uint64
  integer, intent(in) :: idx
  integer(c_short), dimension(8) :: longFictId

  longIdMatch = .false.
  if( (FictKNPR_uint64(idx)%bytes(1) .eq. longFictId(1)) .and. &
      (FictKNPR_uint64(idx)%bytes(2) .eq. longFictId(2)) .and. &
      ...
      (FictKNPR_uint64(idx)%bytes(8) .eq. longFictId(8)) ) then
    longIdMatch = .true.
  end if
end function longIdMatch
```

Given a DOF index `idx` and a particle's byte-representation `longFictId`, it
compares `FictKNPR_uint64(idx)%bytes` against `longFictId` byte-by-byte.
This is the canonical way to ask "does DOF `idx` belong to the particle
identified by `longFictId`?".

**Usage in the force routine:** `ForcesLocalParticlesSerial_KVEL` evaluates
the ownership of each local DOF of an element once per element and reuses it
for the boundary-element test (`NJALFA`/`NIALFA`) and for the alpha gradient
in the cubature loop:

```fortran
DO I=1,IDFL
  IG=KDFG(I)
  IF (bUseIdx) THEN
   LOWN(I) = (FictKNPR_IP(IG) == IP)        ! integer compare (dense index)
  ELSE
   LOWN(I) = longIdMatch(IG, theParticles(IP)%bytes)
  END IF
  ...  DALPHA_E(I) = 1d0 or 0d0
ENDDO
```

where `theParticles(IP)%bytes` is the 8 × c_short representation of the PE
system ID for particle `IP`, obtained via `getAllParticles()`. `bUseIdx` is set
per particle only if the index is valid and `ParticleVertexCache(IP)%longId`
still equals `theParticles(IP)%bytes`; in that case the integer compare is
exactly equivalent to `longIdMatch`. The brute-force reference
`ForcesLocalParticlesSerial_Standard` keeps calling `longIdMatch` directly.

---

## Why `FictKNPR(i) == IP` Does Not Work in the PE Path

The Wangen FBM path (`QuadScalar_FictKnpr_Wangen`) processes particles in
a loop `DO IPP = 1, myFBM%nParticles` and calls geometry functions that
set `FictKNPR(dof) = IP` for the owning particle index. Cache building can
therefore use `FictKNPR(i) == IP` directly.

The PE path (`QuadScalar_FictKnpr` with `fbm_getFictKnprFC2`) only sets
`FictKNPR(dof) = 1` for any inside DOF — there is no particle-index
discrimination in the integer flag. All particle identity information lives
in `FictKNPR_uint64`. The cache must therefore be keyed by the system id
(`longIdMatch` semantics), which is what the dense index `FictKNPR_IP` encodes:

```fortran
! WRONG for PE path:
if (FictKNPR(i) == IP) then ...

! CORRECT for PE path (equivalent forms):
if (FictKNPR(i) /= 0 .and. longIdMatch(i, cacheParticles(IP)%bytes)) then ...
if (FictKNPR_IP(i) == IP) then ...
```

---

## Cache Building (`QuadScalar_FictKnpr`)

The cache is rebuilt each timestep immediately after the alpha field computation
loop, inside `if (myid /= 0)` (i.e. on all worker ranks), when
`bUseKVEL_Accel` is set and there are particles. The work is O(N_DOF)
(plus O(N_solid · log N_p) for the lookups) instead of the former
2 · N_p · N_DOF `longIdMatch` scans. Steps:

1. **Get particle list** via `numTotalParticles()` + `getAllParticles()` —
   same ordering as `ForcesLocalParticlesSerial` will use later (no PE step
   happens between classification and the force computation).
2. **Sort the particle ids** (`sortParticleIds`, merge sort over the 8 shorts,
   `dem_query.f90`). If two particles share an id the cache is left
   unallocated and the force routine falls back to all elements.
3. **Pass 1 — index and count:** for every DOF with `FictKNPR(i) /= 0`, find
   `FictKNPR_uint64(i)` by binary search (`findParticleId`) and store the
   result in `FictKNPR_IP(i)`; count the DOFs per particle. Then allocate
   `dofIndices(nVertices)`.
4. **Pass 2 — fill:** one more pass over the DOFs in ascending order appends
   each solid DOF to its particle's list. The lists are identical (content
   and order) to the former per-particle scans.

Location: `source/src_quadLS/QuadSc_boundary.f90`, subroutine
`QuadScalar_FictKnpr`, after the totalInside counting loop (inside
`#if defined(HAVE_PE) && defined(ENABLE_FBM_ACCELERATION) && defined(PE_SERIAL_MODE)`).

With `DEBUG_FBM_OPTIMIZATION`, the routine additionally checks
`FictKNPR_IP(i) == IP` against `longIdMatch(i, id_IP)` for every DOF and
particle (O(N_p · N_DOF), debug builds only) and prints
`DEBUG_FBM: Index verification PASSED/FAILED`.

---

## Candidate Element Building (`ForcesLocalParticlesSerial`)

For each particle `IP` in the force loop, the candidate set is built from the
cache using the mesh connectivity arrays:

| DOF type | Connectivity array | Dimension |
|----------|--------------------|-----------|
| Corner vertex (`ivt <= NVT`) | `mg_mesh%level(NLMAX)%kvel(j, ivt)` | `nvel` entries per vertex |
| Edge midpoint (`NVT < ivt <= NVT+NET`) | `mg_mesh%level(NLMAX)%keel(j, iedge)` | `neel` entries per edge |
| Face midpoint (`NVT+NET < ivt <= NVT+NET+NAT`) | `mg_mesh%level(NLMAX)%kaal(j, iface)` | `naal` entries per face |

A boolean flag array `bCandidateElement(NEL)` prevents duplicates. The result
is `CandidateList(1:nCandidates)` containing only elements touching the
particle surface. Both arrays are allocated once per call; after each
particle only the flags listed in `CandidateList` are reset.

If `nCandidates == 0`:

- **cache populated** (the normal case): the particle has no inside DOFs on
  this rank; the rank contributes zero force for it and skips to the next
  particle (`cycle`), without touching any element.
- **no cache** (`bUseKVEL_Accel = .FALSE.` at build time, duplicate ids, or the
  particle ordering changed since the cache was built): all elements are used
  as candidates, i.e. the brute-force loop.

**Known limitation:** element-centre DOFs (`ivt > NVT+NET+NAT`) are cached but
not mapped to candidate elements. An element whose *only* inside DOF is its
centre is therefore skipped by KVEL while `ForcesLocalParticlesSerial_Standard`
integrates it. This is rare for resolved particles but means KVEL and Standard
can differ by more than round-off in such configurations.

Location: `source/src_quadLS/QuadSc_force_serial.f90`, subroutine
`ForcesLocalParticlesSerial_KVEL`.

---

## MPI Reporting

The cache build itself uses no MPI communication (`QuadScalar_FictKnpr` is
reached with rank-asymmetric control flow; see the deadlock note there), and
the cache size is not reported. The candidate statistics are aggregated:

- **Candidate elements:** `COMM_SUMMN` on `myKVEL_Stats%nCandidateElements`
  inside `ForcesLocalParticlesSerial_KVEL`. Printed by rank 0 as
  `KVEL: <candidates> candidates vs <NEL*N_p> brute-force (<ratio>x speedup)`.

---

## Runtime Control

```fortran
LOGICAL :: bUseKVEL_Accel = .TRUE.   ! source/src_quadLS/QuadSc_var.f90
```

Set from `SimPar@UseKVELAccel = Yes|No` in `q2p1_param.dat` (default `Yes`).
When `No`, no cache is built and `ForcesLocalParticlesSerial` calls the
brute-force `ForcesLocalParticlesSerial_Standard`. Both variants integrate the
same boundary elements (up to the element-centre limitation above), but in a
different element order, so the forces agree to round-off, not bitwise.
With `DEBUG_FBM_OPTIMIZATION=ON` both variants run every step and are
compared (tolerance 1e-10); the Standard result is used.

---

## Data Flow Summary

```
QuadScalar_FictKnpr (each timestep)
  │
  ├─ fbm_updateFBMGeom() for all DOFs
  │    └─ fbm_getFictKnprFC2()
  │         ├─ sets FictKNPR(i) = 0 (fluid) or 1 (solid)
  │         └─ sets FictKNPR_uint64(i)%bytes = PE system ID (or -1 if fluid)
  │
  └─ [serial PE + ENABLE_FBM_ACCELERATION] Build ParticleVertexCache
       ├─ getAllParticles(cacheParticles)       ← same order as force routine
       ├─ sort ids, binary search per solid DOF → FictKNPR_IP(i)
       └─ ParticleVertexCache(IP)%dofIndices(:) = inside DOFs for particle IP

ForcesLocalParticlesSerial (each timestep)
  │
  ├─ getAllParticles(theParticles)              ← same order as cache
  │
  └─ for each particle IP:
       ├─ check ParticleVertexCache(IP)%longId == theParticles(IP)%bytes
       ├─ use ParticleVertexCache(IP)%dofIndices to look up inside DOFs
       ├─ use kvel/keel/kaal to find adjacent candidate elements
       └─ integrate force only over candidate elements, ownership via
            FictKNPR_IP(IG) == IP
            (zero force if the particle has no DOFs on this rank;
             full element scan with longIdMatch if there is no usable cache)
```
