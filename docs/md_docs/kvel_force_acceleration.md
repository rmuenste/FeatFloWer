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
END TYPE tVertexCache

TYPE(tVertexCache), ALLOCATABLE :: ParticleVertexCache(:)
```

Allocated as `ParticleVertexCache(1:N_particles)`. For each particle it holds
the list of local DOF indices (corner vertices nvt, edge midpoints net, face
midpoints nat) that are inside that particle. This allows the force routine to
look up "which DOFs belong to particle IP" in O(1) rather than scanning the
full DOF array.

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

**Usage in the force routine:** the element boundary check at each cubature
point uses:

```fortran
IF (longIdMatch(IG, theParticles(IP)%bytes)) THEN
  DALPHA = 1d0
ELSE
  DALPHA = 0d0
END IF
```

where `theParticles(IP)%bytes` is the 8 × c_short representation of the PE
system ID for particle `IP`, obtained via `getAllParticles()`.

---

## Why `FictKNPR(i) == IP` Does Not Work in the PE Path

The Wangen FBM path (`QuadScalar_FictKnpr_Wangen`) processes particles in
a loop `DO IPP = 1, myFBM%nParticles` and calls geometry functions that
set `FictKNPR(dof) = IP` for the owning particle index. Cache building can
therefore use `FictKNPR(i) == IP` directly.

The PE path (`QuadScalar_FictKnpr` with `fbm_getFictKnprFC2`) only sets
`FictKNPR(dof) = 1` for any inside DOF — there is no particle-index
discrimination in the integer flag. All particle identity information lives
in `FictKNPR_uint64`. The cache must therefore use `longIdMatch`:

```fortran
! WRONG for PE path:
if (FictKNPR(i) == IP) then ...

! CORRECT for PE path:
if (FictKNPR(i) /= 0 .and. longIdMatch(i, cacheParticles(IP)%bytes)) then ...
```

---

## Cache Building (`QuadScalar_FictKnpr`)

The cache is rebuilt each timestep immediately after the alpha field computation
loop, inside `if (myid /= 0)` (i.e. on all worker ranks). Steps:

1. **Get particle list** via `numTotalParticles()` + `getAllParticles()` —
   same ordering as `ForcesLocalParticlesSerial` will use later.
2. **Allocate** `ParticleVertexCache(1:numCacheParticles)`.
3. **Pass 1 — count:** for each particle `IP`, scan all corner, edge and face
   DOFs. For each DOF where `FictKNPR(i) /= 0 .and. longIdMatch(i, cacheParticles(IP)%bytes)`,
   increment `nVertices`. Then allocate `dofIndices(nVertices)`.
4. **Pass 2 — fill:** repeat the scan and store the DOF indices.

The `FictKNPR(i) /= 0` pre-check is a cheap early-out; `longIdMatch` is only
called for DOFs that are inside *any* particle.

Location: `source/src_quadLS/QuadSc_boundary.f90`, subroutine
`QuadScalar_FictKnpr`, after the totalInside counting loop (inside
`#ifdef HAVE_PE`).

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
particle surface.

If `nCandidates == 0` (particle not on this subdomain), all elements are used
as fallback — this is correct since such ranks have no inside DOFs and the
element loop will skip everything via the `NJALFA == 27` check anyway.

Location: `source/src_quadLS/QuadSc_force_serial.f90`, lines ~201–269.

---

## MPI Reporting

Because the print `if (myid == 1)` only shows rank 1's local cache (which may
be empty if the particle is in a different subdomain), the cache count and force
stats are aggregated across ranks:

- **Cache size:** `MPI_Reduce(k, reducedVal, 1, MPI_INT, MPI_SUM, 0, ...)` in
  `QuadScalar_FictKnpr` after `end if ! myid /= 0`. Printed by rank 0.
- **Candidate elements:** `COMM_SUMMN` on `myKVEL_Stats%nCandidateElements`
  after `ForcesLocalParticlesSerial` returns. Printed by rank 0.

---

## Runtime Control

```fortran
LOGICAL :: bUseKVEL_Accel = .FALSE.   ! source/src_quadLS/QuadSc_var.f90
```

Set to `.TRUE.` to enable the optimization. When `.FALSE.` the candidate
building block is skipped and the force loop falls through to the full
element scan — identical physics, useful for correctness comparison.

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
  └─ [HAVE_PE] Build ParticleVertexCache
       ├─ getAllParticles(cacheParticles)       ← same order as force routine
       ├─ for each particle: longIdMatch(i, cacheParticles(IP)%bytes)
       └─ ParticleVertexCache(IP)%dofIndices(:) = inside DOFs for particle IP

ForcesLocalParticlesSerial (each timestep)
  │
  ├─ getAllParticles(theParticles)              ← same order as cache
  │
  └─ for each particle IP:
       ├─ use ParticleVertexCache(IP)%dofIndices to look up inside DOFs
       ├─ use kvel/keel/kaal to find adjacent candidate elements
       └─ integrate force only over candidate elements
            (full element scan if cache is empty on this rank)
```
