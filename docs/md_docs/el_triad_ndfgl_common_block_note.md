# `/TRIAD/`, `NDFGL`, and Hidden Mesh-Level State

This note documents a bug pattern that appeared in the Euler-Lagrange transfer
work and is easy to reintroduce when using old FEAT routines manually.

## What `/TRIAD/` Is

`/TRIAD/` is the old FEAT global "current mesh topology dimensions" COMMON block:

```fortran
COMMON /TRIAD/ NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED, ...
```

The most important entries are:

- `NEL`: number of elements on the current level
- `NVT`: number of vertices on the current level
- `NET`: number of edges on the current level
- `NAT`: number of faces/areas on the current level
- `NVE`, `NEE`, `NAE`: per-element topology counts

In this codebase, `SETLEV` updates `/TRIAD/` from `mg_mesh%level(ILEV)`:

```fortran
NVT = mg_mesh%level(ILEV)%nvt
NET = mg_mesh%level(ILEV)%net
NAT = mg_mesh%level(ILEV)%nat
NEL = mg_mesh%level(ILEV)%nel
...
```

So `/TRIAD/` means "the mesh dimensions for whatever multigrid level was last
selected with `ILEV; CALL SETLEV(...)`".

## What It Was Good For

Historically, `/TRIAD/` avoided passing a large argument list through old
Fortran assembly routines. Legacy routines such as element assembly, `NDFGL`,
basis evaluation, surface tension, diffusion, and mass matrices could all read
the current mesh dimensions from global state.

That was convenient in a Fortran 77 style codebase where:

- mesh data lived in global work arrays,
- multigrid levels were selected globally,
- element routines had narrow, fixed signatures,
- performance and minimal argument passing mattered,
- there was one implicit "current level" at a time.

For Q2/E013 velocity unknowns, the global DOF layout is assembled from mesh
entity indices:

```text
vertex DOFs:  1 ... NVT
edge DOFs:    NVT+1 ... NVT+NET
face DOFs:    NVT+NET+1 ... NVT+NET+NAT
cell DOFs:    NVT+NET+NAT+1 ... NVT+NET+NAT+NEL
```

For example, the cell-center DOF for element `IEL` is often computed as:

```fortran
NVT + NET + NAT + IEL
```

That is why routines such as `NDFGL` need `/TRIAD/`: even if the connectivity
arrays are passed explicitly, the current global offsets `NVT`, `NET`, and
`NAT` are still needed to construct global Q2 DOF numbers.

## What `NDFGL` Does

`NDFGL` maps an element-local finite-element DOF numbering to global DOF
indices. Conceptually, for a given element `IEL` and element type `IELTYP`, it
fills:

```fortran
KDFG(i) = global DOF index of local basis function i
KDFL(i) = local FE basis index / local DOF id used by E013/dbas
```

Typical assembly code then uses:

```fortran
CALL NDFGL(IEL, 1, IELTYP, KVERT, KEDGE, KAREA, KDFG, KDFL)

DO i = 1, IDFL
  ig = KDFG(i)   ! global row/DOF index
  il = KDFL(i)   ! local basis slot
  ...
END DO
```

For Q2/E013 hex elements, `NDFGL` combines the explicit element connectivity
with implicit `/TRIAD/` offsets:

```text
vertices: KVERT(local_vertex, IEL)
edges:    NVT + KEDGE(local_edge, IEL)
faces:    NVT + NET + KAREA(local_face, IEL)
cells:    NVT + NET + NAT + IEL
```

The exact local ordering is FEAT/E013-specific and must match what `E013`
stores in `DBAS`.

## Why This Causes Headaches

`/TRIAD/` is dangerous because it is hidden state.

A routine can look like it is using explicit mesh arrays:

```fortran
CALL NDFGL(iel, 1, ieltyp, kvert, kedge, karea, kdfg, kdfl)
```

but internally it also depends on global `/TRIAD/`. Therefore correctness
depends on something outside the function signature.

The real dependency is:

```text
correct NDFGL result =
  explicit kvert/kedge/karea arguments
  + implicit /TRIAD/ offsets
  + implicit current finite element type / local ordering
```

If `/TRIAD/` belongs to level 1 but the arrays belong to level 2, the call still
compiles and may run, but it computes wrong global DOF indices.

## The Euler-Lagrange Bug Pattern

The Euler-Lagrange transfer code correctly passed explicit connectivity from the
coupling mesh level:

```fortran
mesh%level(ilev)%kvert
mesh%level(ilev)%kedge
mesh%level(ilev)%karea
```

It also sized its EL force field from the same level:

```fortran
ndof = mesh%level(ilev)%nvt + mesh%level(ilev)%net + &
       mesh%level(ilev)%nat + mesh%level(ilev)%nel
```

However, `/TRIAD/` still contained stale values from a previous solver phase.
`NDFGL` then built Q2 global indices with the wrong offsets.

Example failure mode:

```text
force_rhs size = mesh%level(ilev)%nvt
               + mesh%level(ilev)%net
               + mesh%level(ilev)%nat
               + mesh%level(ilev)%nel

but NDFGL computes center DOF as:

  stale_NVT + stale_NET + stale_NAT + iel
```

If the stale offsets are from a finer or different level, `KDFG` can exceed
`SIZE(force_rhs,2)`. Then force spreading writes outside the intended EL RHS
layout or effectively loses the feedback contribution. That broke the Newton
pair and showed up as a momentum conservation error.

## Implemented Fix

The current EL fix saves the current `/TRIAD/`, pins it to the EL coupling
level for the particle-mesh pass, then restores it:

```fortran
saved_tr_nel = el_tr_nel
saved_tr_nvt = el_tr_nvt
saved_tr_net = el_tr_net
saved_tr_nat = el_tr_nat

el_tr_nel = mesh%level(ilev)%nel
el_tr_nvt = mesh%level(ilev)%nvt
el_tr_net = mesh%level(ilev)%net
el_tr_nat = mesh%level(ilev)%nat

! EL_INTEGRATE_PARTICLE / EL_DEPOSIT_PARTICLE / NDFGL calls

el_tr_nel = saved_tr_nel
el_tr_nvt = saved_tr_nvt
el_tr_net = saved_tr_net
el_tr_nat = saved_tr_nat
```

This is the right kind of fix when legacy `NDFGL` must be used.

## Recommended Rule Going Forward

Do not use `/TRIAD/` directly in new code unless forced by old FEAT routines.

If code calls `NDFGL`, `E013`, or old assembly routines, assume they read COMMON
state. Set or pin that state explicitly near the call.

For new code, prefer:

- pass `mesh%level(ilev)` explicitly,
- compute `ndof`, offsets, and connectivity from the passed mesh,
- avoid global `NVT`, `NET`, `NAT` values,
- assert `MINVAL(KDFG) >= 1` and `MAXVAL(KDFG) <= ndof` after legacy mappings,
- wrap legacy calls in a small helper that saves, pins, validates, and restores
  `/TRIAD/` and any other required level state.

## Short-Term Safer Wrapper

The low-risk path is a guarded wrapper around `NDFGL`. It keeps FEAT's local
ordering behavior intact while making the hidden `/TRIAD/` dependency explicit:

```fortran
SUBROUTINE EL_NDFGL_SAFE(mesh, ilev, iel, ieltyp, kdfg, kdfl, idfl)
  USE def_FEAT, ONLY: tri_nel => NEL, tri_nvt => NVT, &
                      tri_net => NET, tri_nat => NAT
  TYPE(tMultiMesh), INTENT(IN) :: mesh
  INTEGER, INTENT(IN) :: ilev, iel, ieltyp
  INTEGER, INTENT(OUT) :: kdfg(:), kdfl(:), idfl

  INTEGER :: saved_nel, saved_nvt, saved_net, saved_nat
  INTEGER :: ndof

  saved_nel = tri_nel
  saved_nvt = tri_nvt
  saved_net = tri_net
  saved_nat = tri_nat

  tri_nel = mesh%level(ilev)%nel
  tri_nvt = mesh%level(ilev)%nvt
  tri_net = mesh%level(ilev)%net
  tri_nat = mesh%level(ilev)%nat

  idfl = NDFL(ieltyp)
  CALL NDFGL(iel, 1, ieltyp, mesh%level(ilev)%kvert, &
       mesh%level(ilev)%kedge, mesh%level(ilev)%karea, kdfg, kdfl)

  ndof = tri_nvt + tri_net + tri_nat + tri_nel
  IF (MINVAL(kdfg(1:idfl)) < 1 .OR. MAXVAL(kdfg(1:idfl)) > ndof) THEN
    ERROR STOP 'NDFGL returned out-of-range global DOF index'
  END IF

  tri_nel = saved_nel
  tri_nvt = saved_nvt
  tri_net = saved_net
  tri_nat = saved_nat
END SUBROUTINE
```

Production details to add before using this broadly:

- restore `/TRIAD/` on every normal return path,
- include `iel`, `ilev`, `ieltyp`, `ndof`, `MINVAL(KDFG)`, and `MAXVAL(KDFG)`
  in error diagnostics,
- save and restore any other legacy globals required by the called routine,
- avoid putting a heavy wrapper in the hottest legacy assembly loops unless
  profiling says the cost is acceptable.

## Long-Term Explicit Mapper

A new explicit mapper can be written for the element types used by EL, but it
must be tested against `NDFGL` before replacing it. For Q2/E013 hexes the global
index formulas are straightforward:

```fortran
vertex = kvert(iv,iel)
edge   = nvt + kedge(ie,iel)
face   = nvt + net + karea(ia,iel)
cell   = nvt + net + nat + iel
```

The risk is local ordering. `KDFL` must match what `E013` puts in `DBAS`. If the
ordering is wrong, assembly can still run but couple the wrong basis functions.

Recommended path:

1. Use a guarded `NDFGL` wrapper now.
2. Add a test that compares any explicit mapper against `NDFGL` on known meshes,
   all relevant levels, and all EL-used element types.
3. Only switch EL to explicit mapping after the test proves `KDFG/KDFL` identity.

In short: `/TRIAD/` was useful as a compact global current-level mesh descriptor
for old FEAT assembly. In modern mixed-level/multiphysics code it is a hidden
global dependency and should be isolated behind wrappers or replaced by explicit
mesh-level arguments.
