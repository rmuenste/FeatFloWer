# Review v2 of `chimera-integration-design.md`

Recommendation: proceed with Phases 0-1, but do not start the submesh
assembly in Phase 2 yet. The revision resolves the original eight findings,
but two correctness blockers remain.

## Remaining blockers

### 1. The claimed COMMON-free assembly is not possible with the selected legacy kernels

The plan says Layer M is transitively COMMON/`var_QuadScalar`-free, then
calls existing F77 assembly kernels (`chimera-integration-design.md:169-217`).
Those kernels still depend heavily on global state:

- `LAPLACE` obtains `NEL` and FE workspaces from COMMON blocks
  (`source/assemblies/QuadSc_laplace.f:20-29`).
- `CONVQ2` uses COMMON and `var_QuadScalar` configuration
  (`source/assemblies/QuadSc_conv.f:8-47`).
- `Build_BTMatP1` directly imports `transform` and uses COMMON
  (`source/assemblies/QuadSc_BMatrix.f:9-51`).

Passing mesh arrays explicitly is therefore insufficient. Without switching
global state, these routines can use background-mesh dimensions and
configuration; switching it would violate the reentrant, instance-based
architecture.

For clean design, implement new reentrant Chimera assembly kernels using
`chi_geometry`/`chi_fem_eval`, with `nel`, quadrature, basis data, physical
parameters, and CSR passed explicitly. A legacy save/restore adapter is
possible, but it should be considered a temporary bridge in Layer H, not a
COMMON-free Layer-M solution.

### 2. The signed marker representation is incompatible with `E013Max_SUPER`

The state uses `+k` for holes and `-k` for fringes, then proposes an
`E013Max_SUPER` synchronization (`chimera-integration-design.md:261-266`).
That routine performs a numeric `MAX`
(`source/src_quadLS/QuadSc_mpi.f90:3830-3836`). A shared node classified as
fringe (`-k`) on one partition and free (`0`) on another becomes free.

Use separate arrays, for example:

- `marker_kind`: `0=free`, `1=fringe`, `2=hole`;
- `particle_id`: separate identity.

Synchronize kind by precedence, then synchronize identity among ranks with
the selected kind. Also specify the exact paper definition of fringe nodes:
outside-particle nodes belonging to cells cut by the particle boundary. Test
a cut cell spanning two partitions.

## Smaller changes to make before their respective phases

- The plan says the master has no runtime state, but H5 initializes without a
  master guard (`chimera-integration-design.md:298-300,387`). Prefer making
  every facade operation internally rank-safe instead of relying on scattered
  caller guards.
- H6 finalization is application-local while H5 initialization is in shared
  solver initialization. Either prohibit Chimera in all applications except
  `q2p1_chimera`, or pair initialization/finalization in a common lifecycle.
- The distributed correction CG needs an explicit shared-DOF scalar-product
  and preconditioner contract. Achieving bit identity for `D=0` also requires
  a direct fast path executing the existing diagonal loop; a CG solve will not
  generally be bit-identical.
- Write the complete Robin boundary residual, including the sign of the
  Picard-linearized alpha matrix contribution. The strong-form equation is
  now correct, but "enters the matrix" remains ambiguous.

## Verdict

With the assembly and marker decisions corrected, proceed with the full
implementation. As it stands, Phases 0-1 are safe and useful; Phase 2 should
remain gated.
