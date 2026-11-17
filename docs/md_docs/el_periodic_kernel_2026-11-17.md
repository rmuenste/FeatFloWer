# Periodic images of the Euler-Lagrange particle-mesh kernel (2026-11-17)

Design note for `SimPar@ELPeriodicKernel = off | wrap` (branch `feature/el-periodic-kernel`).
Results of the validation runs are not recorded here.

## The defect

With periodic fluid boundaries (`SimPar@Periodic`, `PeriodicAxis` / `PeriodicLength`) the
transfer kernel ended at the periodic planes as it does at a wall. A particle within one
kernel radius of a plane deposited and sampled on its own side only (the normalisation
squeezed its whole volume and reaction force into the truncated support), and when the
particle engine wrapped it across the plane the footprint jumped by one period within a
step: a pair of opposite d(eps)/dt spikes in the two cell layers at the seam, which the
volume-averaged continuity equation turns into a jet normal to the plane.

## The kernel on a lattice

On a periodic axis the fluid domain is one cell of a lattice with period L. A particle at x
stands for the row of images x + n L, and its kernel is the sum over the images,

    W_per(y - x) = sum_n W(y - x - n L),      n over the lattice shifts of the periodic axes.

Every coupling operation is linear in the kernel, so each is the single-kernel operation
summed over the images that reach the elements of a rank:

| operation | single-kernel routine | record level (sum over images) |
|---|---|---|
| sample: int W dV, int W u dV, grad u, p, grad p, eps | `EL_INTEGRATE_PARTICLE` | `EL_SAMPLE_RECORD` |
| viscous average of models Bfull / A (A, B, int W grad u, N) | `EL_INTEGRATE_PARTICLE_VISCOUS` | `EL_SAMPLE_RECORD_VISCOUS` |
| deposit: alpha_p, grad alpha_p, reaction force, drag coefficient | `EL_DEPOSIT_PARTICLE` | `EL_DEPOSIT_RECORD` |

The normalisation is slot 1 of the sample, int W_per dV summed over images and ranks at the
owner (`EL_REDUCE_TO_OWNERS`), and sent back for the deposit (`EL_BROADCAST_FROM_OWNERS`).
The weights of a particle therefore sum to one over the real plus image support exactly as
they do in the interior; deposited volume and reaction force are conserved by the same
argument as before. Everything downstream of the deposit (void fraction, d(eps)/dt, the
transport porosity of the conservative continuity form, the fluid source, the audits that
integrate these fields) reads the deposited fields and needs no change.

Translation invariance: moving the particle and the fields by a whole number of cells along
a periodic axis moves every result by the same number of cells, to rounding (the image sums
are added in a different order). In particular nothing changes when the particle engine
replaces x by x - L.

## Where the logic lives

- `source/src_el/el_periodic.f90`, module `EL_PERIODIC`: the only place that knows about
  images. Holds the lattice (flag and period per axis), enumerates the lattice shifts whose
  kernel sphere reaches a box (`EL_IMAGE_SHIFTS`), gives the two periodic distance tests of
  the halo and the predicates behind the refusals below. No MPI, no stop (a refusal is
  returned to the caller), no state besides the lattice.
- `source/src_el/el_image_transfer.f90`, module `EL_IMAGE_TRANSFER`: the three record-level
  operations of the table. With the lattice inactive each is one call of the single-kernel
  routine with the record as it is.
- `source/src_el/el_halo.f90`: the neighbour graph and the send list use
  `EL_HALO_BOXES_COUPLED` / `EL_HALO_PARTICLE_REACHES_BOX`, which switch to the periodic
  tests. `EL_HALO_LOCAL_BOX` hands the rank box to the transfer pass, `EL_HALO_DOMAIN_BOX`
  the bounding box of the whole mesh.
- `source/src_el/el_transfer.f90`: `EL_PARTICLE_MESH_PASS` calls the record-level
  operations and verifies the boundary-face mask; nothing else changed in the pass.
- `source/src_el/el_geometry.f90`: the wall search of the wall-corrected lift leaves out
  the faces of a periodic axis.
- `source/src_el/el_config.f90`, `applications/q2p1_el_pipeflow/app_init.f90`: the key,
  its validation and echo; the application sets the lattice from `dPeriodicity` and checks
  it against the fluid domain. `EL_KERNEL_SUPPORT_RADIUS` (EL_CONFIG) is the one definition
  of the kernel radius of a particle, used by the quadrature, the halo and the images.

## Conventions

- Only the period enters, never an origin: image = position + n L. The treatment does not
  depend on where the particle engine keeps a wrapped particle. (pe returns positions in
  [0, L) for the box [0, L]: read off the particle logs of P15_probe_cross, min 5.5e-9, max
  0.0199999855 for L = 0.02; a position exactly on a plane, or outside the cell, is handled
  the same way.)
- A record travels once per rank and carries the owner's position, unshifted. The receiving
  rank enumerates the images of every record, owned or received, against its own bounding
  box. One record per rank and particle (`EL_ASSERT_UNIQUE_RECORDS`) stays true.
- Images of a record on a rank: per periodic axis the image nearest to the centre of the
  rank box and its two lattice neighbours are candidates, the sphere-box test decides
  (2 images at a face, 4 at an edge, 8 at a corner of the lattice cell). With a rank box
  shorter than a period and a kernel radius below a period no other image can reach the
  box, so 27 is a bound by construction. The nearest image is found in real arithmetic: a
  position far outside the cell needs no integer conversion; a position that is not finite,
  or more than 2^52 periods away, has no image and ends in the existing error 901 (particle
  without fluid support).
- A rank is never its own halo neighbour.
- Neighbour graph: ranks whose boxes are within the largest kernel radius of each other,
  the gap measured to the nearest image of the other box. Every rank declares its own
  neighbours, so the test must give two ranks the same verdict on their pair: it is
  symmetric bit for bit (per axis `d1 = a_lo - b_hi`, `d2 = b_lo - a_hi`, gap to the image
  shifted by s = `MAX(0, d1 - s, d2 + s)`, both signs of s; exchanging the boxes exchanges
  d1 and d2 and maps s onto -s) and is evaluated with the box of the lower rank first on
  both sides. A kernel wider than a subdomain only makes the graph larger (the existing
  wide-kernel warning is unchanged).
- Send list: only the owner decides which records exist; the receiver learns the count from
  it and enumerates the images with the same routine, radius and gathered box, so the two
  sides cannot disagree. A received record of which no image touches an element contributes
  nothing.
- Walls: the faces of the domain box normal to a periodic axis are not walls for the
  wall-corrected lift (`ELLiftModel = saffman_mei_wall`); with no wall left the closure
  gives the unbounded `saffman_mei` lift. Without the key all six faces are walls as before.
- The particle log (`el_particles_rank*.log`) holds the positions as the particle engine
  keeps them, wrapped into the cell: a trajectory jumps by one period at a crossing, with
  either setting of the key. Post-processing has to unwrap.

## The key

`SimPar@ELPeriodicKernel = off` (default): the code of before, no image list, no additional
arithmetic or collective; every existing log is reproduced bit for bit (gates: the ten
certified EL gates and a periodic case, with a binary built on the host of the certified
build).

`SimPar@ELPeriodicKernel = wrap`: images on every periodic axis of the fluid. The echo
line `EL periodic kernel = wrap` and one `EL_PERIODIC_KERNEL` line per halo-graph build
(periods of the fluid lattice with 0 = not periodic, kernel radius, largest neighbour count)
appear only then.

Refused under `wrap`:

| condition | where | how |
|---|---|---|
| no periodic fluid axis (a silent no-op would hide a deck whose periodicity keys are missing) | deck parser, every rank | `STOP 1` |
| a period that is not positive and finite | deck parser, every rank | `STOP 1` |
| `ELWallStressDiag` with z periodic (the diagnostic integrates the fluid shear on the z planes as walls) | deck parser, every rank | `STOP 1` |
| period of an axis differs from the extent of the fluid domain by more than 1e-8 periods (`SimPar@Periodic = Yes` sets the period 1, `PeriodicLength` defaults to 1) | after the mesh is read, global extents, every rank | abort 913 |
| kernel radius (2 x ELKernelWidthFactor x largest particle radius) reaches a period | every halo-graph build, all workers | abort 911 |
| a subdomain spans a whole period (extent >= L (1 - 1e-8)): the rank would hold both planes of the axis | every halo-graph build, gathered boxes, all workers | abort 914 |
| particle viscous force (models Bfull, A): a face on a periodic plane is in the boundary-face mask, or the mask is longer than the dof coordinates | every transfer pass, on the mask as it is then, rank-local | abort 912 |

The periodic planes of check 912 are those of the lattice on the bounding box of the whole
mesh (gathered rank boxes), independent of the wall geometry of the closures (box or
cylinder). The mask counts the elements adjacent to a face through the parallel sum of the
fluid, which pairs the two faces of a periodic plane pair across ranks; check 914 keeps
both faces off one rank.

Without the key nothing of this is checked and nothing changes: in particular the
wall-stress diagnostic with a periodic z axis and the wall lift next to a periodic plane
behave as they did.

## Limits and open points

- Particle engine lattice not cross-checked. The periodicity of pe (`periodicX_/Y_/Z_` of
  its json; pe takes the periods from the extents of the fluid domain it is handed) is a
  separate input. No function of the pe interface returns it to Fortran
  (`libs/pe/pe/interface/cfd_queries.h` has no periodic query), so the coupling cannot
  compare: a json without the periodic flags leaves the particles bounded by walls while
  the coupling wraps. The run log shows the fluid lattice (`EL_PERIODIC_KERNEL` line) and
  pe prints its own set-up. Closing the gap needs a getter in pe, e.g.
  `extern "C" void getElPeriodicLattice(int periodic[3], double period[3])` filled by the
  EL set-up (`setupELTerminalVelocity`), and a comparison at start-up under `wrap`.
- A subdomain that spans a period is refused, not supported (914). The particle engine
  needs three processes per periodic axis, but its process grid and the partition of the
  fluid mesh are separate inputs, so that rule does not protect the coupling. The image
  enumeration itself handles such a rank (unit test: one rank holding the whole cell);
  what is not established is the fluid side, the pairing of the dofs of the two planes
  inside one rank, on which the boundary-face mask and the assembled deposit depend.
- Restart. `EL_WRITE_RESTART` / `EL_READ_RESTART` (EL_FIELDS) have no caller in the
  application or the solver (unit tests only): the porosity history of the coupling is not
  restored after a restart but rebuilt by the first transfer passes, and a restart
  re-inserts the particles. A restart is therefore not affected by the key, and a run may
  be continued with the other setting. If a coupled restart is wired later, the file needs
  the mode and the lattice, and a mismatch must reseed the history (`history_valid`,
  `rate_valid`).
- The images of a record are enumerated once per operation (sample, viscous average,
  deposit) and once per neighbour in the send list. Keeping a list per pass would need
  state in the halo; not done.
- Cylinder domains: the wall search keeps the cylinder wall; a lattice transverse to the
  axis of a cylinder is not meaningful and not checked.
- The image enumeration assumes a rectangular lattice. A sheared lattice (Lees-Edwards)
  would add an offset to the shift and a velocity offset to the image record;
  `EL_IMAGE_SHIFTS` and `EL_RECORD_IMAGES` are the two places.
- Binaries are compiled with `-march=native`: rounding-level diagnostics (conservation
  residuals) of a run are bitwise reproducible only with a binary built on the same host
  type. Default-path gates against certified logs need a build on the host of the
  certified build.
- Not part of this change: a coupled EL restart contract; unit tests with the other
  kernels (the tests use `deen_poly`); the stagnation of the pressure multigrid under
  lateral periodicity; the default of the key.

## Tests

`source/src_el/tests/test_el_periodic.f90` (ctest `el-periodic-serial`, one rank holding the
whole cell): image enumeration, halo predicates, default path (bit patterns of sample,
deposit and viscous record against the single-kernel routines), translation invariance at
a face, an edge and a corner (wrapped and unwrapped positions, a particle exactly on a
plane, next to a wall), normalisation and conservation, continuity of alpha_p across a
plane, the viscous average, a kernel wider than half a period, bitwise symmetry of the
neighbour test on tie distances (named cases and a scan that the previous expression
fails), the wall lift next to a periodic plane, the refusal predicates, positions far
outside the cell and positions that are not finite.

`source/src_el/tests/test_el_periodic_mpi.f90` (ctest `el-periodic-mpi-4`, `-8`): the cell
partitioned into 2 x 2 and 4 x 2 blocks; the halo graph build, record exchange, owner
reduction, normalisation broadcast and the three record operations of the application
against one rank holding the cell (owner samples, viscous records, deposited fields element
by element and dof by dof). Contains rank pairs that are direct and wrap neighbours at
once, a corner particle, a particle on a plane, an unwrapped position, particle-free
workers, and on eight ranks the tie distance of the neighbour test. Not covered by unit
tests: `EL_PARTICLE_MESH_PASS` itself and the dof pairing of the fluid across a periodic
plane (application runs).
