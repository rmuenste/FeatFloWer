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
  kernel sphere reaches a box (`EL_IMAGE_SHIFTS`), and gives the two periodic distance tests
  of the halo. No MPI, no state besides the lattice.
- `source/src_el/el_image_transfer.f90`, module `EL_IMAGE_TRANSFER`: the three record-level
  operations of the table. With the lattice inactive each is one call of the single-kernel
  routine with the record as it is.
- `source/src_el/el_halo.f90`: the neighbour graph and the send list use
  `EL_HALO_BOXES_COUPLED` / `EL_HALO_PARTICLE_REACHES_BOX`, which switch to the periodic
  tests. `EL_HALO_LOCAL_BOX` hands the rank box to the transfer pass.
- `source/src_el/el_transfer.f90`: `EL_PARTICLE_MESH_PASS` calls the record-level
  operations; nothing else changed in the pass.
- `source/src_el/el_config.f90`, `applications/q2p1_el_pipeflow/app_init.f90`: the key,
  its validation and echo; the application sets the lattice from `dPeriodicity`.

## Conventions

- Only the period enters, never an origin: image = position + n L. The treatment does not
  depend on where the particle engine keeps a wrapped particle. (pe returns positions in
  [0, L) for the box [0, L]: read off the particle logs of P15_probe_cross, min 5.5e-9, max
  0.0199999855 for L = 0.02; a position exactly on a plane, or outside the cell, is handled
  the same way.)
- A record travels once per rank and carries the owner's position, unshifted. The receiving
  rank enumerates the images of every record, owned or received, against its own bounding
  box. One record per rank and particle (`EL_ASSERT_UNIQUE_RECORDS`) stays true.
- A rank is never its own halo neighbour. A rank that touches both planes of a periodic axis
  finds the second image of its own particle by the same local enumeration (2 images at a
  face, 4 at an edge, 8 at a corner of the lattice cell; up to 3 per axis when the kernel is
  wider than half a period).
- Neighbour graph: ranks whose boxes are within the largest kernel radius of each other,
  the gap measured to the nearest image of the other box. A kernel wider than a subdomain
  only makes the graph larger (the existing wide-kernel warning is unchanged).
- The periodicity of the particle engine (`periodicX_/Y_/Z_` of its json) is a separate
  input and must describe the same axes; it is not cross-checked.

## The key

`SimPar@ELPeriodicKernel = off` (default): the code of before, no image list, no additional
arithmetic or collective; every existing log is reproduced bit for bit.

`SimPar@ELPeriodicKernel = wrap`: images on every periodic axis of the fluid. The echo
line `EL periodic kernel = wrap` and one `EL_PERIODIC_KERNEL` line per halo-graph build
(periods, kernel radius, largest neighbour count) appear only then.

Refused:
- `wrap` without a periodic fluid axis (start-up, every rank). A silent no-op would hide a
  deck whose periodicity keys are missing.
- a kernel radius (2 x ELKernelWidthFactor x largest particle radius) that reaches a period
  (first halo-graph build, abort 911): the image bound of 27 per particle and rank needs it.
- with a particle viscous force (models Bfull, A): a rank whose boundary-face mask contains
  faces on a periodic plane (abort 912). The mask counts the elements adjacent to a face
  through the parallel sum of the fluid, which identifies a periodic pair of faces only
  across ranks; a rank holding both planes of an axis would keep a spurious boundary term.

## Limits and open points

- Whether the fluid solver itself supports a rank that holds both periodic planes of an axis
  is outside this change; the coupling does (unit test: one rank holding the whole cell).
- Wall-distance based closures (`ELLiftModel = saffman_mei_wall`) still take the faces of
  the domain box as walls on a periodic axis; unchanged, independent of this key.
- Cylinder domains: the boundary-mask check is made for box domains only.
- The image enumeration assumes rank bounding boxes inside one lattice cell and a
  rectangular lattice. A sheared lattice (Lees-Edwards) would add an offset to the shift
  and a velocity offset to the image record; `EL_IMAGE_SHIFTS` and `EL_RECORD_IMAGES` are
  the two places.

## Tests

`source/src_el/tests/test_el_periodic.f90` (ctest `el-periodic-serial`): image enumeration,
halo predicates, bitwise default path, translation invariance at a face, an edge and a
corner (wrapped and unwrapped positions, a particle exactly on a plane, next to a wall),
normalisation and conservation, continuity of alpha_p across a plane, the viscous average,
and a kernel wider than half a period.
