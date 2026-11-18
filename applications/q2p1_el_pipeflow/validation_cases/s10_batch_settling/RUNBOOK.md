# s10_batch_settling - batch settling with two concentration fronts (cross-code case 10)

Record of 2026-11-16 ... 2026-11-18. Spec and reference: the peer's case spec 10
(`liggghts-dem/LIGGGHTS-PUBLIC/md_docs/case_specs/10_batch_settling/SPEC.md`, derived from MFiX
PIC04; `tools/batch_reference.py` is a copy of its reference script). Marks: **[V]** verified by a
run or output read here, **[P]** reported by the peer session (CFDEM, MFiX, DPMFoam, waLBerla; not
recomputed here), **[C]** from the consultation of GPT-6 Astra (high) on 2026-11-18, **[I]** inferred.

## Status in one paragraph

A wall-bounded bed on a laterally resolved fluid mesh settles 16 % faster than the one-dimensional
reference because the bed is depleted next to the walls and a column-scale convection cell forms;
three other Euler-Lagrange codes show the same cell, with excesses between 0.93 and 1.45 depending
on code and particle-to-mesh mapping. A laterally periodic bed (no wall depletion; needs the
wrapping kernel `SimPar@ELPeriodicKernel = wrap` of this branch) settles 4-5 % above the reference,
in agreement with MFiX's periodic run (1.053). **The homogeneous-reference gates B1 and B2 are
exceeded as written; whether the one-dimensional reference is the right yardstick for a laterally
resolved suspension is not established (owner's note, 2026-11-18). No pass is claimed and no
solver defect is diagnosed; physical validation is inconclusive** until a trustworthy
particle-resolved result exists. The mesh-refinement run is blocked by the pressure solver on
periodic meshes (section 7).

## 1. Case

Closed column 20 x 20 x 200 mm, g 9.81 m/s2 in -z. Fluid rho 1000 kg/m3, mu 0.05 Pa s. Spheres
d 1 mm, rho_p 2700 kg/m3, N = 18 335 (phi0 0.15 in z < 160 mm), released from rest. Restitution
0.9, friction 0.1. Single-particle Re 0.33, relaxation time 3 ms.

Reference (homogeneous settling, zero net volume flux, Di Felice drag with DallaValle C_d, the
drag law of these runs): hindered settling speed and upper-front speed `eps0 u_r0` = 8.172 mm/s;
lower front `f0 / (phi_sed - phi0)` with the measured sediment fraction. Proposed gates of the
spec: B1 core particle velocity +-3 %, B2 upper front +-5 %, B3 lower front +-10 %.

What the reference is: the value the drag closure gives for a perfectly uniform suspension. A
one-dimensional fluid mesh reproduces it by construction (CFDEM 1 x 1 x 100: core 1.003 / 0.971 /
1.003 at phi0 0.10 / 0.15 / 0.20 [P]). This application cannot run that mesh: the particle-to-mesh
transfer finds no fluid support when elements are four times wider than the kernel
("particle has no positive fluid support", 1 x 1 x 40 attempt, 2026-11-16 [V]).

## 2. Set-up

Model settings common to all runs: `ELMomentumModel = A` with `ELMMPressureOperator = weighted`,
conservative continuity form, drag `difelice_std`, `ELFluidViscosity = dynamic` (default),
kernel `deen_poly` with `ELKernelWidthFactor = 2.5` (2.5 mm), dt 0.5 ms, particle log every
200 steps (0.1 s). Decks and json in this directory.

| Variant | Deck / json | Mesh (level 2) | Sides | Partition, ranks | Bed |
|---|---|---|---|---|---|
| wall-bounded | `q2p1_param_wall.dat`, `example_wall.json` | 8 x 8 x 80, 2.5 mm | free slip (Symmetry), particle box walls | 1 x 1 x 20, 21 | centres in [r, L - r], sha256 5b04d0ef277e3a28 |
| periodic | `q2p1_param_periodic.dat`, `example_periodic.json` | 6 x 6 x 60, 3.33 mm | `PeriodicAxis = xy`, json `periodicX_/Y_`, z walls | 3 x 3 x 3, 28 | centres in [0, L), seed 1 sha256 4f4ae6ea1c79014b, seed 2 af72dbe59dd15297 |

- Mesh: the peer's `make_box_mesh.py bat_box --size 0.02 0.02 0.2 --cells 4 4 40` (wall-bounded)
  or `--cells 3 3 30` (periodic; side `.par` headers `Periodic`), not committed.
- Partition: wall-bounded `featflower-partition 20 -3 20 BAT20 bat_box/file.prj`; periodic
  `PyPartitioner.py 27 axis_uniform x3-y3-z3 BAT27P bat_box/file.prj` from the pe_partpy branch
  `fix/axis-uniform-boundary-precision` (pe needs >= 3 processes on each periodic axis, and the
  fluid partition must own the same slabs).
- Beds: the peer's `reference/make_bed.py 0.15 OUT 1` (wall-bounded) and
  `PERIODIC=1 ... make_bed.py 0.15 OUT <seed>` (periodic); the `.xyz` is used unchanged as
  `particles.xyz`. Not committed (0.9 MB each); hashes above.
- Run directories: `/data/warehouse17/rmuenste/runs/batch10` (`A15`, `K5_15`, `K8_15`,
  `P15_*`, `PK_*`, `PL15_s1`, `PL15_s2`, `PR15_m12_probe`); provenance in `runs_versions.txt`.

Tools (`tools/`): `batch_eval_ff.py RUN 0.15 T0 T1` (the peer's definitions: phi(z) in 2 mm bins,
upper front at phi0/2, lower front at (phi0 + phi_sed)/2, core = particles more than 10 mm from
both fronts, least-squares front speeds); `batch_diag_ff.py` (count, max |v|, phi(z) noise, wall
ring / interior); `batch_lateral_map.py` (column maps of solid fraction and v_z over 40-100 mm);
`batch_seam_stats.py`, `batch_vtk_lateral.py` (periodic planes); `pk_probe_report.py`,
`pk_particle_diff.py` (crossing probe, run-to-run differences); `batch_movie_ff.py` (particle
movie of the wall-bounded run).

## 3. Wall-bounded bed (run A15, job 160893, integration/el-pweight-visc d6b65cae) [V]

Ended on its 5 h limit at 4.38 s (8756 of 12 000 steps), no NaN, 18 335 particles in every frame,
max |v| 28-61 mm/s, deposited volume to 1e-14, boundary flux 1e-21, no burst from the sediment.

| Quantity | Window | / reference |
|---|---|---|
| core particle velocity | 1.0-4.3 s | 1.156 |
| upper front | 1.0-4.3 s | 1.260 |
| core particle velocity | 3.0-4.3 s | 1.141 |
| upper front | 3.0-4.3 s | 1.232 |
| lower front (phi_sed 0.565, sediment 9-14 mm) | 3.0-4.3 s | 1.098 |

Cause: lateral convection driven by the bed. Solid fraction 0.139 in the outer ring of 2.5 mm
columns against 0.161 inside. In 40-100 mm height the interior has mean v_z -15 to -16 mm/s
(1-9 % of its particles rising), the wall ring has mean -0.4 to -1.2 mm/s with std 8-12 mm/s and
40-44 % rising: a cancellation, not a quiet layer. Height-averaged column means at 3 s from
-29 mm/s (sinking core, off-centre) to +14 mm/s (corners, one side). Established by 0.4 s, steady
to 4.3 s. The upper front (10-90 % thickness 2.0-2.5 mm throughout) runs ahead of the core mean
because the plume reaches it first [I].

Kernel width, same mesh and bed, short runs: 5 d (K5_15, job 160939, to 0.74 s) core 1.061,
wall ring -5.9 / interior -10.5 mm/s; 8 d (K8_15, 160940, to 0.47 s) core 1.042, -7.2 / -9.4.
Wider kernels smear the lateral density difference; this is not a convergence statement.

Other codes on this bed, 8 x 8 x 80 [P], each against the reference of its own drag law:

| Run | core 1.0-4.3 s | upper front | column means at 3 s (mm/s) | ring / interior phi |
|---|---|---|---|---|
| FeatFloWer (this record) | 1.156 | 1.260 | -29 .. +14 | 0.139 / 0.161 |
| MFiX, cell-centre mapping | 1.432 | 1.404 | -30.0 .. +11.7 | 0.145 / 0.153 |
| MFiX, SQUARE_DPVM 1 mm | 1.316 | 1.323 | -25.6 .. +25.2 | 0.143 / 0.158 |
| MFiX, GARG_2012 | 1.208 | 1.214 | -13.5 .. -6.6 | 0.137 / 0.160 |
| DPMFoam (Wen-Yu, ref. 9.762 mm/s) | 0.925 | n/a | -33.8 .. +12.4 | 0.150 / 0.154 |
| CFDEM, 3 mm smoothing (to 3.2 s) | 1.45 (1.25 late) | 1.92 | -99 .. +55 | 0.177 / 0.155 |

Reading: a sinking core with a slow or rising wall ring in every code; the size of the excess
depends on the particle-to-mesh mapping in every code (1.04-1.45). The one-dimensional reference
does not apply to this configuration. The intrinsic-convection literature (Geigenmueller and
Mazur; Bruneau et al. 1996) describes the mechanism for no-slip walls and gives no number for
free-slip walls at phi 0.15 [C]. Movie: `runs/batch10/movie/batch10_A15_particles.mp4`.

## 4. Periodic bed

### 4.1 Defect found and fixed: the kernel did not wrap [V]

Smoke with the integration binary (P15_smoke 161466, P15_shift 161470, P15_probe_cross 161469):
each particle crossing of a periodic plane moved its whole footprint by 20 mm in one step: d(eps)/dt
+-23 1/s across the seam (0.5 before), pressure initial defect 4e-7 -> 2.7e-3, a fluid jet normal
to the seam, repeated re-crossings; seam-plane fluid rms 9-32 mm/s against 4-7 elsewhere;
continuity residual 1e-7. Fixed on this branch: `SimPar@ELPeriodicKernel = off | wrap` (default
`off`, bitwise), design note `docs/md_docs/el_periodic_kernel_2026-11-17.md`. With `wrap`
(PK_probe_cross_wrap2 161601, PK_smoke_wrap2 161603, PK_shift_wrap2 161607): no signature of a
crossing in any quantity, seam-plane fluid rms 1.1-1.5 mm/s like elsewhere, continuity residual
3e-11 ... 9e-11, largest fluid velocity 0.022 m/s (0.74 before); a bed and its copy shifted by L/2
agree at 0.5 s to rms 1.2e-7 m in position.

### 4.2 Long runs with `wrap` (PL15_s1 job 161617, PL15_s2 job 161618; e01e143b) [V]

Both COMPLETED, 12 000 steps (6 s), no NaN, 18 335 particles throughout, volume conservation
3e-15. Ratio to the reference, seed 1 / seed 2:

| Window | core particle velocity | upper front |
|---|---|---|
| 0.5-1.0 s | 1.035 / 1.037 | 1.047 / 1.036 |
| 1.0-2.0 s | 1.039 / 1.041 | 1.067 / 1.091 |
| 1.0-4.3 s | 1.047 / 1.054 | 1.085 / 1.084 |
| 3.0-6.0 s | 1.062 / 1.082 | 1.089 / 1.079 |

Lower front, 3.0-6.0 s: 1.011 / 1.061 (phi_sed 0.582). The lower-front fit must start after the
sediment plateau exists; a 1.0-4.3 s fit includes zero positions and is meaningless.

- By displacement between 0.1 s frames the settling speed is flat to 0.1 mm/s from 0.2 s to 2 s
  (seed 1: 8.33-8.36 mm/s in the band 20-40 mm, 8.47-8.55 mm/s in 60-120 mm); the core then drifts
  up by about one point per second.
- Lateral structure, 40-100 mm, 3.33 mm columns: rms of column-mean v_z 1.3-2.3 mm/s up to 4 s,
  3.0-3.7 mm/s at 5 s (wall-bounded: about 10); one wavelength across the box, located where the
  bed is densest; column means between about -12 and -5 mm/s up to 4 s, no rising column until
  then. After about 4.1 s a box-scale cell develops (seed 1): the mean over the fixed band
  40-100 mm climbs from 1.05 to 1.10 of the reference at 6 s (sub-bands 1.06 at 40-60 mm to 1.17
  at 90-100 mm) and particles rise in one corner of the 80 mm slab (movie
  `runs/batch10/movie/batch10_PL15_s1_particles.mp4`). MFiX ends at 4.4 s. Whether this is a
  property of the model equations or of the poorly converged pressure solve of these runs is open
  (section 7, item 1).
- Mean solid fraction between the fronts drifts 0.150 -> 0.153-0.154 by 6 s.
- The pressure multigrid converges poorly under xy-periodicity. In PL15_s1 63 % of the solves hit
  the 20-cycle cap and the mean defect reduction per solve is 6.5e-3 (wall-bounded A15: 2 cycles,
  1e-5). A 200-cycle twin before the fix gave the same ratios (P15_smoke_mg200, 161467: 7.7 cycles,
  reduction 3e-4). Not fixed; see section 7.

Other codes, periodic [P]: MFiX (CYCLIC, GARG_2012, 8 x 8 x 80, same seed-1 bed, reference
8.336 mm/s): core 1.053 (1.0-4.3 s), 1.052 (3.0-4.3 s), upper front 1.034, column means
-10.9 .. -7.6 mm/s. CFDEM periodic: core mean oscillates between +2 and -2 times the reference
(undiagnosed). waLBerla particle-resolved (20 x 20 x 80 mm, 7 237 spheres, 10 cells per diameter,
2 s): 0.74-0.75 of the reference, core still accelerating, whole-bed bounce; its authors do not
count it as an arbiter yet.

### 4.3 Verdict against the written gates

| Gate | Tolerance | Measured (1.0-4.3 s; lower front 3.0-6.0 s) | As written |
|---|---|---|---|
| B1 core | +-3 % | +4.7 % / +5.4 % | exceeded |
| B2 upper front | +-5 % | +8.5 % / +8.4 % | exceeded |
| B3 lower front | +-10 % | +1.1 % / +6.1 % | met |

## 5. Is the one-dimensional reference the right yardstick? (owner's note and consultation)

Owner, 2026-11-18: the comparison may not be against the right reference; the 1-D reference may
not be the correct one for a three-dimensional, laterally resolved case. Consultation of GPT-6
Astra (high) the same day (answer kept with the session records; it recomputed the ratios of 4.2
from the particle logs) [C]:

- The homogeneous state is a solution of the model equations, but that does not make a finite bed
  on a laterally resolved mesh stay on it. With `<q_z> = 0` over a cross-section the
  particle-weighted mean velocity is `<phi q_z>/phi + <phi eps s>/phi`; the reference sets the first
  term to zero and evaluates the second at uniform concentration. Dense regions sinking
  preferentially raise the mean. There is no universal size of the correction.
- Three numbers need not coincide: the homogeneous value (right for checking the closure balance),
  the converged resolved result (right for the model equations as posed), the physical settling
  speed (needs an experiment or a trustworthy particle-resolved run).
- Checked on PL15_s1 / PL15_s2: column concentration and settling speed correlate (0.38-0.57) and
  the column velocity variance grows 3.9 / 5.2 times from 1 s to 5 s, but this weighting explains
  only 0.3-0.6 % of the reference at 1-3 s. The concentration drift and the non-linear averaging of
  the flux over 3.33 mm cubes both have the opposite sign (-2 %). So the early +3.5 % is
  unexplained; the leading candidate is the self-induced fluid disturbance of two-way point-particle
  coupling (each particle samples fluid that carries its own drag reaction). An earlier
  single-particle case of this application showed a two-way bias of about +6 % of that kind [V,
  campaign record].
- The wall-bounded spread 1.04-1.45 shows that prediction is not robust to mapping; widening the
  kernel until the mean approaches the homogeneous value is not validation.
- Recommended record: keep the exceedance, claim no revised pass, mark physical validation
  inconclusive. Recommended next run: mesh refinement at fixed physical kernel (section 7).

## 6. What is and is not established

Established [V]: the runs are numerically sound (conservation, steady, no burst, two beds agree to
about one point); wall depletion drives the wall-bounded excess (removing it cuts the excess from
16 % to 4-5 % and the lateral convection five-fold); the periodic result agrees with MFiX within
one point.
Not established: which of 1.00 and about 1.05 is physically right for a laterally resolved uniform
suspension; the origin of the early +3.5 %; convergence in mesh and kernel; the wall-bounded
number.

## 7. Open

1. Mesh refinement, periodic seed 1, 12 x 12 x 120 elements, same 2.5 mm kernel, bed and dt:
   not obtained with the certified binary. The 20-step probe (`PR15_m12_probe`, job 162042, clean
   otherwise) shows the pressure solve effectively unconverged on that mesh: every solve hits the
   20-cycle cap at 0.994 per cycle, defect reduction 0.5-0.7 per solve, continuity residual after
   the corrector 1e-7 (coarse mesh 1.5e-9 at the same step). The full run (job 162043) was
   cancelled after 31 s.
   Cause (diagnosis 2026-11-18, probes `MG_*`, jobs 162050-162073, scratch branch
   `scratch/mg-periodic-diag` 2e081314 in worktree FeatFloWer-mgdiag): the master's coarse
   pressure matrix `B^T M^-1 B` is assembled on the unpartitioned coarse mesh without identifying
   the velocity dofs of opposite periodic planes; those dofs are free one-sided boundary dofs
   there (like p = 0 on the planes), while the workers' operator is periodic. 1200 of 2160 master
   rows (the coarse elements on a lateral plane) have a non-zero row sum; none on the wall-bounded
   mesh. Not specific to this application (single-phase clone `MG_sp` stagnates too). No deck
   remedy: coarse solver type 2 unchanged, 64 smoothing steps 0.10, 200 cycles 0.015. Prototype
   repair on the master only (periodic dof pairs identified, existing kernels reused): 2.0 cycles,
   reduction 1.8e-5 (fine), 6.2e-5 (6 x 6 x 60); wall-bounded control identical solver lines.
   Proper fix (deck key, default off, gates, review) not written yet.
   Provisional reruns with the prototype binary (unreviewed, not gated; `MGDIAG_MODE=2`):
   `PR15_m12_mgw` (fine mesh, to 4.4 s, job 162080) and `PL15_s1_mgw` (6 x 6 x 60 repeat of
   PL15_s1, 6 s, job 162081); results to be appended.
2. Particle-resolved arbiter (peer, waLBerla; a FeatFloWer resolved run at D/h 4-6 is under
   consideration by the owner).
3. phi0 0.10 and 0.20, models B and Bfull: not run.
4. Pressure multigrid under xy-periodicity; pe-side getter for the periodic lattice cross-check
   (see the design note); default of `ELPeriodicKernel`.
