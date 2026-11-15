# s06_startup_inertia - suspension starting from rest (inertia weighting of the momentum models)

Purpose: the one quantity of a simple case that depends on how the fluid inertia is weighted.
A homogeneous suspension in a closed box starts from rest under gravity; with zero net flux
(eps u + phi v = 0) the slip s = u - v relaxes to the same terminal value in every formulation,
but with a different time constant (case proposed by the LIGGGHTS/CFDEM peer session, 2026-11-15;
derivation re-done here, `tools/startup_reference.py`):

| Formulation | Equation of the slip | tau / tau_A at phi 0.4, rho_p/rho_f 1.1 |
|---|---|---|
| A, Bfull | (rho_f phi + rho_p eps) ds/dt = drho g - n F(s)/(eps phi) | 1 |
| `none` + pressure force + feedback (unweighted inertia) | (rho_f phi/eps + rho_p eps) ds/dt = same | 1.25 |
| B | eps phi rho_p ds/dt = phi drho g - n F(s)/eps | 0.62 |

B is not expected to agree with A and Bfull here (it is equivalent for steady uniform flow only).
The Ergun column (s05) cannot see any of this: at rho_p/rho_f = 200 the effect is 0.2 %.

## Set-up

d = 1 mm, rho_p 1100, fluid rho 1000, mu 1e-3 (water), g 9.81, `difelice_std`, kernel `deen_poly`
2.5 d, nominal particle radius for porosity and drag. Closed box, coarse mesh 4 x 4 x 12 from
the case-04 tool `make_box_mesh.py`, level 2 (cells 2.6 d), six z slabs, 7 ranks, dt 0.25 ms,
200 steps (50 ms), pe 10 sub-steps. Decks `q2p1_param_{none_pf_fb,B,Bfull,A}.dat`, `example.json`.
Terminal Reynolds number about 4, so the drag is not linear in the slip: runs and reference are
compared through the same exponential fit v = A (1 - exp(-t/tau)) on t <= 20 ms
(`tools/startup_post.py`, `tools/startup_reference.py`).

Two boxes:
- **S** (first): 20 x 20 x 60 mm, `tools/make_suspension.py --size 0.02 0.02 0.06 --d 1e-3
  --phi 0.4` (jittered FCC lattice centred in the box, N 16 456). The lattice does not reach the
  side walls: solid fraction 0.374 over the cross-section of the middle third, 0.40 in the core.
- **T** (wall to wall): `tools/make_suspension_tight.py` (24 lattice planes, first plane 0.55 d
  from the wall, L = 21.069355 mm, N 19 584, solid fraction 0.400 in the middle third).

Staging of a run directory: mesh `make_box_mesh.py sus_box --size L L 0.06 --cells 4 4 12`,
`featflower-partition 6 -3 6 SUS6 sus_box/file.prj`, deck to `_data/q2p1_param.dat`, `_data/MG.dat`,
`example.json`, `particles.xyz`; submit `tools/el_slurm_submit_np.sh <dir> <name> 00:30:00 short 7 1 25G`.
Run directories: /data/warehouse17/rmuenste/runs/startup06/.

## Results

Box T on the final commit 08f7d420 (binary build-el-e, clean tree; `runs_versions.txt`), all
6 528 particles of the middle third; reference at phi 0.4:

| Run | job | wall | tau fit (ms) | reference | ratio | amplitude (mm/s) | reference |
|---|---|---|---|---|---|---|---|
| F_none_pf_fb | 160551 | 3 min 31 s | 8.844 | 9.009 | 0.982 | -4.110 | -4.110 |
| F_B | 160552 | 3 min 08 s | 4.493 | 4.351 | 1.033 | -4.227 | -4.051 |
| F_Bfull | 160553 | 4 min 55 s | 7.057 | 7.106 | 0.993 | -4.120 | -4.083 |
| F_A | 160554 | 5 min 11 s | 7.221 | 7.106 | 1.016 | -4.165 | -4.083 |

- Inertia weighting: tau(none)/tau(Bfull) = 1.253 (predicted 1.268), tau(B)/tau(Bfull) = 0.637
  (predicted 0.612), tau(A)/tau(Bfull) = 1.023 (predicted 1).
- Every formulation is within 3.3 % of its own prediction in tau; B is 4.3 % high in the
  amplitude (not explained).
- Particles more than 4 mm from the side walls only (2 314 of them) relax faster and settle
  slower (tau 8.16 / 4.00 / 6.51 / 6.75 ms, amplitude -3.82 .. -3.92 mm/s): a lateral flow
  structure remains in a box of 21 d; the ratios between the formulations are the same within 2 %.
- The same four runs on the pre-review working tree (T0..T3, jobs 160485-160488) gave 8.844 /
  4.483 / 7.057 / 7.221 ms: the porosity priming of the first step moves model B by 0.2 %.
- Cost: the viscous sampling of Bfull and A adds about 55 % to the wall time here (19 584
  particles on 1 536 elements); the per-step mass rebuild is not visible in B against `none`.

Box S (pre-review working-tree binaries, three different builds; kept as the record of why box T
exists, not as certified numbers), particles of the middle third more than 4 mm from the walls
(2 314):

| Run | job | tau fit (ms) | ratio to S2 | predicted (phi 0.374) |
|---|---|---|---|---|
| S0 `none` + pf + fb | 160457 | 8.643 | 1.208 | 1.225 |
| S1 B | 160458 | 4.483 | 0.627 | 0.634 |
| S2 Bfull | 160461 | 7.155 | 1 | 1 |
| S3 A | 160463 | 7.441 | 1.040 | 1 |

- The loose wall region of box S separates A from Bfull by 4 % (2.3 % in box T): the two
  discretisations of -eps grad p (continuous transport porosity in the operator against the
  kernel-deposited particle force) differ where the porosity varies across one kernel width.
- Absolute tau in box S is about 11 % below the reference at the cross-section mean phi 0.374 and
  fits phi of about 0.39: the reference porosity of a non-uniform cross-section is not defined
  sharply, which is why box T exists.

Other observations: Newton pair (particle momentum against expected impulse, viscous force
included) closes to 1e-19 in every run; the conservative residual res1_L1 falls from 1e-9 at step
1 to 1e-13 at step 200 with one corrector in all four; the particle viscous force is 1e-6 of the
drag here (uniform flow, and see the viscosity-convention note in the design note, section 6).

Not done: a periodic variant (pe needs three processes per periodic direction and the partitioner
call for a 3 x 3 column split was not available); dt and level convergence; the CFDEM side (the
peer has the case definition; nothing has run there yet).

## Certification of the final commit 08f7d420 (2026-11-15)

- Default path: the ten certified EL configurations (tier-2 twin, semi-implicit twin, VB_L2
  defaults, T9a, s01 two-way, P4 / P5 controls, conservative T10a4) bitwise against
  runs/coupling_p4p5/ergun_gates2 with the strict comparison, jobs 160541-160550
  (runs/coupling_p4p5/mm_final; versions.txt HEAD 08f7d420, dirty 0). Earlier rounds on the way:
  key only (160443-160452), combined code before the reviews (160466-160475), review fixes
  697d9f49 (160491-160500): all bitwise.
- Negative tests (jobs 160555-160558), each refused at start-up with exit 1: B with
  ELPressureForce = Yes and ELDragCoupling = explicit (both conflicts listed); Bfull with
  ELKernel = gaussian; A with MatrixRenewal C2; Bfull with ELDragModel = difelice.
- DNS twins (build-dns5, certified DNS recipe, HEAD 08f7d420; rundirs
  q2p1_dns_rundir_twin_mmf_{e4l3,dkt} in the DNS clone): e4_l3 job 160562 and DKT job 160563
  bitwise against the certified heads (particle_force.log / bulk_flow.log md5 ca98e7af225b /
  89abce88b675 and 82411f64f661 / cfb8b42d0969, all log streams identical, no NaN).
- Unit tests: test_el_momentum_model (37 checks), test_el_kernel_forces, test_el_va_flux,
  test_el_transfer pass.
- Found by this round and fixed before it passed: the boundary-face mask made a workers-only
  MPI call on the master (bdc70e3e) and was released by the first field allocation (08f7d420);
  both stopped every Bfull / A run at step 1 and neither touched the default path.

## Gates proposed for the datasheet

| Quantity | Expected | Tolerance |
|---|---|---|
| tau(none + pf + fb) / tau(Bfull), box T | 1.268 | 5 % |
| tau(B) / tau(Bfull), box T | 0.612 | 5 % |
| tau(A) / tau(Bfull), box T | 1 | 3 % |
| tau of each formulation / own reference, box T | 1 | 5 % |
