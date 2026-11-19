# Periodic coarse pressure operator (2026-11-19)

Design note for `Pres@PeriodicCoarseOperator = legacy | gathered | wrapped` and
`Pres@PeriodicCoarseCheck = n` (branch `feature/pressure-periodic-coarse`, worktree
`FeatFloWer-pcoarse`, base `feature/el-periodic-kernel` 16ffcbcc). Gate and probe runs of the
reviewed state (review round 2): `/data/warehouse17/rmuenste/runs/pcoarse/round4` (EL) and
`round4_dns` (DNS); round 3 (review round 1) in `round3`, `round3_dns`; rounds 1 and 2 in the
parent directory. The datasheet and the s10 RUNBOOK are not edited here.

## The defect

The pressure multigrid solves its coarsest level on the master (rank 0) with the Schur
matrix C = B^T M^-1 B that the master assembles on its unpartitioned coarse mesh
(`Create_CMat`, `Setup_UMFPACK_CoarseSolver`). On that mesh the two images of a periodic
velocity dof are two separate dofs: the master's C has no coupling across a periodic plane
and treats the plane as a free boundary, while the workers apply the periodic operator (their
C plus the parallel part CP, whose ghost columns include the periodic neighbours). The coarse
correction belongs to another operator. On the s10 boxes (xy periodic, z walls) the pressure
multigrid stagnated: 12 x 12 x 120, every solve at the cap of 20 cycles, DefFinal/DefInit
0.54 (geometric mean); single phase 0.82. The legacy master matrix carries entries the
periodic operator does not have, e.g. 1.37e-2 where the periodic operator has 0 on the
triply periodic QBOX9 (check line of `B_qbox9_gath_chk`).

Three further defects of the coarse solve, fixed on the new paths only (the legacy path keeps
its arithmetic, it replays the certified periodic logs):

- UMFPACK orientation: `myUmfPack_Factorize` hands the row-compressed matrix to UMFPACK as
  column-compressed storage, so types 2/4 factorise and solve the transpose (equal for the
  symmetric matrices of the certified runs; shown on a non-symmetric 4 x 4 matrix in
  `test_pco_coarse`); no status was checked after symbolic, numeric or solve.
- Gauge under `NoOutflow`: setup zeroes the off-diagonals of row 1 of the live matrix
  (`QuadSc_solver_coarse.f90`), the right-hand side entry is not zeroed; the smoother and the
  defect then see the modified operator.
- Types 3/4 average the P0 solution to the vertices (`INTPVB`) with the element stars of the
  unpartitioned mesh: the two images of a periodic vertex get different values.

## The two operators

**gathered.** What the workers apply on the coarse level is, for the rows of their own
elements, (C x)_r + (CP xp)_r with xp the ghost values of the pressure exchange
`GetParPressure`. Both matrices are assembled on the workers with the actual coefficients
(weighted EL mass, slip projection, constraint masks, curved-boundary geometry). The master
collects them into its numbering: the workers exchange the global numbers of their own unknowns
with `GetParPressure` (as the MUMPS exporter `MUMPS_SetUpP1_SLAVE` does) to label the ghost
columns, the element map is the one of the coarse vector transfers (`E012DISTR_L1/2/3`: coarse
element link, then child 1 keeps the parent number, children 2..8 of parent p of a level with
N parents follow at N + 7(p-1) + 1..7). A plan (`PCO_GATHER_BUILD`) maps every entry of the
C and CP storage of the workers to a slot of the master layout; every renewal the workers send
their value arrays as they are (`MPI_Gatherv`) and the master sums them into the slots.
Exact by construction; the check confirms it against the worker product to 1e-15.

**wrapped.** The master assembly with the velocity dofs of each periodic class identified:
C = (B^T Q) M_cls^-1 (Q^T B), Q copying a class value to its images. The B rows and B^T columns
of the images are merged per (element, mode) before the product (no first-match loss of
terms), all images count (edges 4, triply periodic corners 8), the class mask is the common
Dirichlet mask (equal over the class, checked component by component), the class mass is the
sum of the raw image masses, under the weighted EL operator U and W summed separately and
transformed once, m = (sum U)^2 / (sum W), as the workers do (`QuadSc_el_mm.f90`). It is the
master geometry. On the s05 pipe (curved wall) the wrapped entries differ from the gathered
ones by up to 18 % (row-scaled) next to the curved wall, where the legacy master matrix has the
wrapped value; the master does not project the coordinates of its finest level (the Q2 dofs of
the coarse level; `ParametrizeBndr` of level MedLev+1 runs on the workers only, app_init). The
discrepancy is consistent with that missing projection; it is not proven to be the only cause
(a worker-consistent projection private to the wrap was not implemented). Therefore wrapped is
compared with the gathered operator at every renewal and stops above 1e-10.

**Default.** Inside one table (`default_scope` in the coordinator): the absent key selects
gathered only for configurations covered by gates (table below); everywhere else the absent key
keeps legacy with the reason logged, and gathered or wrapped are available on explicit request
after validation. Whether the default should become wider is an open owner decision.

| configuration | default (key absent) | covering gates (round 4, job ids in the gate table) |
|---|---|---|
| q2p1_el_pipeflow, validated periodicity, coarse level 1 (MinLev = MedLev = 1), solver type 1-4 | gathered | operator identity A, gathered vs wrapped B, solvers C (types 1-4), probes D, legacy bitwise E |
| q2p1_el_pipeflow, coarse level 2 or 3 | legacy (logged) | explicit gathered: X_lev2_gath_chk, X_lev3_gath_chk (operator identity through the level-2/3 transfers) |
| other applications of the common driver (q2p1_dns_drag, q2p1_bench_sedimentation, q2p1_dkt, ...) | legacy (logged) | DNS: d11 and QBOX2 (periodic, explicit gathered vs legacy: drag and bulk flow; QBOX2 the census on two cells per axis), e4_l3 twin bitwise under absent / legacy / gathered / wrapped |
| q2p1_cc | legacy assembly, the renewal is not entered | not supported (own coupled coarse solve) |
| no validated periodicity (no finite period, or no accepted periodic dof link) | legacy | E (non-periodic gates), F (walled boxes with declared periods, 0.6-cell periods, master hierarchy) |

## Conventions

- Numbering: unknown 4(e-1)+m is mode m of the P1 pressure of master element e on the coarse
  level: m = 1 the element mean, 2..4 the coefficients of the physical linear basis x-x0,
  y-y0, z-z0 (B assembly `QuadSc_BMatrix.f` and `IntQ1toP1` use this physical basis; the
  transferred coefficients therefore do not depend on the local vertex order of an element,
  only on the correspondence of the physical elements and their centres).
- Layout (`PCO_INSTALL_LAYOUT`): the block row of e lists the own block first, then the blocks
  of the neighbour elements ascending; row 4(e-1)+m starts with its diagonal, then the other
  modes of e ascending, then four entries per neighbour (the legacy `Get_CMatStruct` layout, so
  the mean-mode system is read with stride 4). Neighbours: the element graph through vertex
  classes (elements sharing a vertex or one of its periodic images).
- Raw operator (`pco%a`): ungauged, never modified by a solve; it serves the type-3/4 smoothing
  (`PCO_RAW_SOR`, the sweep of `crs_E012_SOR_Smoother`; a row with a zero diagonal keeps its
  value) and the checks. A non-finite entry stops the run.
- Gauge, decided once per layout at the first refresh: components from every nonzero coupling
  (elements coupled when any entry of their 4x4 block is nonzero; a weak but real coupling, e.g.
  1e-12 against order-one diagonals, joins them, so one pressure constant is pinned once; two
  regions touching only through constrained dofs share an all-zero structural block and are two
  components). Each component indicator e_c is checked on its own: closed when
  max_i |(C e_c)_i| / sum_j |C_ij| (j over the mean columns) < 1e-8 with i over EVERY row, inside
  and outside the component (measured on the gates: closed 1.7e-15 to 2.0e-14, open 0.9; the
  tolerance lies six orders from both); e_c^T C likewise over every column. One reference per
  closed component, the mean mode of its lowest element, set to zero on a solver-owned copy
  (row and column cleared except the diagonal, right-hand-side entry zeroed: symmetric
  elimination of a zero value); a zero diagonal of a reference row (a closed component whose
  mean row vanishes) becomes the largest diagonal modulus of the copy, so the equation reads
  d x = 0. For the SSOR types 1/3 a zero diagonal elsewhere in the copy is refused at the
  refresh. Later refreshes keep the decision (no flipping) and
  report drift (closed above 1e-6 or open below 1e-10). The decision is taken from the operator,
  not from `SimPar@NoOutflow`; a disagreement is logged (`PCO gauge note`). The left measure
  e^T A of the closed components is reported; the right-hand-side compatibility measure assumes
  the left null vector e (symmetric operator). A right-hand side whose sum over a closed
  component does not vanish is reported (component, net value, relative imbalance; above 0.1
  relative: late-cycle defects of compatible single-phase runs reach 2e-2 at defects of
  1e-14), not changed.
- Solver copy: full system (types 1, 2) or mean modes (types 3, 4), pattern fixed (explicit
  zeros kept). Types 2/4: UMFPACK through `tUmfCSC`, the compressed-column storage of the copy
  (UMFPACK solves A x = b), symbolic analysis once per layout, numeric factorisation at every
  refresh, every status checked (factorisation: collective stop with the reason; solve, on the
  master inside the multigrid: message and `MPI_Abort`). Types 1/3: BiCGStab of the solver
  module with an SSOR preconditioner that scans every off-diagonal entry with a lower/upper
  predicate (the legacy `PARID117_MASTER` stops at the first upper column, which in this
  diagonal-first layout skips lower neighbours); a zero right-hand side returns zero,
  convergence is tested before every division, a breakdown aborts with a message.
- Types 3/4 (`crs_P0_Correction`, shared by both master coarse solves): P0 solution x3, vertex
  averaging, `IntQ1toP1` on the physical elements, SOR smoothing. With a periodic operator the
  averaging sums numerators and denominators over each periodic vertex class, divides once,
  gives the value to every image (`PCO_P0_TO_VERTEX`, in its own nvt-sized arrays), and the
  smoother uses the raw operator; the legacy branch is `INTPVB` and `crsSmoother` as before.

## Code

| module (file) | responsibility |
|---|---|
| `pressure_coarse_key` (src_util/pressure_coarse_key.f90) | key values, names, parsing |
| `pco_periodic_topology` (QuadSc_pco_topology.f90) | link census on the workers (shared vertices, `PCO_PAIR_SHIFT`); dof classes on the master (coordinate hash, all combined shifts, union-find, completeness / ambiguity / kind checks, axes realised); constraint conflicts (`PCO_MASKS_DIFFER`); element graph through vertex classes; master coordinate check |
| `pco_gather` (QuadSc_pco_gather.f90) | global element/unknown ids, gather plan with validation, value collection, immutability check |
| `pco_wrap` (QuadSc_pco_wrap.f90) | wrapped structure and assembly |
| `pco_coarse_solver` (QuadSc_pco_solver.f90) | master storage: layout, gauge decision, solver copy, UMFPACK / BiCGStab solves, raw SOR, solver type change, diagnostic solves |
| `pco_p0_reconstruct` (QuadSc_pco_reconstruct.f90) | periodic vertex classes and the periodic P0 -> vertex averaging |
| `pco_diagnostics` (QuadSc_pco_check.f90) | operator comparison, the check (worker product, closure, solver check, imbalance probe) |
| `pco_coarse_operator` (QuadSc_pco_operator.f90) | `Renew_Pressure_Operator`, `PCO_REGISTER_APPLICATION`, `PCO_PRESSURE_RELOAD_ADMIT`, `PCO_PRESSURE_PARAMETERS_RELOADED`; mode decision, default scope, validation and refusal, pattern check before the assembly, renewal, wrapped validation, logging |
| `UMFPackSolver` (UmfpackSolver.f90) | `tUmfCSC` adapter; legacy handles tracked one by one, the legacy solve stops without a live factorisation |
| `PP3D_MPI` (src_mpi/pp3d_mpi.f90), `E013_CreateComm*` (QuadSc_mpi.f90) | `MGE013(lev)%nPerLink`: the accepted velocity dof links through the period of each axis, counted where the links are made (`COUNT_PERIODIC_LINK`; the exchanges are unchanged) |

Call sites: `OperatorRegenaration` (C branch), `EL_MM_Pressure_Operator` (every step, hands W
and rho; the wrap forms U = rho Ml on the master when it needs it; the release of the previous
legacy UMFPACK factorisation for types 2/4 stays there, as before), `ProcessControl`
Reload_Pres (`PCO_PRESSURE_RELOAD_ADMIT` before `ValidateSolverTypes`,
`PCO_PRESSURE_PARAMETERS_RELOADED` after it), q2p1_el_pipeflow `General_init_ext`
(`PCO_REGISTER_APPLICATION`). In the multigrid one branch in `mgCoarseGridSolver_P`
(`PCO_ACTIVE`) and the shared type-3/4 tail.

q2p1_cc is not a user: its coarse solve is the coupled saddle-point system with its own
parameters (`CCuvwp`), so it keeps its legacy assembly and set-up (`Create_CMat`,
`Setup_UMFPACK_CoarseSolver` when the master assembles, `Create_ParCMat`); a periodic coarse
operator for q2p1_cc is unsupported pending a dedicated integration.

Invariants. The layout, the classes, the plan and the gauge belong to one coarse pattern. Every
renewal checks, BEFORE `Create_CMat` / `Create_ParCMat` write into their kept allocations, the
pattern, the ghost numbering and the ownership of the workers against the plan and the extent
of every C / CP value array against its pattern; after the assembly the absence of pressure
Dirichlet elements; a change stops the run (restart). A runtime reload of the pressure
parameters: the file must be visible on every rank (agreed by an allreduce in
`ReloadParameterFile`, so every rank enters the reload or none); the ranks agree on the
proposed (solver type, MinLev, MedLev) (allreduce min/max; disagreement: old values
everywhere); a change of the master levels or a type outside 1-4 is undone before the generic
validation (which would abort on it); a new type 1-4 rebuilds the solver copy (same layout
and gauge); the key and the check value are read once. The solves check the dimensions of
the installed copy.

Collective operations per periodic renewal: the pattern and extent check before the assembly
(one pressure exchange, one allreduce); `Create_ParCMat` (as always); the knprP allreduce; the
value gather (gathered mode, and in wrapped mode for its comparison: one `MPI_Gatherv`); in
wrapped mode the broadcast of the comparison result; the status broadcast (and the message
when it fails); when due, the check (one gather, three distribute/exchange/collect cycles). At
the first renewal in addition: the census (six pressure exchanges and two allreduces: the
structure flag; the consistency and link counts), the plan (one exchange, gathers,
broadcasts) and the broadcasts of the master validation. Per reload under a periodic
operator: two allreduces (the agreed tuple) and, for a new solver type, one status
broadcast.

Extension points: a further operator is a key value and a branch in `refresh_raw`; a further
default configuration is an entry of `default_scope` with its gate; a master hierarchy
(MinLev < MedLev) needs an operator on every master level and a tested master cycle; another
coarse solver is a branch in `solve_copy` on the existing copy.

## Deck keys and activation

`Pres@PeriodicCoarseOperator` (absent, `legacy`, `gathered`, `wrapped`; invalid values stop
every rank in the parser) and `Pres@PeriodicCoarseCheck = n` (0 off; n > 0: the check at
renewals 1, 2 and every n-th). Decision order at the first renewal: (1) no finite period ->
legacy, silently when the key is absent (nothing is entered and nothing logged, so the logs
of walled runs stay byte-identical; an explicit gathered/wrapped gets one line); (2) key legacy
-> legacy; (3) census: no accepted periodic dof link, or no pressure structure on the coarse
level -> legacy (any key; logged); (4) default scope (key absent); (5) validation.

Census. A finite period alone does not establish periodicity (q2p1_dns_drag declares unit
periods on every box). The communication set-up (`E013_CreateComm*`) accepts a velocity dof
link between two ranks when the coordinates coincide or lie one period apart, axis by axis;
the links that matched through the period of axis a are counted there (`MGE013(lev)%nPerLink`).
Every periodic coupling of the worker operator goes through such links (the pressure ghost
elements of CP are the elements at the linked dofs, `ParPresComm_Init`), so axis a is active
exactly when the global count is positive. Round 3 re-derived this from the geometry of
coupled element pairs and lost the axis when a pair shares an ordinary vertex AND a periodic
image (two coarse cells along a periodic axis, one per rank); the provenance count does not.
A walled mesh has no such link (adjacent cells match coincident dofs; a period that is not a
distance between interface dofs matches nothing). In addition every coupled element pair of CP
must share a vertex unshifted or in some periodic image (`PCO_PAIR_SHIFT`), otherwise the
periods do not fit the mesh (nbad: fallback or refusal).

| situation | key absent | `legacy` | `gathered` / `wrapped` |
|---|---|---|---|
| no finite period on any axis | legacy, nothing entered | legacy | legacy, one log line |
| finite periods, no accepted periodic dof link (walled box, also with cells of 0.6 period, MinLev < MedLev, any solver type) | legacy, logged | legacy, logged | legacy, logged |
| validated periodicity, inside the default scope | gathered | legacy (logged) | the requested one |
| validated periodicity, outside the default scope | legacy, logged with the reason | legacy | the requested one, after validation |
| couplings sharing no vertex in any image | fallback, logged | legacy | refused |
| coarse solver type outside 1-4, MinLev /= MedLev, coarse level > 3, master not solving | fallback, logged | legacy | refused |
| pressure Dirichlet elements (knprP) on the coarse level | fallback, logged | legacy | refused |
| master without coordinates of level MedLev+1 | fallback, logged | legacy | refused |
| ambiguous or incomplete classes, mixed dof kinds, unequal Dirichlet flags over a class, an active axis without master images | fallback, logged (counts, first dof and its coordinates) | legacy | refused |
| a worker coupling outside the master layout, a master row with two or no owners | fallback, logged | legacy | refused |
| wrapped differs from gathered (any renewal) by more than 1e-10, or a non-finite entry | - | - | stop |

Refusal: every rank reaches the same decision (counts broadcast from the master), the master
and the protocol rank print `PCO REFUSED: ...`, every rank stops. Output: `PCO periodic coarse
operator = ...` (mode, level, active axes, accepted periodic dof links per axis, summed over the
ranks, "shifted links"), `PCO classes ...`,
`PCO gauge ...` on the terminal and the protocol; `PCO_TIMING`, `PCO wrapped ...` and
`PCO_CHECK` lines on the terminal of the master.

The comparison (`PCO_COMPARE_OPERATORS`, also the wrapped validation): both operands are
checked with `ieee_is_finite` first (a NaN or Inf in either is an infinite difference); rows
are scaled by the largest modulus of the gathered row; a row that is zero in the gathered
operator must be zero in the wrapped one (a nonzero entry there is a difference of 1).
The check (`pco_diagnostics`): both operators; the entry-by-entry comparison gathered -
wrapped (with the location and the legacy master value); for three vectors
(mean modes = 1; a four-mode impulse on the element at the largest vertex class; a
reproducible vector exciting all modes) the worker product with the E012_DAX arithmetic
(distributed and collected with the coarse vector transfers) against the master products of
the gathered, the wrapped and, under the weighted EL mass, the per-image transformed wrap (the
prototype seam mass, to show that the check discriminates), every row, scaled by
sum_j |A_ij| max|v|; closure, left null vector, a compatible solve with the ungauged residual,
the reference value, and the imbalance of an incompatible right-hand side. Diagnostic solves
and probes do not touch the production statistics.

## Tests and gates

Binaries, all built on longhouse1a from the clean tree at 34556706 (the reviewed state of review
round 2): `q2p1_el_pipeflow.34556706` (EL recipe of build-perk-lh, md5 ab25d782),
`q2p1_bench_sedimentation.34556706` and `q2p1_dns_drag.34556706` (the certified DNS twin recipe:
gcc 13.2 wrappers, serial PE, HardContactAndFluid, FBM acceleration), the same DNS targets of the
base 16ffcbcc (`*.base`, worktree FeatFloWer-pcoarse-base) and, to show the census defect, the
round-3 DNS binary `q2p1_dns_drag.1004a4de`. The tree builds with USE_PE=OFF (q2p1_fac3d,
q2p1_fc_ext, build-nope) and the changed q2p1_cc source compiles (its application does not link
in this tree, before this change). Runs: `/data/warehouse17/rmuenste/runs/pcoarse/round4` (EL,
report.txt) and `round4_dns` (`tools/dns_report4.py`); every directory has versions.txt. Probe
metrics from the `XxX` lines: cycles per solve, geometric mean of DefFinal/DefInit, `res1_L1`
of the last `EL_VA_FLUX` line. Two intermediate binaries of review round 2 were discarded after
their first runs (round4a_c36b65d3, round4b_989936c9): the new 2 x 2 x 2 gate exposed the
closure scale (normwise now) and the left measure of the check.

Unit tests (`test_pco_coarse`, ctest `pco-coarse-serial`, job 163070; 234 assertions):
UMFPACK orientation, status, double release; classes 2/4/8 images, incomplete/ambiguous
refusals, neighbours, periodic averaging; gauge on closed and open rings, types 1-4; the 64
Dirichlet mask pairs; two closed blocks sharing a structural zero block (two references) and an
open/closed pair; BiCGStab zero right-hand side and exact convergence at the half step;
averaging on a 1 x 1 x 5 strip; solver type change 2 -> 4 -> 1; census pair consistency;
operator comparison with NaN/Inf in either operand and a small entry in a zero oracle row; two
closed rings joined by a conservative edge of 1e-12 (one component, one reference, the
compatible solution with the pressure jump 1 across the bridge, types 2 and 4); a closed
singleton with a vanishing mean row through full (types 1, 2) and reduced (3, 4) solves, the
raw SOR on its zero row, a zero diagonal outside a reference refused for SSOR; the E013
acceptance test on two cells along a periodic axis, one per rank (periodic links on x counted
while the pair classification says unshifted) and on 0.6-period cells (none); a closed ring
with rounding-noise slope-mean couplings stays closed.

**A. Operator identity** (worker product with the E012_DAX arithmetic against the master
product, all rows, normwise row scale; checks at renewals 1, 2, 10, 20, maxima over them)

| case | run (job) | dax - gathered | dax - wrapped | per-image mass variant | verdict |
|---|---|---|---|---|---|
| 12x12x120 (coarse 6x6x60), weighted EL | D_m12_gath_chk (163079), D_m12_wrap_chk (163081) | 6.8e-16 | 1.4e-14 | 4.1e-04 | PASS, 28 lines, 0 FAIL |
| 6x6x60 (coarse 3x3x30), weighted EL | D_s6_gath_chk (163083), D_s6_wrap_chk (163085) | 7.7e-16 | 9.0e-15 | 5.9e-05 | PASS |
| 12x12x120 single phase | D_sp_gath_chk (163089), D_sp_wrap_chk (163092) | 5.3e-16 | 1.4e-14 | = wrapped | PASS |
| coarse level 2, explicit | X_lev2_gath_chk (163157) | 4.9e-16 | 2.0e-14 | - | PASS |
| coarse level 3, explicit | X_lev3_gath_chk (163159) | 5.7e-16 | 4.8e-14 | - | PASS (7 423 488 entries) |
| DNS d11, QBOX9 triply periodic | d11__gathered (163086), d11__wrapped (163087) | 3.9e-16 | 2.3e-15 | - | PASS |
| DNS QBOX2, 2 x 2 x 2 coarse cells, one per worker | q2__gathered (163094), q2__wrapped (163096) | 1.6e-16 | 8.0e-16 | - | PASS, 0 FAIL |

**B. Gathered vs wrapped entry by entry** (difference / largest modulus of the row)

| case | run (job) | entries | max difference | verdict |
|---|---|---|---|---|
| xy periodic, z walls | D_m12_gath_chk (163079) | 922 752 | 1.1e-13 | PASS |
| z planes Symmetry001 | B_sym_gath_chk (163109), B_sym_wrap_chk (163110) | 922 752 | 1.1e-13 | PASS |
| 3x3x6 partition (55 ranks) | B_p336_gath_chk (163111), B_p336_wrap_chk (163112) | 114 048 | 6.6e-14 | PASS |
| triply periodic QBOX9 (corner class of 8) | B_qbox9_gath_chk (163104), B_qbox9_wrap_chk (163105) | 314 928 | 1.6e-14 | PASS |
| QBOX2 (every element pair coupled ordinarily and periodically) | q2__gathered (163094) | 1 024 | 1.2e-15 | PASS |
| z-periodic pipe, s05 E1 (curved wall) | B_pipe_gath_chk (163106) | 482 112 | 1.8e-01 | wrapped FAIL (dax - gathered 5.5e-16 PASS); wrapped refused: B_pipe_wrap_chk (163107, exit 1, `PCO REFUSED`) |

Wrapped runs compare with the gathered operator at every renewal (1.1e-13 at renewals 1, 2, 20
of D_m12_wrap, 163075).

**C. Gauge and solvers**

| item | run (job) | result |
|---|---|---|
| closed periodic cases | all `_chk` runs | one component, closed, closure 7.7e-16 (QBOX2) to 2.0e-14 (coarse level 3), e^T A the same, ungauged residual below 1e-15, reference 0, incompatible right-hand side reported |
| open periodic case (zmax Outflow) | C_open_gath_chk (163103) | closure 0.90: no reference; gauge note; 2.00 cycles, 0 FAIL |
| type 1 (BiCGStab, full) | C_m12_gath_t1 (163098) | 2.05 cycles, 3.94e-05, res1 3.3445e-11 |
| type 2 (UMFPACK, full) | C_m12_gath_t2 (163099), C_m12_wrap_t2 (163101) | 2.00, 1.38e-05, 3.3442e-11; gathered and wrapped identical |
| type 3 (BiCGStab, P0) | C_m12_gath_t3 (163100) | 2.00, 1.82e-05 |
| type 4 (UMFPACK, P0) | D_m12_gath (163072) | 2.00, 1.81e-05, 3.3442e-11 |

**D. End-to-end probes**

| case | legacy (job) | gathered (job) | wrapped (job) |
|---|---|---|---|
| PR15_m12_probe, weighted EL | 20.00 cycles (40/40 at cap), 0.539, res1 9.81e-08, bitwise vs 162042 (163077) | 2.00, 1.81e-05, 3.3442e-11 (163072) | identical XxX lines (163075) |
| PK_smoke_wrap2, 6x6x60 | 20.00, 0.411, prefix-bitwise (163088) | 2.05, 6.28e-05, 3.97e-12 (163083) | same (163085) |
| MG_sp single phase | 20.00, 0.822, bitwise (163095) | 2.00, 4.97e-06 (163089) | 2.00, 6.46e-06 (163092) |
| MG_crs2 (type 2) | 20.00, 0.541, bitwise (163102) | see C | see C |
| coarse level 2 (MG_L3 deck) | key absent = legacy: 20.00, 0.462 (X_lev2_unset, 163158) | explicit: 2.00, 1.27e-04 (163157) | - |
| DNS d11, 40 steps | 9.25 cycles, 1.77e-06 (d11__legacy 163084) | 3.33, 6.00e-07 (163086) | 3.33, 6.00e-07 (163087) |
| DNS QBOX2, 10 steps | 3.60, 4.60e-07 (q2__legacy 163093) | 3.60, 4.60e-07 (163094) | 3.60, 4.60e-07 (163096) |

DNS d11, gathered (and wrapped) against legacy over 40 steps: axial drag (particle_force.log fz)
at most 6.8e-07 relative (fz(0.4) 8.82883971e-03 vs 8.82884008e-03), superficial velocity
(bulk_flow.log) 3.0e-08: within the solver tolerance of the legacy run. QBOX2: the printed drag,
bulk flow and `XxX` lines are identical to legacy (the coarse correction of 8 elements under a
level-4 hierarchy does not change the printed digits); its role is the activation (P below).

**P. Census on two coarse cells per periodic axis** (QBOX2: unit cube, 2 x 2 x 2 coarse cells,
one per worker, 8 workers + master; every element pair shares an ordinary vertex and periodic
images; mesh written by `tools/make_box_mesh.py`, which reproduces the certified QBOX9 mesh
byte for byte)

| binary, key | run (job) | result |
|---|---|---|
| round 3 (1004a4de), gathered explicit | q2__r3gathered (163097) | `PCO periods declared ... no worker coupling crosses a periodic plane: periodicity inactive, legacy operator` (the defect: an explicit request silently legacy) |
| 34556706, gathered explicit | q2__gathered (163094) | `axes TTT shifted links 200 200 200`, 37 classes (largest 8), operator identity PASS |
| 34556706, absent / legacy | q2__unset (163091), q2__legacy (163093) | bitwise vs the base binary (163090) |

**E. Bitwise** (EL: `el_stream_compare.py --strict` on every EL stream, `XxX` lines as text,
particle logs byte for byte; DNS: particle_force.log, bulk_flow.log, the `SED_BENCH_*`,
`DNS_PART_FORCE` and `XxX` lines)

| case | key | runs (jobs) | result |
|---|---|---|---|
| 4 non-periodic EL gates (SCOL4, VFT10 x 3) | absent / gathered / wrapped | E_gates (163119-163130) | 12/12 bitwise |
| 6 periodic QBOX9 EL gates | legacy | 163131, 163133, 163135, 163137, 163139, 163141 | 6/6 bitwise |
| 6 periodic QBOX9 EL gates | absent (= gathered) | 163132, ..., 163142 | not bitwise (expected); 2.03-2.28 cycles |
| A15 wall-bounded batch | absent / gathered / wrapped | E_A15__* (163143-163145) | prefix-bitwise vs A15, bitwise vs the base binary |
| P15_probe_cross twin | legacy | E_cross_off__legacy (163146) | bitwise |
| s05 pipe | legacy | B_pipe_legacy (163108) | bitwise vs the base binary |
| DNS e4_l3 twin (q2p1_bench_sedimentation, 32 ranks; certified deck, not periodic) | absent / legacy / gathered / wrapped vs base binary | e4l3__* (163071 base, 163073, 163074, 163076, 163078) | 4/4 bitwise on every stream |
| DNS d11 (triply periodic, unregistered application) | absent / legacy vs base binary | d11__unset (163082), d11__legacy (163084), base 163080 | bitwise (key absent is legacy outside the default scope, logged) |
| DNS QBOX2 | absent / legacy vs base binary | q2__unset (163091), q2__legacy (163093), base 163090 | bitwise |

**F. Refusals and inactive periodicity**

| case | run (job) | result |
|---|---|---|
| walled A15 box, unit periods | F_wall_per1 (163113), F_wall_per1_gath (163114) | periodicity inactive, logged; bitwise vs base |
| walled A15 box, periods 0.6 of a coarse cell (x, y) | F_walled06__{unset,gathered,wrapped} (163148-163150), base 163151 | periodicity inactive (no accepted periodic dof link), bitwise vs base under every key |
| walled A15 box, unit periods, MinLev 1 < MedLev 2 | F_walledML__{unset,legacy,gathered,wrapped} (163152-163155), base 163156 | periodicity inactive, bitwise vs base under every key (no refusal) |
| x pair Wall / Periodic | F_conflict_default (163115), F_conflict_gath (163116) | fallback logged (1428 conflicts) / `PCO REFUSED`, exit 1 |
| periodic, MinLev 1 < MedLev 2 | F_minmed_default (163117), F_minmed_gath (163118) | legacy logged (outside the default scope) / refused |
| wrapped on a curved boundary | B_pipe_wrap_chk (163107) | refused |
| explicit key with coarse solver type 5, 7-10 | - | not reachable in this build (ValidateSolverTypes stops first) |

**G. Restart** (round 3; DNS d11, q2p1_dns_drag; a continuous run of 20 steps against a run of
10 steps, dump, restart with `StartingProc = 1`; the restarted run went on to t = 0.30, the
continuous control ended at t = 0.20, so the overlap is t = 0.11 to 0.20, ten rows)

| key | runs (jobs) | result over the overlap |
|---|---|---|
| gathered | d11r_gathered_cont (162918), _seg (162919), _rst (162922) | the restart works (the operator, classes, plan and gauge are rebuilt from the restarted state); the axial drag fz of particle_force.log and the complete bulk_flow.log rows are identical in all ten rows; the complete particle_force.log rows are not (0 of 10 byte-identical: the transverse components, of order 1e-12 against fz 5e-3, differ at the level 2e-16 of fz, e.g. fx 1.61080951e-12 continuous, 1.61080942e-12 restarted at t = 0.11) |
| legacy (control) | d11r_legacy_cont (162920), _seg (162921), _rst (162923) | the same: fz and bulk flow identical, the transverse components differ the same way (0 of 10 complete rows identical, 2.2e-16 of fz) |

What is shown: the restart works under gathered, and the named observables (axial drag, bulk
flow) are identical over the overlap. A bitwise identical restart of the complete state is
not shown; the differences of the transverse force components are present under legacy as
well (a pre-existing property of the restart, not of the operator).

The EL application cannot show this: q2p1_el_pipeflow writes no dumps (the dump writer of
FullC0ntact skips on the rank whose parInfo id is 0, the run logs "cannot write format.dmp"),
so an EL restart is not available in this tree at all, with or without this change.

## Cost (12 x 12 x 120 probe, 28 ranks, 20 steps, weighted EL, type 4)

Master side, per renewal (the weighted EL operator renews every step; `PCO_TIMING` lines,
seconds summed over 20 renewals): gather of the values 0.056 s (2.8 ms each), wrapped
assembly 0.43 s (21.6 ms each), refresh with the numeric factorisation of the P0 copy 0.146 s
(7.3 ms each); type 2 (full 8640-unknown copy, 922 752 entries) 4.35 s (218 ms each); one-time
set-up (census, classes, layout, plan) 0.071 s (gathered) and 0.077 s (wrapped, including the
start-up comparison with the gathered operator). The gather of 922 752 values is about 7.4 MB
per renewal into the master.

Whole run (overall time of the protocol rank, 19 steps after start-up; round 2): legacy
60.9 s (pressure multigrid 32.1 s), gathered 32.1 s (3.64 s), wrapped 32.6 s (3.80 s); 3.20,
1.69 and 1.72 s per step (jobs 162769, 162767, 162768, all on nx nodes shared with other gate
jobs; indicative). Round 3: legacy 61.8 s, gathered 32.4 s, wrapped 33.1 s (162838, 162836,
162837); wrapped gathers and compares at every renewal (0.43 s for 20 renewals including the
assembly, comparison 0.11 s). Round 4: legacy 75.0 s, gathered 31.8 s, wrapped 32.5 s (163077,
163072, 163075; the legacy run on a busier node); the pattern check before the assembly
replaces the one after it (same cost). The worker-side assembly of C is the same in all three (1.16 s).

## Not covered

- Master multigrid hierarchy (MinLev < MedLev): refused / fallback. The master F cycle
  (`crs_F_cycle`, loop from MedLev+1 to MedLev-1) is a V cycle, the W recursion stops at a
  hard-coded level 2, the transfers exist for levels 1-3 only.
- Coarse solver types 5 (MUMPS, refused for periodic runs by `ValidateSolverTypes`), 7, 8, 10
  (HYPRE) and 9 (PARDISO): not built in this tree; they stop in `ValidateSolverTypes` before
  the renewal, so the refusal of an explicit key with these types could not be exercised.
- DNS applications are outside the default scope (key absent = legacy, bitwise). Explicit
  gathered/wrapped is shown on d11 and QBOX2 (fixed sphere, triply periodic); moving particles
  in a periodic DNS under the new operator are not gated.
- Runtime reload of the pressure parameters (Reload_Pres): the solver type change is covered by
  the unit test; the collective agreement (file, parameter tuple) and the admission order are
  not exercised by an application run.
- Explicit admission extends to every application of the common driver whose configuration
  passes the validation (not only the gated ones); the default does not.
- A weighted EL restart: q2p1_el_pipeflow writes no dumps (see G).
- E013_CreateComm_read (not called by any application in this tree) does not count periodic
  links: an application using it would see periodicity inactive (legacy, logged).
- Rotational or non-translational periodicity; meshes whose two periodic faces do not match
  node by node within DEpsPrec (refused as incomplete classes or couplings outside the layout).
- The legacy path keeps its defects (transposed UMFPACK solve, row gauge without the
  right-hand side, one-sided type-3/4 averaging, master geometry on curved boundaries).
- q2p1_cc: keeps its legacy assembly and coarse set-up (restored); its source compiles, the
  application itself does not link in this tree (`app_init.f90`, `PartitionReader.f90`:
  `ctemp` without type, before this change). A periodic coarse operator for its coupled
  saddle-point coarse solve would need a dedicated integration.

## Follow-ups

- Master hierarchy: operators on every master level (gathered per level is direct), a correct
  master F/W cycle, a test with MinLev < MedLev.
- HYPRE/MUMPS/PARDISO coarse paths with a periodic operator.
- Default for DNS apps and re-certification of the periodic certified decks under the
  default (the six QBOX9 EL gates change their iteration history under the default).
- The master geometry on curved boundaries (project level MedLev+1 on the master): changes the
  legacy arithmetic; would make wrapped usable on curved periodic meshes.
- Wider default (other applications, coarse levels 2/3): an owner decision; each needs its
  gate entry in `default_scope`.
