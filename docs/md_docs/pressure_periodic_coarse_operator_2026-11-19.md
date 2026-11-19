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
- Gauge, decided once per layout at the first refresh, from the structure of the operator;
  the matrix is a consistency check.
  - Components: elements coupled when any entry of their 4x4 block is nonzero (a weak but real
    coupling, e.g. 1e-12 against order-one diagonals, joins them, so one pressure constant is
    pinned once; two regions touching only through constrained dofs share an all-zero
    structural block and are two components).
  - Structure (authoritative): does a constant pressure on the component produce a gradient
    at a free velocity dof? With e = 1 on the mean modes, g_d(i) = sum_j B_d(i,j) e_j is that
    gradient at velocity dof i, component d (B rows of the workers: own and ghost elements,
    the periodic images included). Per coarse element and component d the workers form
    max over the FREE dofs of the element of |g_d(i)| divided by max over all its dofs of
    sum_j |B_d(i,j)| (mean-mode columns; free: the masks of `Get_CMat`; slip enters through the
    projected B of `SetSlipOnBandBT`; `PCO_FREE_FLUX`). Where the support of every free dof
    closes (interior, periodic seam) the contributions cancel to rounding; at a free boundary
    dof they do not (order one: 0.5 to 1 measured). The scale is the size of the B entries of
    component d around the element, not the row of the dof itself: element and face bubbles
    have rows whose single entries vanish analytically, so a row-relative ratio would be
    noise / noise (measured 1.0 on closed boxes in a first attempt); the scale is of the same
    component, so the ratio does not depend on the cell aspect ratio (open 0.5 at stretch 1,
    10, 30, 1e4, measured). The largest value per coarse element is collected to the master at
    the first renewal; per component s_c = the largest over its elements. s_c <= 1e-10:
    closed; s_c >= 1e-6: open; in between: undecidable. In wrapped mode the structure comes
    from the workers as well and the matrix check below is applied to the wrapped operator, so
    the two must agree (and the wrapped operator is compared with the gathered one entry by
    entry at every renewal).
  - Matrix (consistency): right(c) = max over the rows i of |(C e_c)_i| / sum_j |C_ij|, j over
    ALL columns (the normwise scale), and e_c^T C over every column scaled by the column sums
    of moduli; one O(nnz) traversal (a nonzero entry between two components would contradict
    the decomposition and counts as an inconsistency). A closed component must have
    right(c) <= 1e-8; an open one right(c) > 1e-13 (otherwise the constant is numerically in
    the kernel of a structurally regular operator). A weak outlet (cells stretched across it)
    gives an open component with a small matrix measure (2.5e-10 in the unit test, below the
    former matrix-only threshold 1e-8): it stays open, no reference, and is solved exactly.
  - Undecidable or inconsistent: nothing is pinned; at the first renewal the run is refused
    (explicit key) or falls back to legacy with the reason logged (key absent).
  - Types 3/4 with a weakly open component: the over-relaxed (x3) P0 correction diverges when
    an open component is only weakly coupled to its outlet. Measured on the 12 x 12 x 120 box
    with a z outlet (single phase): matrix measure of the open component 0.90 (cubic cells):
    2 cycles per solve; 1.2e-2 (cells stretched 10 x in z) and 1.2e-3 (30 x): types 3 and 4
    diverge, type 2 converges (5e-4 and 0.16 per 20 cycles), legacy stagnates (0.86, 0.97).
    A variant that did not over-relax the constant of an open component slowed the cubic case
    to 13.4 cycles and still diverged at 30 x, so it was not kept. Types 3/4 are refused for an
    open component with a matrix measure below 0.1 (`PCO_P0_WEAK_OPEN`): fallback to legacy
    (key absent, logged), refusal (explicit key: use type 2), previous type kept (reload).
    The threshold lies between the measured working (0.90; unit ring 0.11) and failing
    (1.2e-2) cases; the cases between are not measured.
  - One reference per closed component, the mean mode of its lowest element, set to zero on a
    solver-owned copy (row and column cleared except the diagonal, right-hand-side entry
    zeroed: symmetric elimination of a zero value); a zero diagonal of a reference row (a
    closed component whose mean row vanishes) becomes the largest diagonal modulus of the copy,
    so the equation reads d x = 0. For the SSOR types 1/3 a zero diagonal elsewhere in the copy
    is refused at the refresh.
  - Later refreshes keep the decision (no flipping) and report drift of the matrix measure
    (closed above 1e-6, open at or below 1e-13). The decision does not use
    `SimPar@NoOutflow`; a disagreement is logged (`PCO gauge note`). The right-hand-side
    compatibility measure assumes the left null vector e (symmetric operator); a right-hand
    side whose sum over a closed component does not vanish is reported (component, net value,
    relative imbalance; above 0.1 relative: late-cycle defects of compatible single-phase runs
    reach 2e-2 at defects of 1e-14), not changed.
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
| `pco_periodic_topology` (QuadSc_pco_topology.f90) | structural flux of a constant pressure on the workers (`PCO_FREE_FLUX`); link census on the workers (shared vertices, `PCO_PAIR_SHIFT`); dof classes on the master (coordinate hash, all combined shifts, union-find, completeness / ambiguity / kind checks, axes realised); constraint conflicts (`PCO_MASKS_DIFFER`); element graph through vertex classes; master coordinate check |
| `pco_gather` (QuadSc_pco_gather.f90) | global element/unknown ids, gather plan with validation, value collection, immutability check |
| `pco_wrap` (QuadSc_pco_wrap.f90) | wrapped structure and assembly |
| `pco_coarse_solver` (QuadSc_pco_solver.f90) | master storage: layout, gauge decision, solver copy, UMFPACK / BiCGStab solves, raw SOR, solver type change, diagnostic solves |
| `pco_p0_reconstruct` (QuadSc_pco_reconstruct.f90) | periodic vertex classes and the periodic P0 -> vertex averaging |
| `pco_diagnostics` (QuadSc_pco_check.f90) | operator comparison, the check (worker product, closure, solver check, imbalance probe) |
| `pco_coarse_operator` (QuadSc_pco_operator.f90) | `Renew_Pressure_Operator`, `PCO_REGISTER_APPLICATION`, `PCO_PRESSURE_RELOAD_ADMIT`, `PCO_PRESSURE_PARAMETERS_RELOADED`; mode decision, default scope, validation and refusal, pattern check before the assembly, renewal, wrapped validation, logging |
| `UMFPackSolver` (UmfpackSolver.f90) | `tUmfCSC` adapter; legacy handles tracked one by one, the legacy solve stops without a live factorisation |
| `pressure_param_agreement` (src_util/pressure_param_agreement.f90) | collective agreement on the reloaded pressure settings (every mode) |
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
parameters (every mode, also legacy and non-periodic runs): the file must be visible on every
rank (allreduce in `ReloadParameterFile`, so every rank enters the reload or none); every rank
parses it, then the ranks agree on EVERY parsed pressure setting (`PRESSURE_PARAMETERS_AGREE`,
src_util/pressure_param_agreement.f90: MinLev, MedLev, Min/MaxIterCycle, SmoothSteps,
IterCoarse, CrsSolverType, the two keys, the cycle type, DefImprCoarse, CrsRelaxPrm,
Criterion1/2, RLX; allreduce min/max of the integers and of the reals, equal only when bitwise
equal); on any difference every rank keeps the previous values (logged), so no rank can run a
different number of smoothing steps or cycles, i.e. of exchanges. A consistent reload behaves
as before. Under a periodic operator a change of the master levels or a type outside 1-4 is
then undone before the generic validation (which would abort on it); a new type 1-4 rebuilds
the solver copy (same layout and gauge); the key and the check value are read once. The solves check the dimensions of
the installed copy.

Collective operations per periodic renewal: the pattern and extent check before the assembly
(one pressure exchange, one allreduce); `Create_ParCMat` (as always); the knprP allreduce; the
value gather (gathered mode, and in wrapped mode for its comparison: one `MPI_Gatherv`); in
wrapped mode the broadcast of the comparison result; the status broadcast (and the message
when it fails); when due, the check (one gather, three distribute/exchange/collect cycles). At
the first renewal in addition: the census (six pressure exchanges and two allreduces: the
structure flag; the consistency and link counts), the plan (one exchange, gathers,
broadcasts), the broadcasts of the master validation and the structural flux of the gauge
(one coarse vector collection). Per reload (every mode): one allreduce for the file and four
for the agreement of the parsed settings; under a periodic operator, for a new solver type,
one status broadcast.

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
| gauge undecidable at the first renewal (structural flux between 1e-10 and 1e-6, or structure and matrix inconsistent) | fallback, logged | legacy | refused |
| coarse solver type 3/4 with an open component whose matrix measure is below 0.1 | fallback, logged | legacy | refused (use type 2) |
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

Binaries, all built on longhouse1a from the clean tree at 2777989b (the state of review round 3):
`q2p1_el_pipeflow.2777989b` (EL recipe of build-perk-lh, md5 1dda0d02),
`q2p1_bench_sedimentation.2777989b` and `q2p1_dns_drag.2777989b` (the certified DNS twin recipe:
gcc 13.2 wrappers, serial PE, HardContactAndFluid, FBM acceleration), the same DNS targets of the
base 16ffcbcc (`*.base`, worktree FeatFloWer-pcoarse-base) and the round-3 DNS binary
`q2p1_dns_drag.1004a4de` (census defect). The tree builds with USE_PE=OFF (q2p1_fac3d,
q2p1_fc_ext, build-nope) and the q2p1_cc source compiles (its application does not link in this
tree, before this change). Runs: `/data/warehouse17/rmuenste/runs/pcoarse/round5` (EL; report.txt,
report_extra.txt) and `round5_dns` (report.txt, `tools/dns_report4.py`); every directory has
versions.txt. Probe metrics from the `XxX` lines: cycles per solve, geometric mean of
DefFinal/DefInit, `res1_L1` of the last `EL_VA_FLUX` line. Development probes of this round
(stretched boxes, the first forms of the structural measure, the over-relaxation variants) are in
`runs/pcoarse/dev5`; a first full gate run with an intermediate binary is kept in
round5a_5cab5b57 (it showed the over-relaxation variant slowing C_open, see Conventions).

Unit and MPI tests (job 163377, measured counts of `[ ok ]` records): `pco-coarse-serial` 115
(the serial checks of the earlier rounds and: operator comparison with NaN/Inf in either operand
and a small entry in a zero oracle row; the weak bridge; the zero reference row through full and
reduced solves; E013 link provenance on two cells per axis; the closure scale with noise; the
structural gauge - a weak outlet with matrix measure 2.5e-10 stays open, has no reference and is
solved exactly (type 2), the P0 types are refused for it, a well open ring (0.11) is accepted,
structure closed / matrix open, structure open / matrix in the numerical kernel and a structural
value between the thresholds are undecided without a reference; the zero-row worker-product
comparison; the per-coordinate tolerance of the pair check and of the master classes),
`pco-reload-mpi-3` 5 (3 ranks; test hook: the parsed parameter set perturbed on one rank -
smoothing steps, a real in the last bit, the cycle type, a NaN - restores the previous values on
every rank; equal values are kept). The job runs six ctest entries in all (252 records).

**A. Operator identity** (worker product with the E012_DAX arithmetic against the master
product, all rows - zero rows of the gathered operator compared exactly -, checks at renewals 1,
2, 10, 20, maxima over them)

| case | run (job) | dax - gathered | dax - wrapped | per-image mass variant | verdict |
|---|---|---|---|---|---|
| 12x12x120 (coarse 6x6x60), weighted EL | D_m12_gath_chk (163386), D_m12_wrap_chk (163388) | 6.8e-16 | 1.4e-14 | 4.1e-04 | PASS, 28 lines, 0 FAIL |
| 6x6x60 (coarse 3x3x30), weighted EL | D_s6_gath_chk (163391), D_s6_wrap_chk (163393) | 7.7e-16 | 9.0e-15 | 5.9e-05 | PASS |
| 12x12x120 single phase | D_sp_gath_chk (163396), D_sp_wrap_chk (163399) | 5.3e-16 | 1.4e-14 | = wrapped | PASS |
| coarse level 2, explicit | X_lev2_gath_chk (163466) | 4.9e-16 | 2.0e-14 | - | PASS |
| coarse level 3, explicit | X_lev3_gath_chk (163468) | 5.7e-16 | 4.8e-14 | - | PASS |
| DNS d11 (QBOX9) | d11__gathered (163392), d11__wrapped (163394) | 3.9e-16 | 2.3e-15 | - | PASS |
| DNS QBOX2 | q2__gathered (163401), q2__wrapped (163402) | 1.6e-16 | 8.0e-16 | - | PASS |

**B. Gathered vs wrapped entry by entry**

| case | run (job) | entries | max difference | verdict |
|---|---|---|---|---|
| xy periodic, z walls | D_m12_gath_chk (163386) | 922 752 | 1.1e-13 | PASS |
| z planes Symmetry001 | B_sym_gath_chk (163418), B_sym_wrap_chk (163419) | 922 752 | 1.1e-13 | PASS |
| 3x3x6 partition (55 ranks) | B_p336_gath_chk (163420), B_p336_wrap_chk (163421) | 114 048 | 6.6e-14 | PASS |
| triply periodic QBOX9 | B_qbox9_gath_chk (163413), B_qbox9_wrap_chk (163414) | 314 928 | 1.6e-14 | PASS |
| QBOX2 | q2__gathered (163401) | 1 024 | 1.2e-15 | PASS |
| z-periodic pipe, s05 E1 (curved wall) | B_pipe_gath_chk (163415) | 482 112 | 1.8e-01 | wrapped FAIL (dax - gathered 5.5e-16 PASS); wrapped refused: B_pipe_wrap_chk (163416, exit 1) |

**C. Gauge and solvers** (structural flux s and matrix measure, as logged in `PCO gauge`)

| item | run (job) | result |
|---|---|---|
| closed periodic cases | all `_chk` runs, D_m12_gath | one component, closed: structural flux 1.1e-15 to 1.4e-14, matrix 1.3e-15 to 2.0e-14 (open: 0.50; four orders or more from the structure thresholds 1e-10 / 1e-6); ungauged residual below 1e-15, reference 0 |
| open periodic case (zmax Outflow), cubic cells | C_open_gath_chk (163412) | structural flux 0.50, matrix 0.90: open, no reference; 2.00 cycles, checks PASS |
| open, cells stretched 10 x in z, type 4 | S_open_z10_t4 (163469) | structure open (0.50), matrix 1.2e-2: type 4 refused for the weakly open component, fallback to legacy logged (then as legacy: 0.86 per 20 cycles) |
| open, stretched 10 x, type 2 | S_open_z10_t2 (163470) | gathered: open, no reference, compatible solve PASS; converges (5e-4 per 20 cycles; legacy S_open_z10_legacy 163488: 0.86) |
| open, stretched 30 x, types 4 / 2 | S_open_z30_t4 (163478), S_open_z30_t2 (163486) | type 4 falls back to legacy (matrix 1.2e-3, logged); type 2 gathered: open, no reference, solve PASS, 0.16 per 20 cycles |
| closed, stretched 10 x, type 4 | S_closed_z10_t4 (163487) | closed (structure 5.1e-15, matrix 8.6e-15), one reference; 5.1e-3 per 20 cycles (legacy at this stretch, dev5/X_sp_z10_closed4leg: 0.81) |
| open, stretched 1e4 x (DNS d11, z images unmatched: free z faces), explicit gathered, type 4 | d11z__gathered (163406), d11z__legacy (163408) | structure open (0.50), matrix 3.2e-11 - below the former matrix-only threshold 1e-8, which would have pinned it: no reference; type 4 with this weakly open component is refused (`PCO REFUSED ... use type 2`). The legacy run of the same deck diverges in step 1 (the aspect ratio, not the coarse operator). In the development probe dev5/U_d11_z10000 (before the P0 refusal) the same decision was taken and the compatible coarse solve had the ungauged residual 3e-16 |
| type 1 / 2 / 3 / 4 (m12) | 163405 / 163407 / 163409 / 163379 | 2.05 / 2.00 / 2.00 / 2.00 cycles |

**D. End-to-end probes**

| case | legacy (job) | gathered (job) | wrapped (job) |
|---|---|---|---|
| PR15_m12_probe, weighted EL | 20.00 cycles, 0.539, bitwise vs 162042 (163384) | 2.00, 1.81e-05, 3.3442e-11 (163379) | identical XxX lines (163382) |
| PK_smoke_wrap2, 6x6x60 | 20.00, 0.411, prefix-bitwise (163395) | 2.05, 6.28e-05 (163391) | same (163393) |
| MG_sp single phase | 20.00, 0.822, bitwise (163403) | 2.00, 4.97e-06 (163396) | 2.00, 6.46e-06 (163399) |
| MG_crs2 (type 2) | 20.00, 0.541, bitwise (163411) | see C | see C |
| coarse level 2 | key absent = legacy: 20.00, 0.462 (163467) | explicit: 2.00, 1.27e-04 (163466) | - |
| DNS d11, 40 steps | 9.25 cycles (163390) | 3.33 (163392); drag within 6.8e-7, bulk flow 3.0e-8 of legacy | 3.33 (163394) |
| DNS QBOX2, 10 steps | 3.60 (163400) | 3.60, printed drag and XxX identical (163401) | same (163402) |

**P. Census on two coarse cells per periodic axis** (QBOX2, mesh by `tools/make_box_mesh.py`,
which reproduces QBOX9 byte for byte): round-3 binary, explicit gathered: periodicity inactive,
legacy (q2__r3gathered, 163404); 2777989b: axes TTT, 200 accepted periodic dof links per axis,
identity PASS (163401); key absent and legacy bitwise vs base (163398, 163400; base 163397).

**R. Runtime reload** (`_data/ProcCtrl.txt`: `Reload_Pres = _data/reload_pres.dat`, executed
after step 1)

| case | runs (jobs) | result |
|---|---|---|
| walled A15 deck (legacy, not periodic), the same values | R_A15_same__new (163489) vs __base (163490) | every stream bitwise: a consistent reload behaves as before |
| walled A15 deck, MGSmoothSteps 16 -> 17 | R_A15_change__new (163491) vs __base (163493) | the reload is taken in both (protocol: 17), every stream bitwise |
| periodic gathered (m12), coarse solver type 4 -> 2 | R_m12_type (163494) | `coarse solver type changed by reload from 4 to 2`, 2.00 cycles continue |
| periodic gathered (m12), MinLev/MedLev 1 -> 2 | R_m12_minlev (163495) | not taken (logged), previous values, 2.00 cycles continue |
| ranks parsing different values | pco-reload-mpi-3 (163377) | test hook (a rank-dependent file cannot be staged on the shared file system): previous values on every rank |

**E. Bitwise**

| case | key | runs (jobs) | result |
|---|---|---|---|
| 4 non-periodic EL gates (SCOL4, VFT10 x 3) | absent / gathered / wrapped | 163428-163439 | 12/12 bitwise |
| 6 periodic QBOX9 EL gates | legacy | 163440, 163442, 163444, 163446, 163448, 163450 | 6/6 bitwise |
| 6 periodic QBOX9 EL gates | absent (= gathered) | 163441, ..., 163451 | not bitwise (expected); 2.03-2.28 cycles |
| A15 wall-bounded batch | absent / gathered / wrapped | 163452-163454 | prefix-bitwise vs A15, bitwise vs base |
| P15_probe_cross twin | legacy | 163455 | bitwise |
| s05 pipe | legacy | 163417 | bitwise vs base |
| DNS e4_l3 twin | absent / legacy / gathered / wrapped vs base | 163378 (base), 163380, 163381, 163383, 163385 | 4/4 bitwise on every stream |
| DNS d11 | absent / legacy vs base | 163389, 163390 (base 163387) | bitwise |
| DNS QBOX2 | absent / legacy vs base | 163398, 163400 (base 163397) | bitwise |

**F. Refusals and inactive periodicity**

| case | run (job) | result |
|---|---|---|
| walled A15 box (coarse cell 0.0025), unit periods | F_wall_per1 (163422), F_wall_per1_gath (163423) | periodicity inactive, logged; bitwise vs base |
| walled A15 box, periods 0.008333 on x, y (3.33 coarse cells: cell width 0.3 period) | F_walled06__{unset,gathered,wrapped} (163457-163459), base 163460 | periodicity inactive (no accepted periodic dof link), bitwise vs base under the three key settings |
| walled A15 box, periods 0.0041667 on x, y (cell width 0.6 period) | F_wallc06__{unset,gathered} (163496, 163497), base 163498 | periodicity inactive, bitwise vs base |
| walled A15 box, unit periods, MinLev 1 < MedLev 2 | F_walledML__{unset,legacy,gathered,wrapped} (163461-163464), base 163465 | periodicity inactive, bitwise vs base under every key |
| x pair Wall / Periodic | F_conflict_default (163424), F_conflict_gath (163425) | fallback logged / `PCO REFUSED`, exit 1 |
| periodic, MinLev 1 < MedLev 2 | F_minmed_default (163426), F_minmed_gath (163427) | legacy logged / refused |
| wrapped on a curved boundary | B_pipe_wrap_chk (163416) | refused |
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
replaces the one after it (same cost). Round 5 (a loaded cluster): legacy 100.6 s, gathered
52.7 s, wrapped 53.4 s (163384, 163379, 163382); the structural flux is evaluated once. The worker-side assembly of C is the same in all three (1.16 s).

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
- Reload_Velo parses rank-locally like Reload_Pres did (no agreement); outside this branch.
- A parse error of a reloaded pressure file (an invalid `Pres@PeriodicCoarseOperator`) stops
  every rank that reads it; a file that differs between ranks in that way would stop some ranks
  only (the parser has no error status to agree on).
