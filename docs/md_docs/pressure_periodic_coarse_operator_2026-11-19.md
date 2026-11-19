# Periodic coarse pressure operator (2026-11-19)

Design note for `Pres@PeriodicCoarseOperator = legacy | gathered | wrapped` and
`Pres@PeriodicCoarseCheck = n` (branch `feature/pressure-periodic-coarse`, worktree
`FeatFloWer-pcoarse`, base `feature/el-periodic-kernel` 16ffcbcc). Gate and probe runs:
`/data/warehouse17/rmuenste/runs/pcoarse/round2` (final binary), round 1 in the parent
directory. The datasheet and the s10 RUNBOOK are not edited here.

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

**gathered (default).** What the workers apply on the coarse level is, for the rows of their
own elements, (C x)_r + (CP xp)_r with xp the ghost values of the pressure exchange
`GetParPressure`. Both matrices are assembled on the workers with the actual coefficients
(weighted EL mass, slip projection, constraint masks, curved-boundary geometry). The master
collects them into its numbering: the workers exchange the global numbers of their own
unknowns with `GetParPressure` (as the MUMPS exporter `MUMPS_SetUpP1_SLAVE` does) to label the
ghost columns, the element map is the one of the coarse vector transfers (`E012DISTR_L1/2/3`:
coarse element link, then child 1 keeps the parent number, children 2..8 of parent p of a
level with N parents follow at N + 7(p-1) + 1..7). A plan (`PCO_GATHER_BUILD`) built once per
pattern maps every entry of the workers' C and CP storage to a slot of the master layout;
every renewal the workers send their value arrays as they are (`MPI_Gatherv`) and the master
sums them into the slots (`PCO_GATHER_VALUES`). Exact by construction; the check below
confirms it against the worker product to 1e-15.

**wrapped.** The master assembly with the velocity dofs of each periodic class identified:
C = (B^T Q) M_cls^-1 (Q^T B), Q copying a class value to its images. The B rows and B^T columns
of the images are merged per (element, mode) before the product (no first-match loss of
terms), all images count (edges 4, triply periodic corners 8), the class mask is the common
Dirichlet mask (equal over the class, checked), the class mass is the sum of the raw image
masses, under the weighted EL operator U and W summed separately and transformed once,
m = (sum U)^2 / (sum W), as the workers do (`QuadSc_el_mm.f90`). It is the master geometry:
on a curved boundary the master does not project the coordinates of its finest level (the Q2
dofs of the coarse level; `ParametrizeBndr` of level MedLev+1 runs on the workers only,
app_init), so its B and lumped mass differ from the workers' there (s05 pipe: 18 % of an
entry; the legacy master matrix has the same values as the wrapped one at that entry).
Therefore wrapped is compared with the gathered operator at its first renewal and refused
above 1e-10 (`validate_wrapped`).

**Why gathered is the default.** It is the operator the fine levels see, whatever the
coefficients, the boundary treatment or the geometry; it needs no second implementation of the
assembly semantics that would have to follow every change of B, the masses or the boundary
conditions. Its cost is one gather of the coarse values per renewal (2.8 ms on the 12 x 12 x
120 probe, master layout 922 752 entries). The wrapped assembly is the alternative for
configurations whose master geometry and coefficients equal the workers' (planar periodic
boxes, ordinary and weighted EL mass, symmetry planes): there it agrees with the gathered
operator to 1e-13 entry by entry and the multigrid histories agree to the printed digits.

## Conventions

- Numbering: unknown 4(e-1)+m is mode m of the P1 pressure of master element e on the coarse
  level (m = 1 the element mean, 2..4 the reference-cube slopes of E012).
- Layout (`PCO_INSTALL_LAYOUT`): the block row of e lists the own block first, then the blocks
  of the neighbour elements ascending; row 4(e-1)+m starts with its diagonal, then the other
  modes of e ascending, then four entries per neighbour (the legacy `Get_CMatStruct` layout, so
  the mean-mode system is read with stride 4). Neighbours: the element graph through vertex
  classes (elements sharing a vertex or one of its periodic images).
- Raw operator (`pco%a`): ungauged, never modified by a solve; it serves the type-3/4 smoothing
  (`PCO_RAW_SOR`) and the checks.
- Gauge: on every connected component whose operator annihilates the constant (closure measure
  max_i |(C e)_i| / sum_j |C_ij| e_j below 1e-8, e = 1 on the mean modes) the mean mode of the
  lowest element is the reference: on a solver-owned copy its row and column are cleared except
  the diagonal and its right-hand-side entry is zeroed (symmetric elimination of a zero value).
  The decision is taken from the operator at every refresh, not from `SimPar@NoOutflow`; a
  disagreement is logged (`PCO gauge note`). An open component gets no reference. A right-hand
  side whose sum over a closed component does not vanish is reported (`PCO_COARSE_RHS`, above
  1e-2 relative; the largest value on the `PCO_TIMING` lines), not changed. The compatibility
  measure assumes the left null vector e (symmetric operator); the check prints e^T A.
- Solver copy: full system (types 1, 2) or mean modes (types 3, 4), pattern fixed (explicit
  zeros kept). Types 2/4: UMFPACK through `tUmfCSC`, the compressed-column storage of the copy
  (UMFPACK solves A x = b), symbolic analysis once per pattern, numeric factorisation at every
  refresh, every status checked (factorisation: collective stop with the reason; solve, on the
  master inside the multigrid: message and `MPI_Abort`). Types 1/3: the master BiCGStab with
  the SSOR preconditioner of the legacy path on the gauged copy (work vectors zeroed before
  each call; a NaN aborts with a message).
- Types 3/4 (`crs_P0_Correction`, shared by both master coarse solves): P0 solution x3, vertex
  averaging, `IntQ1toP1` on the physical elements, SOR smoothing. With a periodic operator the
  averaging sums numerators and denominators over each periodic vertex class, divides once and
  gives the value to every image (`PCO_P0_TO_VERTEX`), and the smoother uses the raw operator;
  the legacy branch is `INTPVB` and `crsSmoother` as before.

## Code

| module (file) | responsibility |
|---|---|
| `pressure_coarse_key` (src_util/pressure_coarse_key.f90) | key values, names, parsing |
| `pco_periodic_topology` (QuadSc_pco_topology.f90) | link census on the workers; dof classes on the master (coordinate hash, all combined shifts, union-find, completeness / ambiguity / kind checks); constraint conflicts; element graph through vertex classes |
| `pco_gather` (QuadSc_pco_gather.f90) | global element/unknown ids, gather plan with validation, value collection |
| `pco_wrap` (QuadSc_pco_wrap.f90) | wrapped structure (elements and block positions per class, B^T column index) and assembly |
| `pco_coarse_solver` (QuadSc_pco_solver.f90) | master storage: layout, components, gauge, solver copy, UMFPACK / BiCGStab solves, raw SOR, solver check |
| `pco_p0_reconstruct` (QuadSc_pco_reconstruct.f90) | periodic vertex classes and the periodic P0 -> vertex averaging |
| `pco_coarse_operator` (QuadSc_pco_operator.f90) | `Renew_Pressure_Operator` (the one renewal entry point), mode decision, validation and refusal, renewal, wrapped validation, check, logging |
| `UMFPackSolver` (UmfpackSolver.f90) | `tUmfCSC` adapter; legacy release only of live factors |

Call sites: `OperatorRegenaration` (C branch), `EL_MM_Pressure_Operator` (every step, hands
the raw coarse masses U = rho Ml and W to the renewal), `q2p1_cc` (`OperatorRegenaration_iso`).
`Create_CMat` only assembles; the release of the previous legacy UMFPACK factorisation moved
from the EL adapter into the legacy branch of the renewal (flag-protected, so a release
without a live factorisation is a no-op). In the multigrid (`QuadSc_mg.f90`) one branch in
`mgCoarseGridSolver_P` (`PCO_ACTIVE`) and the shared type-3/4 tail.

Cached: classes, element graph, layout, gather plan, wrap structure, symbolic analysis.
Refreshed per renewal: values (gather or wrap), closure, gauge, numeric factors. A changed
worker pattern (sizes of C/CP) rebuilds the plan (collective test every renewal), and a plan
that no longer fits the layout stops the run.

Extension points: a further operator is a key value and a branch in `refresh_raw`; a master
hierarchy (MinLev < MedLev) needs an operator on every master level and a tested master
cycle; another coarse solver is a branch in `solve_copy` on the existing copy.

## Deck keys and activation

`Pres@PeriodicCoarseOperator` (absent, `legacy`, `gathered`, `wrapped`; invalid values stop
every rank in the parser) and `Pres@PeriodicCoarseCheck = n` (0 off; n > 0: the check at
renewals 1, 2 and every n-th). The key is read once, at the first renewal (a later
`Reload_SimPar` does not switch the operator).

| situation | key absent | `legacy` | `gathered` / `wrapped` |
|---|---|---|---|
| no finite period on any axis | legacy, nothing entered | legacy | legacy, one log line |
| finite periods, no worker coupling crosses a periodic plane (census) | legacy, logged | legacy, logged | legacy, logged (periodicity inactive, e.g. q2p1_dns_drag unit periods on a walled box) |
| validated periodicity | gathered | legacy (logged) | the requested one |
| coarse solver type outside 1-4, MinLev /= MedLev, coarse level > 3, master not solving | fallback to legacy, logged | legacy | refused |
| couplings inconsistent with the periods (offset off the lattice by more than two element sizes) | fallback, logged | legacy | refused |
| pressure Dirichlet elements (knprP) on the coarse level | fallback, logged | legacy | refused |
| ambiguous or incomplete classes, mixed dof kinds, unequal Dirichlet flags over a class | fallback, logged (counts, first dof and its coordinates) | legacy | refused |
| a worker coupling outside the master layout, a master row with two or no owners | fallback, logged | legacy | refused |
| wrapped differs from gathered at the first renewal by more than 1e-10 | - | - | wrapped refused |

Refusal: every rank reaches the same decision (counts broadcast from the master), the master
and the protocol rank print `PCO REFUSED: ...`, every rank stops. Output: `PCO periodic coarse
operator = ...` (mode, level, active axes, shifted links per axis), `PCO classes ...`,
`PCO gauge ...` on the terminal and the protocol; `PCO_TIMING` and `PCO_CHECK` lines on the
terminal of the master.

The check (`run_check`): both operators; entry-by-entry difference gathered - wrapped (scaled
by the largest modulus of the row, with the location and the legacy master value); for three
vectors (mean modes = 1; a four-mode impulse on the element at the largest vertex class; a
reproducible vector exciting all modes) the worker product with the E012_DAX arithmetic
(distributed and collected with the coarse vector transfers) against the master products of
the gathered, the wrapped and, under the weighted EL mass, the per-image transformed wrap
(the prototype's seam mass, to show that the check discriminates), every row, scaled by
sum_j |A_ij| max|v|; closure, left null vector, a compatible solve with the ungauged
residual, the reference value, an incompatible right-hand side that must be reported.

## Tests and gates

Final binary `q2p1_el_pipeflow.0e9d2575` (built on longhouse1a with the recipe of
build-perk-lh; md5 5935f5ec, identical to a rebuild at 793f4f2b, which changes comments
only). Run directories under `/data/warehouse17/rmuenste/runs/pcoarse/round2` unless noted;
every directory has versions.txt.

Unit tests (`test_pco_coarse`, ctest `pco-coarse-serial`, Slurm job 162766 on nx; the
serial parts without MPI topology) and application probes (Slurm, `Pres@PeriodicCoarseCheck`
where named `_chk`). Probe metrics from the `XxX` lines (40 pressure solves in 20 steps for
the EL probes): cycles per solve, geometric mean of DefFinal/DefInit, post-corrector
continuity residual `res1_L1` of the last `EL_VA_FLUX` line.

**A. Operator identity** (worker product with the E012_DAX arithmetic against the master
product, all rows, normwise row scale; checks at renewals 1, 2, 10, 20, i.e. after the
weights changed; weighted EL under the bed porosity, which is unequal across the seams)

| case | run (job) | dax - gathered | dax - wrapped | per-image mass variant | verdict |
|---|---|---|---|---|---|
| 12x12x120 (coarse 6x6x60), 3x3x3, weighted EL | D_m12_gath_chk (162770), D_m12_wrap_chk (162771) | 6.8e-16 | 1.4e-14 | 4.1e-04 (impulse 9.0e-06) | PASS, 28 check lines, 0 FAIL |
| 6x6x60 (coarse 3x3x30), 3x3x3, weighted EL | D_s6_gath_chk (162772), D_s6_wrap_chk (162773) | 7.7e-16 | 9.0e-15 | 5.3e-04 | PASS |
| 12x12x120, single phase (ordinary mass) | D_sp_gath_chk (162775), D_sp_wrap_chk (162776) | 5.3e-16 | 1.4e-14 | = wrapped | PASS |

The per-image variant (seam mass sum U_i^2/W_i of the prototype) passes the constant vector
and fails the others by 1e-5 to 5e-4: the identity test discriminates the seam mass, the
constant (row-sum) test does not. At renewal 1 (ordinary start-up mass) both variants agree.

**B. Gathered vs wrapped entry by entry** (difference / largest modulus of the row)

| case | run (job) | entries | max difference | verdict |
|---|---|---|---|---|
| xy periodic, z walls, 3x3x3 | D_m12_gath_chk (162770) | 922 752 | 1.1e-13 | PASS |
| same, z planes Symmetry001 (periodic + symmetry intersection) | B_sym_gath_chk (162789) | 922 752 | 1.1e-13 | PASS |
| xy periodic, 3x3x6 partition (55 ranks) | B_p336_gath_chk (162791) | 114 048 | 6.9e-14 | PASS |
| fully periodic QBOX9 3x3x3 (corner class of 8) | B_qbox9_gath_chk (162784) | 314 928 | 1.6e-14 | PASS |
| z-periodic pipe, s05 E1 mesh (curved wall) | B_pipe_gath_chk (162786) | 482 112 | 1.8e-01 | FAIL (wrapped): the largest difference is at an element next to the curved wall, not at the seam, and the legacy master matrix has the wrapped value there (7.0763e-7 vs gathered 6.0095e-7): master geometry. dax - gathered 5.5e-16. Wrapped is refused on this mesh: B_pipe_wrap_chk (162787), `PCO REFUSED ... differs ... by 1.775E-01` |

Classes found: 2783 (largest 4) on the 6x6x60 coarse box, 671 (4) on 3x3x30, 919 (8) on
QBOX9, 209 (2) on the pipe; no worker coupling outside the class layout, no uncovered slot.

**C. Gauge and solvers**

| item | run (job) | result |
|---|---|---|
| unit: orientation | test_pco_coarse | tUmfCSC solves A x = b (1e-13), the legacy wrapper solves A^T x = b, singular matrix status 1 reported, double legacy release safe: PASS |
| unit: classes, graph, averaging | test_pco_coarse | 2/4/8 images for x/xy/xyz, refusals for a moved corner (incomplete) and a duplicated dof (ambiguous), conflicts; neighbours 17/26; averaging equal on all images and = sum(num)/sum(den), INTPVB differs between images: PASS |
| unit: gauge, types 1-4 | test_pco_coarse | closed ring: C e = 0 detected, ungauged residual 3e-16 (direct) / 1e-9 (BiCGStab), reference 0, incompatible rhs reported; open ring: no reference, exact solve: PASS |
| closed periodic case: C e, e^T A, compatible solve | D_m12_gath_chk and all `_chk` | C e 8.9e-15, e^T A 8.8e-15, ungauged residual 5e-16, reference value 0, incompatible rhs reported: PASS |
| open periodic case (zmax Outflow) | C_open_gath_chk (162783) | closure 0.9: no reference; `PCO gauge note` (NoOutflow forced by PeriodicAxis disagrees); 2.00 cycles, 0 FAIL |
| type 1 (BiCGStab, full) | C_m12_gath_t1 (162778) | 2.05 cycles, 3.66e-05, res1 3.3444e-11, checks PASS |
| type 2 (UMFPACK, full) | C_m12_gath_t2 (162779), C_m12_wrap_t2 (162781) | 2.00 cycles, 1.38e-05, res1 3.3442e-11, gathered and wrapped identical XxX lines |
| type 3 (BiCGStab, P0) | C_m12_gath_t3 (162780) | 2.00 cycles, 1.82e-05 |
| type 4 (UMFPACK, P0) | D_m12_gath (162767) | 2.00 cycles, 1.81e-05; type 2 vs 4: mean slip of the bed equal to 2e-6 relative |

**D. End-to-end probes** (reference values of the prototype: 2.0 cycles, 1.8e-5)

| case | legacy (job) | gathered (job) | wrapped (job) |
|---|---|---|---|
| PR15_m12_probe, 12x12x120 weighted EL | 20.00 cycles (40/40 at cap), 0.539, res1 9.81e-08, bitwise vs 162042 (162769) | 2.00, 1.81e-05, 3.3442e-11 (162767) | 2.00, 1.81e-05, 3.3442e-11, XxX lines identical to gathered (162768) |
| PK_smoke_wrap2, 6x6x60 weighted EL, 20 steps | 20.00 (40 at cap), 0.411, 1.53e-09, prefix-bitwise vs the 1000-step log (162774) | 2.05, 6.28e-05, 3.97e-12 (162772) | same as gathered (162773) |
| MG_sp, 12x12x120 single phase | 20.00, 0.822, bitwise vs MG_sp (162777) | 2.00, 4.97e-06 (162775) | 2.00, 6.46e-06 (162776) |
| MG_crs2 (type 2) | 20.00, 0.541, bitwise vs MG_crs2 (162782) | see C | see C |

Gathered against wrapped in the EL audit lines: bitwise or rounding-level quantities only
(conservation residuals of order 1e-15 to 1e-21). The check does not perturb a run:
D_m12_gath and D_m12_gath_chk are bitwise identical.

**E. Bitwise** (`el_stream_compare.py --strict` on every EL stream, the `XxX` lines as text,
every particle log byte for byte; references: PK_gates_lh2 (jobs 161588-161596, 161599, bitwise
against coupling_p4p5/int_gates), batch10 logs)

| case | key | runs (jobs) | result |
|---|---|---|---|
| 4 non-periodic gates (SCOL4, VFT10 x 3) | absent / gathered / wrapped | E_gates/<gate>__{unset,gathered,wrapped} (162799-162810) | 12/12 bitwise (also in round 1, binary 11f1e044); with the key given, one log line "requested, but no axis is periodic" |
| 6 periodic QBOX9 gates | legacy | E_gates/<gate>__legacy (162811, 162813, 162815, 162817, 162819, 162821) | 6/6 bitwise |
| 6 periodic QBOX9 gates | absent (= gathered) | E_gates/<gate>__unset (162812, ..., 162822) | not bitwise (expected: the coarse operator changes); cycles per solve 3.15-5.57 (legacy logs) -> 2.03-2.28, geometric mean DefFinal/DefInit 2.7e-5-1.4e-4 -> 4.6e-6-9.7e-6 |
| A15 wall-bounded batch, 20 steps | absent / gathered / wrapped | E_A15__* (162823-162825) | prefix-bitwise vs the 12000-step A15 log (time averages excluded), bitwise vs the 20-step run of the base binary (F_wall_base, round 1, job 162743) |
| P15_probe_cross periodic twin | legacy | E_cross_off__legacy (162826) | bitwise vs PK_probe_cross_off2 |
| PR15_m12_probe, MG_sp, MG_crs2, PK_smoke_wrap2 | legacy | see D | bitwise / prefix-bitwise |
| s05 E1 smoke (z-periodic pipe) | legacy | B_pipe_legacy (162788) | bitwise vs the base binary (B_pipe_base, round 1, job 162735); the E1_smoke log itself is from an older tree (FeatFloWer-va, 369d512b) and is not reproduced by the base binary either |
| baseline of the build | - | B0_base_PR15 (round 1, job 162697, base 16ffcbcc built here) | bitwise vs PR15_m12_probe |

**F. Refusals and inactive periodicity**

| case | run (job) | result |
|---|---|---|
| walled A15 box with unit periods (`SimPar@Periodic = Yes`), key absent / gathered | F_wall_per1 (162793), F_wall_per1_gath (162794) | `PCO periods declared ... no worker coupling crosses a periodic plane ... legacy`; bitwise vs the base binary |
| x pair with Wall on one side, Periodic on the other (inconsistent pair) | F_conflict_default (162795), F_conflict_gath (162796) | default: fallback to legacy, logged (1428 constraint conflicts, first dof with coordinates); explicit gathered: `PCO REFUSED`, every rank stops (exit 1) |
| MinLev 1 < MedLev 2 | F_minmed_default (162797), F_minmed_gath (162798) | default: fallback, logged; explicit: refused |
| wrapped on a curved boundary | B_pipe_wrap_chk (162787) | refused (see B) |
| explicit key with coarse solver type 5, 7-10 | - | not reachable in this build (ValidateSolverTypes stops first) |

## Cost (12 x 12 x 120 probe, 28 ranks, 20 steps, weighted EL, type 4)

Master side, per renewal (the weighted EL operator renews every step; `PCO_TIMING` lines,
seconds summed over 20 renewals): gather of the values 0.056 s (2.8 ms each), wrapped
assembly 0.43 s (21.6 ms each), refresh with the numeric factorisation of the P0 copy 0.146 s
(7.3 ms each); type 2 (full 8640-unknown copy, 922 752 entries) 4.35 s (218 ms each); one-time
set-up (census, classes, layout, plan) 0.071 s (gathered) and 0.077 s (wrapped, including the
start-up comparison with the gathered operator). The gather of 922 752 values is about 7.4 MB
per renewal into the master.

Whole run (overall time of the protocol rank, 19 steps after start-up): legacy 60.9 s
(pressure multigrid 32.1 s), gathered 32.1 s (3.64 s), wrapped 32.6 s (3.80 s); 3.20, 1.69 and
1.72 s per step (jobs 162769, 162767, 162768, all on nx nodes shared with other gate jobs;
indicative). The worker-side assembly of C is the same in all three (1.16 s).

## Not covered

- Master multigrid hierarchy (MinLev < MedLev): refused / fallback. The master F cycle
  (`crs_F_cycle`, loop from MedLev+1 to MedLev-1) is a V cycle, the W recursion stops at a
  hard-coded level 2, the transfers exist for levels 1-3 only.
- Coarse solver types 5 (MUMPS, refused for periodic runs by `ValidateSolverTypes`), 7, 8, 10
  (HYPRE) and 9 (PARDISO): not built in this tree; they stop in `ValidateSolverTypes` before
  the renewal, so the refusal of an explicit key with these types could not be exercised.
- DNS applications (q2p1_dns_drag, fac3d, ...) were not built or run. They renew through
  `OperatorRegenaration` and get the default; with FBM pressure Dirichlet elements (knprP)
  the default falls back to legacy (logged).
- Rotational or non-translational periodicity; meshes whose two periodic faces do not match
  node by node within DEpsPrec (refused as incomplete classes or couplings outside the layout).
- The legacy path keeps its defects (transposed UMFPACK solve, row gauge without the
  right-hand side, one-sided type-3/4 averaging, master geometry on curved boundaries).
- q2p1_cc: the changed call site compiles; the application itself does not link in this tree
  (`app_init.f90`, `PartitionReader.f90`: `ctemp` without type, before this change).

## Follow-ups

- Master hierarchy: operators on every master level (gathered per level is direct), a correct
  master F/W cycle, a test with MinLev < MedLev.
- HYPRE/MUMPS/PARDISO coarse paths with a periodic operator.
- Default for DNS apps and re-certification of the periodic certified decks under the
  default (the six QBOX9 EL gates change their iteration history under the default).
- The master geometry on curved boundaries (project level MedLev+1 on the master): changes the
  legacy arithmetic; would make wrapped usable on curved periodic meshes.
- Legacy type-3/4 tail: the transferred P1 slopes are in the master element orientation;
  the partitions checked here (BAT27P, PIPEZ27) have the same local vertex order on master and
  workers, which is what keeps that consistent.
