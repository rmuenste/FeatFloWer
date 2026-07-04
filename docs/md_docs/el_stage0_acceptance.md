# Euler-Lagrange Stage-0 Acceptance Battery

Date: 2026-07-04
Branch: `feature/euler-lagrange-phase1`
Head used by harness: `6908d4f8264078e0aacd69218c56f384cc6a1490`

## Tier-2 Harness Sweep

All six enabled E-L tier2 harness definitions passed with the in-tree
`tools/featflower_test` package on `PYTHONPATH`.

| Case | Run ID | Metric | Measured | Tolerance | Verdict |
| --- | --- | --- | ---: | ---: | --- |
| `q2p1_el_pipeflow_tier2_terminal` | `20260704-101520-914763e3` | `vz`, `speed` | `-3.556384e-02`, `3.556384e-02` | `2.0e-3` | PASS |
| `q2p1_el_pipeflow_tier2_straddling` | `20260704-102118-3d1c1be3` | `rel_error`, `residual` | `0.0`, `5.170037e-26` | `1.0e-10` | PASS |
| `q2p1_el_pipeflow_tier2_momentum` | `20260704-102157-e4eeb118` | `drift_rel`, `residual` | `4.290711e-07`, `1.694066e-21` | `1.0e-5`, `1.0e-10` | PASS |
| `q2p1_el_pipeflow_tier2_momentum_semi` | `20260704-102400-85702e8f` | `drift_rel`, `residual` | `1.510406e-03`, `1.185846e-20` | `3.1e-3`, `1.0e-10` | PASS |
| `q2p1_el_pipeflow_tier2_momentum_long` | `20260704-104433-19bc5f15` | `drift_rel`, `residual` | `1.000013e-05`, `1.229441e-27` | `2.0e-5`, `1.0e-10` | PASS |
| `q2p1_el_pipeflow_tier2_saffman` | `20260704-121702-8ff39325` | `drag_x`, `lift_z` | `2.356194e-05`, `1.276770e-06` | `1.0e-9` | PASS |

Notes:

- The installed `featflower-test` entry point imported a different checkout.
  The sweep above used `PYTHONPATH=$PWD/tools/featflower_test/src`.
- The long case was run through Slurm. The previous `constraint: nx` was
  unusable on this cluster state because all `nx` nodes were allocated or
  drained; removing that stale constraint let the job run on `worldgames`.

## Random Seeding Reproducibility

The PE setup now prints `EL_SEED_POS id= x= y= z=` at startup. Position lists
were sorted by global seed ID and compared byte-for-byte.

| Seed | Comparison | Count | Achieved phi | Diff | Verdict |
| --- | --- | ---: | ---: | --- | --- |
| `12345` | PE `1x1x3` vs `3x3x3` | `1492` | `0.049997399884330371` | clean | PASS |
| `99991` | PE `1x1x3` vs `3x3x3` | `1492` | `0.049997399884330371` | clean | PASS |

The two seeds produce different sorted position lists, so the comparison is not
passing accidentally on a seed-independent ordering.

True one-PE-rank execution could not be run from the staged QBOX9 case:
`mpirun -np 2` aborts in FeatFloWer MPI mesh setup before PE initialization with
`Attempt to DEALLOCATE unallocated 'mg_mpi'`, even when `SubMeshNumber=1`.
The reproducibility gate was therefore run between the smallest valid staged
decomposition, PE `1x1x3`, and the production PE `3x3x3` layout.

## File-Mode Regression

The single-entry `particles.xyz` path was rerun from the unmodified terminal
case at the tier2 rank count. Extracted `EL_TERMINAL_VEL` series were compared
against the Part-A terminal harness run.

| Comparison | Lines | Diff | Verdict |
| --- | ---: | --- | --- |
| Part-A terminal run vs file-mode rerun | `120` vs `120` | clean | PASS |

