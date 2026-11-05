# Extracts from dns-validation-campaign-plan.md (FeatFloWer repo root, branch feature/dns-validation, state 2026-10-02)

Verbatim line-range extracts; line numbers refer to the source file. Only the ten Cate-relevant parts are reproduced.

## Lines 170-195: §D0.4 ten Cate truth chain, §D0.5 Reproducibility floor

```
### D0.4 ten Cate truth chain (recertification)
- ~~Digitize the ten Cate PIV trajectories~~ **DONE 2026-07-31**: digitized
  E1–E4 curves provided by the user and committed as `tc-ref/`
  (`case_E*_h.csv` = t vs gap h/d; `ref_E*.dat` = t vs v_z [m/s]; format
  and sanity checks in `tc-ref/README.md`). This repairs the dangling
  references in
  `pipemesh_v1/handoff_euler_lagrange_drag_validation_ten_cate.md`.
  Sharpened gate: the confined-PIV E4 peak is −0.1230 m/s, so the
  documented serial DNS −0.1329 is **+8.1% vs PIV** (not +3.8% vs the
  unbounded u_∞ = 0.128) — the honest baseline the D1 resolution matrix
  must explain or close.
- Commit the mesh fixtures this campaign owns (`git add -f`):
  `benchSym/` quarter-box + `benchSym/mesh12/NEWFAC` (1×1×12 Cartesian),
  and the full-box `ten_cate_mesh_v1/` (already mirrored in the EL v1b
  case). Record explicitly that quarter-box and full-box are different
  discretizations of the same experiment.
- Re-run E4 serial and parallel with the unified diagnostics; regression
  baseline moves from "match an undocumented L2 run" to "match the pinned
  L3 trajectory + PIV within stated tolerance".

### D0.5 Reproducibility floor
- Re-green the two existing featflower_test cases on the campaign branch.
- One bitwise twin: 10-step sedimentation smoke, campaign binary vs current
  branch head, byte-identical `SED_BENCH_VEL`/`DNS_PART_STATE` lines (after
  D0.3 this becomes the walls-off-style regression anchor).

```

## Lines 359-390: §D1.3 ten Cate E1-E4 resolution x dt matrix, §D1.4, D1 deliverables and exit gate

```
### D1.3 ten Cate E1–E4 resolution × dt matrix (the core study)
- Full-box mesh (`ten_cate_mesh_v1`), all four regimes (Re 1.5 / 4.1 /
  11.6 / 31.9), levels L2/L3/L4, dt ladder ×{2, 1, 1/2, 1/4} around the
  guide value at the middle level.
- Observables per run: full v(t) trajectory vs PIV (peak velocity, time of
  peak, approach-to-wall deceleration), plus `DNS_CFL` statistics.
- Gates: L-convergence monotone toward PIV; quantify Δu_t(D/h, dt, Re).
  Bottom-wall approach phase analyzed separately (it is lubrication/contact
  dominated — belongs to D2's crossover map, recorded not gated here).
- Cost estimate before launch (probe-first): one L4 probe run timed before
  committing to the full matrix; expected O(20–30) jobs, minutes (L2) to
  ~day (L4 E1, viscous time scale) each.

### D1.4 Free-fall cross-check at moderate Re (optional, decision point)
- Uhlmann & Dušek (2014) sphere-settling regimes give published DNS with
  documented resolution requirements (D/h up to ~24) — one matched case
  (e.g. Ga ≈ 144, vertical regime) would anchor the guidelines against an
  independent DNS. Include only if the D1.3 cost model says an adequate
  mesh fits the cluster (decide at stage midpoint).

**Deliverables D1**: `dns_practitioners_guide.md` v1 with (a) recommended
D/h per target accuracy and Re band, (b) dt rule (particle CFL bound +
grid-crossing rate), (c) measured FBM force convergence order and noise
floor, (d) cost-per-(D/h, level) table for sizing future runs. Datasheet
rows for every matrix cell (PASS against convergence-model prediction, or
RECORDED).

**Exit gate D1**: E4 at recommended resolution within stated tolerance of
PIV with an error bar that the D1.1/D1.2 analysis explains; guidelines
published.

---
```

## Lines 711-740: §11 Risks and standing decision points (item 5 = the two ten Cate meshes)

```
## 11. Risks and standing decision points

1. **Cost blow-up** is the main risk: resolved DNS scales as (D/h)³ per
   particle per level. Mitigation: D1's cost table is a *gate* for D3.2/D3.3
   scoping; every stage sizes N and resolution from measured cost, not
   hope. Descopes are recorded, not hidden.
2. **`q2p1_dns_drag` is serial-PE-only** (parallel path aborts by design).
   Fine for fixed arrays (no PE dynamics needed), but rank-count-limited
   via the CFD side; if L4 arrays are needed, either extend the setup to
   parallel PE or accept the resolution ceiling — decision at D3 design.
3. **Grid-crossing noise may dominate tolerances** at affordable D/h; if
   D1.2 shows a large noise floor, gates shift from instantaneous values to
   time-averaged/filtered observables (state the filter in the RUNBOOK).
4. **DKT chaos**: post-contact trajectories are not gateable; the design
   pre-commits to pre-contact observables so a FAIL cannot be argued away.
5. **Two ten Cate meshes** (quarter benchSym vs full box) are different
   discretizations; every cross-comparison states which is used. The
   quarter-box symmetry assumption itself gets one full-box vs quarter-box
   twin row in D1.
6. **Parallel-PE Cartesian partitioning** depends on the unmaintained
   METIS-4 `tools/partpy` path — D0.1 must produce a supported recipe
   before any parallel-PE campaign runs.
7. **FullC0ntact backend is out of scope** except D4's OBJ-mesh twin;
   recorded so the campaign's claims are clearly scoped to the PE backend.
8. **Existing DNS regression baselines were produced under
   `HardContactEulerLagrange` builds** (including the committed
   featflower_test sedimentation baseline). The D0.6 solver split
   therefore needs its own equivalence twin, and baselines are re-pinned
   under the DNS solver once that twin is green — never silently reused
   across solvers.
```

## Lines 759-771 and 805-845: §13 Literature references (header; D0/D1 single-sphere block incl. ten Cate 2002, Uhlmann & Dusek 2014, Causin 2005)

```
## 13. Literature references

Full citations, grouped by campaign role. Status tags: **[in repo]** = PDF
available locally — DNS-campaign papers live in `literature/` (untracked;
index in `literature/README.md`), EL-era papers remain at repo root until
the dashboard's path check is updated; **[wanted]** = please provide the
PDF (the EL campaign showed having the originals on hand pays off during
case design); **[optional]** = useful context, not gate-defining;
*(verify)* marks bibliographic details recalled from memory that should be
checked against the actual paper. First batch of 12 PDFs received
2026-07-31; the three *(verify)*-tagged entries among them were checked
against the title pages and confirmed.

[...]
### D0/D1 — single sphere, resolution & dt
- A. ten Cate, C.H. Nieuwstad, J.J. Derksen, H.E.A. Van den Akker,
  "Particle imaging velocimetry experiments and lattice-Boltzmann
  simulations on a single sphere settling under gravity", *Phys. Fluids*
  14 (2002) 4012–4025. **[in repo]** (`ten_cate_piv.pdf`) — E1–E4
  definitions and the PIV trajectories to digitize in D0.4.
- M. Schäfer, S. Turek, "Benchmark computations of laminar flow around a
  cylinder", in *Flow Simulation with High-Performance Computers II*,
  Notes Numer. Fluid Mech. 48, Vieweg (1996) 547–566. **[optional]** — DFG
  C_D/C_L reference already regression-pinned.
- H. Hasimoto, "On the periodic fundamental solutions of the Stokes
  equations and their application to viscous flow past a cubic array of
  spheres", *J. Fluid Mech.* 5 (1959) 317–328. **[in repo]**
  (`literature/hasimoto_1959.pdf`) — D1-ARR-CONV dilute-array drag
  expansion.
- A.A. Zick, G.M. Homsy, "Stokes flow through periodic arrays of spheres",
  *J. Fluid Mech.* 115 (1982) 13–26. **[in repo]**
  (`literature/zick_homsy_1982.pdf`) — D1-ARR-CONV tabulated drag at
  general φ.
- A.S. Sangani, A. Acrivos, "Slow flow through a periodic array of
  spheres", *Int. J. Multiphase Flow* 8 (1982) 343–360. **[optional]** —
  cross-check on Zick–Homsy.
- M. Uhlmann, J. Dušek, "The motion of a single heavy sphere in ambient
  fluid: a benchmark for interface-resolved particulate flow simulations
  with significant relative velocities", *Int. J. Multiphase Flow* 59
  (2014) 221–243. **[in repo]** (`literature/uhlmann_dusek_2014.pdf`) —
  includes documented resolution requirements (D/h up to ~24), directly
  comparable to our guideline study; D1.4 go/no-go now unblocked on the
  reference side.
- N. Mordant, J.-F. Pinton, "Velocity measurement of a settling sphere",
  *Eur. Phys. J. B* 18 (2000) 343–352. **[optional]** — experimental v(t)
  time series, alternative free-fall anchor.
- P. Causin, J.-F. Gerbeau, F. Nobile, "Added-mass effect in the design of
  partitioned algorithms for fluid–structure problems", *Comput. Methods
  Appl. Mech. Eng.* 194 (2005) 4506–4527. **[in repo]**
  (`literature/causin_2005.pdf`, title page verified) — theory for the
  empirically observed lower dt stability bound of the loose FBM↔PE
  coupling (E4/L3 dt=0.25 ms blow-up, datasheet `dt_stability`): proves
  decreasing dt aggravates the added-mass instability in partitioned
  schemes at density ratios near 1.

```
