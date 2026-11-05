# Conservative continuity form for the EL volume-averaged constraint (`ELVAContinuityForm = conservative`)

Branch `feature/el-va-conservative-form` (worktree `FeatFloWer-va-formc`, build `build-formc`),
from the VA head e7ae87fc. Design consultation before coding: GPT-6 Astra high,
`el_va_conservative_form_design_review_2026-11-05.md` (brief and answer verbatim).

## Why

The volume-averaged continuity constraint of the EL branch enforces `div u = g` per element with
`g = -(d_cont + adv)/eps`, where `adv = u8 . grad eps_e` is the 8-corner mean velocity times the
cell-averaged kernel gradient (`el_volume_averaged.f90`). That product is not a conservative flux:
the discrete form of `int div(eps u) dV`, exactly zero in a closed box, does not vanish, and the s04
VB return-flow twin keeps +0.55 % (level 2) and +0.21 % (level 3) in G5 even with BDF2, one
corrector and the alpha sink (RUNBOOK §10-§11). The peer's spec 04 §12.4 and its Codex report
ranked the "conservative form C" (shared face flux, `D_eps(1) = 0` and `D_eps(z) = -int eps u_z`
exactly) as the durable fix, staged as a conservative residual iterated through the existing
projection. This note records the discrete design as implemented.

## Discrete design

**Transport porosity (a reconstruction, not a new storage variable).**
`epsilon_transport_q2(ndof)` holds the trilinear (Q1) interpolant of vertex values in the Q2 dofs:
the vertex value is the `dvol`-weighted mean of the P0 `epsilon_f` of the elements sharing the
vertex, summed across ranks with `E013Sum2` (the certified Q2 nodal sum, global level pinned as
for `E013Sum3`), and the edge, face and centre dofs are the means of their 2, 4 and 8 defining
vertices (`EL_VA_FLUX`: `EL_VA_RECONSTRUCT_EPS_LOCAL`, `EL_VA_EMBED_Q1_IN_Q2`; rebuilt in every
transfer pass right after `EL_FINALIZE_VOID_FRACTION`, never restarted). Continuous across ranks,
bounded by the element values, degree Q1. The P0 `epsilon_f` stays the conserved storage variable
(deposition, rate `d_cont`, sink weights, restart, `EL_MEAN_SLIP`).

**Flux moments** (`EL_VA_ELEMENT_FLUX_MOMENTS`; 3x3 Gauss per face, 3x3x3 in the volume, the Q2
isoparametric geometry of the B matrices with all 27 nodes; the outward normal
`n dS = sign_f (dx/dxi_a x dx/dxi_b) dxi_a dxi_b` from the cyclic tangent pair of the fixed
reference coordinate of face f, no centre test):

```
F_1(u) = sum_faces int_f eps_h u_h.n dS
F_m(u) = sum_faces int_f (x_m - xc_m) eps_h u_h.n dS - int_e eps_h u_h,m dV     m = 2..4
```
with `xc` the element-centre Q2 node (the P1 basis convention of `QuadSc_BMatrix.f`), plus
`int_e eps_h u_h dV`, `int_e u_h dV`, `int_e 1 dV`, `int_e (x_m - xc_m) dV`. Both elements adjacent
to an interior face evaluate the same nodal data at the same nine points with opposite normals, so
`sum_e F_1 = Phi` (boundary flux) is an algebraic identity on any conforming mesh, partition faces
included; no neighbour lookup, no face exchange, `kadj` unused. No-slip faces carry nine Dirichlet
nodes and contribute exactly zero. The coordinate identities
`sum_e [xc_m F_1 + F_1+m] + int eps_h u_m dV = boundary term` hold the same way.

**Storage moments** with the deposition volume `|e| = dvol` (stored single precision by SETARE):
`Q_1 = |e|`, `Q_m = |e| int_e (x_m - xc_m) dV / |e|_geom`; the gap `sum dvol - sum |e|_geom` is
logged (`dvol_gap`).

**Residual and projection right-hand side** (per source evaluation with the current velocity
`u^k`; the assembled `defp` is `-(1/dt) B^T u^k`, so `D_m = -dt defp_m` is `int_e q_m div u^k`
in the quadrature of B itself; sign chain verified: `Build_BTMatP1` positive, `BT_Mul_U_mod`
factor -1/dt, `Velocity_Correction` gives `D du = dt rhsP`):

```
r_m   = F_m(u^k) + d_cont Q_m            unsinked conservative residual
eps_e = MAX(epsilon_f, 1e-3)             preconditioning coefficient (as before)
gV_m  = D_m - r_m/eps_e                  source volume integral, all four rows
int_g = sum_e gV_1  ->  sink s_e as before (uniform gbar, alpha weight), row 1 only
defp_1 := -r_1/(eps_e dt) - s_e |e|/dt ;   defp_m := -r_m/(eps_e dt)   (m = 2..4)
```
i.e. `defp_m + D_m/dt = gV_m/dt`, the direct residual form. Iterated through
`EL_VA_Projection_Correctors` it is form C at the fixed point: `F_m = -d_cont Q_m - delta_m1
eps_e |e| s_e`, and summing row 1 gives `sum_e eps_e |e| s_e = -(S + Phi)` with
`S = sum_e d_cont |e|`: the sink vanishes exactly when storage and boundary flux are compatible
(closed box, no clipping, no relaxation). Row 1 enters the campaign slots as `g = gV_1/|e|`, and
the audit "advective" term becomes `adv = -(eps_e gV_1 + d_cont |e|)/|e|`.

**What the identity controls.** With all four rows the fixed point gives
`int eps_h u_z dV = sum_e d_cont (xc_z |e| + Q_z)` (the storage first-moment rate `m1rate`; the
`Q_z` part vanishes on bricks). That is the reconstructed-eps integral, not the P0-weighted
integral of `EL_MEAN_SLIP`; both are logged (`uf_eps_z`; `uf_p0_z` = P0 eps times the Q2-rule
velocity integral of the same loop, which is close to but not the `EL_MEAN_SLIP` quadrature) and
only like quantities are compared.

## Code

| File | Change |
|---|---|
| `source/src_el/el_va_flux.f90` (new, `EL_VA_FLUX`) | MPI-free: Q2 basis in the FEAT numbering, reference tables, per-face and per-element kernels on plain arrays, driver over the local elements, vertex reconstruction and Q1 embedding |
| `source/src_el/el_volume_averaged.f90` | form id (`EL_VA_FORM`, table `el_va_form_names`), one dispatch outside the element loop (legacy loop verbatim), `EL_VA_CONSERVATIVE_RAW_SOURCE`, the four-row `defp` block, flux block `F_*` and its max (audit steps), `EL_VA_FLUX` / `EL_VA_FLUXMOM` lines, header |
| `source/src_el/el_fields.f90` | `epsilon_transport_q2`, `EL_ENSURE_TRANSPORT_EPS` (separate from the size-gated `EL_ENSURE_FIELDS`), release |
| `source/src_el/el_transfer.f90` | reconstruction after the void-fraction finalisation, deck-only guard, `E013Sum2` under the pinned level |
| `source/src_el/el_config.f90` | key `SimPar@ELVAContinuityForm = advective | conservative`, validation (table; needs `continuity`; `conservative` needs `ELVASourceCorrectors >= 1` and `ELMesoFilter = No`), echo, `EL_VA_CONSERVATIVE_FORM()` |
| `applications/q2p1_el_pipeflow/app_init.f90` | parser case (same file as the guard above) |
| `applications/q2p1_el_pipeflow/app_init.f90` | communication-scheme guard at startup after `EL_VALIDATE_CONFIG` (`E013Sum2` exists for `iCommSwitch` 3 and 4), before the first transfer pass |
| `source/src_el/tests/test_el_va_flux.f90` | unit test (below); CMake target `test_el_va_flux`, ctest `el-va-flux-serial` |
| `tools/el_stream_compare.py` | tags `EL_VA_FLUX t=`, `EL_VA_FLUXMOM` |
| sample deck | `validation_cases/s04_voidfraction_cfdem_twin/q2p1_param_VB_bdf2_alpha_cons.dat` |

Collectives under the form, all keyed on deck values and the pass/step labels: the readiness
reduction at the top of every source evaluation (`EL_VA_CONSERVATIVE_READINESS`: the inherited
rank-local early returns are validated collectively first, an unready worker stops the run
coherently), the layout-check reduction and the reconstruction sum in the post-advance transfer
pass (every worker, with or without particles; the pre-advance pass does not rebuild the field
nobody reads), the flux block and its max on audit steps after the campaign/temporal/weight
reductions, the eps-weighted sink integral on audit steps. The default path has no new statement
in its arithmetic or its collective pattern (gates below).

Audit lines (audit steps, with the `corr=` tag): `EL_VA_FLUX S= Phi= sink_eps= res1_L1= res1_L2=
res1_max= resL_L1=` (storage rate, boundary flux, `sum eps_e s_e |e|`, unsinked row-1 residual
norms, linear rows per element length) and `EL_VA_FLUXMOM id_x= id_y= id_z= uf_eps_z= uf_p0_z=
m1rate= dvol_gap=` (coordinate identities, both superficial velocities over the box volume, the
storage first-moment rate, the volume gap). `EL_VA_RESIDUAL` (dt times the assembled row 1) is
unchanged and is the sinked projection residual; it does not certify conservation.

## Unit test (`test_el_va_flux`, synthetic mesh of `el_test_helpers`, serial)

1. The module basis reproduces FEAT's `E013A` values and reference derivatives (the local
   numbering assumed by `EL_VA_ELEMENT_DOFS`).
2. One cell: face areas and outward orientation from uniform velocities; the divergence theorem
   for `u = x` (`F_1 = 3 |e|_geom`) on the brick and on a warped cell; volume moments.
3. Shared face: the flux through the common face of two adjacent elements cancels pairwise
   (brick and warped, random nodal eps, wall-tangent Q2 velocity).
4. Closed box 3x3x3 (brick and warped interior vertices): `sum_e F_1 = 0` and the three
   coordinate-moment identities to 1e-12 of their scale; volumes sum to the box.
5. Constant eps, linear velocity: `F_1 = eps |e|_geom tr(A)` per element.
6. Reconstruction: a linear element field is reproduced at every dof whose defining vertices are
   interior; across a jump every nodal value stays within the element range.

All pass (2026-11-05). `test_el_kernel_forces` and `test_el_transfer` unchanged and passing.

## Reviews (`el_va_conservative_form_reviews_2026-11-05.md`, brief and answers verbatim)

Codex gpt-6-astra high and Claude Fable 5.1 xhigh on dc2221e0. Both: no blocker, the legacy
loop and the legacy `defp` block byte-identical to the parent, the algebra, sign chain, kernel
tables, signed normal, weights and reconstruction verified by reading. Applied in the fix commit:
collective readiness validation under the conservative key (Codex 1, Fable 2); the
communication-scheme guard moved to startup (Codex 2, Fable 1); `m1rate` completed with the
linear-row storage term `d_cont Q_z` (Codex 3); the layout check made a function returning the
mismatch count, run once in the first post-advance transfer pass with a rank-uniform stop, its
trilinear restriction documented, no hidden module state (Codex 4, Fable 3); the pre-advance
reconstruction skipped (both); one table-driven form lookup `EL_VA_FORM_ID` in `EL_CONFIG`
(Fable 10); the row algebra factored into the pure `EL_VA_CONSERVATIVE_ROWS` with a balanced and
an incompatible-storage test, plus a negative layout test (Fable 8); comment and wording fixes
(sink integral `= -(S + Phi)`, constant-row-only sink as a choice, two centres, log length,
verbatim-block markers, `uf_p0_z` wording). Deferred, recorded as follow-ups: deduplicating the
fifteen campaign-slot statements (needs its own bitwise gate), a curved-Q2 shared-face test, a
real-partition `E013Sum2` test with a particle-free worker, face-trace restriction to the nine
on-face nodes (performance), tolerance-based stopping.

## Runs

<!-- filled in as the smoke, the gates and the T10 ladder complete -->

## Adopted from the design consultation (and what was not)

Adopted: the direct four-row residual right-hand side in its own routine instead of an `adv`
formula; the Q1 vertex reconstruction (Q2 nodal averaging is not bounded); the separately
audited storage rate, boundary flux, eps-weighted sink integral and unsinked residual; logging
both the reconstructed-eps and the P0-eps superficial velocity; an explicitly signed
reference-face normal with one weight convention instead of the `ExtractBoundaryNormals_Q2`
template (inward normal, factor-of-four `SurfDet`); the storage-moment convention with the
deposition volume and the logged volume gap; the separate allocation of the transport porosity;
the `E013Sum2` communication-scheme guard; the meso-filter exclusion; certification on the closed
box first, GAMMA = 0; fixed corrector counts with the convergence read from the unsinked
residual. Not done here (follow-ups): a tolerance-based stopping key with a master-communicated
decision; open and periodic boundaries (the z-moment seam term); double-precision storage
volumes; velocity-coefficient caches; the full operator `F H` in the pressure matrix.
