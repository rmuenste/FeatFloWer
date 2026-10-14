# EL drag inputs: undisturbed-field and voidage-sampling correction — design note

Date 2026-10-14. Branch `feature/el-rz-postfix`. Owner sequencing: after the
v2w level-3 rung (done, null), "look into the undisturbed-field correction".

## 1. What the column twin established

`v2w_column_settling` (RUNBOOK sections 7–10) measured the EL hindrance law in
the DNS D3.2 column with byte-identical clouds:

| rung | u_t | n (through origin) | note |
|---|---:|---:|---|
| baseline δ = 2.5 d, level 2 | 0.3573 | 1.72 | 12 clouds |
| kernel δ = 1.5 d | 0.3189 | 1.15 | deficit grows when the kernel narrows |
| level 3 (h = d/8) | 0.3571 | 1.73 | identical to level 2 to 3 digits |
| DNS, same column | 0.4061 | 4.57 | resolved spheres |

Kernel smearing and mesh resolution are eliminated. What remains are the
two inputs the Di Felice closure receives from the kernel sample
(`EL_INTEGRATE_PARTICLE` → `EL_COMPUTE_PARTICLE_FORCES`):

1. **Voidage ε seen by the particle.** `sample(18)` is the kernel-weighted
   average of the element field `epsilon_f`, which is itself `1 − α_p` with
   α_p deposited by the same kernel: the voidage is smoothed TWICE with
   radius δ. In a cloud only 5–6 d across, particles see a suspension about
   a third as dense as the one they form (N = 120: ε seen 0.927 vs
   1 − φ_cloud 0.888). Di Felice's own homogeneous-suspension prediction
   U/U₀ ≈ ε^(χ+1), χ(Re≈80) = 3.1, gives 0.73 with the sampled ε but 0.61
   with the true cloud voidage; DNS is 0.588. **The voidage dilution alone
   carries most of the exponent deficit.**
2. **Carrier velocity seen by the particle.** `sample(2:4)` is the
   kernel-averaged fluid velocity, which contains the particle's own
   disturbance (the deposited feedback force drives a co-moving fluid
   column). This is the undisturbed-field problem proper. It lowers the
   slip, hence the drag: V1b two-way E4 settles +6% faster than one-way; the
   k15 rung's u_t fell 11% when the kernel narrowed because the
   self-disturbance concentrates. For a Stokes regularised point force of
   Lucy profile the self-induced velocity at the particle is
   u_self = F/(6πμ δ)·c with c = ⟨1/r⟩_W·δ (see §3), i.e. a fraction
   c·a/δ of the slip at Re → 0, reduced by finite-Re (Oseen) effects. This
   bias is the same fraction for every particle in a dilute cloud, so it
   mostly DIVIDES OUT of U/u_t; it governs the single-particle offsets
   (EL u_t −12% vs DNS) more than the exponent.

Consequence for ordering: the voidage-sampling fix is cheaper, affects the
exponent directly, and is a prerequisite for interpreting the velocity
correction; do it first.

## 2. Stage U0 — decoupled voidage sampling (recommended first)

- New deck key `ELVoidageKernelWidthFactor` (default = `ELKernelWidthFactor`
  → bitwise no-op). In `EL_INTEGRATE_PARTICLE`, accumulate `epsilon_f` with
  a second kernel of radius δ_ε and its own normalisation (extend
  `EL_SAMPLE_SIZE` by one slot `EL_SAMPLE_EPS_NORM`); deposition and velocity
  sampling keep δ. Candidate δ_ε = 1.0 d (the voidage of the immediate
  neighbourhood, as in cell-based CFD-DEM) and 0.75 d.
- Alternative U0b (particle-centric): ε_local = 1 − Σ_neighbours V_p·W_ε(r_ij)
  computed from the halo records, no mesh field at all; exact for the
  cloud edge and independent of h. Slightly more code (pair loop over
  records, already present for the lubrication virials), same key.
- Gates: tier-2 twin bitwise with the key unset; lone-particle ε_eff moves
  toward 1 as δ_ε shrinks (predictable, §3); v2w ladder rerun (12 × 1.5 h):
  expected n to rise toward Di Felice's homogeneous ≈ χ+1 ≈ 4; the
  remaining gap to DNS 4.57 is then the velocity-input part.

## 3. Stage U1 — undisturbed carrier velocity

Options, in order of fidelity vs cost:

- **U1a, in-house calibrated correction.** u_undist = u_sampled − K(δ/h, Re)·F_fb^(n−1)/(μ δ),
  with K measured once by a unit test: one particle held fixed in a
  quiescent periodic box, constant feedback force F deposited with the
  production kernel, solve to steady state, sample u at the particle →
  K = u_sampled·μδ/F. Repeat at the production Re (Oseen reduction) and
  near walls (a 1-D table in wall distance, using `EL_WALL_DISTANCE`).
  One-step-lagged F keeps the scheme explicit; the drag's semi-implicit
  path (`drag_B`) is unaffected. This is Horwitz & Mani (2016) in spirit,
  with the discrete transfer function measured instead of modelled.
- **U1b, closed form.** Balachandar, Liu & Lakhote (2019) give the
  self-induced velocity for a Gaussian kernel with Re dependence; our
  kernel is Lucy (`deen_poly`), so the Stokes coefficient must be
  re-derived (§3 numbers) and the Re dependence taken from their fit —
  usable as a cross-check of U1a, not as the primary.
- **U1c, Horwitz & Mani discrete Stokes solve** (per-particle local Stokes
  problem on the kernel stencil): most faithful, heaviest; only if U1a's
  table proves insufficient (strong wall or Re dependence).

Gates: tier-2 twin bitwise with the key off; V1b two-way E4 must land on
the one-way value (the +6% co-flow row `u_t_E4_twoway` becomes a PASS);
k15-style kernel sensitivity of u_t must shrink (0.357 vs 0.319 today);
v2w ladder (with U0 on) must not regress the exponent.

## 4. Numbers (computed 2026-10-14, scratch script siv_coeff)

- Lucy kernel: <1/r>*delta = 2.6250  -> Stokes self-induced velocity / slip = c*a/delta = 0.525 at delta=2.5d, 0.875 at delta=1.5d
- delta=2.5d: eps_eff(lone particle) = 0.9724  (measured V1b one-way ~0.975 at 2.5d)
- delta=1.5d: eps_eff(lone particle) = 0.8723  (measured V1b one-way ~0.975 at 2.5d)
- N=120 sampled eps=0.927: Di Felice homogeneous U/U0 ~ eps^(chi+1) = 0.733   (EL measured 0.805, DNS 0.588)
- N=120 true 1-phi_cloud eps=0.888: Di Felice homogeneous U/U0 ~ eps^(chi+1) = 0.614   (EL measured 0.805, DNS 0.588)

The lone-particle ε_eff predicted from the double smoothing matches the
measured V1b value, confirming the voidage-dilution mechanism is the
kernel bookkeeping and not a solver effect.

## 5. Decision requested

Implement U0 (decoupled voidage kernel, both variants behind one key) and
rerun the twelve-cloud ladder (~18 node-hours), then decide U1 on the
residual. Estimated code: ~60 lines Fortran for U0a, ~120 for U0b.
