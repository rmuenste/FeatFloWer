# Euler–Lagrange framework: implemented methods and their origins

Inventory of every physical/numerical component in the unresolved
Euler–Lagrange (CFD-DEM) extension on `feature/euler-lagrange-phase1`,
with the commonly used name and an author/origin hint for proper citation.
Config keys and source locations are given so each entry can be traced to
code. Citation hints are from memory — verify year/venue before quoting.

## 1. Coupling framework

| Component | As implemented | Common name / origin hint |
|---|---|---|
| Governing model | Volume-averaged Navier–Stokes with fluid volume fraction ε_f; point particles with closure forces (`source/src_el/`) | Unresolved Euler–Lagrange / CFD-DEM. Volume-averaged equations: **Anderson & Jackson (1967)**, Ind. Eng. Chem. Fundam. First DEM-CFD coupling: **Tsuji, Kawaguchi & Tanaka (1993)**, Powder Technol. (fluidized bed). |
| Two-way coupling convention | Feedback to fluid = drag + lift only; pressure-gradient, gravity, buoyancy never spread (`result%feedback_force`, el_forces.f90) | "Model A / set II" formulation convention; classification per **Zhou, Kuang, Chu & Yu (2010)**, J. Fluid Mech. 661 (also **Feng & Yu (2004)**). |
| Fluid solver | Q2/P1 FEM, multigrid, MPI domain decomposition (FeatFloWer core) | In-house lineage: **Turek (1999)** monograph (FEATFLOW); cite the FeatFloWer code papers of the group. |
| Periodic-suspension body-force treatment | `ELFluidGravity=No`: fluid gravity absorbed into modified pressure; uniform counter-force cancels mean feedback (`AddGravForce` guard + `ConstantForcing`) | Standard periodic-box suspension treatment ("backflow correction" / zero-mean-pressure-gradient frame): cf. **Ladd (1994)**, J. Fluid Mech. 271 (Parts 1–2). |

## 2. Grid–particle transfer (G2P / P2G)

| Component | As implemented | Common name / origin hint |
|---|---|---|
| Interpolation/deposition kernel `deen_poly` | Quartic polynomial kernel with compact support δ (el_kernel_functions) | Polynomial "template/filter function" of **Deen, Van Sint Annaland & Kuipers** (Chem. Eng. Sci., mid-2000s CFD-DEM papers; also Link/Deen et al.). Functional family = **Lucy (1977)** kernel (Astron. J., SPH origin). |
| Alternative kernel `gaussian` | Gaussian coarse-graining kernel | Standard statistical coarse-graining kernel; usage in EL two-way coupling per **Capecelatro & Desjardins (2013)**, J. Comput. Phys. |
| Kernel bandwidth | δ = `ELKernelWidthFactor`·d_p, default 2.5 | Coarse-graining bandwidth practice (several diameters), e.g. **Xiao & Sun (2011)**; δ/h ≈ 4–6 grid-sampling rule is in-house practice. |
| Boundary treatment | Shepard (partition-of-unity) renormalization of truncated kernel support near walls | **Shepard (1968)** interpolation, ACM proceedings. |
| Void-fraction floor / relaxation | ε clipped at `ELEpsFMin` (default 0.4); optional under-relaxation `ELEpsFRelax` | Common CFD-DEM regularization (packing-limit clip near random close packing); in-house parameterization. |

## 3. Hydrodynamic force closures on the particle

| Component | As implemented | Common name / origin hint |
|---|---|---|
| Stokes drag (`ELDragModel=stokes`) | F = 3πμd·u_slip | **Stokes (1851)**. |
| Schiller–Naumann drag (`schiller_naumann`) | C_D = (24/Re)(1+0.15·Re^0.687), Re<1000 | **Schiller & Naumann (1933)**, Z. Ver. Dtsch. Ing. |
| Di Felice drag (`difelice`, campaign default) | One consistent package: interstitial slip, Re_p = ρ ε d|u_slip|/μ, C_D = (0.63+4.8/√Re)², voidage factor ε^(−χ), χ(Re) = 3.7−0.65·exp(−(1.5−log₁₀Re)²/2); zero-slip limit 2.88πμd·ε^(−χ−1) | Voidage function: **Di Felice (1994)**, Int. J. Multiphase Flow 20. C_D correlation: **Dallavalle (1948)**. Low-Re limit of Dallavalle = 0.96×Stokes (documented in code). |
| Pressure-gradient force (`ELPressureForce`) | F_p = −V_p ∇p (kernel-sampled) | "Pressure gradient / generalized buoyancy" term of the particle equation of motion: **Maxey & Riley (1983)**, Phys. Fluids 26 (also Anderson–Jackson form). |
| Gravity/buoyancy | Analytic (ρ_p−ρ_f)·V_p·g (never spread to fluid) | Elementary; convention per Model A above. |
| Saffman shear lift (`ELLiftModel=saffman`) | F = 1.615·d²·√(ρμ·|ω|)·(u_slip×ω)/|ω| | **Saffman (1965)**, J. Fluid Mech. 22 (+ 1968 corrigendum). |
| Mei finite-Re lift correction (`saffman_mei`) | Multiplicative factor f(Re_p, β), β = ½Re_ω/Re_p clamped [0.005, 0.4] | **Mei (1992)**, Int. J. Multiphase Flow 18, fitting **McLaughlin (1991)**. |
| Zeng wall lift (`saffman_mei_wall`, near-wall arms) | Translation-induced wall lift C_L = 3.663/(Re²+0.1173)^0.22 (→5.87 as Re→0) and shear-induced composite with wall-distance coefficients; validity Re ≤ 200 with clamp warning | **Zeng, Najjar, Balachandar & Fischer (2009)**, Phys. Fluids 21, 033302 (Eq. 19 and Eqs. 28–29). |
| Slip-lift blending ("option A") | Re-switch: Zeng arms for Re ≥ 2, Saffman–Mei for Re ≤ 0.5, log-Re blend between | In-house blend (document as such; no external citation). |
| Neutrally buoyant inertial lift (`ELInertialLift=matas_asmolov`) | Tabulated ĝ(r/R) profile digitized from the Rc = 30 matched-asymptotics curve; F_r = ĝ(s)·ρU_max²a⁴/(8√2·R²); zero crossing s_eq = 0.675 | "Inertial (Segré–Silberberg) migration force", matched asymptotic expansions: theory lineage **Ho & Leal (1974)** (channel, regular perturbation), **Schonberg & Hinch (1989)**, **Asmolov (1999)** J. Fluid Mech. 381; the curve used is Figure 14 of **Matas, Morris & Guazzelli (2004)**, J. Fluid Mech. 515 (their computation by Asmolov's method). Phenomenon: **Segré & Silberberg (1961/1962)**, Nature 189 / J. Fluid Mech. 14. |
| NOT implemented (documented gaps) | Added mass (C_M = 0.5), Basset history force, Magnus lift (flag exists, aborts), lubrication (lives only in PE's mutually exclusive HardContactLubricated solver) | Added mass/Basset: **Maxey & Riley (1983)**, added-mass coefficient discussion **Auton, Hunt & Prud'homme (1988)**. Magnus: **Rubinow & Keller (1961)**. |

## 4. Time integration & stiff-drag coupling

| Component | As implemented | Common name / origin hint |
|---|---|---|
| Fluid time stepping | Backward Euler (θ-scheme infrastructure), Q2/P1 saddle-point solve per step | FEATFLOW-lineage (Turek). |
| Particle drag update | PE `elSemiImplicitVelocity`: exact relaxation of the linearized drag (dragB, carrier velocity, other forces) — unconditionally stable for dt ≫ τ_p | "Point-implicit / semi-implicit drag integration", standard stiff-drag treatment in CFD-DEM & MP-PIC codes; general context: **Balachandar & Eaton (2010)** Annu. Rev. Fluid Mech. (review). |
| Fluid-side implicit drag sink (`ELDragCoupling=semi_implicit`) | +θΔt·B on the momentum diagonal (fine level) + old-time defect term; characterized first-order O(3·dt) momentum drift | Semi-implicit two-way momentum exchange (fluid side); common in CFD-DEM literature, e.g. **Xu & Yu (1997)** Chem. Eng. Sci. lineage. In-house characterization (dt-halving study, tier2). |

## 5. Rigid-body / contact side (PE engine)

| Component | As implemented | Common name / origin hint |
|---|---|---|
| Rigid-body engine | `pe` physics engine, MPI domain decomposition with shadow copies | **Iglberger & Rüde (2009/2010)** ("pe" rigid body physics engine papers, Comput. Sci. Eng. / Multibody Syst. Dyn.). |
| Contact resolution | HardContactEulerLagrange solver: non-smooth hard contacts, relaxed projected Gauss–Seidel iteration, Coulomb friction, restitution (campaign: dry contacts, e = 0) | Non-smooth contact dynamics with PGS relaxation: **Preclik & Rüde (2015)**, Comput. Part. Mech. (pe's HardContactSemiImplicitTimesteppingSolvers); NSCD lineage **Moreau (1988)** / **Jean (1999)**. |
| Periodic wrap-around decomposition | `decomposePeriodicX3D/Z3D/XY3D/Periodic3D` half-space process connections with ±L body offsets; per-axis json keys | pe domain-decomposition machinery (Iglberger & Rüde); Z-variant and per-axis dispatch in-house (this campaign). |
| Random seeding | Deterministic lattice-candidate shuffle (mt19937_64, rank-0 broadcast) with min-gap and wall-inset contract | Variant of random sequential addition (RSA, cf. **Widom (1966)**) on a lattice; implementation in-house. |

## 6. Validation references used as ground truth

| Benchmark | Origin |
|---|---|
| Terminal velocity, intermediate Re (E1–E4) | **ten Cate, Nieuwstad, Derksen & Van den Akker (2002)**, Phys. Fluids 14, 4012 (PIV + lattice-Boltzmann). u_∞ defined via **Abraham (1970)** drag correlation. |
| Hindered settling exponent | **Richardson & Zaki (1954)**, Trans. Inst. Chem. Eng. 32 (n ≈ 4.65 at low Re). |
| Tubular pinch / annulus at 0.6R | **Segré & Silberberg (1961, 1962)**; Re-dependence of the equilibrium radius: **Matas, Morris & Guazzelli (2004)**. |
| Laminar pipe profile | Hagen–Poiseuille (classical). |
| Suspension viscosity (V4, planned) | **Einstein (1906/1911)** dilute limit; **Krieger & Dougherty (1959)**, Trans. Soc. Rheol. (in-repo `mu_eff` model). |
| Channel suspension profiles (V5 context) | **Lyon & Leal (1998)**, J. Fluid Mech. 363. |

## 7. Quantified in-house scheme properties (cite as "this work")

- Kernel self-voidage bias: ε_eff ≈ 0.974–0.978 at the particle from its own
  deposit → single-particle drag ×ε^(−χ−1); consistent across ten Cate E1–E4
  and the SS migration-rate check (v1b / v3_ss_frozen RUNBOOKs).
- Two-way self-induced co-flow: +6% on single-particle settling (ten Cate
  E4) — the disturbed-vs-undisturbed carrier velocity problem; mitigation
  (not implemented): undisturbed-field reconstruction, cf. **Horwitz & Mani
  (2016)**, J. Comput. Phys. / **Gualtieri et al. (2015)**.
- Semi-implicit fluid-drag momentum drift ≈ 3·dt (first order, dt-halving
  characterized).

## 8. Component ↔ validation cross-reference

Which concrete FeatFloWer test targets which component. Test name = ctest
target, tier2 harness case (`tools/featflower_test/testcases/definitions/`),
or `validation_cases/<case>` RUNBOOK. Status as of 2026-07-05.

| EL component | Dedicated / primary tests | Supporting evidence | Status |
|---|---|---|---|
| Kernel G2P/P2G transfer + ε_f field | `test_el_kernel_forces` (kernel positivity/support), `el-transfer-mpi-2/8`, `el-convergence-serial`; tier2 `straddling_conservation` (cross-rank deposit, residual ~5e-26) | `EL_VOLUME_CONSERVATION` machine-zero in every campaign run (incl. ~8000 periodic crossings) | PASS |
| Two-way feedback (Model A) | tier2 `momentum_conservation` (+`_long`: 1e-5 per 10⁴ steps) | v2 smoke mixture-momentum balance 0.4% at ρ_p≠ρ_f; v1b E4 two-way co-flow quantified | PASS |
| Di Felice drag (voidage package) | v1b_tencate_settling E1–E4 (EXTERNAL truth: +0.4…−6.5% algebraic vs u_∞) | tier2 `terminal_velocity`; unit Di Felice/Stokes ratio checks; **v2_rz_settling production sweep = ε^χ arbiter (in flight)** | PASS / in flight |
| Stokes & Schiller–Naumann drag | unit analytic checks (`test_el_kernel_forces`) | tier2 terminal (stokes); v1b handoff SN cross-table | PASS |
| Semi-implicit drag coupling (fluid + particle side) | ctest `pe-el-semi-implicit-drag`; tier2 `momentum_conservation_semi` (gate 3.1e-3) | dt-halving characterization (drift ≈ 3·dt, first order); v3 convergence run at dt = 1.8·τ_p stable | PASS |
| Saffman–Mei slip lift | tier2 `saffman_lift` (analytic drag 2.356e-5 / lift 1.277e-6 to 1e-9) | v3_ss_frozen dense variant: sustained migration at 0.73–0.76× pure Saffman = Mei trim, both fan extremes | PASS |
| Zeng wall lift + option-A blend | unit coefficient anchors (Re-switch, contact/shear arms) | exercised in v3_ss_frozen ON runs (near-wall fan members); no dedicated quantitative wall-lift gate yet | PASS (unit) / gate open |
| Matas–Asmolov inertial lift | unit table anchors (s=0.4/0.675/0.8/0.42); v3_ss_frozen rate check (profile shape 1–3%, rates 0.84–0.90× ε=1) + convergence demo (annulus 0.6714 vs 0.675) | **v3_ss_coupled lift ON/OFF pair (in flight)** | PASS / in flight |
| Pressure-gradient force | — (flag `ELPressureForce`, off in all current cases to avoid hydrostatic double-count) | no dedicated test; candidate: prescribed-∇p unit check | gate open |
| Gravity/buoyancy (grav_buoy) | tier2 terminal; v1b (SI, ρ_p≠ρ_f) | v2 settling; force-budget balance 1e-7 at plateau | PASS |
| Self-voidage / disturbed-field bias (scheme property) | quantified: v1b E1–E4 (ε_eff 0.974–0.978) + v3 rate check (0.974) cross-consistent | documented as "this work"; mitigation (Horwitz–Mani) not implemented | QUANTIFIED |
| PE hard contacts (PGS, dry, e=0) | v2 φ=0.20 smoke (`EL_CONTACT_STATS`: overlap 0.0, Tgran bounded, 5.6k contacts) | **Stage-2 four-way sheared box + periodic-image contact check PENDING** | partial |
| PE periodic wrap (per-axis) | v3_ss_frozen wrap verification (z: radii frozen through 2 wraps, 2.5e-15) ; v2 smoke (fully periodic banner + run) | convergence run: ~500 wraps/particle; z-path bit-identical regression after dispatch refactor | PASS |
| Seeding (file + deterministic random) | Stage-0 battery: byte-identical across PE decompositions (2 seeds), achieved-φ 2% gate, file-mode bit-for-bit regression | v2 smoke: N=3056 at φ=0.2000147 | PASS |
| dt/PE-stepsize contract | Stage-0: mismatch hard-abort verified | enforced at every case start | PASS |
| Pipe infrastructure (CFD z-periodicity, cylinder wall geometry) | `pipe_hp_check` (HP profile −0.06% at L3, 6.2× convergence) | v3_ss_frozen lift-OFF radial freeze to all digits (geometry/transfer clean) | PASS |
| Prescribed frozen fields (`linear_shear`, `poiseuille`) | tier2 saffman (shear); v3_ss_frozen (poiseuille, exact parabola) | — | PASS |
| Periodic-suspension body-force treatment (`ELFluidGravity`) | v2 smoke: broken (fluid free-fall −0.2) vs fixed (momentum 0.4%) documented pair | — | PASS |
| Effective suspension viscosity (Krieger–Dougherty) | — | **V4/V5 pipe pressure-drop campaign PENDING (Stage 5)** | pending |
