● Review summary                                                                                                                                                                                                   
                                                                                                                                                                                                                   
  The work is in genuinely good shape. Phases 1–3 of the plan are implemented and the last five commits closed out a textbook root-cause hunt: the ~9% one-time momentum drift was traced to NDFGL reading stale   
  /TRIAD/ offsets from a finer multigrid level, producing out-of-bounds force_rhs writes that silently dropped 26/27 of the spread stencil. The fix (pin /TRIAD/ to the coupling level for the whole pass,         
  mirroring the earlier ILEV pin) brought the element-integrated momentum drift from 9.2e-2 down to 4.3e-7, and the feedback-conservation residual to machine precision. That was the last blocker for the flagshi$
  Tier-2 conservation test.                                                                                                                                                                                        
                                                                                                                                                                                                                   
  What exists and works:                                                                                                                                                                                           
  - Phase 1 transfer: positive compact kernel (Lucy/Gaussian), gather/scatter with owner-reduce/broadcast halo, void fraction with clipping counter and relaxation hook, ∂ε_f/∂t storage, volume/feedback          
  conservation diagnostics — all matching the plan's §3.                                                                                                                                                           
  - Phase 2 forces: Di Felice/Stokes/Schiller-Naumann drag, −V∇p, driver-side grav_buoy = (ρ_p−ρ_f)Vg (correctly resolved against the HardContactEulerLagrange contract), semi-implicit sub-cycled PE drag         
  integrator.                                                                                                                                                                                                      
  - Phase 3 feedback: drag + lift only (correct Model-A convention), captured pre-E013Sum3 to avoid the shared-DOF double-count, plus the issue-D fix where the fluid is fed the implicit drag impulse PE actually 
  applies rather than the explicit old-slip drag — that's a subtle consistency point many CFD-DEM codes get wrong.                                                                                                 
  - Phase 4 scaffolding: Tier-1 convergence tests, four Tier-2 YAML definitions (momentum, terminal, Saffman, straddling), momentum diagnostics.                                                                   
                                                                                                                                                                                                                   
  Suggested next steps                                                                                                                                                                                             
                                                                                                                                                                                                                   
  1. Hygiene first (small, do now).                                                                                                                                                                                
  - Commit the untracked artifacts: the four planning/notes docs, docs/md_docs/el_triad_ndfgl_common_block_note.md, and testcases/terminal-velocity/. They document decisions the code depends on.                 
  - Promote the /TRIAD/ pin into the guarded EL_NDFGL_SAFE wrapper your own note already designs (docs/md_docs/el_triad_ndfgl_common_block_note.md:218), with the MAXVAL(KDFG) ≤ ndof assertion as a hard stop. The
  pin fixes this call site; the wrapper prevents the pattern from recurring — and this bug class has now bitten you twice (ILEV, then /TRIAD/).                                                                    
  - Fix or drop the EL_MOMENTUM_LUMPED diagnostic (known shared-DOF double-count artifact) — a diagnostic that reports a false 12% error will mislead someone eventually.                                          
  - Consider a PR to main for Phases 1–3 + the conservation fix; the feature branch is 17 commits deep and the coupling core is now verified.                                                                      
                                                                                                                                                                                                                   
  2. Finish Phase 4 V&V (the current phase).                                                                                                                                                                       
  - Run Tier-2 #8 terminal velocity — the case is fully specified in phase4-terminal-velocity.md and the QBOX9 mesh is staged; acceptance is u_t = 0.03556 within 1% plus the exponential transient. This is the   
  gate on the whole Phase-2 chain.                                                                                                                                                                                 
  - Run #9 (Saffman) and #10 (straddling), and extend the momentum-conservation run well past 60 steps (the plan's V2 asks for 10⁵-step stability) to confirm the drift stays flat and the ε_f-clip counter stays  
  at zero.                                                                                                                                                                                                         
                                                                                                                                                                                                                   
  3. Pay down the flagged physics debt before the validation campaign — three items will otherwise invalidate V2/V3:                                                                                               
  - Pin the Di Felice voidage exponent. The code implements ε^(2−χ) with an honest "provisional, must re-check" comment (el_forces.f90:103-126), the plan specifies ε^(−χ) with C_D = (0.63 + 4.8/√Re)², and the   
  implementation also uses Schiller-Naumann C_D and omits the ε_f factor in Re_p (plan: Re_p = ρ ε d|u_rel|/μ). All three choices are invisible at ε=1 and only discriminate in exactly the regime Richardson–Zaki 
  (V2) tests — decide the convention against the V2 reference first, then run V2.                                                                                                                                  
  - Replace the placeholder lift. The current saffman_mei is plain Saffman (constant 1.615, no Mei Re-correction), and saffman_mei_wall uses ε_f as the "wall factor" (el_forces.f90:179-181) — that is not a wall 
  correction at all, just a proxy. Segré–Silberberg (V3) is the lift validation and will quietly produce wrong equilibrium positions with this stand-in. Implement Mei (1992) C_L and a Zeng-family wall factor    
  driven by actual wall distance, unit-tested per the plan.                                                                                                                                                        
  - Semi-implicit drag on the fluid side (plan §5.2: drag_coupling = semi_implicit, reaction term in the momentum matrix). Feedback is currently explicit-only; fine for dilute tests, but the target is φ up to   
  20% with stiff drag, where explicit feedback will constrain Δt or destabilize. This is the main remaining architectural piece.                                                                                   
  4. Then the validation ladder: V2 Richardson–Zaki (φ sweep, periodic box) → four-way coupling (enable PE contacts + lubrication in EL mode, sheared-box stress test) → the pipe campaign (V3/V4/V5) with the     
  validation report under docs/el_validation_report.md.

  One process note: the diagnostics-first, one-finding-per-commit methodology of the last five commits is exactly why the OOB bug was findable — I'd keep the momentum audit (EL_MOMENTUM_ELEMINT + feedback
  residual) always-on at low frequency as a permanent regression tripwire, per the plan's el.momentum_audit_freq.

