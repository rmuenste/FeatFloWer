# Runbook: unit_periodic_contact — periodic-image contact check (Stage 2 gate)

Goal: prove the PE contact system sees collision partners THROUGH a
periodic boundary (periodic images), including cross-rank. The campaign
plan makes this a hard gate: if PE missed periodic images, every dense
periodic result (Stage-1 φ = 0.20, fourway_shear) would carry an asterisk.

## Configuration

- Mesh: QBOX9 unit box (symlink to tier2 momentum_conservation), level 2,
  np = 28. CFD `Periodic = Yes`, PE `periodicX_/Y_/Z_: true`.
- Particles: two spheres d_p = 0.05 seeded straddling the periodic z-plane
  IN CONTACT: (0.5, 0.5, 0.99) and (0.5, 0.5, 0.03) — surface gap through
  the wrap is negative (overlap ≈ 0.01 = 20% d_p). `seedAllowContact_:
  true` (seeder min-gap assertion relaxed for this case only).
- Field: `linear_shear` with rate 0.0 (quiescent frozen field — isolates
  the contact response), gravity off, 20 steps, dt = 0.002.

## Result (2026-07-08, local np = 28, level 2, 20 steps) — PASS

- The contact between the pair fires at step 1 (EL_CONTACT_STATS
  ncontacts ≥ 1, overlap detected) — the partner is only visible as a
  periodic image across z, and the pair straddles a rank boundary.
- The hard-contact response separates the pair: contact-free with Tgran
  decaying by step 20 (no sticking, no tunnelling through the plane).
- EL_NEWTON_PAIR worst |mismatch| = 4.5e-20 through the periodic-image,
  cross-rank contact (pe 8b037ae fix regression-guarded here).

This is a permanent fast regression case for the periodic-image contact
path; results recorded in commit 161fec09 and the validation datasheet.
