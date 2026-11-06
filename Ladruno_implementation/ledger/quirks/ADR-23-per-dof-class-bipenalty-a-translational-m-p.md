---
wp: ADR-23
title: "Per-DOF-class bipenalty: a translational m_p CANNOT bound the rotation mode (ADR 23 M1/ES-1)"
legacy_seq: 70
---
### Per-DOF-class bipenalty: a translational `m_p` CANNOT bound the rotation mode (ADR 23 M1/ES-1)
- **Why it bites:** the bipenalty mass penalty `m_p` (lumped on the slave's translational
  DOFs) bounds the explicit `dt_cr` of the TRANSLATIONAL coupling only. The rotation tie's
  penalty `K_r` has different units (moment/rotation), so a translational-only `m_p` leaves
  the rotation mode UNBOUNDED in explicit (`dt_cr → 0`). Fix: give the rotation class its
  OWN inertia `I_p = K_r·(dt/2)²` (the SAME `-dtcr`/`-wcap` budget formula but with `K_r`),
  lumped on the slave's ROTATION DOFs. Then `dt_r = 2√(I_p/K_r) = dt_u` and the `lch²` in
  `K_r` cancels (it's also in `I_p`), so the rotation mode self-bounds at the SAME `dt`.
  `getExplicitCriticalTimeStep` reports the MIN over active DOF classes. (The same pattern
  generalizes to a pressure class if pressure bipenalty is ever added — pressure is
  implicit-recommended for now.) 2026-06-04.
