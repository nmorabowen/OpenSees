---
wp: LEGACY
title: "Finite-strain elastoplastic bending/necking BVPs need KrylovNewton (plain Newton diverges); F-bar needs an unsymmetric solver"
legacy_seq: 63
---
## Finite-strain elastoplastic bending/necking BVPs need KrylovNewton (plain Newton diverges); F-bar needs an unsymmetric solver

From the finite-strain validation Phase P1 (2026-06-02, [[18_finite_strain_validation_report]]).

- **`LogStrain(LadrunoJ2)` + `LadrunoBrick -geom finite` bending *into plasticity* does
  NOT converge under `algorithm Newton` — nor under `NewtonLineSearch`.** The residual
  grows from the very first increment (a `NormDispIncr`/`EnergyIncr` norm that climbs, not
  shrinks). It is not a tangent bug (the consistent tangent passes the FD gate on
  homogeneous states): bending+plasticity on a low-order hex is just a stiff, badly-scaled
  Newton basin for a full step. **`algorithm KrylovNewton` (+ `test EnergyIncr 1e-6`) is
  robust** and converges quadratically-ish; the necking bar (C1) and the isochoric-J2
  locking cantilever (B3) both rely on it. Homogeneous single-element states and elastic
  bending converge fine under plain Newton — the divergence is specific to *inhomogeneous
  plastic* finite-strain BVPs.
- **A 1-element-wide cross-section bends too poorly for stable plastic Newton.** A 1×1×nz
  column under transverse displacement control diverges even with KrylovNewton; a ≥2×2
  cross-section is needed. (Elastic load-control on the 1-wide column is fine — it just
  locks.)
- **F-bar (`-formulation bbar -geom finite`) has an UNSYMMETRIC tangent** (dSNPO eq 15.10)
  ⇒ use `system FullGeneral` (dense) or, much faster for meshed studies, `system UmfPack`
  (unsymmetric sparse). A symmetric solver silently mis-solves. `UmfPack` made the 128–576
  hex necking runs tractable where `FullGeneral` would be dense-O(N³) per step.
- **Plastic finite-strain stress paths are path-dependent** (obvious, but it bites tests):
  a sub-stepped element solve does NOT equal a single-step return-map oracle for
  *non-proportional* loading (simple shear, equibiaxial). Drive ONE increment ref→F when
  comparing to a one-step oracle, or step the oracle incrementally over the same F_k.
