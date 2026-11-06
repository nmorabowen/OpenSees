---
wp: LEGACY
title: "D9 interface material returns FORCE, not stress — no bondScale (unlike the rebar)"
legacy_seq: 71
---
### D9 interface material returns FORCE, not stress — no `bondScale` (unlike the rebar)
- **Why it bites:** `LadrunoEmbeddedRebar`'s axial slot drives its bond material with the
  axial SLIP and the material works in STRESS units (τ–s), so the element multiplies by
  `bondScale = perimeter·L_trib` to get a nodal force. `LadrunoEmbeddedNode`'s D9 interface
  materials are driven by the displacement GAP `g·e_d` (metres) and are expected to RETURN
  FORCE directly (`stress()` in N) — so there is **NO bondScale converter**: `t_d =
  mat_d->getStress()`, `k_d = mat_d->getTangent()` go straight into `t=Σ t_d e_d`,
  `D=Σ k_d e_d⊗e_d`. Pick/define the uniaxial accordingly (an Elastic of "stiffness" K is
  a penalty of force-per-metre K; a cohesive law's peak is a force, not a traction). The
  D9 grammar confines `-mat*` to the TRANSLATIONAL normal/tangent directions, so the
  rotation-unit (M1/ES-1) problem never arises in a material slot. 2026-06-04.
- **v1 uses the REFERENCE frame; `-corot` frame co-rotation is DEFERRED.** A material on a
  specific direction (esp. a unilateral-contact normal) ideally co-rotates with the host;
  v1 keeps the frame fixed (valid for small-rotation interfaces). The v2 corot would reuse
  the Phase-2 continuum-rotation (`skew(∇u)` from host gradients) to rotate the frame — and
  would re-introduce the rebar's dropped `∂e_d/∂u` consistent-tangent caveat (frame-objective
  for explicit, converges under step-halving for implicit, may slow Newton for stiff-normal
  large-rotation contact). Large-rotation RIGOROUS contact is the separate `LadrunoContact`.
