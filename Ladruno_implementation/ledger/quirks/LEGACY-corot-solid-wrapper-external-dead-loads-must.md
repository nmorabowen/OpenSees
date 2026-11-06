---
wp: LEGACY
title: "Corot solid wrapper: external dead loads must NOT pass through globalizeForce; and corot IS objective for kinematic hardening"
legacy_seq: 62
---
## Corot solid wrapper: external dead loads must NOT pass through `globalizeForce`; and corot IS objective for kinematic hardening

Two corot-seam gotchas, from the finite-strain trifecta deep review (2026-06-02,
[[10_solid_corotational_adr]]):

- **External dead loads must stay in the global frame.** `SolidTransformationCorot::
  globalizeForce` pushes the core force forward by R (`f_global = R f_d − …`). That is
  correct ONLY for the *internal* force (`∫ Bᵀσ`, self-equilibrated). A fixed-direction
  body/self-weight load (`-b`, `eleLoad -selfWeight`) is a GLOBAL-frame quantity — if it
  is folded into the core force before `globalizeForce`, corot rotates gravity WITH the
  element (wrong, non-conservative; was the COROT-1 bug). **Fix/pattern:** `LadrunoBrick`
  accumulates the body load in a separate `bodyForce` vector and adds it back AFTER
  `globalizeForce`/`globalizeStiff` (also keeps the spurious body-load term out of the
  corot geometric stiffness). Behavior-neutral under `-geom linear` (identity globalize);
  the `-geom finite` path was already correct (assembles in the spatial frame, no
  globalize). Any new fold-then-globalize site must keep external loads out of the core.

- **Corot is objective for KINEMATIC-hardening materials (unlike the LogStrain finite
  path).** A natural worry is that corot shares the dSNPO §14.11 backstress-frame
  non-objectivity. It does NOT: corot feeds the material `u_d = Rᵀ x_rel − X_rel`
  (REFERENCE frame) with reference-config gradients, so the small-strain material — and
  its backstress α — live in a FIXED reference frame across commits (the element's R
  rotates; the material's frame does not). Since `polar(Q·H) = Q·polar(H)` exactly for
  rigid Q, `u_d` is rigid-rotation-invariant ⇒ identical deformational-strain history ⇒
  exact objectivity (verified, `test_corot_kinematic_hardening_objectivity`). The
  LogStrain path differs precisely because `bᵉ_tr = f_Δ bᵉ_n f_Δᵀ` co-rotates the stress
  `s` into the current frame while α stays fixed — THAT is §14.11. Lesson: "the element's
  R rotates between commits" does NOT imply "the material's frame rotates."
