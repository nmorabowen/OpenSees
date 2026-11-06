---
wp: LEGACY
title: "LadrunoIMKBeam: \"elastic end\" ≠ \"released end\" (pin); fake a release with a tiny Elastic hinge"
legacy_seq: 59
---
## LadrunoIMKBeam: "elastic end" ≠ "released end" (pin); fake a release with a tiny Elastic hinge

A `LadrunoIMKBeam` end with **no hinge material** in that slot is **elastic** — it
still carries moment with the full elastic rotational stiffness (`4EI/L`,
far-end-fixed). That is NOT a structural **release** (pin), where the end moment is
identically zero and the rotation is free (the far end condenses out, dropping the
active end stiffness to `3EI/L`). The element has no `-release` flag (deferred, see
[[14_ladruno_imk_beam]] §8).

- **Workaround:** put a near-zero-stiffness `Elastic` uniaxial in that end/axis
  slot (e.g. `-matZj`). The series flexibility `1/k → ∞` zeroes the end moment.
  Use `k ≈ 1e-5·(4EI/L)` of the bending axis (~1e-5 of the elastic rotational
  stiffness): verified to give `k_i = 3EI/L` to ~5 figures and residual
  `M_j/M_i ≈ 7e-6`. Stay above the element `ktFloor` guard (`1e-8·4EI/L`); a
  factor in `[1e-6, 1e-4]` is the sweet spot (smaller hurts conditioning, larger
  leaves a non-negligible residual moment). Learned 2026-06-02 (LadrunoIMKBeam).
