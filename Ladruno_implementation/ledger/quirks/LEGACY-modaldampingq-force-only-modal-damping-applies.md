---
wp: LEGACY
title: "modalDampingQ (force-only modal damping) applies damping with the WRONG SIGN"
legacy_seq: 16
---
### `modalDampingQ` (force-only modal damping) applies damping with the WRONG SIGN
- **Bites:** `modalDampingQ ζ` on a free-vibration SDOF/MDOF *amplifies* the
  response instead of damping it — measured ζ comes out ≈ **−ζ_target**.
  `modalDamping ζ` (the matrix form) is correct (+ζ).
- **Why (evidence, not yet line-pinned):** the sign error is **Δt-independent** —
  refining the step (steps/period 100→800) converges Q to −0.04996…→−0.04998 and
  the matrix form to +0.05000. A lag/explicit-stability artifact would vanish on
  refinement; this doesn't, so it's a **structural sign inconsistency**, not a
  numerical one. The only live force routine is the M-weighted
  `IncrementalIntegrator::addModalDampingForce` (`IncrementalIntegrator.cpp:502`,
  `setB` at :556); the two earlier variants (:303, :349) are commented out. Both
  `modalDamping` and `modalDampingQ` add that SAME `−Cv` force in `formUnbalance`
  (`TransientIntegrator.cpp:135`); the ONLY difference is `modalDamping` also adds
  `+c·C` to the tangent (`addModalDampingMatrix` :563, `formTangent` :88). So the
  matrix term is somehow compensating a force that, alone, has the wrong net sign.
- **PURE SIGN INVERSION — confirmed:** `modalDampingQ(-0.05)` damps correctly at
  **+0.05003**. So the force-only path injects `+Cv` energy instead of dissipating
  `-Cv`. The residual sign convention itself is fine (`FE_Element::addRIncInertia
  ToResidual` does `theResidual += -1.0·getResistingForceIncInertia`,
  `FE_Element.cpp:517`; the modal `-Cv` from `addModalDampingForce`/`setB` matches
  it). So the inversion is NOT in the force expression — it's that a velocity-
  proportional force applied through the *residual only* (no `+c·C` in the tangent,
  the term `modalDamping` adds and `modalDampingQ` omits) is integrated by Newmark
  with the opposite effective sign. i.e. `modalDampingQ` as a standalone force-only
  mode appears to have never worked; only `modalDamping` (matrix + force together)
  is correct.
- **Workaround/status (2026-06-01):** **use `modalDamping`, never `modalDampingQ`.**
  Reproduced under both `Newton` and `Linear`, Δt-independent. Pinned as a `strict`
  xfail: `tests/test_damping_channels.py::test_modalDampingQ_force_only_matches_matrix`.
  Genuine upstream bug. **DECISION (2026-06-01): DOCUMENT-ONLY** — not patched, not reported
  upstream; `modalDamping` (matrix) covers the use case. The `strict` xfail auto-detects any
  future upstream fix (would start XPASSing). A naive "flip the force sign" would break
  `modalDamping` (shares the same force), so any fix must special-case the `inclMatrix==false`
  path. Full map: [[12_damping_channels]].
