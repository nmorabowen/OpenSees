---
wp: LEGACY
title: "Stage 2 PR-2a (element-side embedded hinge, 2D) SHIPPED: LadrunoDispBeamColumn2d -hinge $matTag (reuses ELE_TAG 33013)…"
section: "history"
legacy_seq: 120
---
- **Stage 2 PR-2a (element-side embedded hinge, 2D) SHIPPED:** `LadrunoDispBeamColumn2d` `-hinge $matTag` (reuses `ELE_TAG 33013`) — single scalar rotation jump `α` carried by any UniaxialMaterial, Armero–Ehrlich strong-discontinuity split (bulk sees `κ_bulk=B·v−α/L` → unloads, no double-count; cohesive `M([[θ]])` carries `Gf`), `α` inner-Newton + **guarded** static condensation to the 3-DOF basic system **before** `crdTransf` (pinned invariant). Gated (`hingeOn`) → no-hinge path bit-identical; `-hinge`+`-nl` rejected. `tests/test_ladrunoDispBeamColumn2d_hinge.py` (8/8): patch test 1e-9, energy gate `∫M d[[θ]]==Gf` (LINEAR ~4e-7), element total-dissipation==Gf (no double-count), tangent-consistency via tight Newton through softening, commit/revert + DB roundtrip. ADR 32 carries a 4-agent adversarial review that corrected the plan (E–B not Timoshenko; no section "freeze"; `K_αα` indefinite not zero-at-peak; Dirac out of quadrature).
