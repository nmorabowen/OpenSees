---
wp: LEGACY
title: "Stage 2 PR-3a (3D strong-axis embedded hinge) SHIPPED (33_ladruno_dispbeamcolumn3d_hinge_adr): LadrunoDispBeamColumn3d…"
section: "history"
legacy_seq: 122
---
- **Stage 2 PR-3a (3D strong-axis embedded hinge) SHIPPED ([[33_ladruno_dispbeamcolumn3d_hinge_adr]]):** `LadrunoDispBeamColumn3d` `-hinge $matTag` (reuses `ELE_TAG 33014`) — the strong-axis (Mz) rotation jump `α_z`, the literal 2D scalar algebra on the Mz row of the 6-DOF 3D basic system. One guarded rank-1 condensation to the basic system **before** `CorotCrdTransf3d` (pinned invariant through the quaternion triad); `hingeKvZ` a 6-vector incl. cross-tangent rows. Gated → no-hinge bit-identical; `-hinge`+`-nl` rejected. `tests/test_ladrunoDispBeamColumn3d_hinge.py` (12/12): patch test 1e-9, `∫Mz d[[θz]]==Gf` (19.999992), total-dissipation==Gf, **finite-rotation invariance under `CorotCrdTransf3d`** (>0.5 rad member rotation still dissipates Gf), solver robustness, nIP-objectivity, DB roundtrip. Full 65/65 (2D+3D). ADR 33 scoped + adversarially reviewed by a 17-agent workflow (killed the two-rank-1 / 2D-twice / det-floor shortcuts). **PR-3b** = weak-axis `α_y` + true coupled 2×2 `K_αα` (eigenvalue-floored inverse). Reserve `MAT_TAG_LadrunoCohesiveHingeBiaxial = 33004` (coupled biaxial cohesive, v2, not yet built).
