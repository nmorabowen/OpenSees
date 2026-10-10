---
wp: LEGACY
title: "ADR scoping doc (now SHIPPED — see the LadrunoDispBeamColumn2d/3d, LadrunoCohesiveHinge + LadrunoCohesiveHingeBiaxial r…"
section: "history"
legacy_seq: 117
---
- **ADR scoping doc (now SHIPPED — see the `LadrunoDispBeamColumn2d/3d`, `LadrunoCohesiveHinge` + `LadrunoCohesiveHingeBiaxial` rows in the table above)** — Regularized displacement-based frame ADR ([[32_ladruno_dispbeamcolumn_regularization_adr]]): `LadrunoDispBeamColumn2d` 33013 / `LadrunoDispBeamColumn3d` 33014 (Element registry; numerically equal to `ND_TAG` 33013/33014 — not a collision) / `LadrunoCohesiveHinge` (uniaxial, MAT_TAG 33003). Two tiers: Tier-1 per-IP `lch` channel (mirrors `ForceBeamColumn`, fixes [[LEDGER_quirks]] §59), Tier-2 embedded strong-discontinuity hinge (Armero–Ehrlich). Large-disp via Corotational with the pinned 3-DOF-basic condensation contract. Scoped via 11-agent workflow (regularization/kinematics/cohesive/plumbing/validation × verify + synthesis). Docs only, no `SRC/` change
