---
wp: PR-877
title: "LadrunoBrick swallowed a bare −1 from setTrialStrain (sentinel-only contract) — a failed material state reached the solver as a success (fixed, C3b)"
pr: "#877"
date: 2026-10-10
---
### LadrunoBrick swallowed a bare −1 from `setTrialStrain` (sentinel-only contract) — a failed material state reached the solver as a success (fixed, C3b)
- **Bites:** a material that fails with the plain OpenSees code −1 (LadrunoRCConcrete's loud crack-band failure, a condensation miss, `StagedStrain -maxStrain`) inside a `LadrunoBrick` returned 0 from `update()`: the step "converged" on the stale state, while the same material inside `LadrunoQuad` or `TenNodeTetrahedron` cut the step.
- **Why:** ADR-86b made every brick call site test `== LADRUNO_MATERIAL_REFUSED` only, to keep ASDConcrete3D's advisory negative codes (−10 IMPL-EX error control, −1000 eigen) from failing steps; a bare −1 fell through with them.
- **Fix/status (2026-09-26, Ladruno C3b):** cut on the sentinel OR −1, nothing else. Still NOT a blanket `< 0`. Materials that mean "best state, carry on" must not return −1.
