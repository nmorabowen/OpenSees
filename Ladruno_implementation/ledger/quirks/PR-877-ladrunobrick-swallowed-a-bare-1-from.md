---
wp: PR-877
title: "A material that fails with a bare −1 inside a LadrunoBrick reached the solver as a success (sentinel-only host contract; producers now return the sentinel, C3b)"
pr: "#877"
date: 2026-10-10
---
### A material that fails with a bare `-1` inside a `LadrunoBrick` reached the solver as a success (sentinel-only host contract; producers now return the sentinel, C3b)
- **Bites:** a material that fails with the plain OpenSees code −1 (LadrunoRCConcrete's loud crack-band failure, a condensation miss, `StagedStrain -maxStrain`) inside a `LadrunoBrick` returned 0 from `update()`: the step "converged" on the stale state, while the same material inside `LadrunoQuad` or `TenNodeTetrahedron` cut the step.
- **Why:** ADR-86b made every brick call site test `== LADRUNO_MATERIAL_REFUSED` only, to keep ASDConcrete3D's advisory negative codes (−10 IMPL-EX error control, −1000 eigen) from failing steps; a bare −1 fell through with them. ADR-94 R5 keeps the host that way on purpose.
- **Fix/status (2026-10-10, C3b as amended by the third review of #877):** the host is NOT changed (six `setTrialStrain(x) == LADRUNO_MATERIAL_REFUSED` sites; `tests/test_adr94_contract.py` pins them). The producers convert: `StagedStrain -maxStrain` and `LadrunoRCConcrete`'s crack-band failure returns now return the sentinel. A material that means "best state, carry on" must not return the sentinel; a material that fails must not return a bare −1 into a brick. Still open outside the pin: `LadrunoQuad` tests `!= 0` and `LadrunoRCFiniteStrain::setTrialF` returns a bare −1.
