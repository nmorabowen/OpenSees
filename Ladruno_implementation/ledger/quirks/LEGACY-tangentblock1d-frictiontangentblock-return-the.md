---
wp: LEGACY
title: "tangentBlock1D/frictionTangentBlock return the pressure-coupling term with kn ALREADY folded in -- reapplying it at the call site gives a kn^2 tangent that onl…"
legacy_seq: 313
---
### `tangentBlock1D`/`frictionTangentBlock` return the pressure-coupling term with `kn` ALREADY folded in -- reapplying it at the call site gives a `kn^2` tangent that only the non-symmetric path can see
- **Bites:** anyone scattering the consistent friction cross term in a new lane, and any review that gates `-mu` but not `-mu -consistanttan`.
- **Why:** `dTN = -(dcap/dN)*kn*sgn` -- the `kn` is INSIDE, matching the shipped 3D `frictionTangentBlock`'s `(-dCap_dN*kn*nh[i])*n[j]`. The 2D scatter initially wrote `dTN_s * kn * th[i] * n[j]`, i.e. `kn` twice; with a realistic `kn ~ 1e6..1e9` that is a catastrophically wrong tangent, not a subtle one. It is INVISIBLE on the default path because `dTN` is identically 0 unless `consistent == true` (the symmetric default, design-gate Q2), so every gate that does not exercise `-consistanttan` passes a build carrying it.
- **Workaround/status (2026-08-18, ADR-85 T2):** caught in the orchestrator's review pass BEFORE merge and fixed at the call site (`LadrunoContactFE::addFrictionTang2D` -- `kn` appears exactly once); both misleading doc comments corrected to say the caller must NOT reapply `kn`. Lesson for future lanes: a term that is zero on the default path needs its own gate, or it ships wrong.
