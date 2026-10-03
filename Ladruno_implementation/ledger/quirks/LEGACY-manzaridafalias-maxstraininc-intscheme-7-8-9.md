---
wp: LEGACY
title: "ManzariDafalias::MaxStrainInc (IntScheme 7, 8, 9) sub-steps with UNINITIALISED nG, nK — the sub-steps can do nothing at all"
legacy_seq: 366
---
## `ManzariDafalias::MaxStrainInc` (`IntScheme 7, 8, 9`) sub-steps with UNINITIALISED `nG, nK` — the sub-steps can do nothing at all

**Found 2026-10-01, WP-158, measured.** Same defect as `MaxEnergyInc`'s (IntScheme 4, entry
"IntScheme 4 (`MaxEnergyInc` -> ForwardEuler) is NON-DETERMINISTIC"), in the sibling function:
when the largest strain component of the increment exceeds `maxStrainInc = 1e-5`,
`MaxStrainInc` declares `double nDGamma, nVoidRatio, nG, nK;` and passes `nG, nK` by reference
as the moduli of every `ForwardEuler` sub-step, which builds `aC = GetStiffness(K, G)` from
them without writing them first. Its loop also never advances `cEStrain` (as in
`MaxEnergyInc`). And its `switch` sends every case — 7 (`MAXSTR_MFE`), 8 (`MAXSTR_RK`), 9 — to
`ForwardEuler`, so the "ModifiedEuler"/"Runge-Kutta" variants are not.

- **Measured:** one step of 1e-4 or 2.5e-5 strain from a plastic, on-surface state under
  IntScheme 7, 8 or 9 leaves the yield function EXACTLY unchanged (Delta f = 0.0) on both
  the d63f49750 and the WP-158 build — consistent with zero garbage moduli (Ce = 0, no
  stress change). The WP-129 decks `ls3d_s7/8/9` pin whatever that garbage gives on the
  capture build; they did not move with WP-158.
- **Workaround/status:** do not use IntScheme 4, 7, 8 or 9. Not fixed (vanilla; it would move
  every scheme-4/7/8/9 deck that sub-steps). Fix shape: initialise `nG = G, nK = K` (or
  re-evaluate the moduli per sub-step) and advance `cEStrain`.
- **✅ FIXED — WP-160, [#914](https://github.com/nmorabowen/OpenSees/pull/914)** (2026-10-02):
  committed moduli for every sub-step of `MaxStrainInc` and `MaxEnergyInc`, `cEStrain` advanced;
  7/8 dispatch left on ForwardEuler (owner decision). Measurements in the FIXED bullet of
  "IntScheme 4 (`MaxEnergyInc` -> ForwardEuler) is NON-DETERMINISTIC".
