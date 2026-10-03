---
wp: WP-108
title: "LadrunoSANISAND::schemeReachesModifiedEuler() returns false for IntScheme 2, so the class prints \"-maxSubsteps has NO EFFECT\" — MEASURED FALSE — fixed (WP-108,…"
legacy_seq: 463
---
### `LadrunoSANISAND::schemeReachesModifiedEuler()` returns false for `IntScheme 2`, so the class prints "`-maxSubsteps` has NO EFFECT" — MEASURED FALSE — fixed (WP-108, #845)
- **Bites:** trust the constructor's own warning and you would conclude `-maxSubsteps`/`-honorTolR`
  are inert on scheme 2 and skip capping it. Measured on the `p -> p_min` floor path (WP-105 / F12):
  a `-maxSubsteps 100` cap on scheme 2 turned a run that completed **40 of 40 uncapped** into one
  that **refuses at step 18**. The seam is very much live.
- **Why (verified on `634824e1f`):** `LadrunoSANISAND::schemeReachesModifiedEuler()`
  (`LadrunoSANISAND.cpp:1035-1048`) returns `false` for `mScheme == 2`, and that false is what
  gates the constructor's warning text at `:1206-1213` ("`-maxSubsteps N has NO EFFECT with
  IntScheme 2`") and the identical claim for `-honorTolR` at `:1187-1199`. But
  `ManzariDafalias::explicit_integrator`'s `switch (mScheme)` (`ManzariDafalias.cpp:1070-1101`)
  does not enumerate `INT_BackwardEuler` among its cases, so it falls to `default:` ->
  `ModifiedEuler`, which is exactly where both seams (`mMaxSubstepsInME`, `mHonorTolR`) are read.
  Scheme 2 DOES route through `ModifiedEuler` whenever the CPPM falls back (the previous entry), so
  both seams are live on it. This directly contradicts
  `LadrunoSANISAND_implex_guide.md` §3, which REQUIRES `-maxSubsteps > 0` on scheme 2 under
  `-implex` and refuses the deck without it (`LadrunoSANISAND.cpp:2047-2058`) — one of the two had
  to be wrong, and the measurement says it is the warning, not the guide.
- **Fixed (WP-108, PR [#845](https://github.com/nmorabowen/OpenSees/pull/845)):** the one-line fix
  this row originally flagged as "not applied" has landed — `schemeReachesModifiedEuler()` now
  returns `true` for `mScheme == 2`, ahead of the existing `s > 9` catch-all. See the earlier row in
  this file ("... — FIXED (WP-108)") for the fix detail and verification. WP-105 / F12 itself made
  no source edits; this row records the finding, and the fix is credited to WP-108.
