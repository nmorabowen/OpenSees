---
wp: WP-108
title: "LadrunoSANISAND::schemeReachesModifiedEuler() returned false for IntScheme 2, so the class printed a false \"-maxSubsteps has NO EFFECT\" warning — FIXED (WP-108)"
legacy_seq: 460
---
### `LadrunoSANISAND::schemeReachesModifiedEuler()` returned false for `IntScheme 2`, so the class printed a false "-maxSubsteps has NO EFFECT" warning — FIXED (WP-108)
- **Note on provenance:** this defect was FOUND by WP-105 (F12, PR [#844](https://github.com/nmorabowen/OpenSees/pull/844), still open/unmerged as of this writing) and its own `LEDGER_quirks.md` row lives only on that branch, not on `ladruno` — this WP-108 branch was cut from `ladruno` before #844 merged, so that row could not be edited here; this is therefore a NEW row recording the same finding plus the fix. **When #844 merges, its own row ("schemeReachesModifiedEuler() returns false for IntScheme 2 ... MEASURED FALSE") should be updated to point at this fix (or merged into this row) rather than left saying "not fixed."**
- **Bites:** trust the constructor's own warning and you would conclude `-maxSubsteps`/`-honorTolR`
  are inert on scheme 2 and skip capping it. WP-105 (F12) measured this false on the `p -> p_min`
  floor path: a `-maxSubsteps 100` cap on scheme 2 turned a run that completed **40 of 40 uncapped**
  into one that **refuses at step 18**. The warning was wrong; the seam was live all along.
- **Why (verified on `634824e1f`):** `LadrunoSANISAND::schemeReachesModifiedEuler()`
  (`LadrunoSANISAND.cpp:1035-1048`) returned `false` for `mScheme == 2`, which gated the
  constructor's warning text (`:1206-1213`, "`-maxSubsteps N has NO EFFECT with IntScheme 2`") and
  the identical claim for `-honorTolR` (`:1187-1199`). But `ManzariDafalias::explicit_integrator`'s
  `switch (mScheme)` (`ManzariDafalias.cpp:1070-1101`) does not enumerate `INT_BackwardEuler`, so a
  call into it falls to `default:` -> `ModifiedEuler`, exactly where both seams
  (`mMaxSubstepsInME`, `mHonorTolR`) are read. Scheme 2 routes there whenever
  `BackwardEuler_CPPM`'s own recursive-halving retry ladder (`ManzariDafalias.cpp` ~2472-2588)
  falls back on non-convergence or ladder exhaustion — conditionally, not on every step the way
  schemes 0/1 do, but "sometimes reaches it" is not "never reaches it", and the warning claimed the
  latter. This directly contradicted `LadrunoSANISAND_implex_guide.md` §3, which REQUIRES
  `-maxSubsteps > 0` on scheme 2 under `-implex` and refuses the deck without it
  (`LadrunoSANISAND.cpp:2047-2058`, the `setLadrunoImplexOptions` `s == 2` branch) — one of the two
  had to be wrong, and WP-105's measurement said it was the warning.
- **Workaround/status (2026-09-16, WP-108):** FIXED. `schemeReachesModifiedEuler()` now returns
  `true` for `mScheme == INT_LSANISAND_BackwardEuler` (2) as well as 0/1/>9-not-45, with the
  function's header comment and a new `tests/test_ladruno_sanisand_intscheme2_maxsubsteps.py`
  (reproducing WP-105's own floor-path control) pinning both that the false warning no longer
  prints and that the cap actually fails a step on scheme 2. Scheme 1's output is unchanged
  (byte-identical) — the fix only adds a branch for `s == 2`, ahead of the existing `s > 9` catch-all.
  Full measurement: `Ladruno_files/testbed/hypo_bearing/adr92_f12/F12_intscheme2_verdict.md`
  section 5.4/6 (WP-105).
