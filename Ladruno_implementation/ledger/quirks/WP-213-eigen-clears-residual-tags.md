---
wp: WP-213
title: "Re-running eigen drops residualModes pseudo-modes and their tags"
date: 2026-10-26
---
### Re-running `eigen` drops the `residualModes` pseudo-modes (and `-kinds` reads all 0 again)

- **Bites:** `eigen`, `residualModes -dir 1 2`, then another `eigen` (a new
  mode count, after a stage change, after `wipeAnalysis`): the domain's mode
  set is replaced, the pseudo-modes are gone and `residualModes -kinds`
  reports every mode as eigen. RSA after that is the truncated basis again.
- **Why:** the mode-kind registry (`LadrunoModeKinds`) is keyed by `Domain*`
  and only valid while the domain's eigenvalue vector is bit-identical to the
  snapshot taken by `residualModes`; any new eigenvalue vector makes it stale.
  This is how "a new eigen clears the tags" without editing the eigen path.
- **Fix/status:** call `residualModes` again after every `eigen`, then
  `modalProperties` (residualModes unsets the modal properties on purpose).
  Run it after any `modalDamping`, too. Learned 2026-10-26 (WP-213).
