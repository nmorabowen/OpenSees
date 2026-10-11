---
wp: WP-182
title: "region ... -rayleigh 0 0 0 0 does nothing: the region keeps an earlier global rayleigh pair (WP-182)"
date: 2026-10-10
pr: "#968"
---
### `region <tag> ... -rayleigh 0 0 0 0` does nothing — the all-zero factors skip `setRayleighDampingFactors`, so the region keeps whatever an earlier global `rayleigh` set (WP-182)
- **Bites:** zeroing the damping of a sub-model (e.g. the interior soil of an SSI box) with a region after a global `rayleigh`: the region silently keeps the global pair.
- **Why:** `TclRegionCommands.cpp` l.342-343 calls `theRegion->setRayleighDampingFactors(...)` only `if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)`.
- **Workaround:** use region-only damping (no global `rayleigh`), or read the factors back in an audit. Not changed (vanilla behaviour).
