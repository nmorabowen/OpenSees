---
wp: WP-129
title: "RungeKutta45 never computes dAlpha3 / dAlpha4: its α update weights sum to 301/336 (WP-129, code-read)"
legacy_seq: 491
---
### `RungeKutta45` never computes `dAlpha3` / `dAlpha4`: its α update weights sum to 301/336 (WP-129, code-read)
- **Bites:** in `ManzariDafalias::RungeKutta45` stages 3 and 4 compute `dSigma`, `dFabric`, `dPStrain` but not `dAlpha`; the function-scope static vectors are zeroed at entry, so `dAlpha3 = dAlpha4 = 0` always. Stage 4's and 5's α arguments, the 5th-order α update `(14 dα1 + 35 dα4 + 162 dα5 + 125 dα6)/336` and the α error estimate all read the zeros: every accepted step under-advances α by ~10 %. On top of WP-128's findings (dT_min 1e-3 hard-coded, Mc-clamp force-accept, no drift correction).
- **Workaround/status:** do not use IntScheme 45 as a reference (WP-129's test (b) uses WP-128's validated α-aware port instead). Not fixed (byte-identity of existing schemes).
