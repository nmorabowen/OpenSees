---
wp: WP-142
title: "A bit-exact replay of a committed uniaxial state needs ONE trial strain per step: many uniaxials ignore a trial change below an ABSOLUTE DBL_EPSILON (WP-142)"
legacy_seq: 523
---
### A bit-exact replay of a committed uniaxial state needs ONE trial strain per step: many uniaxials ignore a trial change below an ABSOLUTE `DBL_EPSILON` (WP-142)
`LadrunoUniaxialJ2::setTrialStrain` returns at once when `fabs(Tstrain - strain) < DBL_EPSILON` (same idiom in `ElasticPP`, `Hardening`, `Hardening2`, `FlagShape`, `Ratchet`, `SteelDRC`). The threshold is absolute: at a strain of 1e-3 it is ~1000 ulps. Under Newton the last iterate often moves the strain by less than that, so the committed bar state belongs to an EARLIER iterate and differs from the final element strain by < 2.2e-16. Physically nothing (Δσ ≤ E·2.2e-16), but an identity oracle that replays the element's final strain on a standalone copy then misses by 1 ulp in stress, and looks like a forwarding or rotation bug. Seen building O6: the 0° bar's forwarded strain was 1 ulp off the plate strain it must equal exactly.
- **Rule:** for a bit-exact replay drive the element with `algorithm Linear` (one trial strain per committed step; with Penalty-prescribed displacements the path is still imposed to O(K/α)), or replay the material's OWN committed strain (`... material strain`) and hold the formula to the dead-band.
- **Workaround/status:** ✅ both drivers in `tests/test_plateRebar_response_forwarding.py` (Linear: formula strain, bit-exact; Newton: own strain, formula within `DBL_EPSILON`). *2026-09-27.*
