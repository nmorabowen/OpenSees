---
wp: LEGACY
title: "AUTO-resolved penalties already carry the element thickness -- never h-scale an auto value (the h^2 class has a second member)"
legacy_seq: 316
---
### AUTO-resolved penalties already carry the element thickness -- never h-scale an auto value (the h^2 class has a second member)
- **Bites:** any lane that composes `-thickness h` (or any per-unit-thickness density) with an `auto` stiffness resolved from `getInitialStiff()` -- the assembled element K already folds the element's real thickness, so multiplying the auto value by h again scales the penalty as t*h ~ h^2 (the same class as the gated lambda double-scale, one call site earlier).
- **Why:** the ADR SS How/7 states it for NTS auto-kn ("absorbs h automatically... no thickness parameter") but the T3 mortar injection initially h-scaled `epsUse` unconditionally; caught by 3 independent review finders. A defaulted `-epsT` must inherit the PROVENANCE of the value it derives from (auto-derived => no h; explicit => h), not a blanket rule.
- **Workaround/status (2026-08-18, ADR-85 T3):** fixed at the single injection site (`epsAuto ? epsUse : epsUse*h`, `epsTFromEpsN` provenance flag); no test targets `-epsN auto -thickness` yet -- flagged for the T4 battery.
