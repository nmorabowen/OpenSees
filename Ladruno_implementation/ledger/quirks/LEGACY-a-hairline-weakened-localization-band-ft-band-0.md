---
wp: LEGACY
title: "A hairline-weakened localization band (ft_band ~ 0.95 ft) drives BULK GPs into the Concrete3D return-map apex — seed bands with >= 20% strength margin"
legacy_seq: 151
---
### A hairline-weakened localization band (ft_band ~ 0.95 ft) drives BULK GPs into the Concrete3D return-map apex — seed bands with >= 20% strength margin
- **Bites:** weakened-band localization tests/models with `LadrunoConcrete3D`: with the band only 5% weaker, the bulk sits at ~95% of its own tensile onset at band peak, and the global Newton's trial excursions through the localization transition push bulk GPs past onset into the deep-tension apex regime — the run drowns in "return map did not converge -> step-cut" (the kernel's safe fallback), grinding to a crawl without ever failing outright.
- **Why:** Newton iterates are not monotone: mid-iteration trial strains overshoot the converged state by far more than the 5% margin; the Concrete3D apex regime is exactly where the return map is trajectory-fragile (documented handoff §6 gap).
- **Workaround/status (2026-07-06, ADR 66 P5.2 G5):** seed localization with ft_band = 0.8*ft (dissipation is Gf-governed, so energy gates are unaffected by the seed strength); the churn disappears entirely.
