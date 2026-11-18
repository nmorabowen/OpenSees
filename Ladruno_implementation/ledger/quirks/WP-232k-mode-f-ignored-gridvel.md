---
wp: WP-232k
title: "MPM mode F accepted `mpm gridvel` and `mpm pforce` and ignored them (#1109)"
date: 2026-11-18
---
### Before WP-232k, the explicit engine silently dropped `gridvel` and `pforce`

- **Bites:** any explicit (mode F) MPM deck that prescribed a grid velocity or a particle point
  force on a build before WP-232k. Both commands were parsed and stored (WP-232f, for mode A), but
  the explicit core had no prescribed-velocity condition (`LadrunoMPMCore.h`, `FaceConditions`
  only), and `LadrunoMPMPattern::analyze` never read `gridVels` or `pforce`.
- **Fixed (WP-232k, #1109):** mode F applies `gridvel` in both projections' grid conditions
  (`Engine::presc`, built per step by `buildPrescribed`), and refuses `pforce` loudly;
  `mpm condForce` is mode A's and is refused in mode F. Test:
  `tests/test_ladruno_mpm_soil_gridvel.py`.
