---
wp: WP-134
title: "ModifiedEuler freezes K and G at the COMMITTED state for the whole increment — an error that does not shrink with TolR and is invisible to the error test (U9,…"
legacy_seq: 493
---
### `ModifiedEuler` freezes K and G at the COMMITTED state for the whole increment — an error that does not shrink with TolR and is invisible to the error test (U9, WP-134 → WP-129)
- **Bites:** `ModifiedEuler` builds `aC = GetStiffness(K, G)` from the K, G it is handed (the committed-state moduli, `commitState()`), and neither stage re-evaluates them at the stage stress, although G ∝ √p. Both Heun stages share the same wrong moduli, so their difference -- the error estimate -- cannot see it. WP-134's exact oracle measured it at 0.6 / 6 / 24 % of the stress increment for strain increments of 1e-5 / 1e-4 / 1e-3, independent of TolR. After the elastic→plastic stage flip `mK/mG` even hold the stage-0 (pressure-independent) moduli until the first commit.
- **Workaround/status:** SAS-ME (IntScheme 129) evaluates K, G at every stage's state (and uses Heun on the moduli for the elastic predictor and the elastic part of an intersected increment). ModifiedEuler unchanged (byte-identity).
