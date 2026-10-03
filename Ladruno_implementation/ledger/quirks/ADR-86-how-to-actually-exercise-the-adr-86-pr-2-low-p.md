---
wp: ADR-86
title: "How to actually exercise the ADR-86 PR-2 low-p clamp diagnostic — it does NOT fire on the ADR §5 deck"
legacy_seq: 339
---
### How to actually exercise the ADR-86 PR-2 low-`p` clamp diagnostic — it does NOT fire on the ADR §5 deck
- **Bites:** commit `6b0a85c1e`'s message claims *"the p0 = 1.01 kPa slow leg emits 8 events"*. **It does not.** Re-measured on a freshly rebuilt engine at HEAD: the `ModifiedEuler` clamp fires **0 times on BOTH legs** of `test_presidual_is_the_low_p_defect` (vanilla constants and `p_r = 0`). The ledger row for that commit makes no firing claim, so the build-control record was never wrong — but the commit message is, and anyone reading it will look for events that are not there.
- **Why it does not fire:** that deck ends near `p ~ 3-5 kPa`, far above either floor (`m_Pmin + m_Presidual` = 1.0201 with vanilla constants, 0.101 at the class defaults).
- **How to exercise it:** raise `-Pmin` above the deck's working pressure. Measured on the same deck with `-Presidual 0`: `-Pmin 1.5` -> **10 events**, `-Pmin 3.0` -> **10 events** (10 is the full throttle budget, after which a suppression notice prints and it goes silent). That is also the honest way to regression-test the diagnostic.
- **The wider lesson, which cost this session real time twice:** "I verified the warning fires" is only meaningful with the deck AND the leg AND the constants named. Two different agents reported "8 events" in this PR about two DIFFERENT diagnostics (this clamp, and `Elastic2Plastic`'s `M_c`-inflation notice, which really does fire 8x per leg on the 3D ramp deck — once per Gauss point). Number collisions between unrelated diagnostics are how a false claim survives review.
- **Learned:** 2026-08-27, closing the discrepancy between PR-2 commits 2 and 5.
