---
wp: LEGACY
title: "Eigen members \"initialized\" with *= 0 keep NaN heap garbage (NaN*0 == NaN)"
legacy_seq: 298
---
### Eigen members "initialized" with `*= 0` keep NaN heap garbage (`NaN*0 == NaN`)
- **Bites:** `ASDPlasticMaterial3D`'s constructor zeroed its Trial/Commit stress/strain Eigen members with `*= 0` rather than `setZero()`. Freshly allocated heap storage is uninitialized, and `*= 0` is a no-op on NaN/Inf bit patterns (`NaN*0 == NaN`) — it only "zeroes" values that already happen to be finite. A fresh OS process gets zero-filled pages, so a standalone probe run in its own interpreter always passed; pytest's long-lived, churned heap recycles dirty blocks, so roughly 40% of full-suite runs picked up NaN in the shear slots of `CommitStress` at construction. Every yf comparison against NaN silently evaluates false — plain MC trips its Newton NaN guard and the analysis fails loudly, but the MCTC escalation chain (ADR-84) would instead classify the NaN-poisoned trial as TC-dominant, land on the apex, and return a CLEAN-LOOKING `T_eff·δ` on what should have been an ordinary compression path — exactly the garbage-into-plausible failure mode this feature exists to prevent.
- **Why:** `*=` on Eigen types is elementwise multiply-in-place; it is not a substitute for `setZero()`/`Zero()` when the storage's initial content is unknown.
- **Workaround/status (2026-08-12, PR #741):** fixed by switching the constructor to `setZero()` — benefits every ASDP material, not just MCTC. Before trusting any other constructor in the tree, grep it for the same `*= 0` idiom; the bug pattern is generic to Eigen-backed members and not specific to this class.
