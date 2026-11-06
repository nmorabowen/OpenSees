---
wp: WP-142
title: "PlateRebarMaterial::recvSelf rebuilt the bar direction with a truncated degree constant — a received off-axis bar (MP/SP, database restore) was not the bar tha…"
legacy_seq: 524
---
### `PlateRebarMaterial::recvSelf` rebuilt the bar direction with a truncated degree constant — a received off-axis bar (MP/SP, database restore) was not the bar that was sent (WP-142)
The constructor sets `rang = ang * 4.0 * asin(1.0)/360.0`; `recvSelf` used `rang = angle * 0.0174532925` (π/180 cut after 10 digits, ~1e-9 relative). Angles 0 and 90 take shortcut branches and never read `c, s`, so only off-axis bars were affected: after a `database` save/restore, or in any OpenSeesSP/MP partition that receives its elements through `recvSelf`, the bar strain `ε11c² + ε22s² + γ12cs` differed ~1e-10 relative from the uninterrupted run from the first step on. Physically negligible, but it breaks every bit-identity oracle across a restart and made SP/MP runs differ from serial ones.
- **Rule:** a class that recomputes derived data in `recvSelf` must call the SAME code the constructor does — factor it into one helper; never retype a constant. Gate any restore path with a save → restore → one more step bit-identity test on a case that actually exercises the derived data (here: off-axis angles).
- **Workaround/status:** ✅ FIXED in vanilla `PlateRebarMaterial.cpp` (WP-142, owner-approved upstream fix): both call `plateRebarCosines`. Gated by `test_database_roundtrip_off_axis_bar_is_bit_identical` (30° and 17°); the old constant fails it. *2026-09-28.*
