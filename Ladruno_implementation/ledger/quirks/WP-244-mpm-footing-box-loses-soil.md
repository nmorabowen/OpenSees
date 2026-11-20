---
wp: WP-244
title: "An MPM rigid footing whose box ends at the footing base loses the soil before each node row joins"
pr: "#TBD"
date: 2026-11-20
---
### MPM footing pressure collapses at s/h ≈ 0.75 (q 868 → 38 kPa), then rises like a fresh contact

**Symptom.** A rough rigid footing written as `mpm gridvel footing box <x0> 0.0 <B/2> 1.0e9 v 0 -v
start t0 track y` (the 232f / 232j L4 / 232k S3 footings) gives an irregular q increment at
s/h ≈ 0.25 and a collapse from s/h ≈ 0.75 until s = h, at every resolution (HS r3: 868 → 38 kPa
at s/h = 0.768; r2: 1690 → 228 kPa at s/h ≈ 0.76–0.97). The location scales with h, the soil
model has no softening and the FE curve is monotone: it is a discretisation artifact (#1112).

**Why.** The grid is fixed. `track y` lowers the box's lower bound to −s, so the footing
prescribes every node row at or above its base, and the next row joins only when s reaches it.
With ppc 2 the soil's top particle row sits at −s − h/4, so the gap between the lowest
prescribed row and the soil grows to a full cell; from s/h ≈ 0.75 the top cell holds no
particle centre, and the soil next to the footing interpolates almost entirely from free nodes.

**Not the fix: `translate`.** WP-232k's `translate y` moves both bounds, but these boxes are
unbounded above (`hi_y = 1e9`), so `translate y` and `track y` select the same rows. Translating
a zero-thickness box (`lo_y = hi_y`) is worse: it holds a node row only when s is a multiple of h.

**Fix.** Measured in WP-244 (#TBD); see the follow-up below.
