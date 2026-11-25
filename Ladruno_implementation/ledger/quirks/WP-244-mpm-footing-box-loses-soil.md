---
wp: WP-244
title: "An MPM rigid footing whose box ends at the footing base loses the soil before each node row joins"
pr: "#1136"
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

**Fix (WP-244, #1136).** Start the box's lower bound half a cell below the footing base
(`box <x0> -h/2 <B/2> 1.0e9 … track y`): a node row is prescribed once the base is nearer to it
than to the row above, so the soil–row gap stays within ±h/2. No C++ change. Before s = h/2 the
rows, and so the results, are bit-identical. The HS footing then runs to s = 2h with no step drop
above 4 % (55–92 % before). Residual: q steps up 2–16 % when a row joins (s/h = k + ½) and dips when a
particle row crosses it (k + ¾): 1.5–4 % for HS, 13–15 % for softening DP (L4). Results: `Ladruno_implementation/244_mpm_footing_track.md`.

**Extreme case (P1).** When the overlap is thin enough, the only prescribed row's nodal mass
falls under the mode-A mass cutoff, the row leaves the active set, and nothing is prescribed for
that step (a free block lags exactly one step; `mpm massCutoff 1e-14` removes it).

**Also learned (a record is tied to its BLAS build).** On the stagnating non-associated
Drucker–Prager ladder of 232j L4 r3, where the run gives up depends on the MKL the module links.
The branch's C++ on esmeralda reproduces the 232j record bit for bit with the default configure
(apt MKL 2020, sequential preload: stop at s = 0.094 m, 502 failed attempts). With spack oneMKL
2024 it stops at s = 0.010 m, on 232j's own merge commit too. The Windows build stops at 0.010 m
as well. Compare such a record only on its own build configuration. Within one build, the thread
count alone moves results by about 1e-13.
