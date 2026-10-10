---
title: "WP-164 — Ladruno recorder write path: open handles, flush cadence, chunk plan, -compress, in-place envelope"
project: Ladruno
type: performance work package (follows WP-163)
status: "in progress — baseline measured; implementation under benchmark"
owner: nmora
related:
  - "[[163_ladruno_recorder_hardening_roadmap]] (the review: findings P1-P8; WP-164 = section 1.3)"
  - "[[03_ladruno_recorder]] (recorder design)"
tags: [recorder, hdf5, performance, scale, wp-164]
updated: 2026-10-03
---

# WP-164 — Ladruno recorder write path

> [!summary] The short version
> The WP-163 review found that the recorder's per-step HDF5 work, not the data volume,
> dominated small and long runs. Measured on the WP-163 build (`f150d4e50`), a 20 000-step
> explicit run with two tiny channels spent **0.9 s in the analysis and 88 s in the
> recorder**. WP-164 keeps the result layout readers depend on (same groups, same
> `[T x nIds x nComp]` DATA, same TIME/STEP) and changes only how it is written.

## 1. What changes

| Finding | Before | After |
|---|---|---|
| P2 | DATA/TIME/STEP reopened every step (chunk cache evicted → the partial chunk re-read, inflated, re-deflated) + `H5Fflush` every step | handles held open for the stage, a counter instead of extent queries, reused memory spaces; flush on a wall-clock cadence `-flush <s>` (default 10 s; 0 = every step) and at stage change / close |
| P3 | chunk `{ct, nIds, nComp}` — one chunk spans every id; one entity's history inflates the whole dataset | chunks ~1 MiB: small slabs keep all ids and stack up to 1024 steps; large slabs tile the id axis and stack ≤ 16 steps; chunk cache sized to a row of chunks so each chunk is deflated once |
| P4 | deflate 4 hard-coded | `-compress <0..9>` (0 = no filter); **default 1** (owner decision, 2026-10-03) |
| P1 | `-envelope` deleted and recreated every envelope group (+ attrs + COLUMN_MAP) every recorded step | datasets created once per stage, overwritten in place; COLUMN_MAP once; envelopes ≤ 8 MiB still rewritten every step (plain `H5Dwrite`, so a deck that exits without `wipe` keeps its latest extremes), larger ones on the `-flush` cadence; the ending stage is finalized before the stamp moves |
| P5 | `Domain::getNode(tag)` (a `std::map` walk) per node per channel per step | `Node*` resolved once per source (sources are rebuilt on every stamp change) |
| P7 | a non-reaction channel reset the reaction flag → repeated `calculateNodalReactions` | the last computed flag is remembered |

Crash-safety trade-off: between flushes, the newest output lives in HDF5's chunk cache. A hard
crash (not a normal exit — HDF5's atexit close still writes it out) loses at most the last
`-flush` seconds; `-flush 0` restores the per-step behaviour.

## 2. Benchmark (`Ladruno_scripts/ladruno_recorder_tests/bench/recorder_bench.py`)

Each case runs in a fresh interpreter with the recorder ON and OFF; *recorder* = wall_on − wall_off.
Same machine (a shared, loaded workstation — treat single runs as ±20 %).

| Case | Model | Recorder request |
|---|---|---|
| small_explicit | 2D 40-quad strip, central difference, 20 000 steps | `-R` 10 nodes `-N displacement -G energy` |
| envelope_medium | 24×24×8 stdBrick (4 608 el.), static, 40 steps | `-N displacement reactionForce -E stresses -envelope` |
| large_slab | 40×40×12 stdBrick (19 200 el.), static, 20 steps | `-N displacement -E stresses` |

### Baseline — WP-163 build `f150d4e50`

| Case | Analysis only | Recorder | File | Read one history |
|---|---|---|---|---|
| small_explicit | 0.89 s | **88.3 s** | 3.5 MB | 0.020 s |
| envelope_medium | 5.64 s | **2.91 s** | 8.4 MB | — |
| large_slab | 27.75 s | **6.97 s** | 120.4 MB | **0.565 s** (DATA [20, 19200, 48], chunk [1, 19200, 48]) |

### WP-164 build (same machine, same day)

| Case | Recorder before | Recorder after | Read one history | File |
|---|---|---|---|---|
| small_explicit | 88.3 s | **0.97 s (91×)** | 0.020 → 0.010 s | 3.5 → 3.5 MB |
| envelope_medium | 2.91 s | 1.27 s (median of 3) | — | 8.4 → 8.4 MB |
| large_slab | 6.97 s (1 run) | not resolvable by wall time (see below) | **0.565 → 0.018 s (31×)**, chunk `[4, 682, 48]` | 120.4 → 126.7 MB |

**large_slab write cost.** The analysis-only wall time of this 30 s solve moved between 23.5 s and
34.1 s from run to run on the shared box, so a few seconds of recorder cost cannot be resolved by wall
time (one WP-164 run even came out negative). It was also never expected to drop much: with a 12 MB
slab the old layout had one step per chunk, so each chunk was already deflated once. The levers for big
slabs are the read side (above) and the deflate level. CPU time (process time, robust to other load),
median of 2, WP-164 build:

| `-compress` | Recorder CPU | File | Read one history |
|---|---|---|---|
| 4 (the old hard-coded level) | 4.16 s | 126.7 MB | 0.018 s |
| **1 (new default)** | 3.38 s (−19 %) | 127.4 MB (+0.6 %) | 0.020 s |
| 0 (no filter) | 1.81 s (−56 %) | 167.6 MB (+32 %) | 0.005 s |

Smooth elastic data compresses about equally at levels 1 and 4; noisy nonlinear fields typically lose
5–10 % more at level 1. **The default is 1** (owner decision, 2026-10-03): files written before WP-164 used 4; `-compress 4` reproduces them, `-compress 0` trades file size for the least CPU.

**The envelope file did not shrink** (8.4 MB both builds): HDF5 reused the freed blocks of the old
delete/recreate within one session, as the review's blue team predicted. The gain there is the work
per step (no group walk, no group/attribute/COLUMN_MAP recreation); the file-growth gate is in
`tests/test_ladruno_recorder_write_path.py`.

## 3. Verification

- `tests/test_ladruno_recorder_write_path.py` (6): `-flush 0` and the 10 s cadence write identical
  DATA/TIME (STEP compared by increment); a > 1 MiB slab tiles the id axis and one history reads back;
  `-compress 0/1/9` sets the filter; the envelope file does not grow 3 → 30 steps (COLUMN_MAP present).
- WP-163's `tests/test_ladruno_recorder_hardening.py` (10) still passes.
- Not changed: group names, dataset names and shapes, attributes, FORMAT_VERSION; readers need nothing.
- Stale mirror noted: `Ladruno_scripts/zfp_benchmark/reencode_bench.py` copies the OLD 256 KiB chunk rule.
