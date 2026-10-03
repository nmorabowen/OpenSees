---
wp: WP-164
title: "A chunk that spans every id makes ONE entity's time history read (and inflate) the whole dataset (WP-164)"
date: 2026-10-03
---
### A chunk that spans every id makes ONE entity's time history read (and inflate) the whole dataset (WP-164)
- **Bites:** `[T x nIds x nComp]` chunked `{ct, nIds, nComp}`: `data[:, k, :]` touches every chunk, i.e.
  decompresses all of DATA (0.565 s for one element of a 19 200-element, 20-step stress dataset; ~38 TB of
  inflate for a 10 M-element, 1e4-step run).
- **Workaround/status:** WP-164 chunk plan: ~1 MiB chunks; a slab above 1 MiB tiles the id axis and stacks
  ≤ 16 steps (`[4, 682, 48]` there → 0.017 s). Readers need no change (chunking is transparent). The
  `Ladruno_scripts/zfp_benchmark/reencode_bench.py` helper mirrors the OLD 256 KiB rule — re-read the chunk
  shape from the file rather than mirroring the writer's rule.
