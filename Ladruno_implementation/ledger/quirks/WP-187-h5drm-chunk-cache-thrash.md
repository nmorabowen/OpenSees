---
wp: WP-187
title: "H5DRM serial / free-field runs dominated by re-inflating gzip chunks (fixed by WP-187)"
date: 2026-10-11
pr: "#976"
status: "fixed-by WP-187"
---
### A serial or few-rank H5DRM run slows 4-5x once the motion starts; logs grow by MBs of `H5DRM(0) : t = ...` lines (fixed by WP-187)
- **Bites:** any `pattern H5DRM` on gzip-chunked files (ShakerMaker writes them chunked + deflate; the Paper 1 Legacy files are 192 x 125 chunks, deflate 4), most visibly a serial free-field run: S5 job 689 went from 0.59 s/step before the motion to 2.9 s/step in it (engstruct PROFILE_4D section 4.5). In a 32-rank 4D run it is hidden behind MUMPS (soil ranks wait anyway).
- **Why:** `drm_direct_read` read four 3x1 hyperslabs per node and step through HDF5's default 1 MiB chunk cache (5 chunks of 192 KB). Every read whose chunk was evicted inflated the whole 192 x 125 chunk again: about 50x the bytes needed. The constant part of the record compresses to nothing, so the waste appears only when the motion arrives. Rank 0 also printed 2-3 lines per step.
- **Fix:** WP-187 (`ledger/vanilla/WP-187-*`): a chunk cache sized per rank and a time block of chunk length served from memory (bit-identical); per-step lines behind `-verbose`. Measured on a synthetic Legacy-layout file: 309 -> 6.5 ms/step for the read itself. `-prefetch 0 -chunkCacheMB 0` reproduces the old path for A/B checks.
- **Check:** with `-verbose` the pattern prints `WP-187 read time ... ms/step` when it is destroyed (`wipe`).
