---
wp: WP-164
title: "Closing an HDF5 chunked dataset (or H5Fflush) every step makes a partially filled compressed chunk be deflated, written, re-read and re-inflated EVERY step (WP…"
date: 2026-10-03
---
### Closing an HDF5 chunked dataset (or `H5Fflush`) every step makes a partially filled compressed chunk be deflated, written, re-read and re-inflated EVERY step (WP-164)
- **Bites:** the Ladruno `StreamingSink` reopened DATA/TIME/STEP on each `accept()` and the recorder called
  `H5Fflush` after every recorded step. The chunk cache belongs to the open DATASET, so `H5Dclose` writes the
  dirty partial chunk (shuffle + deflate) and drops it, and the next `H5Dopen2` + partial write has to read and
  inflate it again. A chunk that stacks `ct` steps (up to 1024 for a small slab) was compressed ~`ct` times, and
  a recompressed chunk grows, so HDF5 also reallocates file space. Measured (WP-164 bench): a 20 000-step
  explicit run with two tiny channels — analysis 0.9 s, recorder **88 s**.
- **Why:** HDF5 caches raw chunks per open dataset (`H5Pset_chunk_cache` on the dataset ACCESS list, default
  1 MiB); `H5Fflush` writes every dirty cached chunk through the filter pipeline even when it is incomplete.
- **Workaround/status:** WP-164: keep the handles open for the stage, size the chunk cache to a row of chunks,
  flush on a wall-clock cadence (`-flush <s>`, default 10 s) and at stage end / close → recorder 0.97 s on the
  same deck. Holding handles means the sinks must close them before `H5Fclose` (the recorder deletes sinks
  first) and silently at process exit (HDF5's atexit may have closed the IDs: `H5E_BEGIN_TRY`).
