---
wp: WP-232
title: "An HDF5 file written without periodic flushes is unreadable after a kill"
date: 2026-11-14
---
### An HDF5 file written without periodic flushes is unreadable after a kill

- **Bites:** a writer that creates an HDF5 file and appends to it without `H5Fflush` until it
  closes. A process killed mid-run (a timeout, a scheduler kill, Ctrl-C hard) leaves a file h5py
  cannot open ("bad object header version number"): the object headers and chunk indexes live in
  HDF5's metadata cache and were never written. Measured on the first `mpm output` writer
  (WP-232c, 742b9981c): killed after 1 s (19 MB) and 3 s (59 MB), both files unreadable; a graceful
  exit was fine.
- **Why:** without SWMR, the on-disk structure is consistent only at a flush or a close. Raw
  chunk data written after the last flush lands beyond the end-of-allocation the file records,
  so it is ignored on reopen, but only if a flush ever happened.
- **Between flushes a kill can CORRUPT, not just truncate.** Without SWMR, HDF5 may evict
  metadata and partly written deflated chunks to disk at any time; a kill then leaves headers that
  open but data that fail with "filter returned failure during read" (measured on build 5 with a
  10 s cadence, killed at 1 s and 3 s). A file is crash-consistent only at a flush.
- **Workaround:** `mpm output` flushes after every output by default (outputs are throttled by
  `-every` / `-dt`), lays the data out so a flushed output is never rewritten (one chunk per output
  for `ID` and `DATA`; uncompressed index datasets), and flushes at the end of every run, failed or
  not. A wall-clock cadence (`-flush s`, as `LadrunoRecorder -flush`, WP-164) trades that safety
  for speed. The full fix is HDF5 SWMR, deferred. Readers then apply the contract's crash
  rule for a partial last output (232a §2.5). `tests/test_ladruno_mpm_contracts_io.py::
  test_a_killed_run_leaves_a_readable_file` kills a writer mid-run.
- **"Last in the API" is not "last on disk".** Without SWMR, an eviction can put a later write
  (the output's `OFFSETS` commit entry) on disk before an earlier one (its `ID` / `DATA` chunks):
  a kill then leaves an output the crash rule counts whose ids are fill (all zero), measured at
  0.5 s and 1.5 s kills. Fix: a write barrier, `H5Fflush` after the output's rows and before its
  `OFFSETS` append (WP-232c follow-up), so a counted entry always has its rows on disk.
- **A flush barrier is not enough: metadata-cache EVICTION corrupts between flushes.** Without
  SWMR, HDF5's metadata cache may evict a dirty object header or chunk-index B-tree node at any
  time. It then reaches disk pointing at a node allocated past the end-of-allocation that the
  superblock recorded at the last flush, and a kill leaves a file whose reads fail with
  "addr overflow, addr = …, size = 2096, eoa = …" (the size of a rank-1 v1 B-tree node). Measured
  on build 7, with the barrier: 1 kill in 137 (at 0.8 s). Fix: the particle file's access property
  list turns eviction off (`H5Pget_mdc_config`; `evictions_enabled = 0`, which HDF5 accepts only
  with `incr_mode`, `flash_incr_mode` and `decr_mode` off; an 8 MiB initial cache; then
  `H5Pset_mdc_config`). Metadata then reaches disk only at the writer's explicit flushes, so a kill
  between flushes leaves the last flushed, consistent state; raw chunks written past it are
  unreferenced. The cost: the cache never shrinks, so it holds the file's metadata (a few KiB per
  output). Only the MPM particle writer does this; `LadrunoRecorder` keeps the default cache.
- Learned 2026-11-14 (WP-232c, the c2 review).
