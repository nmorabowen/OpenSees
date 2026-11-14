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
  for speed. The full fix is HDF5 SWMR (below). Readers then apply the contract's crash
  rule for a partial last output (232a §2.5). `tests/test_ladruno_mpm_contracts_io.py::
  test_a_killed_run_leaves_a_readable_file` kills a writer mid-run.
- **"Last in the API" is not "last on disk".** Without SWMR, an eviction can put a later write
  (the output's `OFFSETS` commit entry) on disk before an earlier one (its `ID` / `DATA` chunks):
  a kill then leaves an output the crash rule counts whose ids are fill (all zero), measured at
  0.5 s and 1.5 s kills. Fix: a write barrier, `H5Fflush` after the output's rows and before its
  `OFFSETS` append (WP-232c follow-up), so a counted entry always has its rows on disk.
- **A flush barrier is not enough: a kill INSIDE a flush tears a chunk-index split.** Behind the
  barrier, kills still left unreadable files: "addr overflow, addr = …, size = 2096, eoa = …"
  (1 kill in 137 on build 7; 3 in 400 on build 8). The cause, pinned by a replay model that
  rebuilds the file at every byte-level kill point
  (`Ladruno_scripts/mpm_crash_tests/hdf5_crash_replay.py`, mirroring the writer's HDF5 call
  pattern in h5py, not the C++ itself): all failing kill points fall
  inside the BARRIER flush at the outputs where a dataset's v1 B-tree chunk index splits (64 and
  121; an `ID` node is 2096 B, a `DATA` node 2616 B). A split of the rightmost node keeps 57 of
  64 entries and moves the newest 7 (committed outputs) into a new node at the end of the file; a
  root split moves every entry (HDF5 1.14.6 `H5B.c:395-418`). `H5Fflush` writes the user ring
  (the in-place parent) before the superblock ring (the new EOA), so a kill between them leaves a
  parent pointing past the stored EOA. Metadata-cache eviction plays no part: the model fails
  identically with eviction on or off, and so did a build with eviction disabled (build 8; that
  change was reverted). Larger chunks only make the splits rarer, never absent.
- **Fix: SWMR.** The particle file uses the latest format (`H5Pset_libver_bounds(LATEST, LATEST)`;
  extensible-array chunk indexes) and the writer arms `H5Fstart_swmr_write` once a stage's layout
  is complete and flushed; it then only appends. A new stage closes the file, reopens it, lays the
  stage out and arms SWMR again (no object or attribute may be created while SWMR is armed). The
  barrier stays: SWMR orders each dataset's writes, not `OFFSETS` against `DATA`. In the replay
  model, SWMR recovers every completed output, bit-exact, at all 17,439 kill points of a
  600-output run. On the real binary (build 9, cbe1b7149) the kill stress
  (`Ladruno_scripts/mpm_crash_tests/kill_stress.py`) left 0 unreadable files in 400 kills
  (build 8: 3 in 400).
- **A killed SWMR file needs an SWMR-read open.** Its consistency flags stay set, so a plain
  read-only open fails ("file is already open for write"). Open with `swmr=True`; never
  `h5clear` first (the stored EOA lags; only an SWMR-read open skips that check). After a kill the
  SWMR-read open needs no `HDF5_USE_FILE_LOCKING=FALSE` on Windows (measured, h5py 3.16 /
  HDF5 2.0.0). Reader rule: 232a §2.5, `tests/mpm_contracts.py::open_particle_file`.
- Learned 2026-11-14 (WP-232c, the c2 review).
