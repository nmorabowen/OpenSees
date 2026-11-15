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
- **The stage layout is a window SWMR does not cover; a kill there leaves a file NOTHING opens.**
  The layout (at the start of the run, and at each `mpm reset`, after the writer closes and
  reopens the file read-write) is written without SWMR. A read-write open of a version-3 file
  always writes WRITE_ACCESS (0x01) to the superblock (HDF5 1.14.6 `H5Fint.c:2164-2179`; the flag
  is set whenever the driver has a lock method, `H5Fint.c:1977-1979`, and sec2 has one). A plain
  open refuses 0x01 ("already open for write", or "truncated file" first when the file is shorter
  than its stored EOA), and an SWMR-read open accepts only 0x05 or 0x00 (`H5Fint.c:2196-2208`:
  "not already open for SWMR writing"). Measured on build 9 with a reset every 100 steps
  (`Ladruno_scripts/mpm_crash_tests/reset_stress.py`): 21 of 400 kills (5.25 %) were refused by
  both opens. Recovery, on a COPY (`tests/mpm_contracts.py::recover_particle_file`): flags 0x00,
  stored EOA (superblock offset 28) = max(stored EOA, file size), zero-pad the copy up to the
  stored EOA when it is shorter (1 of the 21), re-checksum the superblock (lookup3 over bytes
  0–43, at 44), then open plainly. 21 of 21 recovered every earlier stage; in 2 the stage being
  laid out did not open ("bad object header version number"), and it held no committed output.
  Real runs reset rarely, but a hit locks out every earlier stage until the copy is recovered.
- Learned 2026-11-14 (WP-232c, the c2 review).

#### Follow-up (WP-232h, #1085, PR #1088): the reset window can lose the ROOT; a locking reader breaks the writer

- **An unsized root header turns a torn stage into a lost file.** At each `mpm reset` the root
  gains a `MODEL_STAGE` link inside the layout written without SWMR. By default a group goes
  dense after 8 links or grows a continuation chunk, so adding the link is several in-place
  writes in no causal order, and a kill between them leaves `/` unreadable ("addr overflow" or
  "incorrect metadata checksum"): every earlier stage lost, and reader rule step 3 cannot
  recover it. Measured by `hdf5_crash_replay.py --outputs 240 --resets 10 --no-presize` (every
  kill point, 24 stages): 634 of 9,376 kill points lost the root. Fix (as WP-233 did for the
  recorder): `Ladruno_Hdf5.h::plist::presizeLinks(fcpl, 1024, 20)` on the particle file's fcpl
  (links compact for good, header sized up front), so adding a stage is ONE in-place header
  write: 0 of 9,249. No other group of the particle file gains a link after its stage's layout.
  The stored fcpl does not carry the setting (`H5Fget_create_plist` reads back 8/6): check the
  ROOT GROUP's gcpl, or that the root's storage is compact past 8 links.
- **The last stage can be torn after a PLAIN open too.** sec2 writes WRITE_ACCESS into the
  superblock lazily, so a file killed inside a layout can still carry flags 0x00. The writer now
  flushes right after the reset's read-write reopen, but a reader must still treat a torn LAST
  stage as holding no output after a plain open as after step 3, and refuse it after the
  SWMR-read open, which only a file killed with SWMR armed (no layout running) needs
  (`mpm_contracts.particle_stage_report`).
  A checker that opens `f[name]` outside its `try` reports such a file as unrecoverable: that
  was `reset_stress.py` (2 false "RECOVERY_FAILED" in 400 before the fix of the checker).
- **A reader with HDF5 locking on breaks the SWMR writer** (as found for the recorder, #1084):
  measured on the pre-fix pyd (bf7ad4ee8), a reader polling with h5py's default locking made the
  writer's reset reopen fail ("cannot reopen ... for a new stage") and the run abort. Fix: HDF5
  locking off for the particle file (`h5::file::accessList`) and the guard lock
  (`SRC/recorder/Ladruno_GuardLock.h`, shared with the recorder) held from `H5Fcreate` to after
  the final close; the reset's reopen retries for 60 s. A locking reader is refused
  (`WriterBusyError`) and the finished file equals an unread run
  (`test_a_live_reader_does_not_disturb_the_writer`).
- **A created file is never `H5Fcreate`'d again.** Before WP-232h a failed reset reopen left the
  writer with no file, and the next stage called `createFile` (`H5F_ACC_TRUNC`): every earlier
  stage truncated away. The writer now retries the reopen instead.
- Real binary (1d64d529e): `kill_stress.py` 0 bad in 400; `reset_stress.py` 0 bad in 400 with a reset every 100 steps
  (up to 138 stages; 16 refused by both opens, 16 recovered, 2 with the last stage torn) and 0
  in 400 with a reset every 10 (up to 400 stages; 43 refused, 43 recovered).
- Learned 2026-11-15 (WP-232h).
