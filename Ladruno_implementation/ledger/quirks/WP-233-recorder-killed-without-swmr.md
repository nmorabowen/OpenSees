---
wp: WP-233
title: "A killed recorder ladruno run left an unreadable or silently wrong .ladruno"
pr: "#1084"
date: 2026-11-14
---
### A killed `recorder ladruno` run left an unreadable or silently wrong `.ladruno` (the `-flush` promise was false)

- **Bites:** any `recorder ladruno` run that is killed (TerminateProcess, a scheduler
  kill, a hard Ctrl-C) rather than closed. Before WP-233 the file could fail to read
  ("filter returned failure during read"), losing the whole unfinished chunk (up to 1,024
  steps, often the whole run). Or it read **silently wrong**: zeros for the steps since the
  last flush, or another result's bytes. This held even though every earlier flush had
  completed. The documented "a hard crash loses at most the last `-flush` seconds" was
  false. Issue #1082. It is the same failure class as
  `WP-232-hdf5-unflushed-file-unreadable-after-kill` (the MPM particle file), seen a second
  time.
- **Why (three mechanisms in one flush; pinned by a byte-level replay):**
  1. Every flush deflated the partially filled `DATA` chunk again. It grew, so HDF5 moved
     it and handed its old space to the next chunk written, while the on-disk chunk index
     still pointed there. This happened at 72–79 % of the kill points inside a flush with
     the default `-compress 1`.
  2. A new chunk's index entry, the dataset extent and the end-of-allocation land in
     separate writes. The new rows then read as fill zeros, or as "addr overflow". This
     fires even with `-compress 0`.
  3. The v1 B-tree chunk-index split of the MPM quirk.

  Large models also had a window *between* flushes, when a completed row of chunks was
  evicted and recompressed.
- **SWMR alone bloats the file.** Under SWMR, HDF5 never reuses the space a moved chunk
  frees. With the WP-164 layout (deflated partial chunks), a 300-step run at `-flush 0`
  grew from 2.6 MB to **206 MB**. HDF5's own `H5D_CHUNK_DONT_FILTER_PARTIAL_CHUNKS` leaks
  too (1.75–2x), because it deflates the edge chunk on completion and frees the unfiltered
  copy. The fix is a layout whose chunks never move:
  - `-compress N > 0` (default 1): one step per chunk, compressed and written once. The
    file size is the same as WP-164's, but reading one id's history is 2–15x slower;
  - `-compress 0`: unfiltered ~64 KiB chunks, rewritten in place. Fastest, but the file is
    1.13x (real nonlinear runs) to 1.7x (a static slab) larger.

  An unfiltered chunk is allocated in full, so a 1 MiB target made a short run with 200
  element buckets 21x its payload. That is why the unfiltered target is 64 KiB
  (`recorder_swmr_layout_bench.py`).
- **The layout window is not SWMR, and an unsized parent tears.** SWMR forbids creating
  objects, so a new stage, the `KIND` rewrite and a modes snapshot close the file, reopen it
  read-write and lay out without SWMR. A plain `H5Fflush` writes metadata in no causal
  order. In the replay, the root group's header (with the new `MODEL_STAGE` link and a
  continuation message) reached the disk before the chunk it references, and a kill there
  left the **root**, and every earlier stage, unreadable. Likewise, a `KIND` delete +
  create left a stage with no `KIND`, and a modes snapshot rewrote `ON_NODES`' links. Fix:
  - groups that gain links during a run (the root; a modes result's `DATA`) keep compact
    links in a header presized for 256 stages / 1,024 snapshots
    (`Ladruno_Hdf5.h::plist::presizeLinks`), so adding one is a single in-place write;
  - the provisional `KIND` is a fixed-length 9-byte string, rewritten in place;
  - the modes groups are laid out with the stage.

  A kill in a layout step now leaves at most the stage being laid out unreadable, and it
  holds no committed step. Beyond 256 stages a new header chunk is needed again; that is a
  residual, documented in schema §7b.2.
- **The MPM particle writer has the same unsized root.** It creates `MODEL_STAGE[n]` at
  every `mpm reset` without presizing (`LadrunoMPMIO.cpp`). Its 21-kill sample never hit
  the window; the issue is filed separately.
- **Reader:** a killed file needs schema §7b.2:
  - a plain open, else SWMR-read (flags 0x05), else repair a copy (flags 0x01: killed in a
    layout step);
  - then cut each stage to T' = min over its streaming results of min(len `DATA`, `TIME`,
    `STEP`). `DATA`/`TIME`/`STEP` are separate datasets, and a kill leaves them ragged at
    ~40 % of the in-flush kill points.

  `ladruno_format.open_ladruno` / `repair_copy` / `committed_steps`; `validate()` refuses a
  killed file.
- **Evidence:**
  - `recorder_crash_replay.py` checks every byte-level kill point of an h5py mirror of the
    writer;
  - `recorder_binary_kill_stress.py` kills the real binary per lifecycle window;
  - the numbers are on #1082 / PR #1084.
- Learned 2026-11-14 (WP-233).
