---
title: "WP-251 — recorder ladruno: off-thread chunk compression and threaded reactions"
project: Ladruno
tags:
  - recorder
  - performance
---

# WP-251 — recorder ladruno: off-thread compression, threaded reactions

## Why

Started as a study of how FEniCS/DOLFINx records results "without a performance drop".
DOLFINx's writers are plain (VTXWriter copies each field and does a sync ADIOS2 `Put`
per step; XDMF writes one HDF5 dataset per step). What keeps its output cheap is the
data model: a field is one contiguous per-rank array, nothing is derived at write time,
and ADIOS2 does not compress by default.

Profiling `recorder ladruno` with the stack profiler (`lr.collect` / `lr.write` /
`lr.reactions` scopes added here) showed the per-element collect walk is NOT the cost.
8192 stdBricks, Newmark + `Linear -factorOnce`, Pardiso, ~83–103 ms/step without a
recorder; recorder ms/step before WP-251:

| recorded | total | collect | write | reactions | flush |
|---|---|---|---|---|---|
| displacement | 3.7 | 0.5 | 3.0 | — | 0.1 |
| reactionForce | 44 | 0.3 | 3.2 | 41 | 0.1 |
| stresses | 46 | 5.4 | 38 | — | 2.0 |
| disp+vel+acc+reac+stress+strain | 145 | 10 | 89 | 42 | 4 |

Two costs dominate: **shuffle+deflate on the solver thread** (~60 %) and
**`Domain::calculateNodalReactions`**, a full resisting-force sweep per recorded step
(~28 %).

## What

1. **Off-thread compression** (`Ladruno_ChunkCompress.{h,cpp}`, `StreamingSink`).
   The compressed layout already stores one step per chunk (WP-233), so a step is a set
   of complete chunks. Each chunk (one id tile of one step) is a `ChunkJob`: the worker
   narrows to f32 if asked, applies HDF5's byte shuffle and `compress2` at the dataset's
   level — the exact bytes the declared filter pipeline would store. The solver thread
   writes the finished bytes in step order with `H5Dwrite_chunk` plus the row's
   `TIME`/`STEP`. **Only the solver thread calls HDF5** (it is not built thread-safe).
   - `-compressThreads <n>`: `-1` auto (`hw/4`, 1..4), `0` = the old in-line path.
     Serialized in `sendSelf`/`recvSelf`.
   - In flight: at most 8 rows per result; `accept()` writes every finished row first.
   - Drained (blocking) before every `flushFile()`, every layout step (still in SWMR),
     `releaseHandles`, `closeHandles`, `finalize`.
   - **Process exit without `wipe`:** HDF5's atexit close cannot see queued rows, so a
     hook registered after HDF5 initialised (atexit is LIFO) compresses them on the
     exiting thread and writes them. It never waits on the pool: in the `.pyd` Windows
     has already terminated the workers when the CRT runs atexit. The pool is leaked
     and its threads detached for the same reason.
2. **Threaded reactions** (`Domain::ladrunoThreadedReactions`). Under `ladrunoThreads > 1`
   the flag-0/1 loop computes the forces of `Element::ladrunoThreadSafeReaction()`
   elements in parallel (default = `ladrunoThreadSafeForm()`), then scatters EVERY
   element in serial order with the base rule — bit-identical sums. Rayleigh (flag 2),
   MPI builds, partitioned domains and the deep profiler stay serial. Elements that
   override `addResistingForceToNodalReaction` (ASDAbsorbingBoundary2D/3D) opt out; the
   `reaction-scatter` quirk rule fails a future override that forgets to.

## Result (same model, recorder ms/step)

| recorded | before | WP-251 | step before → after (no recorder) |
|---|---|---|---|
| stresses, serial | 38.0 | 5.1 | 121 → 91 (83) |
| all six, serial | 123.5 | 44.0 | 207 → 130 (83) |
| all six, `LADRUNO_THREADS=8` | 98.6 | 18.8 | 141 → 62 (41) |
| reactionForce, 8 threads | 33 | 8.7 | — |

File sizes unchanged (the bytes are identical).

## Verification

- `tests/test_ladruno_recorder_compress_threads.py`: stored chunk bytes + decoded arrays
  equal the in-line path (f64, f32), flush cadence / auto threads, exit without `wipe`,
  a killed async run reads by the WP-233 rule with every committed row equal to a full
  run's, threaded reactions bit-identical, option parse error.
- Existing recorder / threading / absorbing / overlay suites: 503 passed on the WP-251
  build; the SWMR kill and live-reader tests included.

## Left out (measured, not worth it now)

- Threading the collect walk: 7 ms/step here, and many elements' `getResponse` fill
  class-static buffers (stdBrick's stresses), so it needs a per-class audit.
- Reactions at constrained nodes only: exact for a node subset, but `-N reactionForce`
  records every node by default; the threaded loop helps every deck.
