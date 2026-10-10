---
wp: LEGACY
title: "Amendment to the collective-getB() guard: the predicate is rank-uniform because setSize() GLOBALIZES size, which is a stronger guarantee than \"no rank has size…"
legacy_seq: 317
---
### Amendment to the collective-`getB()` guard: the predicate is rank-uniform because `setSize()` GLOBALIZES `size`, which is a stronger guarantee than "no rank has sized yet"
- **Why this needs saying:** the entry `printA` / `printB` KILL the interpreter whenever the
  SOE was never sized (earlier in this file) justifies the early return in the five
  collective `getB()` overrides with "the wrappers are allocated in `setSize()`, reached only
  through the collective `domainChanged()`, so in the unsized state NO rank has them". That is
  true and it covers the reported bug, but it is an argument about the **never-sized** state
  only -- and it is not the whole guarantee. There is a second way to reach a null wrapper
  that the sentence does not reach: a `setSize()` that RAN and still left the wrappers null.
- **The mechanism, verified in all five classes.** `setSize()` is ITSELF a collective -- P0 and
  the subprocess SOEs exchange `sendID`/`recvID` -- and it **globalizes `size`**: P0 reduces
  `maxVertexTag` across every rank and broadcasts the result back, so all ranks leave
  `setSize()` with the same value. The wrapper-creation block is then gated on that
  rank-identical `size != oldSize` (`size` first, `oldSize` captured on entry). So "the
  wrappers are null" is a function of a value **every rank agrees on** -- and the ranks cannot
  diverge on the way in either, because they must rendezvous inside `setSize()` to get there
  at all. That is what makes the predicate rank-uniform, and it is checkable per class rather
  than inferred from the call graph.
- **What the stronger form buys, concretely: the ZERO-FREE-EQUATION door.** On a model whose
  every DOF is fixed or `sp`-prescribed the graph is empty, `size` comes out 0, `size !=
  oldSize` is FALSE, and a **completely successful** `setSize()` leaves the wrappers at their
  constructor nulls. Every rank HAS run `setSize()`, so the "no rank has them yet" phrasing
  does not apply -- but the guard is still symmetric, because `size` is 0 on every rank by
  construction. This is the same door as the `FullGenLinSOE::getX - vectX == 0` entry
  elsewhere in this ledger, arriving at the parallel classes; the guard covers it, and it is
  worth knowing that it does.
- **Rule, restated:** when you bail out of a collective, do not justify it with "this state
  cannot happen on one rank only" -- justify it by naming the **rank-identical value** the
  predicate is a function of. If you cannot name one, the guard is not safe yet.
- **Measured (2026-08-18), separately from the sibling probe:** `mpiexec -n 2` on
  `dist/openseesmp`, asserting on the launcher EXIT CODE. Pre-fix `printB` on an unsized SOE
  under `system Mumps` kills rank 1 with `c0000005` and rank 0 follows at `-1`. Post-fix
  **both** ranks emit the warning and the job exits 0 in ~0.1 s against a 120 s timeout -- the
  early return is demonstrably taken on both ranks, and nothing hangs. Two controls made the
  verdict mean something: a POSITIVE control (a case that actually analyzes -- guard silent,
  collective still merges a real `B`, `len=6` both ranks) so "returns empty" cannot pass
  everything, and a NEGATIVE control (rebuild with a single guarded file reverted -- exactly
  its rows crash `c0000005` again while the others stay green). `system ParallelProfileSPD`
  gives `DistributedProfileSPDLinSOE` a second runtime-gated door from the same build; the
  remaining three are constructed only under `_PARALLEL_PROCESSING` from the Tcl `OpenSeesSP`
  binary and stay compile-checked.
