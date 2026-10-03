---
wp: LEGACY
title: "ID::getLocation (or a full constraint/domain iteration) inside a per-node / per-element / per-DOF-group loop is a per-rank quadratic — the \"scan-in-loop\" family"
legacy_seq: 203
---
### `ID::getLocation` (or a full constraint/domain iteration) inside a per-node / per-element / per-DOF-group loop is a per-rank quadratic — the "scan-in-loop" family
- **Bites:** five independent instances cost real wall: `TransformationConstraintHandler::handle()` (N^1.94, 14.4 s at 2.0 M np8 — element classification scanning the SP list per element-node), the `TransformationDOF_Group` SP-only ctor (swept EVERY domain SP per constrained node — invisible until the handler fix landed), `PlainHandler` per-node `getMPs()` AND `getEQs()` sweeps, and the `-4` fixup full-MP sweep in `DOF_Numberer` + `PlainNumberer` (both variants each). All MP/SP-count-driven: zero cost on unconstrained decks, quadratic on slab meshes (constrained area ∝ N) and tie-heavy decks.
- **Why:** `ID::getLocation` is a linear scan (`ID.cpp`); constraint iterators restart from scratch each call. O(outer) × O(scan) with both ∝ N/P.
- **Workaround/status:** ALL FIXED with one-pass hash/multimap indexes, order-preserving, byte-identity-gated ([#595](https://github.com/nmorabowen/OpenSees/pull/595), [#598](https://github.com/nmorabowen/OpenSees/pull/598); the parallel numberer's own copies in #592). Audit method that found them: `dc.*` profiler brackets + a fixed-np rung sweep (exponent), then a fixed-V np-sweep (per-rank vs global discrimination: per-rank quadratics FALL ~1/P², global-serial ones are np-invariant and mimic an Amdahl fraction in scaling studies). **Fix one, re-measure — the ctor sweep was invisible behind the handler scans.** *2026-07-22.*
