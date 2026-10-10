---
wp: WP-163
title: "A recorder that rebuilds on the domain-change stamp must release its cached Response* BEFORE writing the new stage — a failed rebuild otherwise dangles (WP-163)"
date: 2026-10-03
---
### A recorder that rebuilds on the domain-change stamp must release its cached `Response*` BEFORE writing the new stage — a failed rebuild otherwise dangles (WP-163)
- **Bites:** `LadrunoRecorder::record()` commits the new stamp, then `writeModel()`; any writer that failed
  (e.g. every node of the `-R` region removed → "no nodes to write") returned before `clearSources()`, so the
  old element `Response*` (wrapping deleted elements) survived, and the next commit (same stamp → no rebuild)
  called `getResponse()` on them: access violation 0xC0000005 (reproduced). `Domain::commit()` ignores the
  recorder's -1, so nothing else stops it.
- **Workaround/status:** fixed (WP-163 R2): release sources first, latch `stage_failed` until the next stamp
  change. Note `Domain::removeLoadPattern` bumps the stamp only when the pattern owns SPs, so a recorder
  caching a `LoadPattern*` must look it up by tag every step (the overlay source already does).
