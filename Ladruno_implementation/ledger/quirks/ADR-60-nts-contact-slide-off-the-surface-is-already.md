---
wp: ADR-60
title: "NTS contact slide-off-the-surface is already safe — no explicit detection needed (ADR-60 R5)"
legacy_seq: 141
---
## NTS contact slide-off-the-surface is already safe — no explicit detection needed (ADR-60 R5)

A slave that slides clean off the master (out of every segment's parametric domain) needs **no** special
"slide-off" code: `LadrunoContactProjection::evalSegment` gates on penetrating-AND-in-bounds, so an
out-of-bounds slave yields zero force; `project()` returns the out-of-bounds parametric coords
**UNclamped** (no edge-clamp that could hold spurious traction); the migration trigger + friction-slot GC
drop the now-stale adapter within a bounded window; and D4 fresh-slot re-engagement keeps any later
re-pairing traction-continuous. Empirically (`Ladruno_scripts/_probe_r5_slideoff.py` →
`tests/test_adr60_reemit_p4_slideoff.py`): a frictional slave flung off a finite strip's end departs with
force→0, falls freely, and retains its tangential velocity. CAVEAT for diagnostics: the `ladrunoContactForce`
(B3) snapshot is only refreshed on a re-handle (cleared in `frictionGCBegin`), so on the **frozen**
non-`-reemit` path it can report a STALE force after the slave has left contact; `-reemit` clears it each
re-emit so the readout is live there.
