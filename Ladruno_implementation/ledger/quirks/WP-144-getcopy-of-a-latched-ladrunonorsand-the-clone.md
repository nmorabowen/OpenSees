---
wp: WP-144
title: "getCopy of a LATCHED LadrunoNorSand: the clone must NOT inherit the refusal latch (WP-144 G2; the opposite of the \"work already done\" latch above)"
date: 2026-11-06
---
### getCopy of a LATCHED LadrunoNorSand: the clone must NOT inherit the refusal latch (WP-144 G2; the opposite of the "work already done" latch above)
- **Bites:** two different latches look alike. "nDMaterial that regularizes in place: getCopy must PROPAGATE the latch" is a WORK-DONE flag: a clone that forgets it redoes the work (double regularization). The LadrunoNorSand commit-refusal latch (WP-99 F7) is a REFUSAL flag bound to ONE integration point: a refused update reached `commitState` (the host element discarded the code), so that point refuses every trial and commit until `revertToStart()`.
  - A clone that inherited it would make a FRESH element refuse from its first `update()` with no refusal of its own to explain it.
  - The source's committed state is untouched by a refused trial (the commit is aborted and the trial restored from the committed state), so what a clone copies is a valid state: there is nothing for the clone to be protected from.
- **Rule:** decide per latch what it guards. A work-done latch is copied; a refusal latch is per instance and a clone starts clear. `copyFrom` therefore resets `latched`, `trialRefused`, the one-time warning flags and the refusal record (cause, count, `lastInfo`) and re-freezes the clone's trial at its committed state. `sendSelf/recvSelf` is different: it moves THE SAME point (a restart), so it keeps the latch.
- **Workaround/status:** FIXED (WP-144): `LadrunoNorSand::copyFrom`; `tests/test_ladruno_norsand.py::test_getcopy_of_a_latched_source_is_not_latched` (the prototype is latched through NDTest, a stdBrick is built from it, the cube must analyze and commit and report an all-zero `refusal` response; the same test fails on the previous behaviour).
