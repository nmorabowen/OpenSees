---
wp: ADR-92
title: "Growing a MaterialResponse vector needs TWO edits, and a mismatch is SILENT (ADR-92 P2-9, 2026-09-07)"
date: 2026-09-07
legacy_seq: 426
---
## Growing a `MaterialResponse` vector needs TWO edits, and a mismatch is SILENT (ADR-92 P2-9, 2026-09-07)

The `setResponse` / `getResponse` idiom this codebase uses everywhere allocates a
function-local `static Vector` at BOTH ends — `setResponse` builds the recorder's
probe (`static Vector probe4g(6); return new MaterialResponse(this, ID, probe4g);`)
and `getResponse` builds the value it fills (`static Vector out4g(6); ...
matInformation.setVector(out4g);`). The two live ~130 lines apart in
`LadrunoSANISAND.cpp` and nothing ties their sizes together.

Growing a census by one slot (here `implexGuards` 6 -> 7 for P2-9's "backed off"
counter) means editing BOTH. Change only `getResponse` and the extra slot never
reaches a recorder — `Information::setVector` copies into a `Vector` the
`MaterialResponse` sized from the probe, so the value is dropped with no warning
and a harness reads a short vector or a stale one. Change only `setResponse` and
the recorder allocates a slot that is never written. Neither is a compile error
and neither prints anything.

Rule: when you add a slot, grep the response id and fix every `static Vector`
that mentions it in the same file, then update the slot table in the guide
(`LadrunoSANISAND_implex_guide.md` §6) in the same commit — the table IS the
contract the TIMs harnesses read by index.
