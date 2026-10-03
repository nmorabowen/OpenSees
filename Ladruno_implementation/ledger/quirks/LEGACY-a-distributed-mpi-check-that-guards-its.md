---
wp: LEGACY
title: "A distributed MPI check that guards its collective leg on size == N silently tests NOTHING at other sizes"
legacy_seq: 230
---
### A distributed MPI check that guards its collective leg on `size == N` silently tests NOTHING at other sizes
- **Bites:** ADR-1000's P3 plan prescribed "run the checks at np=2, they must now pass" as the acceptance gate for a 2-rank MUMPS fix. Both distributed legs (`testDistributedMumps`, `checkDistributedFourRankFlow`) opened with `if (size != 4) return;`, so at np=2 they returned instantly and the check printed `passed` having exercised no distributed code at all. A green np=2 run would have "validated" the fix while proving nothing.
- **Workaround/status:** when a check's coverage depends on `MPI_Comm_size`, either make the fixture size-generic or make the skip **loud** (print what was skipped and why). The tiny 2x2 distributed fixture was size-generic all along — only the guard was wrong. Note the second-order trap: a fixture can be the right *size* and still be the wrong *shape* — the 2x2 one runs at np=2 but is far too small for MUMPS to make an ordering decision, so it could never have caught the bug the gate was written for. *2026-07-26 (ADR-1000 Part 0).*
