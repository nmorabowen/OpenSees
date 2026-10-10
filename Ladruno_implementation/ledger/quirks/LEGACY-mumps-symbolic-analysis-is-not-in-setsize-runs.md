---
wp: LEGACY
title: "MUMPS: symbolic analysis is NOT in setSize (runs at first solve), and error −13 in \"substitution\" is plain out-of-memory"
legacy_seq: 207
---
### MUMPS: symbolic analysis is NOT in `setSize` (runs at first solve), and error −13 in "substitution" is plain out-of-memory
- **Bites:** (a) profiling `dc.setSize` on the implicit lane and expecting the MUMPS analysis there — it's in the first `solveCurrentStep` instead (`MumpsParallelSolver::setSize` only sets `needsSetSize`). (b) A 1.0 M-node (3.1 M-eq) LU does not fit a 64 GB box; the failure surfaces as `Error -13 returned in substitution dmumps()` per rank, easily misread as a numerical bug.
- **Workaround/status:** attribution: setSize = graph build + triplet fill (measured linear, [#604](https://github.com/nmorabowen/OpenSees/pull/604)); budget first-solve separately. Local implicit ceiling ≈ the 0.5 M rung. *2026-07-22.*
