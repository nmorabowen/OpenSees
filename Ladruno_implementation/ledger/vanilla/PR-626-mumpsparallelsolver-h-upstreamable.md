---
wp: PR-626
title: "626 -- upstreamable-table row(s)"
pr: "#626"
files: ["`SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.{h,cpp}`, `SRC/interpreter/OpenSeesCommands.cpp`"]
table: "upstreamable"
legacy_seq: [364]
---
| `SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.{h,cpp}`, `SRC/interpreter/OpenSeesCommands.cpp` | `// Ladruno ADR-75 P2b`: `system Mumps -stats` dumps MUMPS `INFOG(9)/(21)/(22)`, `RINFOG(3)` and (under BLR) `RINFOG(14)/(15)` after each numeric factorization on rank 0 — BLR's whole justification is factor MEMORY, previously invisible from OpenSees. Handles MUMPS's negative-INFOG = −value×10⁶ convention; BLR fields printed raw + index-labelled rather than reinterpreted. **Also widened the `system Mumps` option loop from `> 1` to `> 0`:** its "every option takes a value" premise (an ADR-43 fix) broke for the bare `-stats` flag — a trailing `-stats` was silently ignored, which would have produced a false "BLR shows no stats" reading. | [#626](https://github.com/nmorabowen/OpenSees/pull/626) |
