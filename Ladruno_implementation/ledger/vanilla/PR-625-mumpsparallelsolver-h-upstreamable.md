---
wp: PR-625
title: "625 -- upstreamable-table row(s)"
pr: "#625"
files: ["`SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.{h,cpp}`", "`SRC/interpreter/OpenSeesCommands.cpp`"]
table: "upstreamable"
legacy_seq: [356, 357]
---
| `SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.{h,cpp}` | `// Ladruno ADR-75 P2`: expose MUMPS **Block Low-Rank** — new `icntl35`/`cntl7` members + ctor params (default 0/0.0 ⇒ byte-identical to stock). `ICNTL(35)`/`CNTL(7)` are applied at BOTH the analysis (job=1) and factorization paths (MUMPS uses ICNTL(35) when planning the assembly tree, so setting it only at factor time under-delivers). **`sendSelf`/`recvSelf` extended 2→3 ints + a `Vector(1)` for the double**: without shipping `icntl35`, rank 0 would factor BLR while every subordinate factored full-rank — an inconsistent distributed factorization that would produce wrong answers rather than an error. | [#625](https://github.com/nmorabowen/OpenSees/pull/625) |
| `SRC/interpreter/OpenSeesCommands.cpp` | `// Ladruno ADR-75 P2`: `system Mumps -BLR <eps>` (sets ICNTL35=1 + CNTL7=eps) plus raw `-ICNTL35`/`-CNTL7` escape hatches. BLR is an **approximate** factorization — opt-in, off by default, and must stay off byte-identical/oracle lanes. | [#625](https://github.com/nmorabowen/OpenSees/pull/625) |
