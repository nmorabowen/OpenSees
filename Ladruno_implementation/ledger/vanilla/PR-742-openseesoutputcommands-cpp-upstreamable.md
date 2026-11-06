---
wp: PR-742
title: "742 -- upstreamable-table row(s)"
pr: "#742"
files: ["`SRC/interpreter/OpenSeesOutputCommands.cpp`"]
table: "upstreamable"
legacy_seq: [427]
---
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | `// Ladruno` (ADR-78 P4 finding 1 — closes the defect the row above opened): the three contact DECLARATION verbs each become a two-line wrapper (`OPS_LadrunoContactSurface` / `OPS_LadrunoContact` / `OPS_LadrunoContactPlane`) over their unchanged bodies (now `static ladrunoContact*Impl()`), routing any `<0` result through `ladrunoContactFatal()`. The declaration-time refusal above is a rank-local parser `return -1` with no MPI awareness: the refusing rank aborts its script, every peer blocks in the next collective, and the job HANGS — measured on the preserved P0 mutation deck `contact_parallel/mp_noghost.tcl`, which P1 measured at a 1.1 s teardown and which hung >30 s / 45 s (once a nondeterministic 17 s hydra reap) after #737 moved the check out of `handle()`. Wrapping the RESULT rather than editing the ~75 individual `return -1` sites is deliberate: it also covers the refusals that propagate out of `addSurface`/`addContact`/`addMortarContact`/`addRigidPlane` (unknown surface tag, duplicate tag — the other partition-dependent class), and a check added later inherits the teardown instead of silently reopening this defect. Serial and `mpiexec -n 1` byte-identical (INV-5): `ladrunoContactFatal()` returns −1 untouched at np ≤ 1. The QUERY verbs are deliberately NOT wrapped — their −1 is legitimately rank-local (`ladrunoMortarTieResidual` is 0 on the non-owner by design). | [#742](https://github.com/nmorabowen/OpenSees/pull/742) |
