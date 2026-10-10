---
wp: PR-674
title: "674 -- upstreamable-table row(s)"
pr: "#674"
files: ["`SRC/interpreter/OpenSeesCommands.cpp`, `SRC/tcl/commands.cpp`"]
table: "upstreamable"
legacy_seq: [400]
---
| `SRC/interpreter/OpenSeesCommands.cpp`, `SRC/tcl/commands.cpp` | `// Ladruno ADR-75 P1i`: **fix the profiler run-header size normalizers** — one `ops_profiler_fillModelMeta(meta, domain)` call added at each of the **four** `buildMeta()` sites (`profiler report` + `profiler checkpoint`, which exist **twice**: once in the Python ladder and once in the completely separate Tcl one — the banked "wiring one does not wire the other" trap, which is exactly how this drifted). Fills `nElem`/`nNode` (previously a hard 0 that a `Profiler.cpp` comment *claimed* the command layer populated — it only ever filled `nDOF`) and, under `-D_PARDISO`, overrides `threads` with `mkl_get_max_threads()`. Both TUs carry `-D_PARDISO` while `OPS_Utility` does not, which is why the MKL branch lives in a header these two include rather than in the profiler core. Reporting-only; no analysis behaviour changes. | [#674](https://github.com/nmorabowen/OpenSees/pull/674) |
