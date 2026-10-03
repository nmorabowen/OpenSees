---
wp: PR-500
title: "AnalysisLog — per-run logs/log.log convention"
pr: "#500"
status: "shipped"
section: "table"
legacy_seq: 102
---
| **AnalysisLog — per-run `logs/log.log` convention** — every analysis run on the fork leaves a self-contained log: `with AnalysisLog(ops) as log:` creates a `logs/` folder next to the run and tees ALL terminal output into `logs/log.log` — the Python side (prints, tracebacks) via a stdout/stderr tee, the C++ side (`opserr` warnings / convergence failures) via the interpreter's own `logFile <file> -append` on the SAME file (console echo unchanged). Writes a **model-info block** (ndm/ndf, node/element counts, per-type element breakdown via `getEleClassTags`+`eleType`, equation count) and a **run-time summary** (named `log.timed(label)` phases, `log.analyze(...)` timed-passthrough wall time + non-zero-return count, total wall time, OK/FAILED status). An exception inside the block lands in the log with full traceback before propagating — a crashed run leaves a diagnosable log, not a truncated one. Pure stdlib + an `ops` handle, no rebuild; Tcl runs get the opserr half natively with `file mkdir logs; logFile logs/log.log`. In-module self-test (7 checks against the produced log, incl. a provoked opserr warning proving the C++ tee). | Tooling (Python driver, script-only) | — | `Ladruno_scripts/ladruno_logs.py` | shipped | [#500](https://github.com/nmorabowen/OpenSees/pull/500) |
