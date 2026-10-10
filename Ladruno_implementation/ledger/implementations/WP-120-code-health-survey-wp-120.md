---
wp: WP-120
title: "Code-health survey (WP-120)"
pr: "#855"
status: "ready (owner merges)"
section: "table"
legacy_seq: 18
---
| **Code-health survey (WP-120)** ([[120_code_health_survey]]) — read-only Phase-0 measurement of duplication and dead code in fork-authored (stamped) sources; no production code changed. Tooling in `Ladruno_implementation/wp120_code_health/` (dependency-free, reuses the quirk lint's C++ scanner): `clones.py` (CPD-style token-window clone finder, exact + renamed-identifier modes, fork/vanilla/cross scopes), `history.py` (clone families × first-parent PR history × `LEDGER_quirks`), `deadcode.py`, `inventory.py` (fork-added per git vs stamp vs GLOBS). Headline: fork 14.5 % exact / 26.5 % renamed duplicated lines vs vanilla 38.3 % / 62.1 %; 7.2 % of fork lines copied from vanilla; trend flat (~12.6–12.8 % once WP-116's scope change is removed); fork has no dead guards / `#if 0`, 37 never-referenced functions; **31 fork-added files lack the stamp** (invisible to the quirk lint). Refactor candidates only where a replicated defect is proven: coupling/embedded `getDamp` trio (#219→#220), continuum element shells (8 replicated-fix PRs), explicit-integrator guards. Acceptance: the clone finder groups the #562 plane family on the pre-fix tree `4a975edee`. | tooling / survey | — | `Ladruno_implementation/120_code_health_survey.md`, `Ladruno_implementation/wp120_code_health/*` | **ready (owner merges)** | #855 |
