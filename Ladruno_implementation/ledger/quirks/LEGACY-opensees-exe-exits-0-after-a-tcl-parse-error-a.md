---
wp: LEGACY
title: "OpenSees.exe exits 0 after a Tcl parse error — a deck that never ran looks like a deck that passed"
legacy_seq: 235
---
### `OpenSees.exe` exits 0 after a Tcl **parse** error — a deck that never ran looks like a deck that passed
- **Bites:** any CI/harness that gates on the process exit status. A brace mismatch (or any Tcl syntax error) aborts the script with `missing close-brace` on stderr, the deck's own `exit 1` on the failure path is never reached, and the process still exits **0**. Hit live while extending `Ladruno_implementation/lapack_singular_regression/` — the run printed its banner, printed nothing else, and reported success. A harness would have recorded a green run for a deck that executed zero assertions.
- **Why:** the Tcl interpreter reports the error and returns; `tclMain` does not translate a script error into a nonzero process status. Same family as the two exit-status traps already recorded here (`analyze()` returning 0 on a NaN field; `mpiexec` returning 0 when every rank died).
- **Workaround/status (2026-07-25):** never gate solely on the exit code. Gate on a **positive terminal marker** the deck prints only on the success path (e.g. grep for a final `=== ... all checks passed ===` line), or count the expected number of PASS lines. Both the ADR-76 smoke and the LAPACK regression print such a marker; the checker should require it, not merely tolerate its absence.
