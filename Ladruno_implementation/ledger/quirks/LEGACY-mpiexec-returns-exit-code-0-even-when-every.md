---
wp: LEGACY
title: "mpiexec returns exit code 0 even when every rank died in Tcl-init or hit a script error — the process rc is NOT the MP-harness failure signal"
legacy_seq: 210
---
### `mpiexec` returns exit code 0 even when every rank died in Tcl-init or hit a script `error` — the process rc is NOT the MP-harness failure signal
- **Bites:** an MP test/sweep harness that gates on `$LASTEXITCODE` / subprocess rc passes green while the run produced nothing — a missing `TCL_LIBRARY`, a bad `source`, a deck `error`, or a per-rank abort all leave `mpiexec` reporting success. Sibling to the serial `analyze() rc=0 on NaN` quirk above, but worse because it hides TOTAL failure, not just a bad answer.
- **Why:** the fork's Intel-MPI `mpiexec` propagates the launcher's exit status, not the ranks' Tcl interpreter status; OpenSees does not `MPI_Abort` with a nonzero code on a Tcl error.
- **Workaround/status:** MP harnesses must assert on **artifact existence + content** (dump-file count, expected line count, a sentinel `puts` like `TIEGATE_DONE`/`ANALYZE_MS`), never on rc. Build-tree exes additionally need `TCL_LIBRARY` exported (the packaged `openseesmp.sh` sets it; a raw `mpiexec … OpenSeesMP.exe` does not) or every rank dies in init — silently, rc=0. *2026-07-22 (banked across the ADR-74 rung/tie/checkpoint harnesses).*
