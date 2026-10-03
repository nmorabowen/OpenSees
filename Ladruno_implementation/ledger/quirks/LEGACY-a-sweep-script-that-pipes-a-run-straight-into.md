---
wp: LEGACY
title: "A sweep script that pipes a run straight into grep can hide the reason it failed — including from itself"
legacy_seq: 227
---
### A sweep script that pipes a run straight into `grep` can hide the reason it failed — including from itself
- **Bites:** the ADR-75 P2h sweep piped each run's output into `grep -E "P2H_RESULT|...|ERROR|rror"`. When `srun` failed with `srun: fatal: ...` / `command not found`, **none of the filter's patterns matched**, so the job produced a clean-looking log, printed its `..._DONE` banner and exited 0 — with zero results and zero explanation. Two takes were burned before the cause was visible. Compounding it: `mpiexec`/`mpirun`/`srun` wrappers **return 0 even when every rank died**, so `rc` proved nothing either.
- **Workaround/status:** write each run's FULL output to its own file, then `grep` the **file**; and gate on an **artifact** (`if grep -q P2H_RESULT "$LOG"`), dumping the log head when the artifact is missing. This is the banked "never let a bench script grep away its own log" rule — it was violated in the very sweep meant to honour it, which is why it is re-banked here with the concrete symptom. *2026-07-26 (ADR-75 P2h).*
