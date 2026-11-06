---
wp: PR-306
title: "306 -- upstreamable-table row(s)"
pr: "#306"
files: ["`SRC/interpreter/PythonStream.h`"]
table: "upstreamable"
legacy_seq: [315]
---
| `SRC/interpreter/PythonStream.h` | `// Ladruno`: format-string bug in `err_out` — the already-formatted message was passed as the *format* to `PySys_FormatStderr(msg.c_str())`, so any literal `%` in an `opserr` message (e.g. `"% of model mass"`, `"exceeds -maxAddedMass cap 5%"`) was consumed as a bogus printf conversion and silently dropped/garbled under **openseespy** (Tcl `StandardStream` path unaffected). Fix: `PySys_FormatStderr("%s", msg.c_str())`. Surfaced by the SMS cap-warning validation (T-CAP). | [#306](https://github.com/nmorabowen/OpenSees/pull/306) |
