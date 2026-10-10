---
wp: ADR-43
title: "ADR43 P3c-serial -- 1 vanilla row(s)"
files: ["`SRC/interpreter/OpenSeesCommands.cpp`"]
table: "main"
legacy_seq: [265]
---
| `SRC/interpreter/OpenSeesCommands.cpp` | `// Ladruno` ADR43 P3c: `eigen -feast ... -rci` flag — parses the token and stashes it on the SOE (`setRci`); the solver then drives the FEAST solve through the `dfeast_srci` RCI with `LadrunoBlockZKernel` as the inner contour solve instead of the packed `dfeast_scsrgv` driver (see [[LEDGER_implementations]]). | ADR43 P3c-serial |
