---
wp: ADR-43
title: "ADR43 P3b -- 2 vanilla row(s)"
files: ["`SRC/interpreter/OpenSeesCommands.cpp`", "`SRC/system_of_eqn/eigenSOE/CMakeLists.txt`"]
table: "main"
legacy_seq: [263, 264]
---
| `SRC/interpreter/OpenSeesCommands.cpp` | `// Ladruno` ADR43 P3b: `eigen -feast ... -blockZGate` flag — parses the token and stashes it on the SOE (`setBlockZGate`); the solver then self-tests the block-real complex-shift kernel on the SOE's own CSR post-solve (see `LadrunoBlockZKernel` in [[LEDGER_implementations]]). | ADR43 P3b |
| `SRC/system_of_eqn/eigenSOE/CMakeLists.txt` | `// Ladruno` ADR43 P3b: add `LadrunoBlockZKernel.{cpp,h}` to `OPS_SysOfEqn` target_sources (build wiring only). | ADR43 P3b |
