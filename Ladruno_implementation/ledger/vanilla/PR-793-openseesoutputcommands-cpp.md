---
wp: PR-793
title: "793 -- 3 vanilla row(s)"
pr: "#793"
files: ["`SRC/interpreter/OpenSeesOutputCommands.cpp`", "`SRC/tcl/commands.cpp`", "`CMakeLists.txt`"]
table: "main"
legacy_seq: [160, 161, 162]
---
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | `// Ladruno` (ADR 91): add `SHELLMOD` to the `ladrunoMutation` family table AND derive the loop bound from `sizeof(fams)` -- it was a hardcoded `i < 5`, so a 6th family was silently unreported and a genuine mutant answered `none`. Found by running the gate. | 793 |
| `SRC/tcl/commands.cpp` | `// Ladruno` (ADR 91): same fix in the classic-Tcl twin of `ladrunoMutation` (the family table is duplicated across the two interpreters). | 793 |
| `CMakeLists.txt` | `# Ladruno` (ADR 91): add `SHELLMOD` to `_ladruno_valid_families` and the `LADRUNO_MUTATE_FAMILY` cache STRINGS property -- an unlisted family is a `FATAL_ERROR`, so the mutant configure would refuse outright. | 793 |
