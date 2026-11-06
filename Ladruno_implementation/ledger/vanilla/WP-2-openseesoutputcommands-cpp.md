---
wp: WP-2
title: "WP-2 -- 6 vanilla row(s)"
files: ["`SRC/interpreter/OpenSeesOutputCommands.cpp`", "`SRC/interpreter/OpenSeesCommands.h`", "`SRC/interpreter/PythonWrapper.cpp`", "`SRC/interpreter/TclWrapper.cpp`", "`SRC/tcl/commands.cpp`", "`CMakeLists.txt` (root)"]
table: "main"
legacy_seq: [5, 6, 7, 8, 9, 10]
---
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | `// Ladruno: ADR-87 D2`: adds `OPS_LadrunoMutation()` — reports which physics family this binary was deliberately sabotaged for ("none" for a normal build). The mutation gate's wrong-build detector: without it a silently-failed mutant build would re-run the previous binary, every test would pass, and the gate would report the OPPOSITE of the truth. | WP-2 |
| `SRC/interpreter/OpenSeesCommands.h` | `// Ladruno: ADR-87 D2`: declares `OPS_LadrunoMutation()`. | WP-2 |
| `SRC/interpreter/PythonWrapper.cpp` | `// Ladruno ADR-87 D2`: registers the `ladrunoMutation` verb (Python). | WP-2 |
| `SRC/interpreter/TclWrapper.cpp` | `// Ladruno ADR-87 D2`: registers the `ladrunoMutation` verb (Tcl-in-OpenSeesPy interpreter). | WP-2 |
| `SRC/tcl/commands.cpp` | `// Ladruno ADR-87 D2`: classic-Tcl twin of `ladrunoMutation` + registration (the five-site registration rule). | WP-2 |
| `CMakeLists.txt` (root) | `# Ladruno ADR-87 D2`: `LADRUNO_MUTATE_FAMILY`/`LADRUNO_MUTATE_MODE` cache vars → one global `-DLADRUNO_MUTATE_<FAMILY>=<code>`. An unknown family or mode is a FATAL error, never a silently unmutated build. Empty family (default) adds no definition. | WP-2 |
