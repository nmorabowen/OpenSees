---
wp: WP-232
title: "WP-232c: the explicit MPM host -- vanilla rows"
pr: "#1078"
date: 2026-11-13
table: "main"
files: ["`SRC/classTags.h`", "`SRC/domain/pattern/CMakeLists.txt`", "`SRC/interpreter/OpenSeesCommands.h`", "`SRC/tcl/commands.cpp`", "`SRC/interpreter/OpenSeesOutputCommands.cpp`", "`CMakeLists.txt`"]
---
| `SRC/classTags.h` | `PATTERN_TAG_LadrunoMPM 33027` (Ladruno comment on the define): the MPM host's PATTERN-registry tag. | WP-232c [#1078](https://github.com/nmorabowen/OpenSees/pull/1078) |
| `SRC/domain/pattern/CMakeLists.txt` | `# N. Mora-Bowen (Ladruno) WP-232`: `add_subdirectory(ladrunoMPM)`, the overlay's wiring. | WP-232c [#1078](https://github.com/nmorabowen/OpenSees/pull/1078) |
| `SRC/interpreter/OpenSeesCommands.h` | `// Ladruno WP-232`: declares `OPS_LadrunoMPM()` for the generated DL / Python bridges of the `mpm` table row. | WP-232c [#1078](https://github.com/nmorabowen/OpenSees/pull/1078) |
| `SRC/tcl/commands.cpp` | `// Ladruno WP-232`: the classic-Tcl bridge `ladrunoMPM` (OPS_ResetInputNoBuilder + the shared entry, refusing before `model` because `OPS_GetNDM` dereferences the builder); `{"MPM", LADRUNO_MUTATE_MPM}` in the `ladrunoMutation` table. | WP-232c [#1078](https://github.com/nmorabowen/OpenSees/pull/1078) |
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | `// Ladruno WP-232`: `{"MPM", LADRUNO_MUTATE_MPM}` in the `ladrunoMutation` reporting table. | WP-232c [#1078](https://github.com/nmorabowen/OpenSees/pull/1078) |
| `CMakeLists.txt` | `# Ladruno ... MPM by WP-232`: `MPM` in the mutation family STRINGS and `_ladruno_valid_families`. | WP-232c [#1078](https://github.com/nmorabowen/OpenSees/pull/1078) |
