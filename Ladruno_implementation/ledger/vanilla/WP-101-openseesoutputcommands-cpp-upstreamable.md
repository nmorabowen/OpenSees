---
wp: WP-101
title: "839 -- upstreamable-table row(s)"
pr: "#839"
files: ["`SRC/interpreter/OpenSeesOutputCommands.cpp`", "`SRC/interpreter/OpenSeesCommands.cpp`", "`SRC/domain/domain/Domain.cpp`"]
table: "upstreamable"
legacy_seq: [619, 620, 621]
---
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | `// Ladruno (WP-101 r2)`: `OPS_LadrunoBeginAugment()` warns when an augmentation sweep is ALREADY open (a missed `ladrunoEndAugment`). While the flag is on `Domain::commit()` is recorder-silent, so the omission does not fail — it silently voids every later recorder sample. Warning only; the setter is unchanged, and the Tcl `ladrunoBeginAugment` routes through this same function so both interpreters are covered. | [#839](https://github.com/nmorabowen/OpenSees/pull/839) |
| `SRC/interpreter/OpenSeesCommands.cpp` | `// Ladruno (WP-101 r2)`: `OpenSeesCommands::wipeAnalysis()` clears `Domain::contactAugmenting`. An ADR-41 held-load augmentation sweep is analysis-scoped; a deck that forgets `ladrunoEndAugment` and rebuilds its analysis would otherwise carry the recorder silence into the new one. (`wipe` is already covered by `Domain::clearAll()`.) | [#839](https://github.com/nmorabowen/OpenSees/pull/839) |
| `SRC/domain/domain/Domain.cpp` | `// Ladruno (WP-101 r2)`: `Domain::clearAll()` resets `contactAugmenting = false`, next to the existing ADR-39 contact-engine teardown. Leaking an unclosed sweep across a `wipe` would make the NEXT model produce empty recorder files. | [#839](https://github.com/nmorabowen/OpenSees/pull/839) |
