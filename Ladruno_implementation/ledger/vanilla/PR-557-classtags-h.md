---
wp: PR-557
title: "557 -- 5 vanilla row(s)"
pr: "#557"
files: ["`SRC/classTags.h`", "`SRC/interpreter/OpenSeesElementCommands.cpp`", "`SRC/element/TclElementCommands.cpp`", "`SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`", "`SRC/element/CMakeLists.txt`"]
table: "main"
legacy_seq: [286, 287, 288, 289, 290]
---
| `SRC/classTags.h` | `// Ladruno` ADR-71 P1: promote the P0 placeholder comment (line ~943) to the real `#define ELE_TAG_LadrunoUP 33017` with the cross-registry occupants stamped (ND 33017 = LadrunoConcrete3D, ND 33018 = LadrunoRCFiniteStrain, LADRUNO 33019 = ComplexEigen — per-registry, not collisions); FIX the stale "33016-33019 remain free ELE slots" note inside the LadrunoSolidShell 33020 comment (now: 33016 reserved LadrunoLST/ADR-70, 33017 = LadrunoUP, 33018 = LadrunoBrick20/ADR-72, 33019 pencilled H27 — none free). | [#557](https://github.com/nmorabowen/OpenSees/pull/557) |
| `SRC/interpreter/OpenSeesElementCommands.cpp` | `// Ladruno` ADR-71 P1: `OPS_LadrunoUP` prototype + `functionMap` rows `LadrunoUP`/`ladrunoUP` (python/interpreter dispatch). Strictly additive. | [#557](https://github.com/nmorabowen/OpenSees/pull/557) |
| `SRC/element/TclElementCommands.cpp` | `// Ladruno` ADR-71 P1: `extern OPS_LadrunoUP` + one `ladrunoElementTable[]` row (`LadrunoUP`/`ladrunoUP`) in the fork-band chain slot (classic-Tcl dispatch; same OPS_ routine as the python path — no parallel parser). Strictly additive. | [#557](https://github.com/nmorabowen/OpenSees/pull/557) |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno` ADR-71 P1: `#include "ladrunoUP/LadrunoUP.h"` + `case ELE_TAG_LadrunoUP: return new LadrunoUP();` (MP broker reconstruction). Strictly additive. | [#557](https://github.com/nmorabowen/OpenSees/pull/557) |
| `SRC/element/CMakeLists.txt` | ADR-71 P1: refresh the `add_subdirectory(ladrunoUP)` comment — the P0 "header-only stub" note is stale now that the element + parser compile. Comment-only. | [#557](https://github.com/nmorabowen/OpenSees/pull/557) |
