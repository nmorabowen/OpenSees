---
wp: PR-874
title: "874 -- upstreamable-table row(s)"
pr: "#874"
files: ["`SRC/material/nD/soil/PressureDependMultiYield.cpp`, `SRC/material/nD/soil/PressureDependMultiYield02.cpp`, `SRC/material/nD/soil/PressureDependMultiYield03.cpp`"]
table: "upstreamable"
legacy_seq: [665]
---
| `SRC/material/nD/soil/PressureDependMultiYield.cpp`, `SRC/material/nD/soil/PressureDependMultiYield02.cpp`, `SRC/material/nD/soil/PressureDependMultiYield03.cpp` | `// Ladruno WP-135` — **bound the substep count.** `setSubStrainRate()` sized `getStress()`'s substep loop as `|Δε_oct|/1e-5` (PDMY01: `1e-4`) or `Δε_v/1e-5` with no cap, converted to `int` (undefined past `INT_MAX`); a wild Newton iterate asked for ~1e6–1e9 substeps per Gauss point per call and `analyze()` never returned (see [[LEDGER_quirks]]). Adds, per file: `#include <LadrunoMaterialStatus.h>`, a file-static `ladrunoPdmyTooManySubIncre()` that evaluates the SAME two expressions against a cap of 1e5 substeps; `setTrialStrain`/`setTrialStrainIncr` return `LADRUNO_MATERIAL_REFUSED` (with an `opserr` line) when it trips in stage 1, so a propagating host fails the step and the analysis can cut it; `setSubStrainRate` caps the loop at 1e5 so a host that swallows the code (`SSPquad`, `Brick`, …) still returns in bounded time. Strictly additive (29 lines per file, 0 deleted); below the cap nothing on the vanilla path changes: byte-identity gates `tests/test_wp135_pdmy_substep_cap.py` B1 (PDMY01/02, baseline from the unmodified build) and WP-133 G1 (PDMY03). | [#874](https://github.com/nmorabowen/OpenSees/pull/874) |
