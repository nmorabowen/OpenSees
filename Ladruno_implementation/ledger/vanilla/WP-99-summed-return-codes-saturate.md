---
wp: WP-99
title: "Summed material return codes saturate instead of wrapping past INT_MIN (apeGmsh #2619)"
date: 2026-10-30
table: "main"
files: ["`SRC/element/triangle/Tri31.cpp`", "`SRC/element/fourNodeQuad/FourNodeQuad.cpp`", "`SRC/domain/domain/Domain.cpp`"]
---
| `SRC/element/triangle/Tri31.cpp` | `// Ladruno (apeGmsh #2619)`: `update()` keeps the most negative `setTrialStrain` code instead of `ret +=`. With `LADRUNO_MATERIAL_REFUSED` = -33086 the model-wide sum wraps past INT_MIN at 64,907 refusing points and the refusal is accepted. | apeGmsh #2619 |
| `SRC/element/fourNodeQuad/FourNodeQuad.cpp` | `// Ladruno (apeGmsh #2619)`: same fold in `update()`, 4 points per element (16,227 fully refusing quads reach the wrap). | apeGmsh #2619 |
| `SRC/domain/domain/Domain.cpp` | `// Ladruno (apeGmsh #2619)`: file-static `ladrunoKeepWorstCode`; the serial `update()` loop folds instead of `ok +=`, and the WP-192 threaded loop keeps each body's code in a per-element slot and folds those, ignoring the integer SUM `ladruno_parallelFor` returns. The sentinel reaches the integrators as itself. | apeGmsh #2619 |
