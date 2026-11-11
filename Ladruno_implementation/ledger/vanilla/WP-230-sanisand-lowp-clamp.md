---
wp: WP-230
title: "WP-230 -- ManzariDafalias.h: low-pressure clamp options and trial measures in the SAS-ME block"
date: 2026-11-11
files: ["`SRC/material/nD/UWmaterials/ManzariDafalias.h`"]
---
| `SRC/material/nD/UWmaterials/ManzariDafalias.h` | `// Ladruno WP-230` — `LadrunoSasOptions` gains `lowPMode`, `lowPFloor`, `lowPAlpha` (default 0 = WP-152's separation). `LadrunoSasState` gains the trial-only `clampFr`, `clampBound`. Two census columns are APPENDED to `LSAS_*` before `LSAS_COUNT`; earlier indices are unchanged. One member declaration, `ladrunoSasSetClamped`, defined in the fork file `LadrunoSANISANDSasME.cpp`. Additive and default-constructed, every added line marked. `ManzariDafalias.cpp` is NOT touched. | WP-230 |
