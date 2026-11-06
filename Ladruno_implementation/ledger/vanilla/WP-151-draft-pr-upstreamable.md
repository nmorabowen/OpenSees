---
wp: WP-151
title: "WP-151 draft PR -- upstreamable-table row(s)"
files: ["`SRC/material/nD/UWmaterials/ManzariDafalias.h`"]
table: "upstreamable"
legacy_seq: [657]
---
| `SRC/material/nD/UWmaterials/ManzariDafalias.h` | `// Ladruno WP-151` — **R1 options in the WP-129 SAS-ME block.** `LadrunoSasOptions` gains `hFloor`, `reseatHyst` and `softCap` (default 0, default-constructed, so no constructor is touched). Three census columns are APPENDED to `LSAS_*` before `LSAS_COUNT`; earlier indices are unchanged. The `ladrunoSasBracketH` declaration takes `b0`, and two new members (`ladrunoSasSoftCapH`, `ladrunoSasReseatDelta`) are defined in the fork file `LadrunoSANISANDSasME.cpp`. Additive; every added line marked. `ManzariDafalias.cpp` is NOT touched. | WP-151 draft PR |
