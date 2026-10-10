---
wp: WP-152
title: "WP-152 #894 -- upstreamable-table row(s)"
pr: "#894"
files: ["`SRC/material/nD/UWmaterials/ManzariDafalias.h`"]
table: "upstreamable"
legacy_seq: [658]
---
| `SRC/material/nD/UWmaterials/ManzariDafalias.h` | `// Ladruno WP-152` — **Tension-cutoff options and state in the WP-129 SAS-ME block.** `LadrunoSasOptions` gains `tcPsep`, `tcPcontact` and `tcP0Max` (default 0 = OFF). `LadrunoSasState` gains `sep`/`sep_n`, `sepTr`/`sepTr_n`, and the trial-only `sepEvent`, `sepCode`, `sepP0` (default false/0). Eight census columns are APPENDED to `LSAS_*` before `LSAS_COUNT`; earlier indices are unchanged. Two new member declarations, `ladrunoSasSetIsotropic` and `ladrunoResetSasSep`, are defined in the fork file `LadrunoSANISANDSasME.cpp`. Additive and default-constructed, with every added line marked. `ManzariDafalias.cpp` is NOT touched. | WP-152 [#894](https://github.com/nmorabowen/OpenSees/pull/894) |
