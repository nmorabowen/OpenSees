---
wp: PR-871
title: "871 -- upstreamable-table row(s)"
pr: "#871"
files: ["`SRC/material/nD/UWmaterials/ManzariDafalias.h`", "`SRC/material/nD/UWmaterials/ManzariDafalias.cpp`", "`SRC/material/nD/CMakeLists.txt`"]
table: "upstreamable"
legacy_seq: [654, 655, 656]
---
| `SRC/material/nD/UWmaterials/ManzariDafalias.h` | `// Ladruno WP-129` — **SAS-ME (IntScheme 129) state block + declarations.** `#define LADRUNO_INT_SAS_ME 129`; structs `LadrunoSasOptions` (errFloor, alphaBoundTol, alphaProject, alphaInMode, errorVars) and `LadrunoSasState` (`allowed` seam, options, `stats[LSAS_COUNT]` census, per-update `refused`), both default-constructed so NO constructor is touched; enum `LSAS_*`; member `mLadrunoSas`; declarations of the scheme's member functions, which are DEFINED in the fork file `SRC/material/nD/LadrunoSANISANDSasME.cpp`. Additive; every added line marked. | [#871](https://github.com/nmorabowen/OpenSees/pull/871) |
| `SRC/material/nD/UWmaterials/ManzariDafalias.cpp` | `// Ladruno WP-129` — `integrate()`: resets `mLadrunoSas.refused` and ONE dispatch branch `else if (mLadrunoSas.allowed && mScheme == LADRUNO_INT_SAS_ME) ladrunoSasIntegrate();` between the CPPM branch and `explicit_integrator`. `allowed` is written only by `LadrunoSANISAND::applyLadrunoConstants()`, so on vanilla ManzariDafalias IntScheme 129 keeps its vanilla meaning (`explicit_integrator` default = ModifiedEuler). Every existing scheme byte-identical on 18 decks (`tests/wp129_sanisand_byteid_baseline.json`, captured on the unmodified WP-127 binary). | [#871](https://github.com/nmorabowen/OpenSees/pull/871) |
| `SRC/material/nD/CMakeLists.txt` | `LadrunoSANISANDSasME.cpp` added to `OPS_Material` sources (WP-129). | [#871](https://github.com/nmorabowen/OpenSees/pull/871) |
