---
wp: WP-231
title: "WP-231 #1073 -- upstreamable-table row(s)"
pr: "#1073"
date: 2026-11-13
files: ["`SRC/material/nD/UWmaterials/ManzariDafalias.h`"]
table: "upstreamable"
---
| `SRC/material/nD/UWmaterials/ManzariDafalias.h` | `// Ladruno WP-231` — **The WP-129 SAS-ME member declarations take `double*` (6-component stack arrays) instead of `Vector&`.** `ladrunoSasStage`, `ladrunoSasSubsteps`, `ladrunoSasDrift`, `ladrunoSasProject`, `ladrunoSasBracketH`, `ladrunoSasIntersect` change signature; `ladrunoSasElastic` returns through an out-array; new declarations `ladrunoSasAlphaRatio(const double*, ...)`, `ladrunoSasF`, `ladrunoSasNormal`, `ladrunoSasModuli`, `ladrunoSasStateDep` (stack-array twins of `GetF`, `GetNormalToYield`, `GetElasticModuli`, `GetStateDependent`). All are defined in the fork file `LadrunoSANISANDSasME.cpp`; `ladrunoSasAlphaRatio` keeps its `Vector` overload for the replay command. Every changed line is a Ladruno line; upstream declarations and `ManzariDafalias.cpp` are NOT touched. | WP-231 [#1073](https://github.com/nmorabowen/OpenSees/pull/1073) |
