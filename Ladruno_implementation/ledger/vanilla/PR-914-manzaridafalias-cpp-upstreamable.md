---
wp: PR-914
title: "914 -- upstreamable-table row(s)"
pr: "#914"
files: ["`SRC/material/nD/UWmaterials/ManzariDafalias.cpp`"]
table: "upstreamable"
legacy_seq: [675]
---
| `SRC/material/nD/UWmaterials/ManzariDafalias.cpp` | `// Ladruno WP-160` — **`MaxStrainInc` (IntScheme 7, 8, 9) and `MaxEnergyInc` (0, 4, 6): the sub-steps' moduli and elastic strain.** Both declared `double nDGamma, nVoidRatio, nG, nK;` and passed the uninitialised `nG, nK` as every sub-step's moduli; both loops never advanced `cEStrain`. Now `nG = G, nK = K` (MaxStrainInc) / the entry `G, K` saved before the full-increment call and restored for each half (MaxEnergyInc: RungeKutta4 writes them), and `cEStrain = nEStrain` in both loops. Dispatch unchanged. Gated by `tests/test_manzari_substep_moduli.py` (6 regression gates fail on d63f49750). WP-129 byte-identity: `ls3d_s4, s6, s7, s8, s9` re-pinned, `NONDETERMINISTIC` emptied. | [#914](https://github.com/nmorabowen/OpenSees/pull/914) |
