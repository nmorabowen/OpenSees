---
wp: WP-229
title: "WP-229 -- ManzariDafalias b0 pressure exponent (LadrunoSANISAND -sasPlasticExp)"
pr: "#1071"
date: 2026-11-11
files: ["`SRC/material/nD/UWmaterials/ManzariDafalias.h`", "`SRC/material/nD/UWmaterials/ManzariDafalias.cpp`"]
---
| `SRC/material/nD/UWmaterials/ManzariDafalias.h` | `// Ladruno WP-229` — `LadrunoSasOptions` gains `plasticExp` (default 0.5, set in the struct's own constructor; no `ManzariDafalias` constructor is touched). Additive, every added line marked. | WP-229 [#1071](https://github.com/nmorabowen/OpenSees/pull/1071) |
| `SRC/material/nD/UWmaterials/ManzariDafalias.cpp` | `// Ladruno WP-229` — `GetStateDependent`: after the vanilla line `b0 = G0 h0 (1 - ch e) / sqrt(p / P_atm)`, which is left as it is, two added lines replace `b0` by `G0 h0 (1 - ch e) (p / P_atm)^-nh` when `mLadrunoSas.opt.plasticExp != 0.5`. At 0.5 (the default, and every vanilla `ManzariDafalias` instance) the added branch is not taken, so the result is bit-identical. Only `LadrunoSANISAND` can set the option, and its parser refuses it for any scheme but IntScheme 129. | WP-229 [#1071](https://github.com/nmorabowen/OpenSees/pull/1071) |
