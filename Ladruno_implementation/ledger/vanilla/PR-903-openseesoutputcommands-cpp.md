---
wp: PR-903
title: "903 -- 1 vanilla row(s)"
pr: "#903"
files: ["`SRC/interpreter/OpenSeesOutputCommands.cpp`"]
table: "main"
legacy_seq: [34]
---
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | `// Ladruno ADR-159`: `OPS_LadrunoContact` gains the `-mortar`-only options `-smoothN <g0>` (> 0) and `-smoothT <r>` (0 < r < 1); refused without `-mortar` and with `-soft`/`-visc` after the loop, applied through `LadrunoContactDomain::setMortarSmoothing` (which refuses `-tie`, an augmenting contact and `-smoothT` without friction) only when given => byte-identical otherwise. | [#903](https://github.com/nmorabowen/OpenSees/pull/903) |
