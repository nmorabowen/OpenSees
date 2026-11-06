---
wp: PR-377
title: "377 -- 1 vanilla row(s)"
pr: "#377"
files: ["`SRC/interpreter/OpenSeesOutputCommands.cpp`"]
table: "main"
legacy_seq: [57]
---
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | `// Ladruno` ADR-41 C3.1: the `-mortar` contact parse gains `-mu <v>` / `-epsT auto\|<v>` / `-cohesion <v>` / `-tauMax <v>` (Coulomb/Tresca friction on the unified cone min(μN+c, τmax); all ≤0 ⇒ the frictionless C2 path) → threaded into `addMortarContact`. `-consistanttan` (already parsed for NTS) is reused for the mortar friction tangent (C3.2). | [#377](https://github.com/nmorabowen/OpenSees/pull/377) |
