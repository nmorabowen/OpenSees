---
wp: PR-358
title: "358 -- 1 vanilla row(s)"
pr: "#358"
files: ["`SRC/interpreter/OpenSeesOutputCommands.cpp`"]
table: "main"
legacy_seq: [55]
---
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | `// Ladruno` ADR-39 P2.5: the `contact` command gains an optional `-cell <frac>` flag (bucket-sort cell = frac·median seg diagonal; default 1.0; a huge value ⇒ 1 bucket = brute force) parsed in the trailing-option loop alongside `-outward`. | [#358](https://github.com/nmorabowen/OpenSees/pull/358) |
