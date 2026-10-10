---
wp: PR-354
title: "354 -- 1 vanilla row(s)"
pr: "#354"
files: ["`SRC/interpreter/OpenSeesOutputCommands.cpp`"]
table: "main"
legacy_seq: [53]
---
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | `// Ladruno` ADR-39 P2b: extend the `contact` command — optional `-outward ox oy oz` (segment-normal orientation direction) parsed after the optional `kn kt mu` triple (peek/un-read via `OPS_ResetCurrentInputArg(-1)` keeps `kn kt mu` positional + flag-disambiguated). | [#354](https://github.com/nmorabowen/OpenSees/pull/354) |
