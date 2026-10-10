---
wp: PR-357
title: "357 -- 1 vanilla row(s)"
pr: "#357"
files: ["`SRC/interpreter/OpenSeesOutputCommands.cpp`"]
table: "main"
legacy_seq: [54]
---
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | `// Ladruno` ADR-39 P2b-2b: the `contact` kn slot also accepts the literal `auto` (`contact tag m s auto [-outward …]`) → sets the `knAuto` flag so the handler auto-sizes kₙ from the master element stiffness; numeric `kn kt mu` path unchanged. | [#357](https://github.com/nmorabowen/OpenSees/pull/357) |
