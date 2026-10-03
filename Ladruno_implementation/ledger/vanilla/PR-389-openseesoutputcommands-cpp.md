---
wp: PR-389
title: "389 -- 1 vanilla row(s)"
pr: "#389"
files: ["`SRC/interpreter/OpenSeesOutputCommands.cpp`"]
table: "main"
legacy_seq: [62]
---
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | `// Ladruno` ADR-39 B3 (P2b-2c): the NTS `contact` parse gains `-geomtan` (opt into the consistent ∂n/∂u geometric NORMAL tangent `kn·gN·∂²gN/∂u²` ⇒ quadratic Newton on curved / large-sliding interfaces) → threaded into `addContact(..., consistentNormal)`. SYMMETRIC ⇒ no special solver (unlike `-consistanttan`). REFUSED with `-mortar` (the mortar geometric tangent is separately deferred). Off (default) ⇒ the shipped `kn·BᵀB` (byte-identical; EXACT for a flat/fixed master). Additive; existing invocations unchanged. | [#389](https://github.com/nmorabowen/OpenSees/pull/389) |
