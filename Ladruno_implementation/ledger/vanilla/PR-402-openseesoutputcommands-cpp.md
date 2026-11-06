---
wp: PR-402
title: "402 -- 1 vanilla row(s)"
pr: "#402"
files: ["`SRC/interpreter/OpenSeesOutputCommands.cpp`"]
table: "main"
legacy_seq: [63]
---
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | `// Ladruno` ADR-39 B1 (P4): the NTS `contact` parse gains `-soft <SOFSCL>` (optional value, default 0.10; the LS-DYNA §26.15 SOFT=1 Courant-stable explicit penalty `k_soft=SOFSCL·4·m_eff/dt²`) → threaded into `addContact(..., softScale)`; `contactPlane` gains the same optional `-soft <SOFSCL>` parse loop → `addRigidPlane(..., softScale)`. Optional value via the safe peek-`OPS_GetString`-then-`OPS_ResetCurrentInputArg(-1)` idiom (string read consumes on both Tcl + Py — dodges the failed-numeric-rewind quirk). REFUSED with `-mortar` (NTS-only) and without a base penalty (a modifier — needs a positional `auto`/`kn`; that base kn is what an implicit run uses). `>1` warns (`ω·dt=2√SOFSCL>2` unstable). Off (default) ⇒ byte-identical; explicit-only ⇒ implicit byte-identical. Additive; existing invocations unchanged. | [#402](https://github.com/nmorabowen/OpenSees/pull/402) |
