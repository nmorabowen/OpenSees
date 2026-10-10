---
wp: PR-385
title: "385 (D2.2 part #387) -- 1 vanilla row(s)"
pr: "#385, #387"
files: ["`SRC/interpreter/OpenSeesOutputCommands.cpp`"]
table: "main"
legacy_seq: [61]
---
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | `// Ladruno` ADR-41 D2: the `contact` parse gains `-visc <μ_c>` (viscous normal-stabilization coefficient; 0 default ⇒ off, byte-identical) → threaded into `addContact(..., muc)` (NTS, D2.1) AND `addMortarContact(..., muc)` (mortar CONTACT, D2.2); REFUSED with `-tie` (a bond has no contact-chatter regime). `contactPlane` gains an optional trailing `-visc <μ_c>` parse loop → `addRigidPlane(..., muc)`. Additive; existing `contact`/`contactPlane` invocations unchanged (muc defaults to 0). | [#385](https://github.com/nmorabowen/OpenSees/pull/385) (D2.2 part [#387](https://github.com/nmorabowen/OpenSees/pull/387)) |
