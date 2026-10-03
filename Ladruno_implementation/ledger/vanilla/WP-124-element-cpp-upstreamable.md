---
wp: WP-124
title: "859 -- upstreamable-table row(s)"
pr: "#859"
files: ["`SRC/element/Element.cpp`"]
table: "upstreamable"
legacy_seq: [666]
---
| `SRC/element/Element.cpp` | `// Ladruno (WP-124)` — `Element::getResponse` case `444444` (`inertialForce`): the one-expression `getResistingForceIncInertia() - getRayleighDampingForces() - getResistingForce()` becomes a fixed-order computation into an owned local copy. The call order of the three operands is unspecified in C++ and the accessors return references into element storage the other calls overwrite, so GCC recorded EXACTLY 0.0 for elements that fall back to the base vocabulary AND share that storage — LadrunoQuad/CST/LST/CSTPair (same `P` from both residual accessors) and LadrunoBrick/Brick20 (GRFII refills the `resid` `getResistingForce` returns); MSVC's order happened to be right. Found by the WP-124 Zone-A run (Ubuntu). Upstream bug — owner approved the vanilla fix 2026-09-27. | [#859](https://github.com/nmorabowen/OpenSees/pull/859) |
