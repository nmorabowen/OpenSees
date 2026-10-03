---
wp: PR-676
title: "676 -- upstreamable-table row(s)"
pr: "#676"
files: ["`SRC/element/Element.cpp`"]
table: "upstreamable"
legacy_seq: [402]
---
| `SRC/element/Element.cpp` | `index == -1` self-heal calls the **virtual** `setRayleighDampingFactors`, so any element overriding it leaves `index == -1` and the next line dereferences `theMatrices[-1]` — access violation in `getRayleighDampingForces`, `getResistingForceIncInertia`, the 5 sensitivity getters and `getGeometricTangentStiff`. One-qualifier fix (`this->Element::`) at 11 sites. **Not fork-specific — upstream's own `Subdomain` overrides the virtual as a forwarder to `Domain::` and hits it too.** | [#676](https://github.com/nmorabowen/OpenSees/pull/676) |
