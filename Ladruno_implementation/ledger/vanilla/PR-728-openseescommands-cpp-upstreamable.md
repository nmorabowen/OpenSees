---
wp: PR-728
title: "728 -- upstreamable-table row(s)"
pr: "#728"
files: ["`SRC/interpreter/OpenSeesCommands.cpp`"]
table: "upstreamable"
legacy_seq: [422]
---
| `SRC/interpreter/OpenSeesCommands.cpp` | `// Ladruno` (#728, extends the #243 block): third `ladrunoDR` subcommand `stabilityMargin` → `LadrunoDynamicRelaxation::getStabilityMargin()`, returning `(ω_max·dt/2)²` for the fictitious mass actually in use, measured against the tangent at the last `M*` rebuild (`≤1` stable, `== massSafety²` on an unchanged tangent, `>1` = marching at/over the explicit boundary where DR can relax to a SILENTLY wrong state, `−1` = no gershgorin bound). Also widens the two usage/`unknown subcommand` strings. Additive — no existing subcommand changed. | [#728](https://github.com/nmorabowen/OpenSees/pull/728) |
