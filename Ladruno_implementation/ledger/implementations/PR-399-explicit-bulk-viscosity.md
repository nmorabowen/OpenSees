---
wp: PR-399
title: "Explicit bulk viscosity"
pr: "#399, #403"
status: "shipped"
section: "table"
legacy_seq: 92
---
| **Explicit bulk viscosity** (`-bulkViscosity b1 b2`) ([[52_ladruno_integrator_strengthening_adr]] W2-E1) — LS-DYNA/Abaqus-style artificial bulk-viscosity pressure (linear `b1` + quadratic `b2` on the volumetric strain rate) on all three Ladruno continuum elements, to damp shock/impact ringing under explicit integration. Off by default ⇒ byte-identical. | Element flag | — (no new tag) | `SRC/element/ladrunoBrick/`, `SRC/element/ladrunoPlane/` (LadrunoQuad + LadrunoCST) | shipped | [#399](https://github.com/nmorabowen/OpenSees/pull/399) (LadrunoBrick), [#403](https://github.com/nmorabowen/OpenSees/pull/403) (LadrunoQuad/CST) |
