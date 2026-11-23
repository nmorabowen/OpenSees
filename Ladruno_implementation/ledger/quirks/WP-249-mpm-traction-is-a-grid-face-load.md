---
wp: WP-249
title: "`mpm traction` loads a grid boundary face row, not the material surface"
pr: "#1143"
date: 2026-11-23
---
### `mpm traction ymax …` on a column puts no load on the soil (or loads the wrong row)

**Symptom.** A surcharge written as `mpm traction <name> ymax <lo> <hi> <t>` on a soil column
whose grid extends above the soil (232k S2: 21 cells for a 10 m column, h 0.5) changes nothing:
the soil sees no surcharge.

**Why.** 232l's traction (interface P10b; `LadrunoMPMTwoPhase.cpp` `prepareStepUP`,
`routeSLoads`) puts consistent nodal forces λ t ∫hat_a on the nodes of a **grid boundary** row,
`row = 0` or `G.n[d]` of the named face. It does not follow a material surface. When the soil top
lies below the grid's top face, that row carries no mass, so in mode A (which loads only active
nodes, `prepareStepUP`) it is not an active node and receives no force. Mode F's `routeSLoads`
loads every node of the face row with no activity test; what a massless row does there was not
checked. The traction is also refused on `bspline2` / `bspline3` ("GIMP / linear
only in phase 6").

**Use.** Put the grid face at the loaded surface (S2: `-cells 2 20`, the grid top at y = l0). The
load is scaled by the host's λ; to ramp it apart from gravity, apply gravity as the unscaled
`mpm particleBodyForce` (B11). On S2 (mode A, gimp) this gives σ_yy max / ρgl0 = 0.108 (HS)
and 0.033 (SANISAND) against S2's 0.03, from 0.27 with the top-row body force. The remaining
error sits in the top cell. WP-249 record: `Ladruno_scripts/mpm_249/s2_traction_probe_results.json`;
the open part (B-spline, a non-grid surface, the top-cell residual) stays on #1110.
