---
wp: WP-226
title: "A new solid element missing from the recorder's getGeometryAndIntRuleByClassTag writes no GP_PARAM; a class with varying node counts used to vanish from MODEL/ELEMENTS (WP-226)"
pr: "#1066"
date: 2026-11-07
---
### A new solid element missing from the recorder's `getGeometryAndIntRuleByClassTag` writes no `GP_PARAM`; a class with varying node counts used to vanish from MODEL/ELEMENTS (WP-226)
- **Bites (GP_PARAM):** `recorder ladruno` derives TOPOLOGY and the standard Gauss rule
  from the class-tag table in `SRC/recorder/Ladruno_ElementResults.h`. A class with no
  arm there gets `TOPOLOGY custom`, `INTEGRATION_RULE 0` and **no
  `QUADRATURE/GP_PARAM`**. Its `stresses` still record (8 × 6 columns), so nothing
  looks wrong until a viewer cannot place the Gauss points. `LadrunoBrick` (33002)
  was missing for its whole life. So were `LadrunoQuad` and `LadrunoCST`.
  `MPCORecorder.cpp` has the same gap (an upstream file, not fixed).
- **Rule:** a new element with standard Gauss points goes into the arm whose
  `getStandardQuadrature` table matches its GP ORDER, not just its count. Check the
  element's own loop against the table comment. An element with its own rule answers
  `basisInfo` / `integrationPoints` instead (BezierTri6 pattern).
- **Bites (variable nodes):** `mapElements` assumed one node count per class tag. RBE3 /
  RBE2 elements (one footprint each) printed one error per odd-sized element and were
  dropped from the file (815 RBE3s → ~4 000 lines, 0 written). Fixed in WP-226 with
  the schema §3.4 CSR layout. A reader of `MODEL/ELEMENTS` must handle a group with
  `CONNECTIVITY_LAYOUT="csr"` and no `CONNECTIVITY`.
