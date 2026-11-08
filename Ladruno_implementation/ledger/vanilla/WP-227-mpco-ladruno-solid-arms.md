---
wp: WP-227
title: "MPCO recorder: LadrunoBrick / LadrunoQuad / LadrunoCST integration-rule arms"
date: 2026-11-07
pr: "TBD"
table: "main"
files: ["`SRC/recorder/MPCORecorder.cpp`"]
---
| `SRC/recorder/MPCORecorder.cpp` | `// Ladruno (WP-227)`: `getGeometryAndIntRuleByClassTag` adds `ELE_TAG_LadrunoBrick` to the Hex GL2 arm (Brick's i,j,k order, every `-formulation`), `ELE_TAG_LadrunoQuad` to the Quad GL2 arm (FourNodeQuad CCW order) and `ELE_TAG_LadrunoCST` to the Tri GL1 arm. Same orders as the Ladruno recorder (WP-226). Before, these classes were written as custom geometry with no rule. Three lines; other elements are byte-identical. | TBD |
