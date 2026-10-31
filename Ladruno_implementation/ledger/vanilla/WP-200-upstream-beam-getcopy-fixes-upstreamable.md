---
wp: WP-200
title: "Upstream beam element fixes: GradientInelasticBeamColumn2d/3d pass 1.0 to getSectionLocations in constructor and refresh in setDomain; DispBeamColumnNL3d fixes 12-DOF localForce in responseID == 2 and un-nests responseID == 9 in getResponse"
date: 2026-10-18
pr: "#1000"
table: "upstreamable"
files: ["`SRC/element/gradientInelasticBeamColumn/GradientInelasticBeamColumn2d.cpp`", "`SRC/element/gradientInelasticBeamColumn/GradientInelasticBeamColumn3d.cpp`", "`SRC/element/dispBeamColumn/DispBeamColumnNL3d.cpp`"]
---
| `SRC/element/gradientInelasticBeamColumn/GradientInelasticBeamColumn2d.cpp` | `// Ladruno`: pass `1.0` instead of uninitialized `L` (`0.0`) to `beamIntegr->getSectionLocations(numSections, 1.0, secX)` in constructor and re-evaluate in `setDomain()`. | #1000 |
| `SRC/element/gradientInelasticBeamColumn/GradientInelasticBeamColumn3d.cpp` | `// Ladruno`: pass `1.0` instead of uninitialized `L` (`0.0`) to `beamIntegr->getSectionLocations(numSections, 1.0, secX)` in constructor, re-evaluate in `setDomain()`, and fix warning message string at line 95. | #1000 |
| `SRC/element/dispBeamColumn/DispBeamColumnNL3d.cpp` | `// Ladruno`: fix 2D copy-paste in `responseID == 2` (`localForce` 12-DOF 3D mapping) and un-nest `else if (responseID == 9)` in `DispBeamColumnNL3d::getResponse()` so it is not unreachable inside `else if (responseID == 10)`. | #1000 |
